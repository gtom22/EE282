import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scvi
import numba
import scipy.sparse as sp
from pynndescent import NNDescent
from scib_metrics.benchmark import Benchmarker


BASE         = "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts"
QUERY_H5AD   = f"{BASE}/fibro_scvi.h5ad"
REF_H5AD     = f"{BASE}/hca_ref.h5ad"   # downloaded from CellxGene portal

REF_MODEL_DIR    = f"{BASE}/scanvi_ref_model"
QUERY_MODEL_DIR  = f"{BASE}/scanvi_query_model"
OUT_QUERY        = f"{BASE}/fibro_query_labeled.h5ad"
OUT_COMBINED_EMB = f"{BASE}/fibro_combined_embedding.h5ad"
FIG_DIR          = f"{BASE}/figures_scarches"
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

LABEL_KEY = "author_cell_type"
#LABEL_KEY        = "cell_type"       # e.g. "cell_type", "cell_state", "Author_Annotation"

# Batch key — must exist in BOTH ref and query
REF_BATCH_KEY    = "donor_id"        # column in ref obs
QUERY_BATCH_KEY  = "donor_id"        # column in query obs

N_HVG            = 4000              # HVGs for the ref model (more = better for fibroblast subtypes)
N_LATENT         = 30
SCVI_EPOCHS      = 100               # ref training
SCANVI_EPOCHS    = 20                # ref fine-tuning (scANVI on top of scVI)
SURGERY_EPOCHS   = 500               # query surgery (with early stopping)
KNN_K            = 30                # neighbors for label transfer (HLCA uses 30)
UNCERTAINTY_THRESH = 0.2             # cells above this threshold → flagged as ambiguous
LEIDEN_RES       = 0.5               # slightly higher than before to get more clusters

@numba.njit
def weighted_prediction(weights, ref_cats):
    """For each query cell, compute weighted-majority-vote label and uncertainty."""
    N = len(weights)
    predictions = np.zeros(N, dtype=ref_cats.dtype)
    uncertainty  = np.zeros(N)
    for i in range(N):
        obs_weights = weights[i]
        obs_cats    = ref_cats[i]
        best_prob   = 0.0
        best_cat    = ref_cats[0, 0]
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            if cand_prob > best_prob:
                best_prob = cand_prob
                best_cat  = c
        predictions[i] = best_cat
        uncertainty[i]  = max(1.0 - best_prob, 0.0)
    return predictions, uncertainty


def transfer_labels(ref_obs, ref_neighbors, ref_distances, label_key):
    """Convert distances → affinities → weighted KNN prediction + uncertainty."""
    stds = np.std(ref_distances, axis=1, keepdims=True)
    stds = (2.0 / (stds + 1e-9)) ** 2
    dist_tilda = np.exp(-ref_distances * stds)
    weights    = dist_tilda / dist_tilda.sum(axis=1, keepdims=True)

    # encode categories as integers for numba
    ref_labels   = ref_obs[label_key].values
    cats, codes  = np.unique(ref_labels, return_inverse=True)
    ref_cats_idx = codes[ref_neighbors].astype(np.int32)

    pred_idx, uncert = weighted_prediction(weights, ref_cats_idx)
    predictions      = cats[pred_idx]
    return predictions, uncert


adata_ref = sc.read_h5ad(REF_H5AD)
adata_ref.var_names = adata_ref.var["feature_name"].astype(str).values  # ENSG → symbol
adata_ref.var_names_make_unique()


if "counts" not in adata_ref.layers:
    if adata_ref.raw is not None:
        raw_adata = adata_ref.raw.to_adata()
        raw_adata.var_names = raw_adata.var["feature_name"].astype(str).values
        raw_adata.var_names_make_unique()
        # subset raw to the genes in adata_ref
        common = raw_adata.var_names.intersection(adata_ref.var_names)
        raw_adata = raw_adata[:, common].copy()
        adata_ref  = adata_ref[:, common].copy()
        adata_ref.layers["counts"] = raw_adata.X.copy()
    else:
        adata_ref.layers["counts"] = adata_ref.X.copy()


if REF_BATCH_KEY not in adata_ref.obs.columns:
    adata_ref.obs[REF_BATCH_KEY] = adata_ref.obs.get(
        "sample", pd.Categorical(["ref_batch"] * adata_ref.n_obs))
    print(f"  WARNING: '{REF_BATCH_KEY}' not found in ref — created placeholder")


if "counts" not in adata_ref.layers:
    if adata_ref.raw is not None:
        adata_ref.layers["counts"] = adata_ref.raw.X.copy()
        # raw may have more genes — subset to adata_ref.var_names
        if adata_ref.raw.var_names.tolist() != adata_ref.var_names.tolist():
            raw_adata = adata_ref.raw.to_adata()
            raw_adata = raw_adata[:, adata_ref.var_names].copy()
            adata_ref.layers["counts"] = raw_adata.X.copy()
        print("  used raw.X for ref counts")
    else:
        # try to detect if X is counts (max value should be large integers)
        X_sample = adata_ref.X[:100]
        if sp.issparse(X_sample): X_sample = X_sample.toarray()
        if X_sample.max() > 20:   # likely counts
            adata_ref.layers["counts"] = adata_ref.X.copy()
            print("  WARNING: using adata_ref.X as counts — verify these are raw integers")
        else:
            raise ValueError("Reference does not appear to contain raw counts. "
                             "Download the raw count h5ad from CellxGene.")

Xc = adata_ref.layers["counts"]
d  = Xc.data if sp.issparse(Xc) else Xc.ravel()
adata_query = sc.read_h5ad(QUERY_H5AD)

common_genes = adata_ref.var_names.intersection(adata_query.var_names)
if len(common_genes) < 1000:
    raise ValueError(
        f"Only {len(common_genes)} shared genes. "
        "Check that both objects use the same gene ID format (symbol vs Ensembl). "
        "If ref uses Ensembl, convert query var_names or vice versa."
    )
adata_ref = adata_ref[:, common_genes].copy()
print(f"Selecting {N_HVG} HVGs on reference...", flush=True)

hvg_per_batch = []
for batch in adata_ref.obs[REF_BATCH_KEY].unique():
    mask = adata_ref.obs[REF_BATCH_KEY] == batch
    sub = adata_ref[mask].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, flavor="seurat")
    hvg_per_batch.append(set(sub.var_names[sub.var["highly_variable"]]))

# take genes that are HVG in at least 2 batches
from collections import Counter
gene_counts = Counter(g for s in hvg_per_batch for g in s)
hvg_genes = [g for g, c in gene_counts.items() if c >= 2]
hvg_genes = hvg_genes[:N_HVG]  # cap at N_HVG

adata_ref_hvg = adata_ref[:, hvg_genes].copy()
print(f"  Using {adata_ref_hvg.n_vars} HVGs")

arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
    n_latent=N_LATENT,
)

print("\nSetting up scVI on reference...", flush=True)
scvi.model.SCVI.setup_anndata(
    adata_ref_hvg,
    layer="counts",
    batch_key=REF_BATCH_KEY,
)

vae_ref = scvi.model.SCVI(adata_ref_hvg, **arches_params)
print(f"Training scVI on reference ({SCVI_EPOCHS} epochs)...", flush=True)
vae_ref.train(max_epochs=SCVI_EPOCHS, gradient_clip_val=10.0)

print("\nFine-tuning scANVI on reference...", flush=True)
vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",   # placeholder; all ref cells are labeled
    labels_key=LABEL_KEY,
)
vae_ref_scan.train(
    max_epochs=SCANVI_EPOCHS,
    n_samples_per_label=100,        # balances rare cell types during training
    gradient_clip_val=10.0,
)

# Store reference latent representation
adata_ref_hvg.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()

# Save reference model — needed for surgery
os.makedirs(REF_MODEL_DIR, exist_ok=True)
vae_ref_scan.save(REF_MODEL_DIR, overwrite=True)

scvi.model.SCANVI.prepare_query_anndata(adata_query, REF_MODEL_DIR)

vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    REF_MODEL_DIR,
)

vae_q.train(
    max_epochs=SURGERY_EPOCHS,
    plan_kwargs={"weight_decay": 0.0},      # CRITICAL: preserves ref latent space
    early_stopping=True,
    early_stopping_monitor="elbo_train",
    early_stopping_patience=10,
    early_stopping_min_delta=0.001,
)

# Get query latent representation (in the reference latent space)
adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation(adata_query)

os.makedirs(QUERY_MODEL_DIR, exist_ok=True)
vae_q.save(QUERY_MODEL_DIR, overwrite=True)
ref_emb   = adata_ref_hvg.obsm["X_scANVI"]
query_emb = adata_query.obsm["X_scANVI"]

# Build approximate NN index on the reference
knn_index = NNDescent(ref_emb, n_neighbors=KNN_K, metric="euclidean", random_state=42)
knn_index.prepare()

# Query: find K nearest reference neighbors for each query cell
ref_neighbors, ref_distances = knn_index.query(query_emb, k=KNN_K)

preds, uncerts = transfer_labels(adata_ref_hvg.obs, ref_neighbors, ref_distances, LABEL_KEY)

adata_query.obs[f"{LABEL_KEY}_transfer"]   = preds
adata_query.obs[f"{LABEL_KEY}_uncertainty"] = uncerts

# Flag ambiguous cells
adata_query.obs["transfer_confident"] = uncerts < UNCERTAINTY_THRESH

n_total      = adata_query.n_obs
n_confident  = (uncerts < UNCERTAINTY_THRESH).sum()
n_ambiguous  = n_total - n_confident
ref_emb_adata = anndata.AnnData(
    X=ref_emb,
    obs=adata_ref_hvg.obs[[REF_BATCH_KEY, LABEL_KEY]].copy()
)
ref_emb_adata.obs["ref_or_query"] = "reference"
ref_emb_adata.obs[f"{LABEL_KEY}_transfer"]   = ref_emb_adata.obs[LABEL_KEY]
ref_emb_adata.obs[f"{LABEL_KEY}_uncertainty"] = 0.0

query_emb_adata = anndata.AnnData(
    X=query_emb,
    obs=adata_query.obs[[QUERY_BATCH_KEY, f"{LABEL_KEY}_transfer",
                          f"{LABEL_KEY}_uncertainty", "transfer_confident"]].copy()
)
query_emb_adata.obs["ref_or_query"] = "query"
query_emb_adata.obs[LABEL_KEY]      = query_emb_adata.obs[f"{LABEL_KEY}_transfer"]

combined = anndata.concat(
    {"reference": ref_emb_adata, "query": query_emb_adata},
    label="dataset_origin", merge="unique"
)
combined.obsm["X_scANVI"] = combined.X
sc.pp.neighbors(combined, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(combined, min_dist=0.3)
sc.tl.leiden(combined, resolution=LEIDEN_RES, key_added=f"leiden_r{LEIDEN_RES}")

# ── UMAP plots
def save_umap(adata, color, title_suffix="", filename_suffix=""):
    sc.pl.umap(adata, color=color, title=f"{color}{title_suffix}",
               frameon=False, show=False,
               save=f"_{color.replace('.','_')}{filename_suffix}.png")

# Overview
sc.pl.umap(combined, color=["ref_or_query", LABEL_KEY],
           ncols=2, frameon=False, show=False, save="_overview.png")

# Reference labels vs transferred labels
sc.pl.umap(combined, color=f"{LABEL_KEY}_transfer",
           title="Transferred labels (all)", frameon=False, show=False,
           save="_transferred_all.png")
sc.pl.umap(combined, color=f"{LABEL_KEY}_uncertainty",
           color_map="RdYlGn_r", vmin=0, vmax=1,
           title="Transfer uncertainty (0=confident)", frameon=False, show=False,
           save="_uncertainty.png")

# High-confidence only (query)
q_hi = combined[(combined.obs["ref_or_query"] == "query") &
                (combined.obs[f"{LABEL_KEY}_uncertainty"] < UNCERTAINTY_THRESH)].copy()
sc.pl.umap(q_hi, color=f"{LABEL_KEY}_transfer", legend_loc="on data",
           title=f"Query — confident transfers (uncertainty < {UNCERTAINTY_THRESH})",
           frameon=False, show=False, save="_query_confident.png")
bm_adata = combined[
    (combined.obs["ref_or_query"] == "reference") |
    (combined.obs[f"{LABEL_KEY}_uncertainty"] < UNCERTAINTY_THRESH)
].copy()

bm = Benchmarker(
    bm_adata,
    batch_key=QUERY_BATCH_KEY,
    label_key=LABEL_KEY,
    embedding_obsm_keys=["X_scANVI"],
    n_jobs=-1,
)
bm.benchmark()
df_results = bm.get_results(min_max_scale=False)
print("\nscIB results:")
print(df_results.to_string())
df_results.to_csv(f"{BASE}/scib_scarches_results.csv")
bm.plot_results_table(min_max_scale=False, save_dir=FIG_DIR)


adata_query.write_h5ad(OUT_QUERY)
combined.write_h5ad(OUT_COMBINED_EMB)

