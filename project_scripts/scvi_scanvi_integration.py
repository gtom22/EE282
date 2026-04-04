import os
import tempfile

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
import scipy.sparse as sp
import scipy.io as sio


BASE = "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts"

MTX_PATH = f"{BASE}/RNA_counts.mtx"
GENES_PATH = f"{BASE}/RNA_genes.txt"
BARCODES_PATH = f"{BASE}/RNA_barcodes.txt"
OBS_PATH = f"{BASE}/RNA_obs.csv"

OUT_H5AD = f"{BASE}/brca_ref_fibroblasts_counts.h5ad"       
OUT_H5AD_SCVI = f"{BASE}/fibro_scvi.h5ad"                   
FIG_DIR = f"{BASE}/figures"

# scVI params
N_LAYERS = 2
N_LATENT = 30
MAX_EPOCHS = 200
LR = 1e-3
GRAD_CLIP = 10.0
NUM_WORKERS = 8  

# Leiden/UMAP params
LEIDEN_RES = 0.3
UMAP_MIN_DIST = 0.3


def sparse_minmax_nonneg(X, name="matrix"):
    if sp.issparse(X):
        d = X.data
        if d.size == 0:
            return 0.0, 0.0, True
        mn = float(d.min())
        mx = float(d.max())
        nonneg = not bool((d < 0).any())
        return mn, mx, nonneg
    else:
        arr = np.asarray(X)
        mn = float(arr.min())
        mx = float(arr.max())
        nonneg = bool((arr >= 0).all())
        return mn, mx, nonneg



def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    sc.settings.figdir = FIG_DIR
    

    try:
        X = sio.mmread(MTX_PATH).tocsr()  
        genes = pd.read_csv(GENES_PATH, header=None)[0].astype(str).tolist()
        barcodes = pd.read_csv(BARCODES_PATH, header=None)[0].astype(str).tolist()
        obs = pd.read_csv(OBS_PATH)
    except:
        scipy.io.mmwrite(f"{BASE}/RNA_counts.mtx", adata.X)
        # genes
        pd.DataFrame(adata.var_names).to_csv(f"{BASE}/RNA_genes.txt", index=False, header=False)

        # barcodes
        pd.DataFrame(adata.obs_names).to_csv(f"{BASE}/RNA_barcodes.txt", index=False, header=False)

        # obs metadata
        adata.obs.to_csv(f"{BASE}/OBS_PATH.csv", index=False)
        X = sio.mmread(MTX_PATH).tocsr()
        genes = pd.read_csv(GENES_PATH, header=None)[0].astype(str).tolist()
        barcodes = pd.read_csv(BARCODES_PATH, header=None)[0].astype(str).tolist()
        obs = pd.read_csv(OBS_PATH)

    # align obs to barcodes
    if "cell_barcode" not in obs.columns:
        raise ValueError(f"'cell_barcode' column not found in {OBS_PATH}. Columns: {list(obs.columns)}")
    obs = obs.set_index("cell_barcode")
    obs = obs.loc[barcodes]

    # shape checks
    if X.shape[0] != len(genes):
        raise ValueError(f"MTX rows (genes) {X.shape[0]} != len(genes) {len(genes)}")
    if X.shape[1] != len(barcodes):
        raise ValueError(f"MTX cols (cells) {X.shape[1]} != len(barcodes) {len(barcodes)}")
    if obs.shape[0] != len(barcodes):
        raise ValueError(f"obs rows {obs.shape[0]} != len(barcodes) {len(barcodes)}")

    # AnnData expects cells x genes
    adata = sc.AnnData(X.T, obs=obs)
    adata.var_names = genes
    adata.obs_names = barcodes

    # store true counts in a layer
    adata.layers["counts"] = adata.X.copy()
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
    if "counts" in adata.layers and sp.issparse(adata.layers["counts"]):
        adata.layers["counts"] = adata.layers["counts"].tocsr()
    # sanity check counts
    mn, mx, nonneg = sparse_minmax_nonneg(adata.layers["counts"], name="counts")
    print(f"counts layer: shape={adata.layers['counts'].shape} min={mn} max={mx} nonneg={nonneg}", flush=True)
    if not nonneg:
        raise ValueError("Counts layer has negative values. Something is wrong with the export.")

    adata.write_h5ad(OUT_H5AD)

    # scVI
    batch_key = "donor_id"
    print(f"Using batch_key='{batch_key}'", flush=True)

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)

    model = scvi.model.SCVI(adata, n_layers=N_LAYERS, n_latent=N_LATENT)

    print("Training scVI...", flush=True)
    if sp.issparse(adata.layers["counts"]):
        adata.layers["counts"] = adata.layers["counts"].tocsr()

    model.train(
        max_epochs=MAX_EPOCHS,
        gradient_clip_val=GRAD_CLIP,
    )
    print("Computing latent + neighbors/UMAP/Leiden...", flush=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, key_added=f"leiden_r{LEIDEN_RES}")

    # Save plots
    sc.pl.umap(adata, color=f"leiden_r{LEIDEN_RES}", show=False, save="_scvi_umap.png")
    sc.pl.umap(
        adata,
        color="donor_id",
        show=False,
        save="_donor_id.png"
    )   
    sc.pl.umap(
        adata,
        color="predicted.id",
        show=False,
        save="_cell_state.png"
    )
    sc.pl.umap(
        adata,
        color="group",
        show=False,
        save="_group.png"
    )
    
    adata_high = adata[adata.obs["prediction.score.max"] >= 0.8].copy()

    sc.pl.umap(
        adata_high,
        color="predicted.id",
        legend_loc="on data",
        title="High-confidence cells (score ≥ 0.8)"
    )
    # Save final
    print(f"Writing scVI AnnData: {OUT_H5AD_SCVI}", flush=True)
    adata.write_h5ad(OUT_H5AD_SCVI)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
