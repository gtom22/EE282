# Analysis Proposal

## Introduction

A common mutation that occurs in cancer is the BRCA1 mutation. However, its role in the tumor microenvironment is poorly misunderstood. Specifically in breast cancer, the presence of a BRCA1 mutation significantly increases a patient's susceptibility to the disease, and recent investigations have documented a notable prevalence of pre-cancer associated fibroblasts (pre-CAFs) within the stroma (Buechler et al. 2023). In terms of the immune infiltration, we see that BRCA1 wild-type patients experienced increased T-cell exhaustion and dendritic cell tolerance (Sun et. al 2025). But, we do not know if there are specific subtypes of T-cells or what other immune response (ie. mast cells, macrophages, etc.) are also involved. 

We propose a single-cell analysis using R tools such as Seurat, MiloR and Monocle on the integrated human breast cancer atlas (iHBCA). Using these, would address some of the main questions in BRCA1 carriers and will be broken down into 2 sections: 1) Fibroblast subclustering and 2) immune landscape classification. The main questions we want to address for the fibroblasts are: Do preCAFs activate markers of the DNA damage response? Do preCAFs activate markers of NF-kB and the NF-kB DDR axis? What cells are actively dividing? Do those correlate with high BRCA1 expression? Do those correlate with preCAF markers or stress markers such as DNA Damage response? For the immune landscape we want to answer: Do T cells express exhaustion markers in BRCA1 mutation carriers? MDSCs? Neutrophils? Macrophages etc? And what is their ISG/IFN expression in BRCA1 mutation carriers vs noncarriers? The dataset is publicly available [here](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637).

## Pertinence to Topic

To address the questions, Seurat is very powerful in zooming into the fibroblasts to analyze. Seurat allows you to create clusters of gene expression that we can correlate to cell types or cell states. We can then further keep breaking them down and make comparisons for specific sets of genes. Where are they located in this projected space? Are there hotspots of expression for specific genes within the clusters? We can then use Monocle to see if there's some sort of expression trajectory from different subtypes of fibroblasts. Once the final Seurat object that has been properly classified is made, we can observe the pseudotime which can tell us about cell division, correlations of fibroblasts going from pre-CAF to CAF and how expression changes to different cell states. 

We will perform a similar pipeline with the immune landscape, but we will also use MiloR to see which population of immune cells are more abundant in certain conditions (BRCA1 carrier vs control). We can then run differential gene expression analysis (DEG) to see which genes are more significant up/downregulated in each condition.

## Specific Analyses

### Seurat Single-cell Clustering

Seurat is an R package that is designed to explore single-cell transcriptomic space that allows for easy integration across datasets (Hao et al. 2024). The iHBCA dataset is already integrated, but it was created using Seurat so we will want to load that object as a starting point. Then we will want to subset the fibroblasts or immune cells depending on the part of the analysis and display the UMAP to observe the gene expression clusters for the cell type of interest. Then using the subset of genes that we are interested in, we can use this paper that classifies fibroblast cell states (Kumar et al. 2023) and create feature plots to annotate those clusters based on the marker genes. I will also create a heatmap showing the top 10 genes that highlight each cluster. Once I determine the pre-CAF cluster, I will subset those and look at the specific genes of interest and where in the cluster it is located while comparing it to other clusters. With the subclusters, we would compare conditions using volcano plots provided by Seurat highlighting the top 10 most differentially expressed genes in the BRCA carriers compared to the controls.

### Monocle Trajectory Analysis

Monocle is an algorithm that learns a cell’s sequence of gene expression changes that can cause it to go from one state to another creating a sort of pseudotime lineage where we can place cells along the developmental process (Trapnell et al. 2014).  We can tailor this to start at the pre-CAF cluster and see where it goes in hopes that there is a trajectory from pre-CAF to CAFs that can be found mostly in BRCA dense clusters. We will be creating a trajectory plot like this linking the different cell states together.

### MiloR

Milo is a package that is used to perform differential abundance testing among various experimental conditions. This will be done using the MiloR framework which takes the processed single cell object and creates “neighborhoods” based on expression patterns which can include different cell types/clusters. This essentially would create a global map within the UMAP connecting clusters together (Dann et al. 2022). We would split by immune cell state per immune cell breaking up the cluster using Seurat into subclusters. Then define the experimental design as BRCA carriers vs Controls and count how many cells are from each sample for each neighborhood to keep track of the variation in cell numbers between replicates for the same experimental condition (BRCA vs Control). We then would use the ‘calcNhoodDistance’ function to store the distances between nearby neighbors in the Milo object. With all of these pieces in place, we can perform DA testing using the ‘testNhoods’ function to test for differences between BRCA carriers and controls. From this test, Milo generates a fold-change and corrected p-value for each neighborhood indicating significant differential abundance. The resulting plots would be the global umap connecting neighborhoods across clusters and a beeswarm plot showing the log-fold change of abundance and significance coloring based on spatial FDR. The beeswarm plot would be generated using the ggplot2 package while the rest of the plots are built into the respective packages.


## Conclusion

This analysis plan would allow us to answer all of the questions proposed and give a comprehensive view of the fibroblasts and the immune landscape. A lot of the methods are repetitive so once I finalize the pipeline for one part, the second part would just involve changing thresholds and marker genes. I plan to subset down to make analysis easier, possibly just using 10k cells from the samples and subset out cell types I'm not interested in. The pipeline is well published so we would not be creating something new rather adapt existing pipelines for our application. It is feasible for computational resources as I have been granted permission to use my lab computational resources as the project is a part of one of my current projects. Large datasets can be accommodated with more CPUs. 


## References

Buechler MB, Fu W, Turley SJ. 2023. Healthy BRCA1/2 mutation carriers exhibit a pre-CAF signature and altered epithelial marker expression in breast tissue Scientific Reports 15(1):32736 doi:10.1038/s41598-025-18171-y

Dann, E., Henderson, N.C., Teichmann, S.A. et al. Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nat Biotechnol 40, 245–253 (2022). https://doi.org/10.1038/s41587-021-01033-z

Hao, Y., Stuart, T., Kowalski, M.H. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol 42, 293–304 (2024). https://doi.org/10.1038/s41587-023-01767-y

Kumar, T., Nee, K., Wei, R. et al. A spatially resolved single-cell genomic atlas of the adult human breast. Nature 620, 181–191 (2023). https://doi.org/10.1038/s41586-023-06252-9

Sun S, Chen S, Li K, Zhang G, Wang N, Xu Y, Wang X, Chi J, Li L and Sun Y (2025) Resolving tumor microenvironment heterogeneity to forecast immunotherapy response in triple-negative breast cancer through multi-scale analysis. Front. Oncol. 15:1538574. doi: 10.3389/fonc.2025.1538574

Trapnell C, Cacchiarelli D, Grimsby J, et al. 2014. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat Biotechnol 32:381–386. doi:10.1038/nbt.2859
