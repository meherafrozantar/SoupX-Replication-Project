# Increase the download timeout limit to prevent cloud server download drops
options(timeout = 600)

# 1. Automated Setup: Create folders and download data directly via code
dir.create("./data", showWarnings = FALSE)
tmp_dir <- tempdir()

message("Downloading raw data directly...")
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmp_dir, "tod.tar.gz"), quiet = TRUE)

message("Downloading filtered data directly...")
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmp_dir, "toc.tar.gz"), quiet = TRUE)

# 2. Extract files into your Posit Cloud workspace
untar(file.path(tmp_dir, "tod.tar.gz"), exdir = "./data")
untar(file.path(tmp_dir, "toc.tar.gz"), exdir = "./data")

# 3. Load mandatory single-cell analysis libraries
library(SoupX)
library(Seurat)

# 4. Load datasets into memory from the extracted path
tod <- Read10X("./data/raw_gene_bc_matrices/GRCh38/")
toc <- Read10X("./data/filtered_gene_bc_matrices/GRCh38/")
sc <- SoupChannel(tod, toc)

# 5. Fast clustering execution via Seurat
srat <- CreateSeuratObject(toc)
srat <- NormalizeData(srat, verbose=FALSE)
srat <- FindVariableFeatures(srat, verbose=FALSE)
srat <- ScaleData(srat, verbose=FALSE)
srat <- RunPCA(srat, verbose=FALSE)
srat <- FindNeighbors(srat, dims=1:10, verbose=FALSE)
srat <- FindClusters(srat, resolution=0.5, verbose=FALSE)

sc <- setClusters(sc, setNames(srat$seurat_clusters, colnames(toc)))

# 6. Apply target lineage baseline markers from the original paper
paper_markers <- c("HBB", "HBA1", "HBA2", "IGKC", "IGHG1", "IGHG2", "IGHG3")
sc <- setEstCont(sc, nonExpressedGeneList = list(PBMC = paper_markers))

# 7. Run final mathematical calculation fraction
sc <- autoEstCont(sc, forceAccept = TRUE)
replicated_rho <- sc$metaData$rho

# 8. Output value check
cat("\n======================================================\n")
cat("FINAL REPLICATED CONTAMINATION VALUE: ", round(replicated_rho * 100, 2), "%\n")
cat("======================================================\n")

