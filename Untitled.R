library(Seurat)
library(SeuratData)
data("pbmc3k")
pbmc <- pbmc3k
pbmc = UpdateSeuratObject(object = pbmc)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top2000 <- head(VariableFeatures(pbmc), 2000)
pbmc <- pbmc[top2000]


print(pbmc) # Seurat object

library(reticulate)
library(sceasy)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

adata <- convertFormat(pbmc, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata) # Note generally in Python, dataset conventions are obs x var

# run setup_anndata
scvi$model$SCVI$setup_anndata(adata)

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

# to specify the number of epochs when training:
# model$train(max_epochs = as.integer(400))