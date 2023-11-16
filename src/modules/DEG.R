library(edgeR)
library(limma)


degs <- function(adata_, filter_low_counts, include_time) {
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group=colData(adata_)$isHuman)

    # filter out genes with low counts
    if (filter_low_counts) {
        print("Dimensions before subsetting:")
        print(dim(y))
        print("")
        keep <- filterByExpr(y)
        y <- y[keep, , keep.lib.sizes=FALSE]
        print("Dimensions after subsetting:")
        print(dim(y))
        print("")
    }

    # normalize
    y <- calcNormFactors(y)
    # create a vector that is concatentation of condition and cell type that we will later use with contrasts
    group <- paste0(colData(adata_)$cluster)
    # replicate <- paste0(colData(adata_)$line, ".", colData(adata_)$replicate)
    isHuman <- paste0(colData(adata_)$isHuman)
    if (include_time) {
        isHuman <- paste0(colData(adata_)$isHuman, ".", colData(adata_)$time)
    }
    
    design <- model.matrix(~ 0 + isHuman + group)
    # estimate dispersion
    y <- estimateDisp(y, design=design)
    # fit the model
    fit <- glmQLFit(y, design)
    return(list("fit"=fit, "design"=design, "y"=y))
}


fit_model <- function(adata_) {
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group=colData(adata_)$isHuman)
    # filter out genes with low counts
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    # normalize
    y <- calcNormFactors(y)
    # create a vector that is concatentation of condition and cell type that we will later use with contrasts
    group <- paste0(colData(adata_)$time, ".", colData(adata_)$cluster)
    replicate <- paste0(colData(adata_)$line, ".", colData(adata_)$replicate)
    isHuman <- paste0(colData(adata_)$isHuman)#, ".", colData(adata_)$time)
    design <- model.matrix(~ 0 + isHuman)
    # estimate dispersion
    y <- estimateDisp(y, design=design)
    # fit the model
    fit <- glmQLFit(y, design)
    return(list("fit"=fit, "design"=design, "y"=y))
}