# The goal of this script is to consider a complete flow that takes an
# abundance matrix and the experiment design matrix (hence supervised) and
# performs differential analysis.
# The flows here can take care of the data modeling, which can include 
# sample-to-sample normalization, variance stabilization (heteroskedasticity),
# and empirical Bayes methodolody for shrinkage (e.g., of fold changes).

###############################################################
# Load required functions
required_libs = c("apeglm","preprocessCore",
                  "NOISeq","edgeR","limma","DESeq2",
                  # For PQC normalization use the Rcpm package
                  # devtools::install_github("ricoderks/Rcpm")
                  # Use the function pqc
                  "Rcpm",
                  # For loess normalization: use normalize.loess
                  # This package can also be used for cubic splines
                  "affy",
                  # For imputation using kNN
                  "impute"
                  )
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}
###############################################################

#'
#'@param design_mat is a model matrix of the experiment with the features/covariates to include 
#'@param x a numeric counts matrix with samples as the columns
edgeR_vool_limma_analysis<-function(x,design_mat){
  dge <- edgeR::DGEList(counts=x) # transform the counts to a dge object
  dge <- calcNormFactors(dge) # use the default method
  v<-limma::voom(dge,design=design_mat,plot=F)
  fit<-limma::lmFit(v,design_mat)
  fit<-limma::eBayes(fit)
  return(fit)
}

#' Run the complete DESeq2 pipeline
#' @count_matrix a matrix, can be before or after filtering lowly expressed genes.
#' @design_mat is a model matrix of the experiment with the features/covariates to include.
#' @ests_to_shrink are the indices of the coefficient for shrinkage, by default we use the first two (assuming there are >2)
deseq2_differential_analysis<-function(count_matrix,design_mat,ests_to_shrink=1:2){
  #if(class(count_matrix)=="matrix"){mode(count_matrix) = "integer"}
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                    colData = as.data.frame(design_mat),design = design_mat)
  dds <- DESeq(dds)
  shr_res = NULL
  if(length(ests_to_shrink)>0){
    # Get shrinked fold changes
    ests_to_shrink = ests_to_shrink+1 # 1 is the intercept, cannot shrink it
    shr_res = get_shrinked_estimates(dds,ests_to_shrink)
    names(shr_res) = colnames(design_mat)[ests_to_shrink]
  }
  cook_matrix = assays(dds)[["cooks"]]
  return(list(original_deseq_res = dds,shr_res = shr_res))
}

deseq2_analyze_cook<-function(count_matrix,design_mat){
  #if(class(count_matrix)=="matrix"){mode(count_matrix) = "integer"}
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = as.data.frame(design_mat),design = design_mat)
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
  cook_matrix = assays(dds)[["cooks"]]
  return(list(original_deseq_res = dds,res = res,cook=cook_matrix))
}


# Auxiliary functions for dealing with DESeq's results
# try(library("apeglm"))
get_shrinked_estimates<-function(dds,coeffs=NULL){
  res = list()
  if(is.null(coeffs)){
    resnames = resultsNames(dds)
    coeffs = 2:length(resnames)
  }
  for(i in coeffs){
    res[[as.character(i)]] = lfcShrink(dds,coef=i,type="apeglm")
  }
  return(res)
}

#' A wrapper for running simple linear regression
#' @param y the dependent variable
#' @param x a matrix with the predictors
#' @form a formula to be used in lm
#' @return a vector with the coefficients and their significance
lm_wrapper_for_diff_abundance_analysis<-function(y,x,form = NULL){
  df = data.frame(y=y,x)
  if(is.null(form)){
    form = y~.
  }
  lm_res = lm(form,data=df)
  lm_sum = summary(lm_res)
  coeffs = lm_sum$coefficients
  x1 = coeffs[,1]
  x2 = coeffs[,4]
  names(x1) = paste("est",rownames(coeffs),sep=";")
  names(x2) = paste("pval",rownames(coeffs),sep=";")
  return(c(x1,x2))
}
