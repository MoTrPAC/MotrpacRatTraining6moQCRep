# The functions in this file implement normalization methods that
# take as input a numeric matrix or a data frame.
# The goal of these functions is to create a matrix after reducing
# sample-to-sample variance, and/or address feature heteroskedasticity.
# Note that methods for differential analysis directly from raw data may 
# use additional steps aimed to reduce error and increase power. One example
# is DESeq2's fold change shrinkage. 
# These approaches are explored in the script differential_analysis_functions.R 
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
                  # For imputation using kNN or missForest for random forests
                  "impute","missForest","parallel","doParallel",
                  "sva" # dealing with batch - combat, sva, fsva
)
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}
###############################################################

# use xmsPANDA for some metabolomics options
# install_github("kuppal2/xmsPANDA")

#' A wrapper for preprocessCore's quantile normalization.
#' Comment: we do not use this by default
#' @param x a numeric matrix with samples as columns
#' @return the numeric matrix after column-based quantile normalization
run_quantile_normalization<-function(x){
  x = as.matrix(x)
  mode(x) = "numeric"
  newx = preprocessCore::normalize.quantiles.robust(x)
  rownames(newx) = rownames(x)
  colnames(newx) = colnames(x)
  return(newx)
}

try(library("NOISeq"))
run_tmm_normalization<-function(x,...){
  return(NOISeq::tmm(x,...))
}

run_median_value_norm<-function(x){
  meds = apply(x,2,median,na.rm=T)
  newx = t(t(x) - meds)
  return(newx)
}

run_median_mad_norm<-function(x,zero_medians=T){
  meds = apply(x,2,median,na.rm=T)
  mads = apply(x,2,mad,na.rm=T)
  mean_mads = mean(mads,na.rm=T)
  mads[mads==0] = 1
  mean_med = mean(meds,na.rm=T)
  for(j in 1:ncol(x)){
    x[,j] = (x[,j]-meds[j])/mads[j]
  }
  x = x*mean_mads
  if(!zero_medians){
    x = x + mean_med
  }
  return(x)
}

#' Takes an FPKM matrix, removes lowly expressed genes and log transform
#' the remaining matrix. This can be used on the FPKM output matrix
#' from the RNA-seq pipeline.
#' @param fpkm_matrix a numeric matrix with samples as the columns
#' @param intensity_threshold a threshold specifying low expression values
#' @param intensity_pct a threshold for excluding genes with low values
#' @return A matrix of log transformed FPKMs
process_fpkm1 <-function(fpkm_matrix,intensity_threshold=1,intensity_pct=0.2){
  lowly_expressed_genes = rowSums(
    fpkm_matrix<=intensity_threshold)/ncol(fpkm_matrix) > intensity_pct
  fpkm_matrix = fpkm_matrix[!lowly_expressed_genes,]
  fpkm_matrix = log(fpkm_matrix+1,base = 2)
  return(fpkm_matrix)
}

#' A simple function to determine low expressions using unsupervised analysis
remove_low_intensity_rows<-function(x,intensity_threshold=1,intensity_pct=0.2){
  lowly_expressed_genes = rowSums(x<=intensity_threshold) /
    ncol(x) > intensity_pct
  return(x[!lowly_expressed_genes,])
}

#' Takes a count matrix compute normalization factors, and
#' return the normalized log cpms.
#' @param x a numeric counts matrix with samples as the columns
edgeR_normalized_log_cpm<-function(x,min_cpm = 0.5,min_num_samples = 2,norm_method="TMM"){
  # raw counts filtered down to samples in this tissue 
  raw_dge = edgeR::DGEList(counts=x) 
  keep = rowSums(cpm(raw_dge) > min_cpm) >= min_num_samples
  filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
  
  # filt --> tmm 
  dge = edgeR::calcNormFactors(filt_dge, method=norm_method)
  tmm = edgeR::cpm(dge,log=TRUE)
  
  return(tmm)
}

# try(library("DESeq2"))
#' Use DESeq2 to estimate sample factors and gene dispersion
#' @param het_norm - character, a method for heteroskedasticity adjustment - can be log2, rlog, vst, or NULL
#' @return a DESeqDataSet
deseq_process_counts<-function(count_matrix,plotFactors=F,
                               het_norm = "rlog"){
  if(class(count_matrix)=="matrix"){mode(count_matrix) = "integer"}
  se <- SummarizedExperiment(count_matrix)
  dds <- DESeqDataSet(se, design = ~ 1 )
  #Estimate size factors
  dds <- DESeq2::estimateSizeFactors(dds)
  if(plotFactors){
    # Plot the size factors
    plot(sizeFactors(dds), colSums(counts(dds)),ylab="Library size",
         xlab = "DESeq estimated size factors")
    abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
  }
  if(is.null(het_norm)){
    newdata = assay(dds)
    return(newdata)
  }
  if(het_norm == "log2"){
    logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
    return(logcounts)
  }
  if(het_norm == "rlog"){
    rld <- rlog(dds)
    newdata = assay(rld)
    return(newdata)
  }
  if(het_norm == "vst"){
    dds <- estimateDispersions(dds) # used by vst
    vsd <- varianceStabilizingTransformation(dds)
    newdata = assay(vsd)
    return(newdata)
  }
  return(NA)
}

#' Plot coefficient of variation as a function of mean intensity
cv_mean_plot<-function(x,...){
  mu = apply(x,1,mean)
  sds = apply(x,1,sd)
  cvs = sds/mu
  ll = lowess(mu,cvs,iter = 10,f=1/50)
  plot(x=mu,y=cvs,...)
  lines(x=ll$x,y=ll$y,col="red",lwd=2)
}


#' MA plot for a pair of vectors with log intensities
ma_plot<-function(x,y,useEdgeR = T,...){
  M = x-y
  A = (x+y)/2
  if(!useEdgeR){plot(y=M,x=A,ylab="M",xlab="A",...)}
  else{
    edgeR::maPlot(logAbundance = A,logFC = M,lowess = T,pch=20,...)
    # MA <- new("MAList")
    # MA$M = M
    # MA$A = A
    # limma::plotMA(MA,...)
  }
  abline(0,0,col="blue",lty=2,lwd=2)
}

#' A hybrid approach for imputing data, mainly for metabolomics.
#' @param x a matrix, rows are metabolites, columns are samples
#' @param alpha a number between 0 and 1, threshold to determine which metabolites are imputed using the half min value 
#' @param ml_method a character: knn or rf (default)
#' @param num_cores a number, the number of cores to use with random forests (default is 6)
#' @description if a metabolite has < alpha*100 percent missing values impute using either kNN or random forest, otherwise use half min value within samples
min_val_knn_hybrid_imputation<-function(x,alpha=0.2,ml_method="rf",num_cores=6,seed=125){
  # check the input first
  if(alpha < 0 || alpha >1){
    stop("alpha must be between 0 and 1")
  }
  if(ml_method!="rf" && ml_method!="knn"){
    stop("ml_method must be knn or rf")
  }
  percent_na = rowSums(is.na(x)) / ncol(x)
  minvals = apply(x,1,min,na.rm=T)
  minvals[is.na(minvals)] = 0
  minvals = pmax(0,minvals)
  for(row in which(percent_na >= alpha)){
    currinds = is.na(x[row,])
    x[row,currinds] = minvals[row]/2
  }
  if(ml_method == "knn"){
    # impute - use capture.output to avoid prints
    capture.output(
      {x_imp = impute::impute.knn(as.matrix(x))$data}
      ,file=NULL)
  }
  if(ml_method == "rf"){
    cl <- makePSOCKcluster(num_cores)
    registerDoParallel(cl)
    set.seed(seed)
    x_imp = missForest(t(x),parallelize = "variables")
    x_imp = t(x_imp$ximp)
    stopCluster(cl)
  }
  return(x_imp)
} 


