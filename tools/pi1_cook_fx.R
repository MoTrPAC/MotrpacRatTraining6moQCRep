#!/bin/R
# Nicole Gay
# 28 July 2020 

# ancillary functions for outlier and covariate analysis 

# load required packages
library(data.table)
library(DESeq2) 
library(qvalue)
library(ggplot2)
library(ggrepel)
library(ggcorrplot)
library(reactome.db)
library(fgsea)
library(org.Rn.eg.db)
library(car) # for cook distances (lm/glm)
library(AnnotationDbi)
library(limma)
library(ggplotify)


#' PARENT FUNCTION
#' ID and remove PCA outliers
#' perform pi1-selection of covariates and ID/remove Cook outliers 
#' summarize results from all interations
#' 
#' @param quant matrix of raw counts (for \code{da_method = 'DESeq'}) or normalized values (for \code{da_method = 'lm'})
#' @param metadata data.frame or data.table with \code{variables_of_interest} \code{outcome_of_interest} and 'viallabel' as column names 
#' @param covariates_of_interest candidate list of technical covariates; can be correlated 
#' @param outcome_of_interest categorical variable on which contrasts are performed, e.g. time+intervention label 
#' @param da_method one of \code{c('DESeq','lm','limma')}. 'lm' and 'limma' perform the same
#' @param min_pc_ve minimum variance a PC must explain in order to be used to call outliers, e.g. 0.05
#' @param pc_iqr_coef flag PC outliers if they are outside of IQR * \code{pc_iqr_coef}
#' @param pi1_cutoff pi1 cutoff to select covariates
#' @param run_cook whether or not to look for outliers in Cook's distance 
#' @param cook_qf_cutoff Cook F distn probability to select genes for colMeans
#' @param cook_median_factor flag an outlier if its average Cook distance in selected genes is greater than this many times the median across samples
#' @param plot bool, whether or not to print plots
#' @param verbose bool, whether or not to print descriptive strings
#' @param ... additional parameters for lm
#' 
#' @return A named list
id_outliers_covariates = function(quant, metadata, covariates_of_interest, outcome_of_interest, 
                                  da_method = 'DESeq',
                                  min_pc_ve=0.075, pc_iqr_coef=3, 
                                  pi1_cutoff=0.1, 
                                  run_cook=T, cook_qf_cutoff=0.5, cook_median_factor=100, 
                                  plot=T, verbose=T, description=NULL, ...){
  
  if(!run_cook){
    cook_qf_cutoff=NULL
    cook_median_factor=NULL
    co=NULL
  }
  
  if(length(covariates_of_interest)==0){
    pi1_cutoff=NULL
  }
  #-- check the supplied data -----------------------------
  
  if(verbose){cat('\nChecking input data...\n')}
  formatted = check_data(quant, metadata, covariates_of_interest, outcome_of_interest, method=da_method)
  # overwrite
  quant = formatted$quant
  metadata = formatted$meta
  full_metadata = formatted$full_meta
  covariates_of_interest = formatted$cov
  
  #-- remove PCA outliers ----------------------------------
  
  if(verbose){cat('\nLooking for outliers in PC space...\n')}
  
  # run PCA
  # correct for size factors when input is counts 
  if(tolower(da_method) == 'deseq'){
    tmm = edgeR_normalized_log_cpm(quant)
    pca_quant = as.data.frame(tmm, check.names=F)
  }else{
    pca_quant = quant 
  }
  if('sex' %in% colnames(meta) & length(unique(meta[,sex]))==2){
    # split by sex 
    pca_outliers = c()
    for(s in unique(meta[,sex])){
      sub_quant = pca_quant[,meta[sex == s,viallabel]]
      pca_outliers = c(pca_outliers, id_pca_outliers(sub_quant, min_pc_ve, plot, verbose, iqr_coef=pc_iqr_coef)$pca_outliers, TITLE=description)
    }
  }else{
    if(length(unique(meta[,sex]))!=2){
      warning("'sex' found in columns of metadata, but there are not 2 unique values. Removing PCA outliers with all samples together.\n")
    }else{
      warning("'sex' not found in columns of metadata. Removing PCA outliers with all samples together.\n")
    }
    # pca with all samples
    pca_outliers = id_pca_outliers(pca_quant, min_pc_ve, plot, verbose, iqr_coef=pc_iqr_coef, TITLE=description)$pca_outliers
  }

  if(length(pca_outliers) > 0){
    if(verbose){cat(sprintf("Removing PCA outliers from data:\n   %s\n", paste(pca_outliers, collapse=', ')))}
    metadata = metadata[!viallabel %in% pca_outliers]
    quant = quant[,!colnames(quant) %in% pca_outliers]
    for (cov in covariates_of_interest){
      if(length(unique(metadata[,get(cov)]))==1){
        warning(sprintf("Variable of interest %s has 0 variance after removing PCA outliers. Removing.\n", cov))
        covariates_of_interest = covariates_of_interest[covariates_of_interest != cov]
      }
    }
  }
  
  #-- calc pi1 ---------------------------------------------

  # calculate pi1 
  if(length(covariates_of_interest) > 0){
    if(verbose){
      cat('\nCalculating pi1 for each covariate of interest...\n')
      cat(sprintf('Dim of quant: %s\n', paste(dim(quant), collapse=',')))
      cat(sprintf('Dim of metadata: %s\n', paste(dim(metadata), collapse=',')))
    }
    if(tolower(da_method) == 'deseq'){
      pi1 = calculate_pi1_deseq(covariates_of_interest, quant, metadata, pi1_cutoff, verbose)
    }else{
      pr_metadata = as.data.frame(metadata[,covariates_of_interest,with=F])
      if(is.null(dim(pr_metadata))){
        pr_metadata = matrix(pr_metadata,ncol=1,
                             dimnames = list(rownames(metadata),covariates_of_interest))
        pr_metadata = data.frame(pr_metadata)
      }
      pi1 = calculate_pi1_lm(quant, pr_metadata, verbose, ...)
    }
    
    pi1_values = pi1$pi1_values
    # plot pi1 values
    g = plot_pi1(pi1_values, pi1_cutoff)
    if(plot){print(g)}
  
    # ID pi1 covariates
    need_correction = names(pi1_values)[pi1_values > pi1_cutoff]
    if(verbose){
      if(length(need_correction) == 0){
        cat(sprintf('No pruned covariates with pi1 > %s.\n', pi1_cutoff))
      }else{
        cat(sprintf('pi1-selected covariates: %s\n', paste(need_correction, collapse=',')))
      }
    }
  }else{
    if(verbose){cat('\nNo covariates supplied. Skipping pi1 selection analysis.\n')}
    need_correction = NULL
    pi1 = list(pvalues = NULL, pi1_values = NULL)
  }
  
  # -- look for outliers in cook distance -----------------
  
  if(run_cook){
    if(verbose){cat("\nLooking for outliers in Cook's distance...\n")}
    if(tolower(da_method) == 'deseq'){
      co = cook_outlier_deseq(need_correction, quant, metadata, outcome_of_interest, cook_qf_cutoff, cook_median_factor, plot, verbose)
    }else{
      co = cook_outlier_lm(need_correction, quant, metadata, outcome_of_interest, cook_qf_cutoff, cook_median_factor, plot, verbose)
    }
    if(length(co$cook_outliers)>0){
      # remove outliers 
      metadata = metadata[!as.character(viallabel) %in% co$cook_outliers]
      quant = quant[!colnames(quant) %in% co$cook_outliers]
      # remove covariates with 0 variance (this can happen if a sample is removed)
      for (cov in covariates_of_interest){
        if(length(unique(metadata[,get(cov)]))==1){
          warning(sprintf("Variable of interest %s has 0 variance after removing Cook outliers. Removing.\n", cov))
          covariates_of_interest = covariates_of_interest[covariates_of_interest != cov]
        }
      }
      
      #-- recalc pi1 if Cook outliers were found ------------
      
      if(length(covariates_of_interest)>0){
        if(verbose){cat("\nRecalculating pi1 values...\n")}
        if(tolower(da_method) == 'deseq'){
          pi1 = calculate_pi1_deseq(covariates_of_interest, quant, metadata, pi1_cutoff, verbose)
        }else{
          pr_metadata = as.data.frame(metadata[,covariates_of_interest,with=F])
          if(is.null(dim(pr_metadata))){
            pr_metadata = matrix(pr_metadata,ncol=1,
                                 dimnames = list(rownames(metadata),covariates_of_interest))
            pr_metadata = data.frame(pr_metadata)
          }
          pi1 = calculate_pi1_lm(quant, pr_metadata, verbose, ...)
        }
        pi1_values = pi1$pi1_values
        # plot pi1 values
        g = plot_pi1(pi1_values, pi1_cutoff)
        if(plot){print(g)}
        # ID pi1 covariates
        need_correction = names(pi1_values)[pi1_values > pi1_cutoff]
        if(verbose){
          if(length(need_correction) == 0){
            cat(sprintf('No pruned covariates with pi1 > %s.\n', pi1_cutoff))
          }else{
            cat(sprintf('pi1-selected covariates: %s\n', paste(need_correction, collapse=',')))
          }
        }
      }else{
        need_correction = NULL
        pi1 = list(pvalues = NULL, pi1_values = NULL)
      }
    }
  }
  
  #-- summarize results -----------------------
  
  final = summarize_iterations(pca_outliers, co, pi1, pi1_cutoff, full_metadata, outcome_of_interest, verbose)
  
  #-- save parameters -------------------------
  
  params = list(covariates_of_interest, outcome_of_interest, 
                min_pc_ve, pc_iqr_coef, pi1_cutoff, 
                cook_qf_cutoff, cook_median_factor, da_method, run_cook)
  names(params) = c('covariates_of_interest', 'outcome_of_interest', 
                    'min_pc_ve', 'pc_iqr_coef', 'pi1_cutoff', 
                    'cook_qf_cutoff', 'cook_median_factor', 'da_method', 'run_cook')
  
  return(list(pca_outliers = final$pca_outliers, 
              cook_outliers = final$cook_outliers,
              pi1_pval_dt = pi1$pvalues,
              pi1_values = pi1$pi1_values,
              final_pi1_covariates = final$pi1_covariates,
              adjusted_metadata = metadata,
              pipeline_params = params))
}


plot_pi1 = function(pi1_values, pi1_cutoff){
  pi1_dt = data.table(variable = names(pi1_values), pi1 = unname(unlist(pi1_values)))
  pi1_dt = pi1_dt[order(pi1, decreasing=T)]
  g = ggplot(pi1_dt, aes(x=variable, y=pi1)) +
    geom_bar(stat='identity',position='dodge',colour='black',fill='black') +
    geom_hline(yintercept=pi1_cutoff, linetype='dashed', colour='red') +
    theme(legend.title=element_blank()) +
    scale_x_discrete(limits=pi1_dt[,variable]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90,hjust=1,colour='black',vjust=0.5),
          axis.title.x=element_blank(),
          title = element_blank()) +
    labs(y='pi1, proportion of true alt. tests')
  return(g)
}


check_data = function(expr, metadata, covariates_of_interest, outcome_of_interest, method){
  
  # method-specific things
  method = tolower(method)
  if(method == 'deseq'){
    # coerce counts to integers, just in case. RSEM does something weird 
    warning("Method 'DESeq' selected. Coercing values in 'quant' to integers.\n")
    expr = as.data.frame(apply(expr, c(1,2), as.integer)) 
  }else if(method %in% c('lm','limma')){
    # expression is normalized
    #
  }else{
    stop(sprintf("DE method '%s' not recognized. Please select from: 'DESeq','lm','limma'.\n", da_method))
  }
  
  # coerce metadata to data.table
  if(!is.data.table(metadata)){metadata = as.data.table(metadata)}
  
  # variables are in meta 
  if(!all(c(outcome_of_interest, covariates_of_interest) %in% colnames(metadata))){
    missing = c(outcome_of_interest, covariates_of_interest)[!c(outcome_of_interest, covariates_of_interest) %in% colnames(metadata)]
    stop(sprintf("The following variables were included as covariates or the outcome of interest but are not columns in the metadata:\n   %s\n",
                 paste(missing, collapse=', ')))
  }
  if(!'viallabel' %in% colnames(metadata)){
    stop("Please define the column of vial labels as 'viallabel' in 'metadata'.\n")
  }
  metadata[,viallabel := as.character(viallabel)] # this makes things easier downstream
  # coerce expr to data.frame
  if(!is.data.frame(expr)){norm = as.data.frame(expr)}
  
  # make sure covariates have more than one level 
  for(cov in covariates_of_interest){
    if(length(unique(metadata[,get(cov)]))==1){
      warning(sprintf("Variable of interest %s has 0 variance. Removing.\n", cov))
      covariates_of_interest = covariates_of_interest[covariates_of_interest != cov]
    }
  }
  
  # make sure there are no missing values in covariates
  # remove missing values 
  new = fix_missing(covariates_of_interest, metadata)
  covariates_of_interest = new$covariates
  metadata = new$meta

  # make sure variables are supplied correctly
  if(length(outcome_of_interest) == 0){
    stop("Please provide an outcome of interest, i.e. intervention/time group.\n")
  }else if(length(outcome_of_interest) > 1){
    warning("More than one outcome of interest was provided, i.e. the variable on which contrasts are made. Was this a mistake?\n")
  }
  if(length(covariates_of_interest) == 0){
    warning("No technical variables of interest were supplied for pi1 analysis.\n")
  }
  if(any(outcome_of_interest %in% covariates_of_interest)){
    warning("Outcome of interest was included in list of covariates of interest. Removing.\n")
    covariates_of_interest = covariates_of_interest[!covariates_of_interest %in% c(outcome_of_interest)]
  }
  
  # check dim of input data 
  if(dim(metadata)[1] != dim(expr)[2]){
    stop("There are different numbers of samples in the rows of 'metadata' and the columns of 'quant'. Please subset both objects to the same set of samples.\n")
  }
  
  # put samples in the same order in both tables (not that this really matters, but it also checks column names)
  expr = expr[,metadata[,viallabel]]
  if(dim(expr)[2] == 0){
    stop(paste0("Vial labels in metadata column 'viallabel' do not match column names 'quant'.\n"))
  }
  
  full_meta = copy(metadata) # save for later 
  
  return(list(quant=expr, 
              meta=metadata,
              full_meta=full_meta,
              cov=covariates_of_interest))
}


#' Internal function to normalize raw counts using TMM
#' 
#' @param x raw counts
edgeR_normalized_log_cpm = function(x,min_cpm = 0.5,min_num_samples = 2,norm_method="TMM"){
  # raw counts filtered down to samples in this tissue 
  raw_dge = edgeR::DGEList(counts=x) 
  keep = rowSums(cpm(raw_dge) > min_cpm) >= min_num_samples
  filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
  
  # filt --> tmm 
  dge = edgeR::calcNormFactors(filt_dge, method=norm_method)
  tmm = edgeR::cpm(dge,log=TRUE)
  
  return(tmm)
}

#' Internal function to plot pairs of PCs; highlights indicated outliers 
#' 
#' @param pcaA string of PC for x-axis, e.g. "PC1"; also a column name in \code{pcax}
#' @param pcaB string of PC for y-axis, e.g. "PC2"; also a column name in \code{pcax}
#' @param pcax e.g. \code{prcomp(data)$x}
#' @param outliers list of viallabels corresponding with PCA outliers 
#' @param pca result returned by \code{prcomp()}
plot_pcs = function(pcA, pcB, pcax, outliers, pca, title=NULL){
  
  pcax = as.data.frame(pcax)
  pcax$viallabel = rownames(pcax)
  
  # get VE 
  explained_var = summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])]
  xlab = sprintf("%s (%s%%)", pcA, round(explained_var[[pcA]]*100, 1))
  ylab = sprintf("%s (%s%%)", pcB, round(explained_var[[pcB]]*100, 1))
  
  g = ggplot(pcax, aes(x=get(pcA), y=get(pcB))) +
    geom_point() +
    theme_classic() +
    labs(x=xlab, y=ylab, title=title)
  
  if(length(outliers) > 0){
    g = g + geom_point(data=pcax[rownames(pcax) %in% outliers,], size=2, colour='red') +
      geom_text_repel(data=pcax[rownames(pcax) %in% outliers,], aes(label=viallabel)) +
      labs(title=title)
  }
  
  return(g)
}


#' ID outliers in PC space
#' 
#' @param norm filtered, normalized values 
#' @param min_pc_ve
#' @param plot bool, whether or not to print plots
#' @param verbose bool, whether or not to print descriptive strings
#' @param iqr_coef flag PC outliers if they are outside of IQR * \code{iqr_coef}
#' @param N select N most variable features
id_pca_outliers = function(norm, min_pc_ve, plot, verbose, iqr_coef=3, N=Inf, TITLE=NULL){
  
  # keep N features with highest CVs
  if(N < nrow(norm)){
    cv = apply(norm, 1, function(x) (sd(x)/mean(x))*100)
    cv = cv[order(cv, decreasing=T)]
    norm = norm[names(cv)[1:N],]
  }
  
  # remove features with 0 variance 
  tnorm = as.data.frame(t(norm))
  novar = names(which(apply(tnorm, 2, stats::var)==0))
  tnorm[,novar] = NULL
  pca = prcomp(x = tnorm,center=T,scale.=T)
  cum_var = summary(pca)[["importance"]][3,1:ncol(summary(pca)[["importance"]])]
  # save for later
  pca_before = pca
  
  # save as many PCs that explain at least X variance
  explain_vars = summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])]
  if(verbose){print(summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])][1:10])}
  explain_var = explain_vars[explain_vars > min_pc_ve]
  num_pcs = length(explain_var)
  if(num_pcs == 0){
    num_pcs = 1
    explain_var[explain_vars[1]]
  }
  
  if(verbose){cat(sprintf("The first %s PCs were selected to identify outliers.\n", num_pcs))}
  pcax = pca$x[,1:max(2, num_pcs)]
  
  # ID outliers 
  # Univariate: use IQRs
  pca_outliers_report = c()
  pca_outliers = c()
  for(j in 1:num_pcs){
    outlier_values = boxplot.stats(pcax[,j],coef=iqr_coef)$out # flag samples beyond IQR*iqr_coef
    for(outlier in names(outlier_values)){
      pca_outliers_report = rbind(pca_outliers_report,
                                  c(paste("PC",j,sep=""),outlier,
                                    unname(format(outlier_values[outlier],digits=5)))
      )
      if(!is.element(outlier,names(pca_outliers))){
        pca_outliers[outlier] = outlier_values[outlier]
      }
    }
  }
  if(verbose){
    if(length(pca_outliers) == 0){
      cat("No PC outliers.\n")
    }else{
      cat("PC outliers:\n")
      colnames(pca_outliers_report) = c('PC','viallabel','PC value')
      print(pca_outliers_report)
    }
  }
  
  if(plot){
    # plot before
    cat("PC plots with any outliers flagged:\n")
    for(i in 2:max(2,num_pcs)){
      print(plot_pcs("PC1", paste0("PC", i), pcax, names(pca_outliers), pca, title=sprintf('%s before outlier removal',TITLE)))
    }
    if(length(pca_outliers) > 0){
      cat("PC plots with outliers removed:\n")
      # plot after 
      filt_norm = norm[,!colnames(norm) %in% names(pca_outliers)]
      # remove features with 0 variance 
      tnorm <- as.data.frame(t(filt_norm))
      novar <- names(which(apply(tnorm, 2, stats::var)==0))
      tnorm[,novar] = NULL
      pca = prcomp(x = tnorm,center=T,scale.=T)
      pcax = pca$x[,1:max(2, num_pcs)]
      for(i in 2:max(2, num_pcs)){
        print(plot_pcs("PC1", paste0("PC", i), pcax, c(), pca, title=sprintf('%s after outlier removal',TITLE)))
      }
    }
  }
  
  # return outlier list 
  return(list(pca_outliers = names(pca_outliers),
              prcomp_obj = pca_before,
              num_pcs = num_pcs,
              pc_outliers_report = pca_outliers_report))
}


#' Get pairwise correlations
get_cor = function(keep, meta_data){
  exprn = paste0('~', paste(keep, collapse = ' + '))
  meta_data_cat2num = model.matrix(object = as.formula(exprn), data = meta_data)[,-1]
  c = cor(meta_data_cat2num, use = "pairwise.complete.obs", method = 'spearman')
  c_melt = data.table(reshape2::melt(c))
  c_melt = c_melt[Var1!=Var2]
  return(c_melt)
}


#' Internally used function to plot a heatmap of pairwise correlations between variables
plot_cor = function(vars, .title, meta_data){
  exprn = paste0('~', paste(vars, collapse = ' + '))
  meta_data_cat2num = model.matrix(object = as.formula(exprn), data = meta_data)[,-1]
  c = cor(meta_data_cat2num, use = "pairwise.complete.obs", method = 'spearman')
  return(ggcorrplot(c,title=.title,lab=F,type='upper'))
}


#' Calculate pi1 for a list of variables using DESeq
#' 
#' @param keep list of technical variables to include as covariates in the regressions 
#' @param counts filtered raw counts 
#' @param meta_data data.table with \code{keep} and "viallabel" as column names 
#' @param .plot whether or not to plot pi1 values 
#' 
#' @return A named list 
calculate_pi1_deseq = function(keep, counts, meta_data, cutoff, verbose){
  
  # center and scale continuous variables (recommended by DESeq)
  sub = meta_data[,keep, with=F]
  for (cov in keep){
    if(is.numeric(sub[,get(cov)])){
      sub[,(cov) := scale(sub[,get(cov)], center = T, scale = T)]
    }
  }
  
  pval_dt = data.table(gene_id=rownames(counts))
  for(cov in keep){
    if(verbose){cat(sprintf('  %s...\n', cov))}
    dds = DESeqDataSetFromMatrix(countData = counts,
                                 colData = sub,
                                 design = eval(parse(text=sprintf('~ %s', cov))))
    dds = DESeq(dds, quiet = T)
    if(is.factor(meta_data[,get(cov)]) | is.character(meta_data[,get(cov)])){
      # get res for all levels; keep most sig p-value
      tmp = data.table(gene_id=rownames(counts))
      tmppi = list()
      for(n in resultsNames(dds)[grepl(cov, resultsNames(dds))]){
        res = results(dds, name=n)
        tmp[,(n) := res$pvalue]
        tmppi[[n]] = 1 - qvalue(p=res$pvalue)$pi0
      }
      print(head(tmp))
      print(tmppi)
      # keep highest pi1
      pval_dt[,(cov) := tmp[,get(names(tmppi)[which.max(tmppi)])]]
      print(head(pval_dt))
    }else{
      res = results(dds, name=cov)
      # get p-values
      pval_dt[,(cov) := res$pvalue]
    }
  }

  # calculate pi1
  pi_values = list()
  for (col in colnames(pval_dt)){
    if(col=='gene_id'){next}
    pi1 = 1 - qvalue(p=pval_dt[,get(col)])$pi0
    pi_values[[col]] = pi1
  }
  
  # pvalues: gene x covariate data.table of beta p-values 
  # pi1_values: named list of pi1 values for each covariate in "keep"
  return(list(pvalues = pval_dt, 
              pi1_values = pi_values))
}


#' Obtain pi1 scores using limma
#' 
#' @param m a matrix or a data.frame with rows as analytes
#' @param covs a data.frame. the set of covariates of interest. Rows are samples.
#' @param ... additional parameters for limma (e.g., model="robust")
calculate_pi1_lm = function(m,covs,verbose,...){
  
  pval_df = data.frame(gene_id = rownames(m))
  pi1s = list()
  for(c in colnames(covs)){
    if(verbose){cat(sprintf('  %s...\n', c))}
    form = sprintf('~ 1 + %s', c)
    modelmat = as.data.table(model.matrix(eval(parse(text=form)), covs))
    limma_fit = lmFit(m,modelmat,method="robust",...)
    e_fit = eBayes(limma_fit)
    cov_res = topTable(e_fit,coef = c,number = nrow(m),sort.by = "none")
    cov_ps = cov_res[,"P.Value"]
    pval_df[,c] = cov_ps
    pi1s[[c]] = 1 - limma::propTrueNull(cov_ps)
  }
  
  # return pval dt as well as pi1s 
  return(list(pvalues = as.data.table(pval_df),
              pi1_values = pi1s))
}


#' Identify outliers in terms of Cook's distance using DESeq
#' 
#' @param need_correction vector of covariate names, e.g. pi1-selected covariates
#' @param counts counts for DESeq2
#' @param meta metadata for DESeq2; columns must include \code{need_correction} and \code{outcome_of_interest}
#' @param outcome_of_interest categorical variable on which contrasts are performed, e.g. time+intervention label 
#' @param cook_qf_cutoff Cook F distn probability to select genes for colMeans
#' @param cook_median_factor flag an outlier if its average Cook distance in selected genes is greater than this many times the median across samples
#' @param plot whether or not to print plots
#' @param verbose whether or not to print descriptive strings
#' 
#' @return A named list 
cook_outlier_deseq = function(need_correction, counts, meta, outcome_of_interest, cook_qf_cutoff, cook_median_factor, plot=T, verbose=T){
  
  # center and scale continuous variables (recommended by DESeq)
  sub = meta[,c(need_correction, outcome_of_interest), with=F]
  for (cov in need_correction){
    if(is.numeric(sub[,get(cov)])){
      sub[,(cov) := scale(sub[,get(cov)], center = T, scale = T)]
    }
  }
  
  # run DEA with pi1 covariates 
  contrast = paste0('~ ', paste(c(outcome_of_interest, need_correction), collapse = ' + '))
  if(verbose){cat(sprintf('Running DESeq2 with the model:\n   %s\n', contrast))}
  dds = DESeqDataSetFromMatrix(countData = counts,
                               colData = sub,
                               design = eval(parse(text=contrast)))
  dds = DESeq(dds, quiet=T)
  if(verbose){
    cat("dds summary:")
    print(summary(results(dds)))
  }

  cook = assays(dds)[["cooks"]]
  cook = as.data.frame(cook)
  
  # select genes that have a Cook's distance above some values 
  matModelMatrix = model.matrix(eval(parse(text=contrast)), data=meta) # same contrast as for DESeq above
  df1 = ncol(matModelMatrix) # num_levels-1 of 'outcome_of_interest' + number of continuous or dummy variables in 'need_correction' + Intercept
  df2 = ncol(counts) - df1 # number of samples - df1
  
  cook_res = process_cook(cook, cook_qf_cutoff, cook_median_factor, df1, df2, plot, verbose)
  return(cook_res)
}


#' Compute linear analysis with cook distances for a single analyte
#' 
#' @param x a numeric vector, the dependent variable
#' @param design_mat design matrix with covariates to adjust for
lm_wrapper_single_analyte = function(x,design_mat){
  df = data.frame(x=as.numeric(x),design_mat)
  model = lm(x~.,data=df)
  cook = cooks.distance(model)
  return(list(model=model,cook=cook))
}


#' Identify outliers in terms of Cook's distance using lm
#' 
#' @param need_correction vector of covariate names, e.g. pi1-selected covariates
#' @param counts counts for DESeq2
#' @param meta metadata for DESeq2; columns must include \code{need_correction} and \code{outcome_of_interest}
#' @param outcome_of_interest categorical variable on which contrasts are performed, e.g. time+intervention label 
#' @param cook_qf_cutoff Cook F distn probability to select genes for colMeans
#' @param cook_median_factor flag an outlier if its average Cook distance in selected genes is greater than this many times the median across samples
#' @param plot whether or not to print plots
#' @param verbose whether or not to print descriptive strings
#' 
#' @return A named list 
cook_outlier_lm = function(need_correction, expr, meta, outcome_of_interest, cook_qf_cutoff, cook_median_factor, plot=T, verbose=T){
  
  # make design mat 
  contrast = paste0(' ~ ', paste0(c(need_correction, outcome_of_interest), collapse = ' + '))
  design_matrix = model.matrix(eval(parse(text=contrast)), data = meta)

  #compute cook
  lm_objects = apply(expr,1,lm_wrapper_single_analyte,
                     design_mat=design_matrix)
  cook_matrix = as.data.frame(t(sapply(lm_objects,function(x)x$cook)))

  # select genes that have a Cook's distance above some values 
  df1 = ncol(design_matrix) # num_levels-1 of 'outcome_of_interest' + number of continuous or dummy variables in 'need_correction' + Intercept
  df2 = ncol(expr) - df1 # number of samples - df1
  
  cook_res = process_cook(cook_matrix, cook_qf_cutoff, cook_median_factor, df1, df2, plot, verbose)
  return(cook_res)
}


process_cook = function(cook_matrix, cook_qf_cutoff, cook_median_factor, df1, df2, plot, verbose){
  
  print(cook_qf_cutoff)
  cutoff = qf(cook_qf_cutoff, df1, df2)
  selected_genes = rownames(cook_matrix)[apply(cook_matrix, 1, function(x) any(x > cutoff))]
  cook_sub = cook_matrix[selected_genes,]
  if(verbose){cat(sprintf("%s out of %s genes have a Cook's distance of at least %s.\n", 
                          length(selected_genes), nrow(cook_matrix), round(cutoff, 2)))}
  
  m = colMeans(cook_sub, na.rm=T)
  m = m[order(m)]
  mdf = data.frame(m)
  mdf$viallabel = rownames(mdf)
  if(median(m) == 0){
    out = cook_median_factor*0.001
  }else{
    out = cook_median_factor*median(m)
  }
  g = ggplot(mdf, aes(x=viallabel, y=m)) +
    geom_point() +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    labs(y=sprintf("Average Cook's distance across genes (M=%s)", length(selected_genes))) +
    scale_x_discrete(limits = names(m))
  # only add the cutoff line if there are points above it
  if(out < max(mdf$m, na.rm=T)){
    g = g + geom_hline(yintercept = out, linetype = 'dashed', colour = 'red') 
  }
  if(plot){print(g)}
  
  cook_outliers = names(m)[m > out]
  if(verbose){
    if(length(cook_outliers) == 0){
      cat('No Cook outliers.\n')
    }else{
      cat(sprintf('Cook outliers: %s\n', paste(cook_outliers, collapse=',')))
    }
  }
  
  return(list(cook_matrix = cook_matrix,
              cook_outliers = cook_outliers))
}


#' Summary pipeline results (pi1 covariates and removed outliers)
#' 
#' @param pca_outliers vector of viallabels IDed as PC outliers
#' @param cook result of \code{cook_outlier_deseq()} or \code{cook_outlier_lm()}
#' @param pi1 result of \code{calculate_pi1_deseq()} or \code{calculate_pi1_lm()}
#' @param meta_data full metadata without outliers removed
#' @param verbose whether or not to print descriptive strings
#' 
#' @return A list with three values: 
#'   a vector of viallabels flagged as Cook outliers 
#'   a vector of viallabels flagged as PCA outliers 
#'   a vector of the final pi1-selected covariates 
summarize_iterations = function(pca_outliers, cook, pi1, pi1_cutoff, meta_data, outcome_of_interest, verbose){
  
  final_pi1 = names(pi1$pi1_values)[pi1$pi1_values > pi1_cutoff]
  
  if(is.null(cook)){
    cook_outliers = NULL
  }else{
    cook_outliers = cook$cook_outliers 
  }
  
  if(verbose){
    cat("\nSummary:\n")
    if(length(pca_outliers)>0){
      cat(sprintf('PCA outliers to exclude: %s\n', paste(pca_outliers,collapse=',')))
      print(meta_data[as.character(viallabel) %in% pca_outliers, c('viallabel', outcome_of_interest), with=F])
    }else{
      cat('No outliers in terms of PCA.\n')
    }
    
    if(!is.null(cook)){
      if(length(cook$cook_outliers)>0){
        cat(sprintf('Cook outliers to exclude: %s\n', paste(cook$cook_outliers,collapse=',')))
        print(meta_data[as.character(viallabel) %in% cook$cook_outliers, c('viallabel', outcome_of_interest), with=F])
      }else{
        cat('No outliers in terms of Cook distance.\n')
      }
    }else{
      cat('Cook distance analysis was not run.\n')
    }
    
    if(length(final_pi1) > 0){
      cat(sprintf('Covariates to include: %s\n', paste(final_pi1,collapse=',')))
    }else{
      cat('No pi1-selected covariates.\n')
    }
  }
  
  return(list(cook_outliers = cook_outliers,
              pca_outliers = pca_outliers, 
              pi1_covariates = final_pi1))
}


#' @param gene_list list of ENSEMBL IDs, ENTREZ IDs, or gene symbols to map to other IDs
#' @param key_type one of "ENSEMBL", "ENTREZID", or "SYMBOL", i.e. what kind of IDs are used in \code{gene_list}
#' 
#' @return data.frame with columns "ENSEMBL", "ENTREZID", and "SYMBOL". Duplicate entries for "ENSEMBL" are removed. 
make_gene_map = function(gene_list, key_type='ENSEMBL'){
  gene_map = AnnotationDbi::select(org.Rn.eg.db, key=gene_list, columns=c('ENSEMBL','ENTREZID','SYMBOL'), keytype=key_type)
  # remove duplicate Ensembl IDs
  gene_map = gene_map[!duplicated(gene_map$ENSEMBL),]
  return(gene_map)
}


## accessory lm/limma functions ############################################################################


#' Extract summary statistics from an lm_wrapper_single_analyte results list
#' 
#' @param obj a list. The result object of running lm_wrapper_single_analyte
#' @param coeff a character. The name of the covariate whose summary statistics are to be extracted
lm_wrapper_get_stats<-function(obj,coeff = "is_case"){
  model = obj$model
  mat = summary(model)$coefficients
  colnames(mat) = c("Est","Std","Tstat","Pvalue")
  res = mat[coeff,]
  return(res)
}


#' A wrapper that runs lm or limma on each row in an input matrix
#' 
#' @param X a matrix or a data.frame with rows as analytes
#' @param design_mat a design matrix with the covariates for the models
#' @param col a character. The name of the covariate of interest
#' @param compute_cook a logical. If TRUE the cook distances matrix is returned
#' @param use_limma a logical. If TRUE limma's lmFit is used with Bayesian adjustment.
#' @param ... additional parameters for limma (e.g., method="robust")
lm_wrapper<-function(X,design_matrix,col="is_case",compute_cook=F,use_limma=F,...){
  if(!compute_cook && !use_limma){
    return(
      t(apply(X,1,function(x,y,z){
        lm_wrapper_get_stats(lm_wrapper_single_analyte(x,y),z)
      },y=design_matrix,z=col))
    )
  }
  limma_res = NULL
  if(use_limma){
    limma_fit = lmFit(X,design_matrix,...)
    e_fit = eBayes(limma_fit)
    limma_res = topTable(e_fit,coef = col,number = nrow(X),sort.by = "none")
  }
  if(!compute_cook && use_limma){
    return(limma_res)
  }
  
  # At this point we know we need to compute cook:
  lm_objects = apply(X,1,lm_wrapper_single_analyte,
                     design_mat=design_matrix)
  cook_matrix = t(sapply(lm_objects,function(x)x$cook))
  if(use_limma){
    return(list(limma_res=limma_res,cook_matrix=cook_matrix))
  }
  lm_res = t(sapply(lm_objects,lm_wrapper_get_stats,coeff = col))
  if(use_limma){
    return(list(lm_res=lm_res,cook_matrix=cook_matrix))
  }
}
# # Test
# y1 = rnorm(500)
# y2 = rnorm(500)
# x = 0.6*y1 + 0.2*y2
# y1[1] = 10 # create one influential point
# df = data.frame(y1=y1,y2=y2)
# lm_res = lm_wrapper(x,df)
# plot(lm_res$cook[1:200],pch=20)
# summary(lm_res$model)


#' Run DESeq
#' 
#' @param counts raw filtered counts. columns will be subset to meta$viallabel, so make sure they match
#' @param meta metadata with columns \code{\param{outcome_of_interest}, \param{covar}, 'viallabel'} at a minimum. \param{counts} are subset to \code{meta$viallabel}
#' @param covar covariates to include in the DESeq model
#' @param outcome_of_interest outcome of interest to include in the model; levels provided in \param{contrasts}
#' @param constrasts list of vectors, where each vector is in the form \code(c(outcome_of_interest, numerator, denominator)), e.g. \code{c('sex_group','female.1w','female.control')}
#' @param shrink bool, whether or not to apply \code{lfcShrink()}
#' @param verbose bool, whether or not to print the design string 
run_deseq = function(counts, meta, covar, outcome_of_interest, contrasts, shrink=T, verbose=F){
  
  meta = as.data.table(meta)
  meta[,(outcome_of_interest) := as.factor(get(outcome_of_interest))]
  counts = counts[,as.character(meta[,viallabel])]
  
  # coerce to counts (RSEM does something weird)
  counts = as.data.frame(apply(counts, c(1,2), as.integer)) 
  
  # remove missing values 
  new = fix_missing(covar, meta)
  covar = new$covariates
  meta = new$meta
  
  # center and scale continuous variables 
  for (cov in covar){
    # remove if constant
    if(length(unique(meta[,get(cov)])) == 1){
      message(sprintf("Covariate %s is constant. Removing.", cov))
      covar = covar[covar != cov]
    }else{
      # center and scale
      if(is.numeric(meta[,get(cov)])){
        meta[,(cov) := scale(meta[,get(cov)], center = T, scale = T)]
      }
    }
  }
  
  # make contrast
  contrast = paste0('~', paste0(c(outcome_of_interest, covar), collapse=' + '))
  if(verbose) message(contrast)
  
  # run DESeq
  dds = DESeqDataSetFromMatrix(countData = counts,
                               colData = meta,
                               design = eval(parse(text=contrast)))
  dds = DESeq(dds, quiet = T)
  
  # get results for each contrast 
  res_list = list()
  for (c in contrasts){
    if(!shrink){
      res = results(dds, contrast = c)
    }else{
      res = lfcShrink(dds, contrast = c, type = 'ashr', quiet = T, 
                      control=list(numiter.em=1000), optmethod='mixSQP') # "apeglm" doesn't work with contrasts in this form 
    }
    res_dt = data.table(gene_id = rownames(counts), 
                        log2FoldChange = res$log2FoldChange,
                        lfcSE = res$lfcSE,
                        pvalue = res$pvalue,
                        numerator = c[2],
                        denominator = c[3],
                        covariates = paste0(covar, collapse=','))
    if(!shrink){
      res_dt[,stat := res$stat]
    }
    res_list[[paste0(c, collapse=' ')]] = res_dt
  }
  all_res = rbindlist(res_list)
  return(list(res=all_res,
              dds=dds))
}

#' @param results data.frame with 2 columns: "feature_ID" and "score"
#' @param gene_map data.frame with columns: "ENSEMBL", "ENTREZID". see \code{make_gene_map()}
reactome_gsea = function(results, gene_map, plot=T, filename=NULL){
  
  set.seed(10)
  gene_map = as.data.table(gene_map)
  
  # get ENTREZ ids
  results = as.data.table(results)
  results[,entrez := gene_map[match(results[,feature_ID], ENSEMBL), ENTREZID]]
  results = results[!is.na(score)]
  
  scores = results[,score]
  names(scores) = results[,entrez]
  
  # get reactome pathways
  pathways = reactomePathways(names(scores))
  # gsea
  fgseaRes = fgsea(pathways, scores, nperm=50000, maxSize=500, minSize=10, BPPARAM=MulticoreParam(workers = 1))
  
  # collapse similar pathways 
  if(plot){
    collapsedPathways = collapsePathways(fgseaRes[order(pval)][padj < 0.05], pathways, scores)
    mainPathways = fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
    h = length(mainPathways)*0.4
    pdf(sprintf('plots/%s.pdf', filename), width=12, height=h)
    plotGseaTable(pathways[mainPathways], scores, fgseaRes, gseaParam = 0.5)
    dev.off()
  }
  return(list(fgsea_res = fgseaRes,
              pathways = pathways,
              scores = scores))
}

#' @param map dl_read_gcp('gs://my-bucket/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt')
#' @param kegg_pathways gmtPathways('/oak/stanford/groups/smontgom/nicolerg/MOTRPAC/PATHWAYS/c2.cp.kegg.v7.2.symbols.gmt')
#' @param results data.frame with 2 columns: "feature_ID" and "score"
kegg_gsea = function(kegg_pathways, results, map, plot=T, filename=NULL){
  
  set.seed(10)
  map = data.table(map)
  
  # get human gene symbols
  results = as.data.table(results)
  results[,human_gene_symbol := map[match(results[,feature_ID], ensembl_gene), human_gene_symbol]]
  results = results[!is.na(score)]
  
  scores = results[,score]
  names(scores) = results[,human_gene_symbol]
  
  fgseaRes = fgsea(kegg_pathways, scores, nperm=100000, maxSize=500, minSize=10, BPPARAM=MulticoreParam(workers = 1))
  
  # collapse similar pathways 
  if(plot){
    collapsedPathways = collapsePathways(fgseaRes[order(pval)][padj < 0.05], pathways, scores)
    mainPathways = fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
    h = length(mainPathways)*0.4
    pdf(sprintf('plots/%s.pdf', filename), width=12, height=h)
    plotGseaTable(pathways[mainPathways], scores, fgseaRes, gseaParam = 0.5)
    dev.off()
  }
  return(list(fgsea_res = fgseaRes,
              pathways = pathways,
              scores = scores))
}

fix_missing = function(covar, meta){
  # make sure there are no missing values in covariates
  for(cov in covar){
    if(any(is.na(meta[,get(cov)]))){
      num_missing = nrow(meta[is.na(get(cov))])
      if(is.character(meta[,get(cov)]) | is.factor(meta[,get(cov)])){
        # don't know how to handle missing values for factor. remove.
        warning(sprintf("Categorical variable of interest %s has %s missing values. Removing.\n", cov, num_missing))
        covar = covar[covar != cov]
      }else{
        # if continuous 
        # remove if missing for >5% of samples
        if(num_missing/nrow(meta) > 0.05){
          warning(sprintf("Numeric variable of interest %s has %s missing values. Removing.\n", cov, num_missing))
          covar = covar[covar != cov]
        }else{
          warning(sprintf("Numeric variable of interest %s has %s missing values. Replacing missing values with mean.\n", cov, num_missing))
          meta[is.na(get(cov)), (cov) := mean(meta[,get(cov)], na.rm=T)]
        }
      }
    }
  }
  return(list(meta=meta, covariates=covar))
}



## NEED TO BE UPDATED ############################################################################
 
#' #' #' Internally used function to run all contrasts with DESeq
#' get_res = function(dds, outcome_of_interest, ref_level, meta){
#'   reslist = list()
#'   for(comp in unique(meta[,get(outcome_of_interest)])){
#'     if(comp == ref_level){next}
#'     res = results(dds, contrast=c(outcome_of_interest, comp, ref_level))
#'     dt = data.table(pvalue=res$pvalue,
#'                     log2FoldChange=res$log2FoldChange, 
#'                     gene=gsub('[.].*','',rownames(res)))
#'     dt[,padj := p.adjust(pvalue, method='BH')]
#'     dt[,numerator := comp]
#'     reslist[[comp]] = dt
#'   }
#'   res = rbindlist(reslist)
#'   return(res)
#' }
#' 
#' 
#' #' Get DESeq results before and after outlier removal and adjustment
#' #' 
#' #' @param pipeline_res output of id_outliers_covariates()
#' #' @param counts raw counts 
#' #' @param outcome_of_interest categorial variable of interest in \code{meta} for DESeq 
#' #' @param ref_level reference level of \code{outcome_of_interest} for DESeq contrasts
#' #' @param verbose whether or not to print descriptive strings
#' deseq_before_after = function(pipeline_res, counts, full_meta, outcome_of_interest='group', ref_level='control', verbose=T){
#'   
#'   # get covariates
#'   adjusted_meta = pipeline_res$adjusted_metadata
#'   pi1_cov = pipeline_res$final_pi1_covariates
#'   keep = c('viallabel', pi1_cov, outcome_of_interest)
#'   adjusted_meta = adjusted_meta[,keep, with=F]
#'   
#'   # before
#'   before_counts = copy(counts)
#'   before_meta = full_meta
#'   before_counts = before_counts[as.character(before_meta[,viallabel])]
#'   stopifnot(all(as.character(before_meta[,viallabel]) == colnames(before_counts)))
#'   
#'   # after 
#'   after_counts = copy(before_counts)
#'   outliers = c(pipeline_res$pca_outliers, pipeline_res$cook_outliers)
#'   if(length(outliers) > 0){
#'     after_counts = after_counts[!colnames(after_counts) %in% outliers]
#'   }
#'   after_meta = adjusted_meta
#'   after_meta = after_meta[as.character(viallabel) %in% colnames(after_counts)]
#'   after_counts = after_counts[as.character(after_meta[,viallabel])]
#'   for (cov in colnames(after_meta)){
#'     if(cov == 'viallabel'){next}
#'     if(is.numeric(after_meta[,get(cov)])){
#'       after_meta[,(cov) := scale(after_meta[,get(cov)], center = T, scale = T)]
#'     }
#'   }
#'   
#'   no_change = F
#'   if(length(outliers) == 0 & length(pi1_cov) == 0){
#'     # only make one set of plots if there are no outliers or covariates
#'     no_change = T
#'   }
#'   
#'   # DESeq - no adjustment or outlier removal 
#'   if(verbose){cat("Running DESeq without adjustment...\n")}
#'   contrast = sprintf('~ %s', outcome_of_interest)
#'   before_dds = DESeqDataSetFromMatrix(countData = before_counts,
#'                                       colData = before_meta,
#'                                       design = eval(parse(text=contrast)))
#'   before_dds = DESeq(before_dds, quiet=T)
#'   before_res = get_res(before_dds, outcome_of_interest, ref_level, before_meta)
#'   
#'   res_list = list()
#'   res_list[['before']] = before_res
#'   
#'   if(!no_change){
#'     if(verbose){cat("Running DESeq with outlier removal and adjustment...\n")}
#'     # DESeq - after outlier removal and with any covariates 
#'     contrast =  paste0('~ ', paste(c(outcome_of_interest, pi1_cov), collapse=' + '))
#'     after_dds = DESeqDataSetFromMatrix(countData = after_counts,
#'                                        colData = after_meta,
#'                                        design = eval(parse(text=contrast)))
#'     after_dds = DESeq(after_dds, quiet=T)
#'     after_res = get_res(after_dds, outcome_of_interest, ref_level, before_meta)
#'     
#'     res_list[['after']] = after_res
#'   }
#'   
#'   return(res_list)
#' }
#' 
#' 
#' #' @param res_list output of deseq_before_after()
#' deseq_qq_plot = function(res_list){
#'   
#'   # one line per item in res_list 
#'   qq_data_list = list()
#'   for(i in 1:length(res_list)){
#'     input_dt = res_list[[i]]
#'     groups = unique(input_dt[,numerator])
#'     dtlist = list()
#'     for(g in groups){
#'       observed = -log10(input_dt[numerator == g,pvalue])
#'       observed = observed[!is.na(observed)] # remove NA
#'       observed = observed[order(observed, decreasing=F)]
#'       expected = -log10(seq(1/length(observed), 1, length.out=length(observed)))
#'       expected = expected[order(expected, decreasing=F)]
#'       dt = data.table(observed=observed, expected=expected)
#'       dt[,numerator := g]
#'       dtlist[[g]] = dt
#'     }
#'     dt = rbindlist(dtlist)
#'     dt[,type := names(res_list)[i]]
#'     qq_data_list[[i]] = dt
#'   }
#'   qq_data = rbindlist(qq_data_list)
#'   
#'   g = ggplot(qq_data, aes(x=expected, y=observed, colour=factor(type))) +
#'     geom_abline(linetype='dashed') + 
#'     geom_point() +
#'     theme_classic() +
#'     theme(legend.title = element_blank()) +
#'     facet_wrap(~numerator, scales='free') +
#'     labs(x='Expected p-value (-log10)', y='Observed p-value (-log10)', title='DESeq p-values')
#'   
#'   return(g)
#' }
