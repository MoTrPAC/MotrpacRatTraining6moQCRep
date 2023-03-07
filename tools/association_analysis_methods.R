###############################################################
# # Load required functions
# required_libs = c("limma")
# for (lib_name in required_libs){
#   tryCatch({ary(lib_name,character.only = T)}, error = function(e) {
#     print(paste("Cannot load",lib_name,", please install"))
#   })
# }
###############################################################

#' The simple time point based linear analysis for PASS1A data
#' 
#' @param x a numeric vector with continuous values (an analyte)
#' @param tps a vector with the time points
#' @param tp the time point for the current analysis
#' @param is_control a binary vector specifying the controls
#' @param covs a data frame with additional covariates for adjustment (over the same samples as x)
#' @param control_tp the control time point to take, defualt is NULL, which means all
#' @param return_model a binary paramater: return the lm object or the summary with a vector with the effect size (log FC), its standard deviation, the t-statistic, and the p-value
#' @return the linear regression results, see the specification of return_model 
pass1a_simple_differential_abundance<-function(x,tps,tp,is_control,
          covs=NULL,control_tp = NULL,return_model=F){
  if(!is.null(control_tp)){
    is_control = is_control & tps == control_tp
  }
  inds = tps == tp | is_control
  inds[is.na(inds)] = F
  labels = tps == tp & !is_control
  labels[is.na(labels)] = F
  # df1 for the analysis
  df1 = data.frame(
     x = x[inds],
     tp = as.factor(labels[inds])
   )
   # add the covariates if not null
   if(!is.null(covs)){
     df1 = cbind(df1,covs[inds,])
   }
   lm1 = lm(x~.,data=df1)
   if(return_model){return(lm1)}
   res1 = summary(lm1)$coefficients[2,]
   names(res1) = c("Est","Std","Tstat","Pvalue")
   return(res1)
}

#' Simple comparison of two time points for a single analyte
#' 
#' @param x a numeric vector with continuous values (an analyte)
#' @param tps a vector with the time points
#' @param tp the time point for the current analysis
#' @param is_control a binary vector specifying the controls
#' @param covs a data frame with additional covariates for adjustment (over the same samples as x)
#' @param control_tp the control time point to take, defualt is NULL, which means all
#' @param return_model a binary paramater: return the lm object or the summary with a vector with the effect size (log FC), its standard deviation, the t-statistic, and the p-value
#' @return the linear regression results, see the specification of return_model 
simple_pairwise_differential_abundance<-function(x,tps,is_control,
  tp1,tp1_iscontrol,tp2,tp2_iscontrol,covs=NULL,return_model=F){
  inds = (tps == tp1 & is_control == tp1_iscontrol) |
    (tps == tp2 & is_control == tp2_iscontrol)
  inds[is.na(inds)] = F
  labels = (tps == tp1 & is_control == tp1_iscontrol)
  labels[is.na(labels)] = F
  labels = as.factor(as.numeric(labels))
  # df1 for the analysis
  df1 = data.frame(
    x = x[inds], y = labels[inds]
  )
  # add the covariates if not null
  if(!is.null(covs)){
    df1 = cbind(df1,covs[inds,])
  }
  lm1 = lm(x~.,data=df1)
  if(return_model){return(lm1)}
  res1 = c(
    "tp1" = tp1,"tp1_isctrl" = tp1_iscontrol,
    "tp2" = tp2, "tp2_isctrl" = tp2_iscontrol,
    summary(lm1)$coefficients[2,]
  )
  names(res1)[5:8] = c("Est","Std","Tstat","Pvalue")
  return(res1)
}
# # test
# a = rnorm(30)
# a_tps = sort(rep(1:3,length(a)/3))
# a_iscontrol = as.logical(rep(0:1,length(a)/2))
# table(a_tps,a_iscontrol)
# simple_pairwise_differential_abundance(a,a_tps,a_iscontrol,1,FALSE,1,TRUE)
# pass1a_simple_differential_abundance(a,a_tps,1,a_iscontrol,control_tp = 1)
# t.test(a[a_tps==1&!a_iscontrol],a[a_tps==1 & a_iscontrol])
# library(limma)
# inds = a_tps==1
# z = a_iscontrol[inds]
# mm = model.matrix(~1+as.numeric(z))
# lf = lmFit(rbind(a[inds],a[inds]),mm)
# topTable(eBayes(lf))

#' Estimate differential abundance between the control sets
#' 
#' @param x a numeric vector with continuous values (an analyte)
#' @param tps a vector with the time points
#' @param tp the time point for the current analysis
#' @param is_control a binary vector specifying the controls
#' @param covs a data frame with additional covariates for adjustment
#' @return a vector with the effect size (log FC), its standard deviation, the t-statistic, and the p-value
pass1a_controls_differential_abundance<-function(x,tps,tp,is_control,covs=NULL){
  df2 = data.frame(
    x = x[is_control],
    tp = as.factor(tps[is_control])
  )
  # add the covariates if not null
  if(!is.null(covs)){
    df2 = cbind(df2,covs[is_control,])
  }
  lm2 = lm(x~.,data=df2)
  res2 = summary(lm2)$coefficients[2,]
  names(res2) = c("Est","Std","Tstat","Pvalue")
  return(res2)
}

#' Regress out covariate information from a single score
#' @param x a numeric vector
#' @param covs a matrix of covariates
#' @return the residual vector
lm_regress_out_vector<-function(x,covs){
  df = data.frame(x=x,covs)
  lm_model = lm(x~.,data=df)
  return(residuals(lm_model))
}

#' Apply regress out to each column in a matrix
#' @param x a numeric vector
#' @param covs a matrix of covariates
#' @return the residual vector
lm_regress_out_matrix<-function(x,covs){
  newx = apply(x,2,lm_regress_out_vector,covs=covs)
  return(newx)
}

#' Compute linear association between a continuous variable y and another
#' variable x, which is transformed into a factor if it is discrete.
#' Z is a dataframe to be conditioned upon: compute y's residuals with it.
#' x,y,z, are all named vectors/matrices
#' @param x - vector (numeric or discrete)
#' @param y - numeric vector
#' @param z - optional - a data frame of variables to adjust for
linear_association_analysis<-function(x,y,z=NULL){
  if(!is.null(z)){
    df1 = data.frame(y=y,z)
    lm1 = lm(y~.,data=df1)
    y = lm1$residuals
  }
  if(!is.null(names(x))&&!is.null(names(y))){
    inds = intersect(names(x),names(y))
    x = x[inds];y=y[inds]
  }
  inds = !is.na(x) & !is.na(y)
  x = x[inds];y=y[inds]
  xcopy = as.numeric(as.character(x))
  if(sum(is.na(x))>length(x)/2){
    x = as.factor(x)
  }
  df = data.frame(y=y,x=x)
  lm_model1 = lm(y~x,data=df)
  lm_model0 = lm(y~1,data=df)
  lm_res = summary(lm_model1)
  lm_an = anova(lm_model1,lm_model0)
  pval = as.matrix(lm_an)[2,6]
  r2 = lm_res$r.squared
  rho = sqrt(r2)
  betas = lm_res$coefficients[-1,1]
  names(betas) = rownames(lm_res$coefficients)[-1]
  if(length(betas)==1 && betas[1]<0){rho=-rho}
  return(c(pval=pval,r2=r2,rho=rho,betas))
}

#' Go over the column pairs in m and run the function func on each one. 
#' Store the value from the field f in the function's output in the result
pairwise_eval_single_mat<-function(m,func,f,...){
  n = ncol(m)
  res = matrix(NA,nrow=n,ncol=n)
  colnames(res) = colnames(m)
  rownames(res) = colnames(m)
  for(i in 1:n){
    for(j in 1:i){
      curr_output = func(m[,i],m[,j],...)
      if(!is.numeric(curr_output[f])){
        curr_output = curr_output[[f]]
      }
      else{
        curr_output = curr_output[f]
      }
      res[i,j] = curr_output
      res[j,i] = res[i,j]
    }
  }
  return(res)
}

pairwise_eval_two_matrices<-function(m1,m2,func,f,...){
  n1 = ncol(m1);n2=ncol(m2)
  res = matrix(NA,nrow=n1,ncol=n2)
  colnames(res) = colnames(m2)
  rownames(res) = colnames(m1)
  for(i in 1:n1){
    for(j in 1:n2){
      curr_output = func(m1[,i],m2[,j],...)
      if(!is.numeric(curr_output[f])){
        curr_output = curr_output[[f]]
      }
      else{
        curr_output = curr_output[f]
      }
      res[i,j] = curr_output
    }
  }
  return(res)
}

pairwise_eval<-function(m,m1=NULL,func,f,...){
  if(is.null(m1)){return(pairwise_eval_single_mat(m,func,f,...))}
  return(pairwise_eval_two_matrices(m,m1,func,f,...))
}

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

#' Compare multiple groups using ANOVA's F-test or Kruskal-Wallis
#' @param x - a continuous vector
#' @param y - a discrete vector
compare_multigroup_means<-function(x,y,f=aov){
  df = data.frame(x=x,y=as.factor(y))
  test_res = f(x~y,data=df)
  if(is.element("aov",set=class(test_res))){
    test_res = summary(test_res)
    return(test_res[[1]]["y","Pr(>F)"])
  }
  return(test_res$p.value)
}

compare_two_groups<-function(x,y,f=t.test,...){
  if(sd(x)==0){return(c(x[1],x[1],1))}
  x1 = x[y==y[1]]
  x2 = x[y!=y[1]]
  if(sd(x1)==0){return(c(x1[1],mean(x2),1))}
  if(sd(x2)==0){return(c(mean(x1),x2[1],1))}
  res = f(x1,x2,...)
  return(c(res$estimate,res$p.value))
}

#' A wrapper for measuring the association between x and y
#' 
#' @description If x and y are discrete use chi-square test. 
#'     If one is discrete and the other is continuous use mean group comparison using t- or f- test. 
#'     If both are continuous use the correlation test.
#' @return the association p-value
pairwise_association_wrapper<-function(x,y,max_num_vals_for_discrete=5){
  is_contin_x = is.numeric(x) && length(unique(x))>max_num_vals_for_discrete
  is_contin_y = is.numeric(y) && length(unique(y))>max_num_vals_for_discrete
  
  if(length(unique(x))<2){return (1)}
  if(length(unique(y))<2){return (1)}
  
  if(is_contin_y && is_contin_x){
    return(cor.test(x,y)$p.value)
  }
  if(!is_contin_y && is_contin_x){
    return(compare_multigroup_means(x,y))
  }
  if(is_contin_y && !is_contin_x){
    return(compare_multigroup_means(y,x))
  }
  tb = table(x,y)
  return(chisq.test(tb)$p.value)
}

# # test
# n=100
# x = rbinom(n,1,0.3)
# y = rnorm(n,sd = 0.2) + x
# z = matrix(rnorm(n*10),ncol=10,nrow=n)
# rownames(z) = 1:n
# names(x) = 1:n
# names(y) = 1:n
# cor(x,y)
# cor.test(x,y,method="spearman")
# lmtest_res = linear_association_analysis(x,y,z)
# abs(lmtest_res["rho"] - cor(x,y))
# corrs = pairwise_eval(z,linear_association_analysis,"rho")
# corrs2 = cor(z)
# hist(corrs - corrs2)
# sum(corrs - corrs2 > 0.00001)
