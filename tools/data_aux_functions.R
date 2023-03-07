#' Aggregate repeated samples
#' @param x a numeric matrix with samples as columns and analytes as rows
#' @param g a vector with the grouping of the samples
#' @param func a function used to merge
#' @param ... additional input to func
aggregate_repeated_samples<-function(x,g,func=mean,...){
  x = t(x)
  df = data.frame(x)
  newdf = aggregate(df,list(g),simplify=T,func,...)
  newx = as.matrix(t(newdf[,-1]))
  colnames(newx) = newdf[,1]
  rownames(newx) = colnames(x)
  return(newx)
}
# # tests
# x = matrix(rnorm(200),20,10)
# rownames(x) = as.character(1:nrow(x))
# colnames(x) = as.character(1:ncol(x))
# newx = aggregate_repeated_samples(cbind(x,x),c(colnames(x),colnames(x)))
# all(newx[,colnames(x)] == x)
# all(rownames(x)==rownames(newx))

#' Aggregate repeated rows (analytes)
#' @param x a numeric matrix with samples as columns and analytes as rows
#' @param g a vector with the grouping of the metabolites
#' @param func a function used to merge
#' @param ... additional input to func
aggregate_repeated_rows<-function(x,g,func=mean,...){
  df = data.frame(x)
  newdf = aggregate(df,list(g),func,...)
  newx = as.matrix(t(newdf[,-1]))
  colnames(newx) = newdf[,1]
  return(newx)
}

#' A function to select the top results within groups.
#' The data set d is assumed to have two column sets: group_cols and score_col.
#' Within each group in d[,group_cols] we take (up to) the top k rows, 
#' ordered by their score.
#' @param group_cols - a vector with column names or inidices to be concatenated for the group definition
get_top_by_group<-function(d,group_cols,score_col,k,decreasing=T){
  newd = c()
  if(length(group_cols)<2){
    groupv = d[,group_cols]
  }
  else{
    groupv = apply(d[,group_cols],1,paste,collapse=",")
  }
  for(g in unique(groupv)){
    currd = d[groupv==g,]
    if(is.null(dim(currd))){
      newd = rbind(newd,currd)
      next
    }
    scoresv = as.numeric(as.character(currd[,score_col]))
    ord = order(scoresv,decreasing = decreasing)
    currk = min(k,length(ord))
    newd = rbind(newd,currd[ord[1:currk],])
  }
  rownames(newd) = NULL
  return(newd)
}
