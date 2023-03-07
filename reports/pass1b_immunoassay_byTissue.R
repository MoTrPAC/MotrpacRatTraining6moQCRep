# load packages
library(dplyr)
library(tidyverse)
library(stringr)

# TODO: replace the path with gcloud bucket
setwd('/Users/mshruti/Desktop/Stanford/MoTrPAC_Luminex/data_freeze/normalized-data')

# extract tissue and panel from a given file name
get_tissue_panel <- function(filename){
  filename_split <- strsplit(filename, "_")
  tissue <- filename_split[[1]][2]
  panel <- filename_split[[1]][4]
  return(c(tissue,panel))
}


# for a given output file (panel and tissue specific), extract bid and add panel name as a column
preprocess <- function(filename){
  #I. fetch first 5 digits from vial_label which should be the first column and always start with 9 for experimental samples
  df <- read.delim2(filename, stringsAsFactors = T)
  bid <- sapply(df$vial_label, function(i) {str_extract(i,"^\\d{5}")})
  new_df <- cbind(bid,df[, !(names(df) %in% c('vial_label'))])
  
  #II. Transpose
  new_df_transpose <- setNames(data.frame(t(new_df[,!(names(new_df) %in% c('bid'))])), new_df[,c('bid')])
  
  #III. add panel name
  file_metadata <- get_tissue_panel(filename)
  panel_name <- file_metadata[2]
  panel <- c(rep(panel_name, nrow(new_df_transpose) ))
  tissue_name <- file_metadata[1]
  tissue <- c(rep(tissue_name, nrow(new_df_transpose) ))
  new_df2 <- cbind(feature_ID=rownames(new_df_transpose),panel, tissue,new_df_transpose, row.names=NULL)
  
  return(new_df2)
}


# list of tissue names, each element contains output file names of the panels that were run on that tissue
filenames_by_tissue <- function(filepath, file_pattern){
  tissue_list <- list()
  
  for(f in list.files(path=filepath, pattern=file_pattern)){
    tissue_name <- get_tissue_panel(f)[1]
    tissue_list[[tissue_name]] <- c(tissue_list[[tissue_name]],f)
  }
  return(tissue_list)
}


merge_outputs_by_tissue <- function(filepath, file_pattern){
  list_of_tissues <- filenames_by_tissue(filepath, file_pattern)

  for(tissues in names(list_of_tissues)){
    merged_tissue_output <- data.frame()
    for(file_name in list_of_tissues[[tissues]]){
      print(file_name)
      df <- preprocess(file_name)
      #use bind_rows instead of rbind to allow merging of dataframes with different number of columns (samples)
      merged_tissue_output <- bind_rows(merged_tissue_output,df)
    }
    output_name <- paste('motrpac', gsub('-','',Sys.Date()), 'pass1b-06', tissues, 'immunoassay_merged', strsplit(file_pattern,"\\*_")[[1]][2], sep="_")
    # TODO: replace the path with gcloud bucket
    write.table(x=merged_tissue_output, file=paste0("../by_tissue/",output_name), sep="\t", row.names=FALSE)
    print(cat("\n"))
  }
}


#file_patterns <- c("pass1b-06_t*_immunoassay_*_mfi-log2-filt-imputed-na-outliers.txt","pass1b-06_t*_immunoassay_*_mfi-log2-filt-imputed.txt")
file_patterns <- c("*_mfi-log2-filt-imputed-na-outliers.txt","*_mfi-log2-filt-imputed.txt")

for(p in file_patterns){
  print(p)
  merge_outputs_by_tissue(filepath=".", file_pattern=p)
  print("##############")
}

# few samples in t61-colon_immunoassay_rat-mag27plex panel were removed due to low quality issues. 
# But these are present in t61-colon_immunoassay_rat-metabolic panel and thus you see some NA in merged data
colon_rat_mag27plex_imputed <- preprocess('pass1b-06_t61-colon_immunoassay_rat-mag27plex_mfi-log2-filt-imputed.txt')
colon_rat_metabolic_imputed <- preprocess("pass1b-06_t61-colon_immunoassay_rat-metabolic_mfi-log2-filt-imputed.txt")
setdiff(colnames(colon_rat_metabolic_imputed), colnames(colon_rat_mag27plex_imputed))
