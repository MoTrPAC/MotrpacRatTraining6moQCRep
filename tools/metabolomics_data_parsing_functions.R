#' Parse a downloaded metabolomics dataset, restricting the input to files with the regex in their names.
#' 
#' @param downloaded_files a list of local files of a single dataset.
#' @param regex a string. Can be "_named" or "_unnamed", 
#' which tells the function which part of the dataset to take.
#' @param platform a string. The name of the current platform.
#' @param site a string. The name of the current site that had generated the dataset.
#' @details
#' Get the sample data and sample metadata files and read their tables. 
#' We assume that the first column in the sample metadata file contains the sample ids.
#' We assume that these are unique for non-control samples.
#' For control samples append the sample order to make their names unique.
#' The method removes data rows that are all NAs.
#' @return a list with the following fields:
#' 1. sample_data a data frame of metabolites \times samples
#' 2. control_data a data frame of metabolites \times control samples
#' 3. sample_meta a data frame of samples \times metadata fields
#' 4. row_annot a data frame of metabolite \times annotation fields
parse_metabolomics_dataset_by_regex<-function(downloaded_files,regex="_named",
                                              platform = "",site=""){
  downloaded_files = downloaded_files[grepl(regex,downloaded_files)]
  # get the files
  data_file = downloaded_files[grepl("-results\\.",downloaded_files)]
  if(is.na(data_file) || length(data_file)!=1){
    data_file = downloaded_files[grepl("_results\\.",downloaded_files)]
  }
  if(is.na(data_file) || length(data_file)!=1){
    data_file = downloaded_files[grepl("results\\.",downloaded_files)]
  }
  if(is.na(data_file) || length(data_file)!=1){
    data_file = downloaded_files[grepl("results.",downloaded_files)][1]
  }
  sample_metadata_file = downloaded_files[grepl("sample",downloaded_files) &
                                            grepl("metadata",downloaded_files)]
  metab_metadata_file = downloaded_files[grepl("metabolite",downloaded_files) &
                                           grepl("metadata",downloaded_files)]
  
  if(length(sample_metadata_file)!=1 || length(metab_metadata_file)!=1 || length(data_file)!=1){
    print("WARNING: sample or meta data files are not unique when analyzing:")
    print(downloaded_files)
    sample_metadata_file = sample_metadata_file[1]
    metab_metadata_file = metab_metadata_file[1]
    data_file = data_file[1]
  }
  
  # read the data in
  # use fread
  data = fread(data_file,data.table = F,stringsAsFactors = F,header = T)
  row_annot = fread(metab_metadata_file,data.table = F,stringsAsFactors = F,header = T)
  rownames(row_annot) = row_annot[[1]]
  sample_meta = fread(sample_metadata_file,data.table = F,stringsAsFactors = F,header = T)
  sample_meta_rownames = as.character(sample_meta[,1])
  # for control samples concatenate the order - will make it unique
  currinds = tolower(sample_meta[,2])!="sample_type"
  sample_meta_rownames[currinds] = paste(sample_meta[currinds,1],"order",sample_meta[currinds,3],sep="_")
  rownames(sample_meta) = sample_meta_rownames
  
  # remove rows that are all NAs
  to_rem = apply(is.na(data[,-1]),1,all)
  data = data[!to_rem,]
  row_annot = row_annot[!to_rem,]
  
  # avoid issues with parsing data from non-standard encodings
  # also, take care of checking the first column
  is_numeric = sapply(data,is.numeric)
  newdata = data[,is_numeric]
  newdata = recast_numeric_data_frame(newdata)
  newdata = as.data.frame(newdata)
  if(is_numeric[1]){
    print("Warning: first column in dataset is numeric, example:")
    print(data[1,1])
    print("Make sure these are metabolite names or ids")
  }
  if(sum(!is_numeric) > 1){
    print(paste("In dataset:",tissue,platform," more than a single non-numeric column"))
    print("Head of these columns:")
    print(head(data[,is_numeric]))
  }
  if(sum(!is_numeric)>0){
    tmp = data
    data = cbind(data[,!is_numeric],newdata)
    # fix the names
    names(data)[1:sum(!is_numeric)] = names(tmp)[!is_numeric]
    rm(tmp)
  }

  # split data into sample and controls
  qc_samples = sample_meta[tolower(sample_meta[,"sample_type"]) != "sample",1]
  # make sure both keep the dirst column
  sample_data = data[,union(1,which(!is.element(colnames(data),set=qc_samples)))]
  control_data = data[,union(1,which(is.element(colnames(data),set=qc_samples)))]
  
  # Add site and platform to metadata
  sample_meta$site = rep(site,nrow(sample_meta))
  sample_meta$platform = rep(platform,nrow(sample_meta))
  
  return(list(sample_data = sample_data, control_data = control_data,
              sample_meta = sample_meta, row_annot = row_annot))
}

#' Parse a single metabolomics dataset
#' 
#' @param download_obj a download object from download_bucket_to_local_dir. 
#'   Provides paths through download_obj$downloaded_files to all the metabolomics targeted or untargeted files in the bucket
#' @param local_path a string. The name of the local directory where the dowloaded data are to be kept
#' @param delete a logical. TRUE: delete the data after downloading; FALSE: keep the local files
#' @param platform a string. The name of the assay platform. If NULL take it from the bucket name.
#' @param verbose a logical. TRUE: print comments to the user.
#' 
#' @details 
#' Search for the experimentaldetails file and use it to obtain the site name. If this field
#' does not exist then an error is printed if verbose=T but we do not stop the function.
#' Given the downloaded data and parsed site and platform information the function parses the 
#' named and unnamed parts of the dataset separately (if unnamed exists) using the
#' parse_metabolomics_dataset_by_regex function.
read_single_metabolomics_dataset_from_bucket<-function(download_obj,tissue,platform,verbose=T){
  downloaded_files = download_obj$downloaded_files
  downloaded_files = downloaded_files[
    grepl(tissue,downloaded_files) & grepl(platform,downloaded_files)
  ]
  # get site from the experimental details file
  if(any(grepl("experimentaldetails",downloaded_files,ignore.case = T))){
    site = get_site(
      read_experimental_details(downloaded_files,verbose = verbose),
      verbose = verbose)
  }
  else{
    if(verbose){
      print("Error: no experimental details file in bucket, assigning empty string to site field")
    }
    site = ""
  }
  
  # read the named data
  named_obj = parse_metabolomics_dataset_by_regex(
    downloaded_files,regex="_named",site=site,platform=platform)
  # if exists, read the unnamed data
  unnamed_obj = NULL
  if(any(grepl("metab-u",downloaded_files))){
    unnamed_obj = parse_metabolomics_dataset_by_regex(
      downloaded_files,regex="_unnamed",site=site,platform=platform)
  }
  return(list(named=named_obj,unnamed=unnamed_obj))
}

#' Auxilary function for reading the experimental details file.
read_experimental_details<-function(downloaded_files,verbose=T){
  exp_details_file = downloaded_files[(
    grepl("experimental",downloaded_files,ignore.case = T) &
      grepl("details",downloaded_files,ignore.case = T))][1]
  if(length(exp_details_file)<1){
    if(verbose){
      print("Error! experimental details file not found")
    }
    return(NULL)
  }
  exp_details = readLines(exp_details_file)
  return(exp_details)
}

#' Auxiliary function for obtaining the site names.
get_site<-function(exp_details,verbose=T){
  st_details = exp_details[grepl("ST:",exp_details)]
  st_details = st_details[grepl("INST",st_details)]
  if(length(st_details)<1){
    if(verbose){
      print("Error: institute name does not appear in the experimental details file, setting to \"\" ")
    }
    return("")
  }
  st_details = strsplit(st_details,split="\t")[[1]][2]
  return(st_details)
}

#' Load all metabolomics datasets from the downloaded files
#' 
#' @param download_obj a download object from download_bucket_to_local_dir. 
#'   Provides paths through download_obj$downloaded_files to all the metabolomics targeted or untargeted files in the bucket
#' @details
#' The assumed structure: root (input bucket)-> tissues -> platforms -> datasets
#' The function goes over the bucket tree and extracts the information of each dataset using
#' the read_single_metabolomics_dataset_from_bucket function.
#' @return a list of datasets, 
#' see the read_single_metabolomics_dataset_from_bucket documentation for the description
#' of each object.
read_metabolomics_datasets_from_download_obj<-function(download_obj,verbose=T){
  downloaded_files = download_obj$downloaded_files
  local_path_reg = gsub("~/","",download_obj$local_path)
  bucket_content = strsplit(downloaded_files,split=local_path_reg)
  bucket_content = lapply(bucket_content,function(x)
    strsplit(x[2],split=.Platform$file.sep)[[1]])
  # this is crucial: we assume that bucket structure is omics/tissue/assay/file, which
  # provides a vector of size five when splitting using "/
  bucket_content = bucket_content[sapply(bucket_content,length)==5]
  tissues = unlist(sapply(bucket_content,function(x)x[length(x)-2]))
  platforms = unlist(sapply(bucket_content,function(x)x[length(x)-1]))
  tissue_plarform_pairs = unique(cbind(tissues,platforms))
  tissue_plarform_pairs = 
    tissue_plarform_pairs[!grepl("qa-qc",tissue_plarform_pairs[,2]),]
  
  metabolomics_parsed_datasets = list()
  failed_buckets = c()
  for(j in 1:nrow(tissue_plarform_pairs)){
    tissue = tissue_plarform_pairs[j,1]
    platform = tissue_plarform_pairs[j,2]
    dataset_name = paste(c(tissue,platform),collapse=",")
    if(verbose){
      print(paste("parsing dataset:",dataset_name))
    }
    # read all datasets in current platform bucket - named and unnamed (if exists)
    # try reading first, if files are missing we will get an error
    datasets = NULL
    try({
      datasets = read_single_metabolomics_dataset_from_bucket(download_obj,tissue,platform)
    })
    if(is.null(datasets)){
      failed_buckets = c(failed_buckets,dataset_name)
      if(verbose){
        print(paste("Error, could not parse the data from:", dataset_name))
        print("Possible reason: missing data or metadata files, please revise")
      }
      next
    }
    if(verbose){print("done")}
    # add new datasets to the list of objects
    for(nn in names(datasets)){
        if(is.null(datasets[[nn]])){next}
        metabolomics_parsed_datasets[[paste(dataset_name,nn,sep=",")]] = datasets[[nn]]
    }
  }
  return(list(metabolomics_parsed_datasets=metabolomics_parsed_datasets,
              failed_buckets=failed_buckets))
}


#' Check for issues in a metabolomics parsed dataset d
#' 
#' @param d a data frame, output of parse_metabolomics_dataset_by_regex
#' @param rownames_ind, the index of the column in d$sample_data that has the metabolite names
#' NULL means that rownames(d$sample_data) can be used
#' @details 
#' The method checks for a few issues in the dataset and resurns a character vector with 
#' the description of the errors (if there are any). Specifically, we check if
#' 1. metabolite names are unique in the metabolite names column (warning if not)
#' 2. all samples appear in the metadata file (error)
#' 3. annotation data and sample data do not have the same nrow (warning)
#' 4. there are metabolites in the data that do not have annotation information (error)
#' 5. there are metabolite names that are NA
check_metabolomics_dataset_issues<-function(d,rownames_ind=1){
  if(is.null(rownames_ind)){
    rnames = rownames(d$sample_data)
  }
  else{
    rnames = d$sample_data[[rownames_ind]]
  }
  res = NULL
  tb = table(rnames)
  if(any(tb>1)){
    res = c(res,paste("Warning:",sum(tb>1), "metab names are not unique in data"))
  }
  if(!all(is.element(colnames(d$sample_data)[-1],
                     set=as.character(d$sample_meta[,1])))){
    res = c(res,"Error: not all samples are in the metadata file")
  }
  if(nrow(d$row_annot) != nrow(d$sample_data)){
    res = c(res,"Warning: annotation and data do not have the same nrow")
  }
  if(!all(rnames %in% d$row_annot[,1])){
    res = c(res,"Error: not all metabolite names have annotation data")
  }
  if(!setequal(rnames,d$row_annot[,1])){
    res = c(res,"Error: metabolite sets in annotation and data do not match")
  }
  na_mets = sum(is.na(rnames))
  if(na_mets>0){
    res = c(res,paste(na_mets,"metabolite names are NAs"))
  }
  return(res)
}

#' A function for merging named and unnamed datasets from the same site and platform.
#' 
#' @param named_data a named list with the named metabolomics dataset
#' @param unnamed_data a named list with the unnamed metabolomics dataset
#' @strict a logical variable, whether to enforce strict perfect agreement between the datasets
#' @details 
#' The function first runs a few tests to make sure that the datasets can be merged, including 4 tests:
#' 1. do the sample sets match?
#' 2. do the control sample sets match?
#' 3. are there discrepancies between the sample metadata tables?
#' If strict is TRUE and there are errors detected then the method stops and returns 
#' a description of the errors.
#' Merging the datasets is done by averaging metabolite rows if they appear in
#' both datasets. Also, the row annotation tables are merged and then reordered to fit
#' the new merged data table.
merge_named_and_unnamed_metabolomics_datasets<-function(named_data,unnamed_data,
                                                        strict = F){
  # three types of error: (1) not the same sample set, (2) not the same control sample set, and
  # (3) not the same metadata tables
  err_message = c()
  if(!all(is.element(colnames(named_data$sample_data),
                     set = colnames(unnamed_data$sample_data)))){
    err_message["1_named"] = 
      "not all samples in named data appear in unnamed data"
  }
  if(!all(is.element(colnames(unnamed_data$sample_data),
                     set = colnames(named_data$sample_data)))){
    err_message["1_unnamed"] = 
      "not all samples in unnamed data appear in named data"
  }
  if(!all(is.element(colnames(named_data$control_data),
                     set = colnames(unnamed_data$control_data)))){
    err_message["2_named"] = 
      "not all control samples in named data appear in unnamed data"
  }
  if(!all(is.element(colnames(unnamed_data$control_data),
                     set = colnames(named_data$control_data)))){
    err_message["2_unnamed"] = 
      "not all control samples in unnamed data appear in named data"
  }
  
  sample_meta1 = named_data$sample_meta
  sample_meta2 = unnamed_data$sample_meta
  if(any(dim(sample_meta1)!=dim(sample_meta2)) ||
     !all(as.matrix(sample_meta1) == as.matrix(sample_meta2))){
    err_message["3"] = "sample metadata is not the same in both datasets, cannot merge"
  }
  
  if(strict && length(err_message)>0){
    stop(paste("Error:",paste(err_message,collapse=",")))
  }
  if(!strict && length(err_message)>0 && any(grepl("1_",names(err_message)))){
    stop(paste("Error:",paste(err_message,collapse=",")))
  }
  
  shared_meta_set = intersect(rownames(sample_meta1),rownames(sample_meta2))
  shared_meta = sample_meta1[shared_meta_set,]
  
  # merge row names by averaging if there are any duplications
  s1 = named_data$sample_data
  s2 = unnamed_data$sample_data
  s1 = s1[,colnames(s2)]
  sample_data = aggregate(rbind(s1,s2),list(c(rownames(s1),rownames(s2))),mean,na.rm=T)
  rownames(sample_data) = sample_data[,1]
  sample_data = sample_data[,-1]
  
  # merge row annotation but make sure we keep the correct order
  r1 = named_data$row_annot
  r2 = unnamed_data$row_annot
  for(colname in setdiff(colnames(r1),colnames(r2))){
    r2[colname] = NA
  }
  row_annot = rbind(r1,r2[,colnames(r1)])
  row_annot = row_annot[rownames(sample_data),]
  
  c1 = named_data$control_data
  c2 = unnamed_data$control_data
  shared_controls = intersect(colnames(c2),colnames(c1))
  control_data = aggregate(rbind(c1[,shared_controls],c2[,shared_controls]),
                           list(c(rownames(c1),rownames(c2))),mean,na.rm=T)
  rownames(control_data) = control_data[,1]
  control_data = control_data[,-1]
  
  res = list(sample_data = sample_data,control_data = control_data,
             row_annot = row_annot,sample_meta = shared_meta)
  
  return(res)
}
