###############################################################
required_libs = c("data.table")
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}

###############################################################
# Set a mapping from tissue codes to their names
tissue_code2name = c(
  "t30"="blood-rna",
  "t31"="plasma",
  "t32"="packed-cells",
  "t33"="hippocampus",
  "t34"="cortex",
  "t35"="hypothalamus",
  "t36"="gastrocnemius",
  "t37"="vastus-lateralis",
  "t38"="tibia",
  "t39"="heart",
  "t40"="kidney",
  "t41"="adrenals",
  "t42"="colon",
  "t43"="spleen",
  "t44"="testes",
  "t45"="ovaries",
  "t46"="brown-adipose",
  "t47"="white-adipose",
  "t48"="aorta",
  "t49"="lung",
  "t50"="small-intestine",
  "t51"="liver",
  "t52"="hippocampus-powder",
  "t53"="cortex-powder",
  "t54"="hypothalamus-powder",
  "t55"="gastrocnemius-powder",
  "t56"="vastus-lateralis-powder",
  "t57"="tibia-powder",
  "t58"="heart-powder",
  "t59"="kidney-powder",
  "t60"="adrenal-powder",
  "t61"="colon-powder",
  "t62"="spleen-powder",
  "t63"="testes-powder",
  "t64"="ovaries-powder",
  "t65"="aorta-powder",
  "t66"="lung-powder",
  "t67"="small-intestine-powder",
  "t68"="liver-powder",
  "t69"="brown-adipose-powder",
  "770"="white-adipose-powder"
)
tissue_code2name = tolower(tissue_code2name)
###############################################################

#' A general function for recasting numeric data frame through a string representation.
#' This is useful as a result of a peculiar case in which doing an as.matrix on a data
#' frame distorts the data.
recast_numeric_data_frame<-function(x){
  newx = sapply(x,as.character)
  if(is.null(dim(newx))){
    newx = matrix(newx,nrow=nrow(x))
  }
  colnames(newx) = colnames(x)
  rownames(newx) = rownames(x)
  newx = as.matrix(newx)
  mode(newx) = "numeric"
  return(as.data.frame(newx))
}

###############################################################

#' Load an RData file from a bucket path
load_from_bucket<-function(file,bucket,delete=T,GSUTIL_PATH="~/google-cloud-sdk/bin/gsutil"){
  system(paste(GSUTIL_PATH,"cp",
               paste(bucket,file,sep=""),
               getwd()))
  objects_before = ls()
  load(file)
  objects_after = ls()
  if(delete){
    system(paste("rm",file))
  }
  added_objects = setdiff(objects_after,objects_before)
  added_objects = setdiff(added_objects,"objects_before")
  res = list()
  for(obj_name in added_objects){
    res[[obj_name]] = get(obj_name)
  }
  return(res)
}

#' Download a single file from a bucket path and read it into R
dl_read_gcp = function(path, OUTDIR='.', sep='\t', GSUTIL_PATH="~/google-cloud-sdk/bin/gsutil", df=FALSE, compressed=FALSE){
  if(grepl("\\.gz$",path)){
    compressed=TRUE 
  }
  # download
  cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, OUTDIR)
  system(cmd,ignore.stdout = T,ignore.stderr = T)
  # read in the data
  new_path = sprintf('%s/%s',OUTDIR,basename(path))
  if(compressed){
    # read with fread even if we want to return a data frame since it's so much faster
    dt = fread(cmd=sprintf('zcat %s',new_path), sep=sep, header=T)
    if(df){
      dt = data.frame(dt)
    }
  }else{
    if(df){
      dt = read.table(new_path,sep=sep,header=T,stringsAsFactors = F)
    }else{
      dt = fread(new_path,sep=sep,header=T)
    }
  }
  return(dt)
}

#' Download all files in bucket to a local dir
download_bucket_to_local_dir<-function(bucket,
                                       local_path = NULL,
                                       remove_prev_files = TRUE,
                                       GSUTIL_PATH = "~/google-cloud-sdk/bin/gsutil",
                                       recursive_download = FALSE){
  
  # make sure the bucket name does not end with "/"
  if(grepl("/$",bucket,perl=T)){
    bucket = substr(bucket,1,nchar(bucket)-1)
  }
  
  if(is.null(local_path)){local_path=paste(getwd(),"gs_download",sep="/")}
  # download files
  # remove old data if exists
  if(remove_prev_files){
    if(file.exists(local_path)){system(paste("rm -r",local_path))}
  }
  # download
  if(!dir.exists(local_path)){
    dir.create(local_path, recursive = TRUE)
  }
  if(recursive_download){
    cmd = paste(GSUTIL_PATH,"cp -r",bucket,local_path)
  }
  else{
    cmd = paste(GSUTIL_PATH,"cp",bucket,local_path)
  }
  system(cmd)
  # get the data
  downloaded_files = list.files(local_path,full.names = T,recursive = T)
  return(list(downloaded_files=downloaded_files,local_path=local_path))
}

#' Download all files in bucket to a local dir
get_single_file_from_bucket_to_local_dir<-function(bucket,
                                                   local_path=NULL,
                                                   remove_prev_files = TRUE,
                                                   GSUTIL_PATH="~/google-cloud-sdk/bin/gsutil"){
  if(is.null(local_path)){local_path=paste(getwd(),"gs_download",sep="/")}
  # download files
  # remove old data if exists
  if(remove_prev_files){
    if(file.exists(local_path)){system(paste("rm -r",local_path))}
  }
  # download
  system(paste("mkdir",local_path))
  cmd = paste(GSUTIL_PATH," -m cp",bucket,local_path)
  system(cmd,ignore.stdout = T,ignore.stderr = T)
  # get the data
  downloaded_files = list.files(local_path,full.names = T)
  return(list(downloaded_files=downloaded_files,local_path=local_path))
}

#' Save R objects into an RData file in a given bucket
save_to_bucket<-function(...,file,bucket,
                         GSUTIL_PATH="~/google-cloud-sdk/bin/gsutil"){
  save(...,file = file)
  system(paste(GSUTIL_PATH,"cp",
               file, bucket))
  system(paste("rm",file))
}

#' Get a list of files in a bucket
get_files_in_bucket<-function(bucket,
             GSUTIL_PATH="~/google-cloud-sdk/bin/gsutil"){
  tmp_name = paste("tmp",runif(1,1,1000),".txt",sep="")
  system(paste(GSUTIL_PATH,"ls",
               bucket,">",tmp_name))
  files = readLines(tmp_name)
  system(paste("rm",tmp_name))
  return(files)
}

#' Write single file to bucket
write_to_bucket<-function(object, file_name, output_bucket, tmpdir='.', sep='\t', col.names=T, row.names=F, quote=F){
  if(grepl("/",file_name)){
    stop("File name should not include a path.")
  }
  local_file = sprintf("%s/%s", gsub("/$", "", tmpdir), file_name)
  gcp_file = sprintf("%s/%s", gsub("/$", "", output_bucket), file_name)
  write.table(object, file=local_file, sep=sep, col.names=col.names, row.names=row.names, quote=quote)
  system(sprintf("gsutil cp %s %s", local_file, gcp_file))
  system(sprintf("rm %s", local_file))
  message(system(sprintf("gsutil ls %s", gcp_file), intern=T))
}

###############################################################
###############################################################

# for GET
# Assumption - suffix is tissue->omic->results->data
GET_get_tissue_from_download_path<-function(p){
  arr = strsplit(p,split = .Platform$file.sep)[[1]]
  return(arr[(length(arr)-3)])
}
GET_get_assay_from_download_path<-function(p){
  arr = strsplit(p,split = .Platform$file.sep)[[1]]
  return(arr[(length(arr)-2)])
}

get_tissue_from_download_path<-function(p){
  arr = strsplit(p,split = .Platform$file.sep)[[1]]
  return(arr[(length(arr)-2)])
}
get_assay_from_download_path<-function(p){
  arr = strsplit(p,split = .Platform$file.sep)[[1]]
  return(arr[(length(arr)-1)])
}


###############################################################
###############################################################
# Pheno
parse_pheno_data<-function(bucket,use_long_field_names=T,...){
  obj = download_bucket_to_local_dir(bucket,...)
  dict_file = obj$downloaded_files[
    grepl("merged",obj$downloaded_files) &
      grepl("dictionary",obj$downloaded_files) &
      grepl(".txt",obj$downloaded_files)
  ]
  if(length(dict_file)==0){
    print("Error in parse_pheno_data: no dict_file file")
    return(NULL)
  }
  if(length(dict_file)>1){
    print("Error in parse_pheno_data: more than one dict_file files")
    return(NULL)
  }
  
  viallabel_data_files = obj$downloaded_files[
    grepl("viallabel",obj$downloaded_files) &
      grepl("data",obj$downloaded_files) &
      grepl(".txt",obj$downloaded_files)
  ]
  if(length(viallabel_data_files)==0){
    print("Error in parse_pheno_data: no viallabel_data_file file")
    return(NULL)
  }
  if(length(viallabel_data_files)>1){
    print("Note: more than one viallabel_data_file, merging using rbind")
  }
  
  dict = fread(dict_file,stringsAsFactors=F,data.table=F,header=T,sep="\t")
  rownames(dict) = dict[,1]
  
  viallabel_data = NULL
  for(viallabel_data_file in viallabel_data_files){
    curr_viallabel_data = fread(
      viallabel_data_file,stringsAsFactors=F,data.table=F,header=T,sep="\t")
    if(is.null(viallabel_data)){
      viallabel_data = curr_viallabel_data
    }
    else{
      viallabel_data = rbind(viallabel_data,curr_viallabel_data)
    }
  }
  
  
  vial_col = rownames(dict)[grepl("viallabel",dict$Field.Name)]
  rownames(viallabel_data) = viallabel_data[[vial_col]]
  if(use_long_field_names){
    tp_field_arrs = strsplit(colnames(viallabel_data),split="_")
    main_field = sapply(tp_field_arrs,function(x)x[1])
    version_field = sapply(tp_field_arrs,function(x)x[length(x)])
    version_field[version_field==main_field]="0"
    version_field = as.numeric(version_field)
    newcolnames = tolower(dict[main_field,"FullName"])
    is_visit_field = version_field>0
    newcolnames[is_visit_field] = 
      paste0(newcolnames[is_visit_field],"_visit",version_field[is_visit_field])
    rownames(dict) = tolower(dict$FullName)
    colnames(viallabel_data) = newcolnames
  }
  return(list(
    dictionary = dict,
    viallabel_data = viallabel_data
  ))
}

###############################################################
###############################################################
# RNA-seq
parse_rnaseq_data<-function(bucket,matrix_type="genes-count",...){
  obj = download_bucket_to_local_dir(bucket,...)
  data_table_paths = obj$downloaded_files[grepl(matrix_type,obj$downloaded_files)]
  data_tissues = sapply(data_table_paths,GET_get_tissue_from_download_path)
  data_tables = lapply(data_table_paths,fread,stringsAsFactors=F,data.table=F,header=T)
  if(all(table(data_tissues)==1)){
    names(data_tables) = data_tissues
  }
  else{
    print("Warning in parse_rnaseq_data: more than a single table per tissue, keeping path as the table names")
  }
  
  for(nn in names(data_tables)){
    x = data_tables[[nn]]
    if(any(table(x[,1])>1)){
      print("Warning in parse_rnaseq_data: gene names are not unique in dataset:")
      print(nn)
    }
    else{
      rownames(x) = x[,1]
      x = x[,-1]
    }
    data_tables[[nn]] = x
  }
  gc()
  
  qc_metrics_file = obj$downloaded_files[grepl("qc-metrics",obj$downloaded_files)]
  if(length(qc_metrics_file)==0){
    print("Warning in parse_rnaseq_data: no qc-metrics file")
  }
  if(length(qc_metrics_file)>1){
    print("Warning in parse_rnaseq_data: more than one qc-metrics files")
  }
  
  sample_metadata = NULL
  if(length(qc_metrics_file)==1){
    sample_metadata = fread(qc_metrics_file,stringsAsFactors=F,data.table=F,header=T)
    if(any(table(sample_metadata[,1])>1)){
      print("Error in parse_rnaseq_data: vial labels (first column) in qc_metrics are not unique")
    }
    rownames(sample_metadata) = sample_metadata[,1]
  }
  
  return(list(
    sample_metadata = sample_metadata,
    data_tables = data_tables,
    data_table_paths = data_table_paths
  ))
}

###############################################################
###############################################################

parse_proteomics_dataset_by_regex<-function(downloaded_files,regex="ratio"){
  
  data_file = downloaded_files[grepl(regex,downloaded_files) &
                                 grepl("results",downloaded_files)]
  metadata_file = downloaded_files[grepl("vial",downloaded_files) &
                                     grepl("metadata",downloaded_files)]
  feature_metadata_file = downloaded_files[!grepl("vial",downloaded_files) &
                                     grepl("meta",downloaded_files)]
  
  if(length(feature_metadata_file)>0){
    print(paste("Identified the feature metadata file:",feature_metadata_file))
  }
  else{
    print(paste("Using the annotation columns in the data file:",data_file))
  }
  
  if(length(metadata_file)==0){
    metadata_file = downloaded_files[grepl("metadata",downloaded_files)]
  }
  
  data = fread(data_file,data.table = F,stringsAsFactors = F,header = T)
  
  #Identify columns that are not data (i.e., feature annotations)
  row_annot_cols <- grep("^9",colnames(data),invert = T)

  
  if(length(feature_metadata_file)>0){
    #print(paste("Reading feature annotation from a metadata file:",feature_metadata_file))
    row_annot = fread(feature_metadata_file[1],data.table = F,stringsAsFactors = F,header = T)
    row_annot_from_data_file = data[,row_annot_cols]
    shared_cols = intersect(colnames(row_annot),colnames(row_annot_from_data_file))
    if(length(shared_cols)==0){
      print("Warning:annotation in data file has no field overlap with the feature metadata file")
    }
    if(length(shared_cols)>0){
      check_meta=all(row_annot[,shared_cols]==row_annot_from_data_file[,shared_cols],na.rm=T)
      if(!check_meta){
        print(paste("Metadata file does not fit the dataset feature annotation columns:",feature_metadata_file))
      }
    }
    if(! "gene_symbol" %in% colnames(row_annot)){
      print("Warning: feature metadata file does not have the gene symbol field")
    }
    if(!all(colnames(row_annot_from_data_file) %in% colnames(row_annot))){
      print("Warning: dataset file has annotation columns that are not in the feature metadata file")
    }
  }
  else{
    row_annot = as.matrix(data[,row_annot_cols])
  }
  
  data = data[,-row_annot_cols]
  
  sample_meta = fread(metadata_file,data.table = F,stringsAsFactors = F,header = T)
  rownames(sample_meta) = sample_meta[,1]
  
  if(!all(is.element(colnames(data)[-c(1:2)],set=rownames(sample_meta)))){
    print("Error in ratio data, not all samples are in the metadata file")
    return(NULL)
  }
  
  print("Column names in the row annotation data:")
  print(colnames(row_annot))
  # set the rownames of the annotation and data
  if("cas_id" %in% colnames(row_annot)){
    print(paste("using cas_id as feature name"))
    rownames(row_annot) = row_annot[,"cas_id"]
  }
  else{
    if("ptm_id" %in% colnames(row_annot)){
      print(paste("using ptm_id as the feature name"))
      rownames(row_annot) = row_annot[,"ptm_id"]
    }
    else{
      print(paste("using protein id as feature name"))
      rownames(row_annot) = row_annot[,"protein_id"]
    }
  }
  
  rownames(data) = rownames(row_annot)
  
  return(list(sample_data=data,sample_meta=sample_meta,
              row_annot=row_annot,data_file=data_file,metadata_file=metadata_file))
}

parse_proteomics_datasets_from_download_obj<-function(obj,matrix_type="ratio",
         platforms=c("prot-ph","prot-pr","prot-ac"),...){
  
  # After splitting the local download we get arrays with the following info, by location:
  # n - data file names
  # n-1 - assay (e.g., prot, ac)
  # n-2 - tissue
  file_paths = strsplit(obj$downloaded_files,split=.Platform$file.sep)
  proteomics_parsed_datasets = list()
  failed_datasets = c()
  tissues_in_bucket = unique(sapply(file_paths,function(x)x[length(x)-2]))
  tissues_in_bucket = tissues_in_bucket[grepl("^t\\d\\d",tissues_in_bucket)]
  
  for(tissue in tissues_in_bucket){
    for(platform in platforms){
      curr_files = obj$downloaded_files[
        grepl(tissue,obj$downloaded_files) & grepl(platform,obj$downloaded_files)
      ]
      if(length(curr_files)==0){
        failed_datasets = rbind(failed_datasets,
             c(tissue,platform,"no files"))
        next
      }
      
      dataset_name = paste(tissue,platform,sep=",")
      print("--------------------------------------------------------------------")
      print(paste("             analyzing dataset:",dataset_name))
      try({
        proteomics_parsed_datasets[[dataset_name]] = 
          parse_proteomics_dataset_by_regex(curr_files,matrix_type)
      })
      if(is.null(proteomics_parsed_datasets[[dataset_name]])){
        print(paste("Error, could not parse the data from:", tissue,platform))
        print("Possible reason: missing data or metadata files, please revise")
        failed_datasets = rbind(failed_datasets,
                c(tissue,platform,"some files are missing"))
        next
      }
    }
  }
  
  return(list(proteomics_parsed_datasets=proteomics_parsed_datasets,
              failed_datasets=failed_datasets))
}






