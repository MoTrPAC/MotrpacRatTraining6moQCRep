#!/bin/R
# Nicole Gay
# custom functions for GET scripts 

#' Read a single file from Google Cloud into a data.table
#'
#' @param path GCP path, i.e. starts with "gs://"
#' @param sep column separator to use with "fread"
#' @param tmpdir scratch path to download files from GCP
#' @param GSUTIL_PATH path to "gsutil" on your computer
#' @param check_first check if file exists before downloading it. read in existing file if it exists. should be set to TRUE if you are running this function in parallel
#'
#' @return A data.table
dl_read_gcp <- function(path,sep='\t',tmpdir='/oak/stanford/groups/smontgom/nicolerg/tmp',GSUTIL_PATH='~/google-cloud-sdk/bin/gsutil',check_first=F){
  system(sprintf('mkdir -p %s',tmpdir))
  # download
  new_path = sprintf('%s/%s',tmpdir,basename(path))
  # only download if it doesn't exist to avoid conflicts when running this script in parallel; clear scratch space when you're done
  if(check_first){
    if(!file.exists(new_path)){
      cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
      system(cmd,ignore.stdout = T,ignore.stderr = T)
    }
  }else{
    cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
    system(cmd,ignore.stdout = T,ignore.stderr = T)
  }
  # read in the data as a data.table
  if(file.exists(new_path)){
    dt = fread(new_path,sep=sep,header=T)
    return(dt)
  }
  warning(sprintf("gsutil file %s does not exist.\n",path))
  return()
}


dl_load_gcp = function(path,object,tmpdir='/oak/stanford/groups/smontgom/nicolerg/tmp',GSUTIL_PATH='~/google-cloud-sdk/bin/gsutil'){
  system(sprintf('mkdir -p %s',tmpdir))
  system(sprintf("%s cp %s %s",GSUTIL_PATH, path, tmpdir))
  load(sprintf("%s/%s", tmpdir, basename(path)))
  system(sprintf("rm %s/%s", tmpdir, basename(path)))
  if(exists(object)){
    return(get(object))
  }else{
    error(sprintf("%s is not an object in %s", object, basename(path)))
  }
}


#' @param gost_object object returned from gprofiler2::gost()
#' @return updated gprofiler2::gost() object with raw p-values in the column "computed_p_value"
calc_gost_pvalue = function(gost_object){
  message('Reminder: Make sure you set "significant" to FALSE when you ran gost.')
  res_dt=data.table(g.res$result)
  setnames(res_dt, old='p_value', new='gost_adj_p_value')
  res_dt[,computed_p_value := phyper(intersection_size-1, 
                                     term_size, 
                                     effective_domain_size-term_size, 
                                     query_size,
                                     lower.tail= FALSE)]
  gost_object$result = as.data.frame(res_dt, stringsAsFactors=F)
  return(gost_object)
}

### Plotting functions ######################################################################

# plot individual time course 
# map = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t')
plot_read_timecourse = function(gene, dea, map, label=NULL){
  symbol = map[ensembl_gene==gene, gene_symbol][1]
  sub = dea[feature_ID == gene, .(sex, comparison_group, comparison_average_intensity, comparison_average_intensity_se,
                                  reference_average_intensity, reference_average_intensity_se)]
  dt = unique(sub[,list(comparison_average_intensity = reference_average_intensity,
                        comparison_average_intensity_se = reference_average_intensity_se,
                        comparison_group = '0w'),
                  by=sex])
  sub = rbindlist(list(sub, dt), fill=T)
  sub[,comparison_group := as.numeric(gsub('w','',comparison_group))]
  if(is.null(label)){
    t = sprintf('%s (%s)', symbol, gene)
  }else{
    t = sprintf('%s %s (%s)', label, symbol, gene)
  }
  g = ggplot(sub, aes(x=comparison_group, y=comparison_average_intensity, group=sex, colour=sex)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=comparison_average_intensity-comparison_average_intensity_se,
                      ymax=comparison_average_intensity+comparison_average_intensity_se), 
                  width=.2,
                  position=position_dodge(0.05)) +
    theme_classic() +
    labs(x='Weeks of training', y='Mean normalized counts', title=t) +
    scale_x_continuous(breaks=c(0,1,2,4,8)) +
    scale_colour_manual(values=sex_cols)
  return(g)
}


# adapted from qqman::qq
# pvector can be a single numeric vector of p-values or an named list, where each list item is a numeric vector of p-values
qq = function(pvector) {
  
  if(!is.list(pvector)){
    pvector_list = list()
    pvector_list[['observed']] = pvector
  }else{
    pvector_list <- pvector
  }
  
  reslist = list()
  for(i in 1:length(pvector_list)){
    ps = pvector_list[[i]]
    if (!is.numeric(ps)) stop("Input must be numeric.")
    # limit to not missing, not nan, not null, not infinite, between 0 and 1
    ps <- ps[!is.na(ps) & !is.nan(ps) & !is.null(ps) & is.finite(ps) & ps<1 & ps>0]
    o = -log10(sort(ps,decreasing=FALSE))
    e = -log10( ppoints(length(ps) ))
    
    df <- data.frame(expected = e,
                     observed = o,
                     name = names(pvector_list)[i])
    colnames(df)[colnames(df)=='observed'] = names(pvector_list)[i] 
    reslist[[i]] = df
  }
  res=rbindlist(reslist)
  
  return(res)
}


#' Plot time course from a single gene using normalized data 
timecourse <- function(gene, tissue, norm, meta, ensembl=F, gene_map=NULL){
  
  tnorm <- data.frame(t(norm))
  tnorm$viallabel <- rownames(tnorm)
  tnorm <- data.table(tnorm)
  
  # pick out a single gene
  if(is.null(gene_map)){
    gene_map = data.table(bitr(geneID = colnames(tnorm)[1:ncol(tnorm)-1], fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = org.Rn.eg.db, drop = TRUE))
    gene_map = gene_map[!duplicated(ENSEMBL)]
  }
  if(ensembl){
    gene_id = gene
    gene_symbol = gene_map[ENSEMBL == gene, SYMBOL][1]
    if(is.na(gene_symbol)){
      gene_symbol = gene_id
    }
  }else{
    gene_symbol = gene
    gene_id = unique(gene_map[SYMBOL == gene_symbol, ENSEMBL])
  }
  
  # check if gene is expressed
  if (length(gene_id) == 0){
    return()
  }
  if (!gene_id %in% colnames(tnorm)){
    return()
  }
  if(length(gene_id)>1){
    warning(sprintf("Multiple Ensembl IDs for %s: %s\nUsing %s.", gene, paste(gene_id, collapse=', '), gene_id[1]))
    gene_id = gene_id[1]
  }
  
  tnorm <- tnorm[,c(gene_id, 'viallabel'),with=F]
  
  Breaks=c('control','1w','2w','4w','8w')
  Limits=c('control','1w','2w','filler','4w',rep('filler',3),'8w')
  Labels=c('0','1','2','4','8')
  xlab='Weeks trained'
  
  # merge with meta
  meta_sub <- meta[,.(viallabel, sex, group)]
  meta_sub[,viallabel := as.character(viallabel)]
  tnorm[,viallabel := as.character(viallabel)]
  m <- merge(meta_sub, tnorm, by='viallabel')
  m <- m[,list(y = mean(get(gene_id)), se = sd(get(gene_id))/sqrt(length(get(gene_id))) ), by=c('group','sex')]
  
  g <- ggplot(m, aes(x=group, y=y, ymin=y-se, ymax=y+se, group=sex, colour=sex)) +
    geom_point() +
    geom_line() +
    geom_errorbar(width=0.3) +
    theme_classic() +
    scale_x_discrete(breaks=Breaks,limits=Limits,labels=Labels) +
    labs(x=xlab,y='Normalized expression',title=sprintf('%s: %s (%s)', tissue, gene_symbol, gene_id)) +
    scale_colour_manual(values=sex_cols,name='Sex') 
  
  return(g)
}

### Functions to preprocess GET data sets ######################################################################

# moved to cluster_viz_fx.R
# dl_format_pheno = function(scratch, gsutil_path, parallel=F){
#   dmaqc_metadata = 'gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/phenotype/pass1b_6m_viallabel_data.txt'
#   dmaqc_dict = 'gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/phenotype/merged_dictionary.txt'
#   
#   # download and format phenotypic data 
#   dmaqc_metadata = dl_read_gcp(dmaqc_metadata, tmpdir = scratch, GSUTIL_PATH = gsutil_path, check_first = parallel)
#   cols = dl_read_gcp(dmaqc_dict, tmpdir = scratch, GSUTIL_PATH = gsutil_path, check_first = parallel)
#   old_cols = colnames(dmaqc_metadata)
#   new_cols = tolower(cols[match(old_cols, BICUniqueID), FullName]) 
#   colnames(dmaqc_metadata) = new_cols # this isn't perfect, but we don't care about the columns it doesn't work for for now 
#   
#   # make some variables human-readable
#   # create new variables "protocol", "agegroup", "intervention", "sacrificetime", "sex" with readable strings 
#   for (var in c('key.protocol','key.agegroup','key.intervention','key.sacrificetime','registration.sex')){
#     d = cols[Field.Name == gsub('.*\\.','',var)]
#     keys=unname(unlist(strsplit(d[,Categorical.Values],'\\|')))
#     values=tolower(unname(unlist(strsplit(d[,Categorical.Definitions],'\\|'))))
#     names(values) = keys
#     # match keys to values; create new column 
#     new_var = gsub(".*\\.","",var)
#     dmaqc_metadata[,(new_var) := unname(values)[match(get(var), names(values))]]
#   }
#   dmaqc_metadata[,time_to_freeze := calculated.variables.frozetime_after_train - calculated.variables.deathtime_after_train]
#   
#   # clean up "sacrificetime"
#   dmaqc_metadata[,sacrificetime := sapply(sacrificetime, function(x) gsub(' week.*','w',x))]
#   
#   # clean up 'intervention'
#   dmaqc_metadata[grepl('training',intervention), intervention := 'training']
#   
#   # make "group" - "1w", "2w", "4w", "8w", "control"
#   dmaqc_metadata[,group := sacrificetime]
#   dmaqc_metadata[intervention == 'control', group := 'control']
#   
#   # # add "time to freeze"
#   # dmaqc_metadata[,time_to_freeze := calculated.variables.frozetime_after_train - calculated.variables.deathtime_after_train]
#   
#   # make tech ID a string
#   dmaqc_metadata[,specimen.processing.techid := paste0('tech',specimen.processing.techid)]
#   
#   # make viallabel char
#   dmaqc_metadata[,viallabel := as.character(viallabel)]
#   
#   return(dmaqc_metadata)
# }


#' Get and fix peak annotations from ChIPseeker
#' 
#' @param counts_dt data.table or data.frame with columns c('chrom','start','end','feature_ID') (or for ATAC only, data.frame with rownames as feature_IDs, format chrom:start-end).
#' 
#' @return data.table with formatted ChIPseeker output and additional columns "custom_annotation" and "relationship_to_gene". 
#'     "relationship_to_gene" is the shortest distance between the feature and the start or end of the closest gene. 
#'     It is 0 if the feature has any overlap with the gene. "custom_annotation" fixes many issues with the ChIPseeker annotation (v1.22.1)
get_peak_annotations = function(counts_dt, species="Rattus norvegicus", release=96, txdb=NULL){
  require("data.table")
  require("ChIPseeker")
  require("GenomicFeatures")
  require("RMariaDB")

  if(!"feature_ID" %in% colnames(counts_dt) & !is.data.table(counts_dt)){
    genomic_peaks = data.table(
      feature_ID = rownames(counts_dt),
      chrom = gsub(":.*","",rownames(counts_dt)),
      start = as.numeric(gsub(".*:|-.*","",rownames(counts_dt))),
      end = as.numeric(gsub(".*-","",rownames(counts_dt)))
    )
  }else if(!"feature_ID" %in% colnames(counts_dt) & is.data.table(counts_dt)){
    counts = counts_dt
    counts[,feature_ID := paste0(chrom,':',start,'-',end)]
    genomic_peaks = counts[,.(chrom, start, end, feature_ID)]
  }else if("feature_ID" %in% colnames(counts_dt) & is.data.table(counts_dt)){
    counts = counts_dt
    genomic_peaks = counts[,.(chrom, start, end, feature_ID)]
  }else if("feature_ID" %in% colnames(counts_dt) & is.data.table(counts_dt)){
    counts = as.data.table(counts_dt)
    genomic_peaks = counts[,.(chrom, start, end, feature_ID)]
  }else{
    stop("Incorrect input format. 'counts_dt' should be a data.frame or data.table with columns 'chrom', 'start', 'end', and 'feature_ID'.")
  }
  
  if(is.null(txdb)){
    # make TxDb object 
    txdb = makeTxDbFromEnsembl(organism=species,
                               release=release)
  }

  accepted_chrom = seqlevels(txdb)
  # remove contigs
  accepted_chrom = accepted_chrom[!grepl("\\.",accepted_chrom)]
  
  # remove contigs from input 
  genomic_peaks = genomic_peaks[!grepl("\\.",chrom)]
  
  genomic_peaks[,chrom := gsub("^chr","",as.character(chrom))]
  if(!all(unique(genomic_peaks[,chrom]) %in% accepted_chrom)){
    stop(sprintf("The following chromosomes are found in the input but not in the txdb object: %s", paste0(unique(!genomic_peaks[,chrom] %in% accepted_chrom), collapse=', ')))
  }
  
  # annotate peaks 
  peak = GRanges(seqnames = genomic_peaks[,chrom], 
                 ranges = IRanges(as.numeric(genomic_peaks[,start]), as.numeric(genomic_peaks[,end])))
  peakAnno = annotatePeak(peak, 
                          level = "gene",
                          tssRegion=c(-2000,1000), 
                          TxDb=txdb,
                          overlap = "all")
  pa = as.data.table(peakAnno@anno)
  
  # add back feature_ID
  if(nrow(pa)==nrow(genomic_peaks)){
    pa[,feature_ID := genomic_peaks[,feature_ID]]
  }else{
    cols=c('seqnames','start','end')
    pa[,(cols) := lapply(.SD, as.character), .SDcols=cols]
    cols=c('chrom','start','end')
    genomic_peaks[,(cols) := lapply(.SD, as.character), .SDcols=cols]
    pa = merge(pa, genomic_peaks, by.x=c('seqnames','start','end'), by.y=c('chrom','start','end'), all.y=T)
  }
  
  # add simpler annotation
  pa[,short_annotation := annotation]
  pa[grepl('Exon', short_annotation), short_annotation := 'Exon']
  pa[grepl('Intron', short_annotation), short_annotation := 'Intron']
  
  pa[,c('geneChr','strand') := NULL]
  
  # add a "relationship_to_gene" (from gene start or end, whichever is closer)
  cols=c('start','end','geneStart','geneEnd','geneStrand')
  pa[,(cols) := lapply(.SD, as.numeric), .SDcols=cols]
  pa[,dist_upstream := ifelse(end-geneStart <=0, end-geneStart, NA_real_)]
  pa[,dist_downstream := ifelse(start-geneEnd >= 0, start-geneEnd, NA_real_)]
  # if feature overlaps with gene body, say dist_to_gene is 0
  pa[end >= geneStart & start <= geneEnd, dist_downstream := 0]
  pa[end >= geneStart & start <= geneEnd, dist_upstream := 0]
  pa[, relationship_to_gene := ifelse(is.na(dist_downstream), dist_upstream, dist_downstream)]
  pa[,c('dist_upstream','dist_downstream') := NULL]
  
  # fix "Downstream" and "Distal Intergenic" annotations
  pa[relationship_to_gene == 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Overlaps Gene"]
  # Downstream
  pa[geneStrand == 1 & relationship_to_gene > 0 & relationship_to_gene < 5000, short_annotation := "Downstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 0 & relationship_to_gene > -5000, short_annotation := "Downstream (<5kb)"]
  # Upstream (promoter excluded)
  pa[geneStrand == 1 & relationship_to_gene > -5000 & relationship_to_gene < 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 5000 & relationship_to_gene > 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[abs(relationship_to_gene) >= 5000, short_annotation := "Distal Intergenic"]
  
  # for(anno in unique(pa[,short_annotation])){
  #   print(anno)
  #   hist(pa[short_annotation==anno & geneStrand ==1, relationship_to_gene], xlim=c(-8000,8000), breaks=10000, main=sprintf("%s, positive strand", anno))
  #   hist(pa[short_annotation==anno & geneStrand ==2, relationship_to_gene], xlim=c(-8000,8000), breaks=10000, main=sprintf("%s, negative strand", anno))
  # }
  setnames(pa, c('short_annotation','annotation','seqnames','geneId'), c('custom_annotation','chipseeker_annotation','chrom','ensembl_gene'))
  
  return(pa)
}


preprocess_pass1b_rnaseq_gcp = function(TISSUE_CODE, SEX, gsutil_path='~/google-cloud-sdk/bin/gsutil', scratch='/oak/stanford/groups/smontgom/nicolerg/tmp', outliers=NULL, parallel=F){

  if(!SEX %in% c('male','female','all')){
    stop('"SEX" must take on one of the following values:\n  "male","female","all"\n')
  }
  
  # fix some things
  if(TISSUE_CODE == 't54-hypothalmus'){
    TISSUE_CODE = 't54-hypothalamus'
  }
  
  dmaqc_metadata = dl_format_pheno(scratch, gsutil_path, parallel=parallel)
  
  # load counts 
  counts = dl_read_gcp(sprintf('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/transcriptomics/%s/transcript-rna-seq/motrpac_pass1b-06_%s_transcript-rna-seq_rsem-genes-count.txt', TISSUE_CODE, TISSUE_CODE), 
                       tmpdir = scratch, 
                       GSUTIL_PATH = gsutil_path,
                       check_first = parallel)
  if(is.null(counts)){return()}
  counts = as.data.frame(counts)
  rownames(counts) = counts$gene_id
  counts$gene_id = NULL
  raw_counts = counts
  
  # filter by genes in normalized data 
  tmm = dl_read_gcp(sprintf('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/transcriptomics/transcript-rna-seq/normalized-data/motrpac_pass1b-06_%s_transcript-rna-seq_normalized-log-cpm.txt', TISSUE_CODE), 
                    tmpdir = scratch, 
                    GSUTIL_PATH = gsutil_path,
                    check_first = parallel)
  # remove pid and bid rows
  tmm = tmm[3:nrow(tmm)]
  counts = counts[tmm[,viallabel],]
  
  # convert tmm to data.frame
  tmm = as.data.frame(tmm)
  rownames(tmm) = tmm$viallabel
  tmm$viallabel = NULL
  
  # subset metadata 
  meta_data = dmaqc_metadata[as.character(viallabel) %in% colnames(tmm)]
  # remove reference standards from counts 
  counts = counts[,as.character(meta_data[,viallabel])]
  
  # coerce counts to integer (RSEM has fractional values for counts sometimes)
  counts_round = as.data.frame(apply(counts, c(1,2), as.integer)) 
  raw_counts_round = as.data.frame(apply(raw_counts, c(1,2), as.integer)) 
  
  # RNA-seq metadata 
  qa_qc = dl_read_gcp('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv', 
                      sep=',', 
                      tmpdir = scratch, 
                      GSUTIL_PATH = gsutil_path,
                      check_first = parallel)
  # adjust column names
  colnames(qa_qc) = tolower(gsub(' .*','',colnames(qa_qc)))
  # remove duplicate columns
  qa_qc[,c('pid','bid') := NULL]
  meta_data[,viallabel := as.character(viallabel)]
  qa_qc[,vial_label := as.character(vial_label)]
  meta = merge(meta_data, qa_qc, by.x = 'viallabel', by.y = 'vial_label')
  
  if(SEX != 'all'){
    # subset to sex 
    meta = meta[sex == SEX]
    tmm = tmm[as.character(meta[,viallabel])]
    counts_round = counts_round[as.character(meta[,viallabel])]
  }
  
  # remove outliers 
  if(!is.null(outliers)){
    meta = meta[,viallabel := as.character(viallabel)]
    meta = meta[!viallabel %in% outliers]
    tmm = tmm[meta[,viallabel]]
    counts_round = counts_round[meta[,viallabel]]
  }
  
  return(list(meta = meta,
              raw_counts = raw_counts_round, 
              counts = counts_round,
              norm = tmm))
}


preprocess_pass1b_atacseq_gcp = function(TISSUE_CODE, # or dataset, if it starts with tissue_code
                                         outliers = NULL, 
                                         gsutil_path='~/google-cloud-sdk/bin/gsutil', 
                                         scratch='/oak/stanford/groups/smontgom/nicolerg/tmp', 
                                         parallel=T){
  
  # ATAC-seq meta 
  wet = dl_read_gcp('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics.csv', sep=',', check_first = T)
  # DMAQC metadata 
  dmaqc_meta = dl_format_pheno(scratch, gsutil_path, parallel=parallel)
  wet[,viallabel := as.character(viallabel)]
  meta = merge(dmaqc_meta, wet, by='viallabel', all.y=T)
  
  colnames(meta) = tolower(colnames(meta))
  
  tissue_number = toupper(gsub("-.*","",TISSUE_CODE))
  tissue = gsub(",.*","",TISSUE_CODE)
  # load raw counts 
  counts = dl_read_gcp(sprintf("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/epigenomics/%s/epigen-atac-seq/motrpac_pass1b-06_%s_epigen-atac-seq_counts.txt.gz",TISSUE_CODE,TISSUE_CODE), check_first = T)
  counts = data.frame(counts, check.names=F)
  rownames(counts) = paste0(counts$chrom, ':', counts$start, '-', counts$end)
  counts[,c('chrom','start','end')] = NULL
  
  # remove outliers
  if(!is.null(outliers)){
    counts[,as.character(outliers)] = NULL
  }
  
  # separate data by site
  data_list = list()
  curr_meta = meta[viallabel %in% colnames(counts) & !grepl("^8", viallabel)]
  for(site in unique(curr_meta[,get_site])){
    
    site_meta = curr_meta[get_site==site]
    site_counts = counts[,site_meta[,viallabel]]
    
    # remove non-auto peaks
    site_counts = site_counts[grepl("^chr[0-9]|^chrY|^chrX", rownames(site_counts)),]
    
    # exclude low count peaks in the current dataset
    # at least 10 counts in N samples
    n_samples = 4
    filt_counts = site_counts[rowSums(data.frame(lapply(site_counts, function(x) as.numeric(x >= 10)), check.names=F)) >= n_samples,]
    
    # quantile normalize
    # this takes a couple of minutes given the size of the peak x sample counts matrix
    norm_counts = voom(filt_counts,normalize.method = "quantile")
    sub_norm = round(norm_counts$E,2)
    
    label = sprintf("%s,%s", tissue, tolower(site))
    data_list[[label]] = list(meta = site_meta, 
                              #voom_obj = norm_counts, # this might break some old code, but voom weights need design
                              norm = sub_norm,
                              unfilt_counts = site_counts,
                              filt_counts = filt_counts)
  }
  return(data_list)
}


### Functions for transcript-rna-seq DEA ######################################################################


transcript_prep_data = function(tissue_code, 
                                meta, 
                                counts, 
                                tmm,
                                covariates=c('pct_globin', 'rin', 'pct_umi_dup', 'median_5_3_bias'), 
                                outliers=NULL, 
                                gsutil_path='~/google-cloud-sdk/bin/gsutil',
                                parallel=F){
  
  if(tissue_code %in% c("t31-plasma", "t57-tibia")){return(NULL)}
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  # load outliers
  if(is.null(outliers)){
    rna_outliers = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/dea/pass1b-06_transcript-rna-seq_removed-outliers_20201028.txt', sep='\t', check_first=parallel)
    outliers = as.character(rna_outliers[,viallabel])
  }
  
  cat(tissue_code, sep = '\n')
  
  # remove outliers 
  curr_outliers = outliers[outliers %in% as.character(meta[,viallabel])]
  if(length(curr_outliers)>0){
    meta = meta[,viallabel := as.character(viallabel)]
    meta = meta[!viallabel %in% outliers]
    counts = counts[meta[,viallabel]]
    tmm = tmm[meta[,viallabel]]
  }
  if(tissue_code == 't65-aorta'){
    # add Ucp1 as a covariate
    ucp1 = data.table(viallabel = colnames(counts), ucp1 = unname(unlist(counts['ENSRNOG00000003580',])))
    meta = merge(meta, ucp1, by = 'viallabel')
    covariates = c(covariates, 'ucp1')
  }
  
  # impute missing values
  new = fix_missing(covariates, meta)
  covariates = new$covariates
  meta = new$meta
  
  # center and scale continuous covariates
  for (cov in covariates){
    if(is.numeric(meta[,get(cov)])){
      meta[,(cov) := scale(meta[,get(cov)], center = T, scale = T)]
    }
  }
  
  meta[,sex_group := paste0(sex, ';', group)]
  
  return(list(fixed_meta = meta, 
              fixed_covariates = covariates, 
              fixed_counts = counts,
              fixed_tmm = tmm, 
              curr_outliers = curr_outliers))
}


transcript_timewise_dea_each_sex = function(tissue_code, meta, 
  counts, covariates, curr_outliers, date, training_groups = c('1w','2w','4w','8w'),
  save_rdata=T, write=T){
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_timewise-dea_%s.txt',
                    tissue_code,date)
  if(file.exists(outfile)){
    dt = fread(outfile, sep='\t', header=T)
    return(dt)
  }

  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
    # subset counts and meta
    curr_samples = meta[sex == SEX, viallabel]
    curr_meta = meta[sex == SEX]
    curr_counts = counts[,curr_samples]
    
    contrasts = list()
    i = 1
    for (tp in training_groups){
      contrasts[[i]] = c('group', tp, 'control')
      i = i+1
    }
    
    # shrunk results 
    # function in pi1_cook_fx.R
    deseq_res_shrunk = run_deseq(curr_counts, # filtered counts
                          curr_meta, # metadata
                          covariates, # covariates
                          'group', # outcome of interest
                          contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                          shrink = T)
    
    # non-shrunk results 
    deseq_res = run_deseq(curr_counts, # filtered counts
                          curr_meta, # metadata
                          covariates, # covariates
                          'group', # outcome of interest
                          contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                          shrink = F)
    
    if(save_rdata){
      save(deseq_res, deseq_res_shrunk, file=sprintf('rdata/%s_%s_timewise-dea_%s.RData', tissue_code, SEX, date))
    }
    
    # collect res
    res_shrunk = data.table(deseq_res_shrunk$res)
    res_nonshrunk = data.table(deseq_res$res)
    setnames(res_shrunk, c("log2FoldChange", "lfcSE"), c("shrunk_logFC","shrunk_logFC_se"))
    setnames(res_nonshrunk, c("log2FoldChange", "lfcSE", "stat"), c("logFC","logFC_se", "zscore"))
    res_shrunk = res_shrunk[,.(gene_id, shrunk_logFC, shrunk_logFC_se, numerator, denominator)]
    res = merge(res_nonshrunk, res_shrunk, by=c("gene_id","numerator","denominator"))
    res[,sex := SEX]
    
    setnames(res, c("numerator","pvalue","gene_id"), c("comparison_group","p_value","feature_ID"))
    res[,denominator := NULL]
    
    # add some columns
    res[,tissue := tissue_code]
    res[,assay := 'transcript-rna-seq']
    res[,removed_samples := paste0(curr_outliers, collapse=',')]
    # res[,covariates := paste0(covariates, collapse=',')] added within run_deseq()
    
    # add average intensities 
    norm_counts = as.data.frame(counts(deseq_res$dds, normalized=T))
    ref_sub = norm_counts[,as.character(curr_meta[group == 'control', viallabel])]
    ref_means = rowMeans(ref_sub, na.rm=T)
    ref_se = apply(ref_sub, 1, function(x) sd(x)/sqrt(sum(!is.na(x))) )
    mlist = list()
    i = 1
    for(tp in unique(res[,comparison_group])){
      # get average values
      counts_sub = norm_counts[,as.character(curr_meta[group == tp, viallabel])]
      counts_means = data.table(sex=SEX,
                                comparison_group=tp,
                                comparison_average_intensity=rowMeans(counts_sub, na.rm=T),
                                comparison_average_intensity_se=apply(counts_sub, 1, 
                                                                      function(x) sd(x)/sqrt(sum(!is.na(x))) ),
                                reference_average_intensity=ref_means,
                                reference_average_intensity_se=ref_se,
                                feature_ID=rownames(counts_sub))
      mlist[[i]] = counts_means
      i = i+1
    }
    cmeans = rbindlist(mlist)
    dt = merge(res, cmeans, by=c('feature_ID', 'sex', 'comparison_group'))
    
    sex_res[[SEX]] = dt
  }
  
  dt = rbindlist(sex_res)
  
  dt = dt[,.(
    feature_ID,
    sex,
    comparison_group,
    assay,
    tissue,
    covariates,
    removed_samples,
    logFC,
    logFC_se,
    shrunk_logFC,
    shrunk_logFC_se,
    zscore,
    p_value,
    comparison_average_intensity,
    comparison_average_intensity_se,
    reference_average_intensity,
    reference_average_intensity_se
  )]
  
  # if aorta, remove 1w, 2w F
  if(tissue_code=="t65-aorta"){
    dt = dt[!(sex=='female' & comparison_group %in% c('1w','2w'))]
  }
  
  if(write){
    write.table(dt, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
  }
  return(dt)
  
}


transcript_timewise_dea = function(tissue_code, meta, counts, covariates, curr_outliers, date, save_rdata=T, write=T){

  # fix some inconsistencies
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }

  outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_timewise-dea_%s.txt',
                    tissue_code,date)
  if(file.exists(outfile)){
    dt = fread(outfile, sep='\t', header=T)
    return(dt)
  }

  # run DESeq with all samples combined
  contrasts = list()
  i = 1
  for (tp in c('1w','2w','4w','8w')){
    for (s in unique(meta[,sex])){
      contrasts[[i]] = c('sex_group', sprintf('%s;%s', s, tp), sprintf('%s;control', s))
      i = i+1
    }
  }
  # function in pi1_cook_fx.R
  deseq_res = run_deseq(counts, # filtered counts
                        meta, # metadata
                        covariates, # covariates
                        'sex_group', # outcome of interest
                        contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                        shrink = T)
  if(save_rdata){
    save(deseq_res, file=sprintf('rdata/%s_timewise-dea_%s.RData', tissue_code, date))
  }

  dt = copy(deseq_res$res)

  # rename some columns
  setnames(dt,
           old=c('gene_id','log2FoldChange','pvalue','lfcSE'),
           new=c('feature_ID','logFC','p_value','logFC_se'))

  # add some columns
  dt[,sex := sapply(numerator, function(x) unlist(unname(strsplit(x, ';')))[1])]
  dt[,comparison_group := sapply(numerator, function(x) unlist(unname(strsplit(x, ';')))[2])]
  dt[,c('numerator','denominator') := NULL]
  dt[,tissue := tissue_code]
  dt[,assay := 'transcript-rna-seq']
  dt[,removed_samples := paste0(curr_outliers, collapse=',')]
  dt[,covariates := paste0(covariates, collapse=',')]

  norm_counts = as.data.frame(counts(deseq_res$dds, normalized=T))
  # add average intensities
  mlist = list()
  i = 1
  for(.sex in unique(dt[,sex])){
    for(tp in unique(dt[,comparison_group])){
      # get average values
      counts_sub = norm_counts[,as.character(meta[group == tp & sex == .sex, viallabel])]
      ref_sub = norm_counts[,as.character(meta[group == 'control' & sex == .sex, viallabel])]
      counts_means = data.table(sex=.sex,
                                comparison_group=tp,
                                comparison_average_intensity=rowMeans(counts_sub, na.rm=T),
                                comparison_average_intensity_se=apply(counts_sub, 1,
                                                                      function(x) sd(x)/sqrt(sum(!is.na(x))) ),
                                reference_average_intensity=rowMeans(ref_sub, na.rm=T),
                                reference_average_intensity_se=apply(ref_sub, 1,
                                                                     function(x) sd(x)/sqrt(sum(!is.na(x))) ),
                                feature_ID=rownames(counts_sub))
      mlist[[i]] = counts_means
      i = i+1
    }
  }
  cmeans = rbindlist(mlist)
  dt = merge(dt, cmeans, by=c('feature_ID', 'sex', 'comparison_group'))

  dt = dt[,.(assay, feature_ID, sex, tissue,
             comparison_group, covariates, removed_samples,
             logFC, logFC_se, p_value,
             comparison_average_intensity, comparison_average_intensity_se,
             reference_average_intensity, reference_average_intensity_se)]

  # get logFC and logFC_se without shrinkage
  dds = deseq_res$dds
  res_list = list()
  for (c in contrasts){
    res = results(dds, contrast = c)
    res_dt = data.table(gene_id = rownames(counts(dds)),
                        log2FoldChange = res$log2FoldChange,
                        lfcSE = res$lfcSE,
                        zscore = res$stat,
                        pvalue = res$pvalue,
                        numerator = c[2],
                        denominator = c[3])
    res_list[[paste0(c, collapse=' ')]] = res_dt
  }
  all_res = rbindlist(res_list)
  # add some columns
  all_res[,sex := sapply(numerator, function(x) unlist(unname(strsplit(x, ';')))[1])]
  all_res[,comparison_group := sapply(numerator, function(x) unlist(unname(strsplit(x, ';')))[2])]
  all_res[,c('numerator','denominator') := NULL]

  # merge
  setnames(dt, old=c("logFC","logFC_se"), new=c("shrunk_logFC","shrunk_logFC_se"))
  merged = merge(dt, all_res[,.(zscore, gene_id, log2FoldChange, lfcSE, sex, comparison_group)],
                 by.x=c('feature_ID','sex','comparison_group'),
                 by.y=c('gene_id','sex','comparison_group'))
  setnames(merged, old=c("log2FoldChange", "lfcSE"), new=c("logFC","logFC_se"))
  merged = merged[,.(
    feature_ID,
    sex,
    comparison_group,
    assay,
    tissue,
    covariates,
    removed_samples,
    logFC,
    logFC_se,
    shrunk_logFC,
    shrunk_logFC_se,
    zscore,
    p_value,
    comparison_average_intensity,
    comparison_average_intensity_se,
    reference_average_intensity,
    reference_average_intensity_se
  )]

  # if aorta, remove 1w, 2w F
  if(tissue_code=="t65-aorta"){
    merged = merged[!(sex=='female' & comparison_group %in% c('1w','2w'))]
  }

  if(write){
    write.table(merged, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
  }
  return(merged)

}


transcript_training_dea_each_sex = function(tissue_code, meta, counts, covariates, curr_outliers, date, write=T){
  
  require("metap")
  require("DESeq2")
  require("data.table")
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  if(write){
    outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_training-dea_%s.txt',
                      tissue_code,date)
    
    if(file.exists(outfile)){
      dt = fread(outfile, sep='\t', header=T)
      return(dt)
    }
  }
  
  # add vena cava outliers
  if(tissue_code=="t65-aorta"){
    curr_outliers = unique(c(curr_outliers, as.character(meta[sex=='female' & group %in% c('1w','2w'), viallabel])))
    meta = meta[,viallabel := as.character(viallabel)]
    meta = meta[!viallabel %in% curr_outliers]
    counts = counts[meta[,viallabel]]
  }
  
  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
    # subset counts and meta
    curr_samples = meta[sex == SEX, viallabel]
    curr_meta = meta[sex == SEX]
    curr_counts = counts[,curr_samples]
    
    # center and scale continuous variables
    curr_cov = covariates
    for (cov in curr_cov){
      # remove if constant
      if(length(unique(curr_meta[,get(cov)])) == 1){
        message(sprintf("Covariate %s is constant for %s. Removing.", cov, SEX))
        curr_cov = curr_cov[curr_cov != cov]
      }else{
        # center and scale
        if(is.numeric(curr_meta[,get(cov)])){
          curr_meta[,(cov) := scale(curr_meta[,get(cov)], center = T, scale = T)]
        }
      }
    }
    
    full = paste0('~', paste0(c(curr_cov, 'group'), collapse=' + '))
    reduced = paste0('~', paste0(curr_cov, collapse=' + '))
    # # custom contrast for vena cava
    # if(tissue_code=="t65-aorta"){
    #   meta[,group := factor(group, levels=c('control','1w','2w','4w','8w'))]
    #   coldata = data.frame(model.matrix(eval(parse(text=contrast)), data=meta))
    #   coldata$`group1w.sexmale` = NULL
    #   coldata$`group2w.sexmale` = NULL
    #   coldata$X.Intercept. = NULL
    #   contrast = paste0('~', paste0(colnames(coldata), collapse=' + '))
    #   meta = coldata
    #   cols = colnames(meta)[grepl('group|sex', colnames(meta))]
    #   meta[cols] = sapply(meta[cols],as.factor)
    #   reduced = paste0('~', paste0(colnames(meta)[!grepl('group', colnames(meta))], collapse=' + '))
    # }
    
    dds = DESeqDataSetFromMatrix(countData = curr_counts,
                                 colData = curr_meta,
                                 design = eval(parse(text=full)))
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds)
    dds = nbinomLRT(dds, reduced=eval(parse(text=reduced)), maxit=500)
    
    res = results(dds)
    res_dt = data.table(feature_ID = rownames(res), 
                        lrt = res$stat,
                        p_value = res$pvalue,
                        tissue = tissue_code,
                        assay = 'transcript-rna-seq',
                        removed_samples = paste0(curr_outliers, collapse=','),
                        full_model=gsub(' ','',full),
                        reduced_model=gsub(' ','',reduced))
    res_dt = res_dt[,.(
      feature_ID,
      assay,
      tissue,
      removed_samples,
      lrt,
      p_value,
      full_model,
      reduced_model
    )]
    
    sex_res[[SEX]] = res_dt
  }
  
  if(length(sex_res) > 1){
    male = sex_res[['male']]
    female = sex_res[['female']]
    merged = merge(male, female, by=c('feature_ID','assay','tissue','removed_samples'),
                   suffixes=c("_male","_female"), all=T)
    missing = merged[is.na(p_value_male) | is.na(p_value_female)]
    complete = merged[!is.na(p_value_male) & !is.na(p_value_female)]
    complete[, p_value := sumlog(c(p_value_male, p_value_female))$p, by=seq(1, nrow(complete))]
    missing[, p_value := ifelse(is.na(p_value_male), p_value_female, p_value_male)]
    res_dt = rbindlist(list(complete, missing))
  }else{
    res_dt = sex_res[[1]]
  }
  
  if(write){
    write.table(res_dt, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
  }
  return(res_dt)
}


# transcript_training_dea = function(tissue_code, meta, counts, covariates, curr_outliers, date, write=T){
#   
#   # fix some inconsistencies 
#   if(tissue_code == 't54-hypothalmus'){
#     tissue_code = 't54-hypothalamus'
#   }
#   
#   if(write){
#     outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_training-dea_%s.txt',
#                       tissue_code,date)
#     
#     if(file.exists(outfile)){
#       dt = fread(outfile, sep='\t', header=T)
#       return(dt)
#     }
#   }
# 
#   # add vena cava outliers
#   if(tissue_code=="t65-aorta"){
#     curr_outliers = unique(c(curr_outliers, as.character(meta[sex=='female' & group %in% c('1w','2w'), viallabel])))
#     meta = meta[,viallabel := as.character(viallabel)]
#     meta = meta[!viallabel %in% curr_outliers]
#     counts = counts[meta[,viallabel]]
#   }
#   
#   if(tissue_code %in% c("t63-testes","t64-ovaries")){
#     contrast = paste0('~', paste0(c(covariates, 'group'), collapse=' + '))
#     reduced = paste0('~', paste0(covariates, collapse=' + '))
#   }else{
#     contrast = paste0('~', paste0(c(covariates, 'group + sex + group:sex'), collapse=' + '))
#     reduced = paste0('~', paste0(c(covariates, 'sex'), collapse=' + '))
#   }
#   
#   # custom contrast for vena cava
#   if(tissue_code=="t65-aorta"){
#     meta[,group := factor(group, levels=c('control','1w','2w','4w','8w'))]
#     coldata = data.frame(model.matrix(eval(parse(text=contrast)), data=meta))
#     coldata$`group1w.sexmale` = NULL
#     coldata$`group2w.sexmale` = NULL
#     coldata$X.Intercept. = NULL
#     contrast = paste0('~', paste0(colnames(coldata), collapse=' + '))
#     meta = coldata
#     cols = colnames(meta)[grepl('group|sex', colnames(meta))]
#     meta[cols] = sapply(meta[cols],as.factor)
#     reduced = paste0('~', paste0(colnames(meta)[!grepl('group', colnames(meta))], collapse=' + '))
#   }
#   
#   dds = DESeqDataSetFromMatrix(countData = counts,
#                                colData = meta,
#                                design = eval(parse(text=contrast)))
#   dds = estimateSizeFactors(dds)
#   dds = estimateDispersions(dds)
#   dds = nbinomLRT(dds, reduced=eval(parse(text=reduced)), maxit=500)
#   
#   if(write){
#     save(dds, file=sprintf('rdata/%s_training-dea_%s.RData', tissue_code, date))
#   }
#   
#   res = results(dds)
#   res_dt = data.table(feature_ID = rownames(res), 
#                       lrt = res$stat,
#                       p_value = res$pvalue,
#                       tissue = tissue_code,
#                       assay = 'transcript-rna-seq',
#                       removed_samples = paste0(curr_outliers, collapse=','),
#                       full_model=gsub(' ','',contrast),
#                       reduced_model=gsub(' ','',reduced))
#   res_dt = res_dt[,.(
#     feature_ID,
#     assay,
#     tissue,
#     removed_samples,
#     lrt,
#     p_value,
#     full_model,
#     reduced_model
#   )]
#   
#   if(write){
#     write.table(res_dt, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
#   }
#   return(res_dt)
# }


# transcript_training_sex_biased = function(tissue_code, meta, counts, covariates, curr_outliers, date, write=T){
#   
#   # fix some inconsistencies 
#   if(tissue_code == 't54-hypothalmus'){
#     tissue_code = 't54-hypothalamus'
#   }
#   
#   outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_training-sex-biased_%s.txt',
#                     tissue_code,date)
#   
#   if(file.exists(outfile)){
#     dt = fread(outfile, sep='\t', header=T)
#     return(dt)
#   }
#   
#   # add vena cava outliers
#   if(tissue_code=="t65-aorta"){
#     curr_outliers = unique(c(curr_outliers, as.character(meta[sex=='female' & group %in% c('1w','2w'), viallabel])))
#     meta = meta[,viallabel := as.character(viallabel)]
#     meta = meta[!viallabel %in% curr_outliers]
#     counts = counts[meta[,viallabel]]
#   }
#   
#   if(tissue_code %in% c("t63-testes","t64-ovaries")){return(NULL)}
#   
#   contrast = paste0('~', paste0(c(covariates, 'group + sex + group:sex'), collapse=' + '))
#   reduced = paste0('~', paste0(c(covariates, 'sex', 'group'), collapse=' + '))
#   
#   # custom contrast for vena cava
#   if(tissue_code=="t65-aorta"){
#     meta[,group := factor(group, levels=c('control','1w','2w','4w','8w'))]
#     coldata = data.frame(model.matrix(eval(parse(text=contrast)), data=meta))
#     coldata$`group1w.sexmale` = NULL
#     coldata$`group2w.sexmale` = NULL
#     coldata$X.Intercept. = NULL
#     contrast = paste0('~', paste0(colnames(coldata), collapse=' + '))
#     meta = coldata
#     cols = colnames(meta)[grepl('group|sex', colnames(meta))]
#     meta[cols] = sapply(meta[cols],as.factor)
#     reduced = paste0('~', paste0(colnames(meta)[!grepl('\\.', colnames(meta))], collapse=' + '))
#   }
#   
#   dds = DESeqDataSetFromMatrix(countData = counts,
#                                colData = meta,
#                                design = eval(parse(text=contrast)))
#   dds = estimateSizeFactors(dds)
#   dds = estimateDispersions(dds)
#   dds = nbinomLRT(dds, reduced=eval(parse(text=reduced)), maxit=500)
#   
#   save(dds, file=sprintf('rdata/%s_training-sex-biased_%s.RData', tissue_code, date))
#   
#   res = results(dds)
#   res_dt = data.table(feature_ID = rownames(res), 
#                       lrt = res$stat,
#                       p_value = res$pvalue,
#                       tissue = tissue_code,
#                       assay = 'transcript-rna-seq',
#                       removed_samples = paste0(curr_outliers, collapse=','),
#                       full_model=gsub(' ','',contrast),
#                       reduced_model=gsub(' ','',reduced))
#   res_dt = res_dt[,.(
#     feature_ID,
#     assay,
#     tissue,
#     removed_samples,
#     lrt,
#     p_value,
#     full_model,
#     reduced_model
#   )]
#   
#   if(write){
#     write.table(res_dt, outfile, sep='\t', col.names=T, row.names=F, quote=F)
#   }
#   
#   return(res_dt)
# }

### Functions for epigen-atac-seq DEA ######################################################################

# atac_training_dea = function(filt_counts, meta, covariates, label=NULL, return_resid=F){
#   
#   full = paste0(c("~1", "group*sex", covariates), collapse=" + ")
#   reduced = paste0(c("~1", "sex", covariates), collapse=" + ")
#   # full = paste0("~ 1 + group*sex + ", paste0(covariates, collapse=" + "))
#   # reduced = paste0("~ 1 + sex + ", paste0(covariates, collapse=" + "))
#   
#   meta_df = as.data.frame(meta, check.names=F)
#   meta_df$Sample_batch = factor(meta_df$Sample_batch)
#   meta_df$group = factor(meta_df$group, levels=c('control','1w','2w','4w','8w')) # IMPORTANT
#   
#   # normalize and get voom weights 
#   design = model.matrix(eval(parse(text=full)), data = meta_df)
#   # check if full rank
#   if(!is.fullrank(design)){
#     if("Sample_batch" %in% covariates){
#       curr_cov = covariates[!covariates == 'Sample_batch']
#       full = paste0("~ 1 + group*sex + ", paste0(curr_cov, collapse=" + "))
#       reduced = paste0("~ 1 + sex + ", paste0(curr_cov, collapse=" + "))
#       design = model.matrix(eval(parse(text=full)), data = meta_df)
#       warning(sprintf("Sample_batch and group or sex are confounded for %s", label))
#     }else{
#       stop(sprintf("Model matrix with design %s is not full rank.", full))
#     }
#   }
#   
#   curr_voom = voom(filt_counts, design=design, normalize.method="quantile")
#   
#   limma_model1 = lmFit(curr_voom, design)
#   eb_Ftest1 = eBayes(limma_model1)
#   res = topTable(eb_Ftest1, n=nrow(eb_Ftest1), coef=colnames(design)[grepl('group',colnames(design))])
#   dt = data.table(assay='epigen-atac-seq',
#                   tissue=label,
#                   feature_ID=rownames(res),
#                   fscore=res$`F`,
#                   p_value = res$P.Value,
#                   adj_p_value = p.adjust(res$P.Value, method='BH'),
#                   full_model=gsub(' ','',full),
#                   reduced_model=gsub(' ','',reduced))
#   
#   if(!return_resid){
#     return(dt)
#   }
#   
#   # residuals
#   residual_mat = residuals(limma_model1, curr_voom)
#   
#   return(list(res=dt, 
#               residuals=residual_mat))
# 
# }


#' @param label t55-gastrocnemius,mssm
atac_training_dea_each_sex = function(filt_counts, meta, covariates, curr_outliers=c(), label=NULL, return_resid=F){
  
  meta_df = as.data.frame(meta, check.names=F)
  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
    meta_df[,which] = factor(meta_df[,which])
  }
  meta_df$group = factor(meta_df$group, levels=c('control','1w','2w','4w','8w')) # IMPORTANT
  
  # remove outliers 
  if(length(curr_outliers) > 0){
    meta_df = meta_df[!as.character(meta_df$viallabel) %in% as.character(curr_outliers),]
  }
  filt_counts = filt_counts[,as.character(meta_df$viallabel)]
  
  # split by sex 
  sex_res = list()
  resids = list()
  for(SEX in c('male','female')){
    
    curr_meta = meta_df[meta_df$sex==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    
    full = paste0(c("~ 1", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 1", covariates), collapse=" + ")
    
    # normalize and get voom weights 
    design = model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 1", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 1", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s", label, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = voom(curr_counts, design=design, normalize.method="quantile")
    
    limma_model1 = lmFit(curr_voom, design)
    eb_Ftest = eBayes(limma_model1)
    res = topTable(eb_Ftest, n=Inf, coef=colnames(design)[grepl('group',colnames(design))], sort.by = "none")
    dt = data.table(assay='epigen-atac-seq',
                    dataset=label, 
                    tissue=gsub(",.*","",label),
                    feature_ID=rownames(res),
                    fscore=res$`F`,
                    p_value = res$P.Value,
                    adj_p_value = p.adjust(res$P.Value, method='BH'),
                    full_model=gsub(' ','',full),
                    reduced_model=gsub(' ','',reduced),
                    removed_samples=paste0(curr_outliers, collapse=','),
                    covariates=paste0(curr_cov, collapse=','))
    sex_res[[SEX]] = dt
    
    # residuals
    residual_mat = residuals(limma_model1, curr_voom)
    resids[[SEX]] = residual_mat
  }
  
  # merge 
  merged = data.table(merge(sex_res[['male']], sex_res[['female']], 
                            by=c("assay","tissue","feature_ID","removed_samples"), 
                            suffixes=c('_male','_female')))
  
  # get a single meta p-value per feature using the male- and female- specific p-values 
  merged[,p_value := sumlog(c(p_value_male, p_value_female))$p, by=seq_len(nrow(merged))]
  merged[,adj_p_value := p.adjust(p_value, method='BH')]
  
  if(!return_resid){
    return(merged)
  }
  
  residual_mat = data.frame(cbind(resids[['male']], resids[['female']]))
  
  return(list(res=merged, 
              residuals=residual_mat))
}


# atac_timewise_dea = function(filt_counts, meta, covariates, label=NULL){
#   
#   meta[,sex_group := paste0(sex,'_',group)]
#   meta_df = as.data.frame(meta, check.names=F)
#   contrast = paste0(c('~ 0', 'sex_group', covariates), collapse=' + ')
#   mod1 = model.matrix(eval(parse(text=contrast)), meta_df)
#   voom_obj = voom(filt_counts, design=mod1, normalize.method="quantile")
#   fit = lmFit(voom_obj,mod1)
#   
#   cont.matrix=makeContrasts(
#     '1W_M'='sex_groupmale_1w - sex_groupmale_control',
#     '2W_M'='sex_groupmale_2w - sex_groupmale_control',
#     '4W_M'='sex_groupmale_4w - sex_groupmale_control',
#     '8W_M'='sex_groupmale_8w - sex_groupmale_control',
#     '1W_F'='sex_groupfemale_1w - sex_groupfemale_control',
#     '2W_F'='sex_groupfemale_2w - sex_groupfemale_control',
#     '4W_F'='sex_groupfemale_4w - sex_groupfemale_control',
#     '8W_F'='sex_groupfemale_8w - sex_groupfemale_control',
#     levels=mod1)
#   
#   fit2=contrasts.fit(fit,cont.matrix)
#   e=eBayes(fit2)
#   
#   dea_list = list()
#   for(tp in c('1W','2W','4W','8W')){
#     for(sex in c('M','F')){
#       comp = paste0(tp,'_',sex)
#       res = topTable(e, number=nrow(e), coef=comp, confint=TRUE)
#       dt = data.table(assay='epigen-atac-seq',
#                       feature_ID=rownames(res),
#                       sex=ifelse(sex=='M','male','female'),
#                       comparison_group=tolower(tp),
#                       logFC=res$logFC,
#                       logFC_se=(res$CI.R - res$CI.L)/3.92,
#                       tscore=res$t,
#                       p_value = res$P.Value,
#                       label = label)
#       dea_list[[comp]] = dt
#     }
#   }
#   
#   dea_ttest = rbindlist(dea_list)
#   dea_ttest[,adj_p_value := p.adjust(p_value, method='BY')]
#   return(dea_ttest)
# }


atac_timewise_dea_each_sex = function(filt_counts, meta, covariates, curr_outliers=c(), label=NULL){
  
  meta_df = as.data.frame(meta, check.names=F)
  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
    meta_df[,which] = factor(meta_df[,which])
  }
  # remove outliers 
  if(length(curr_outliers) > 0){
    meta_df = meta_df[!as.character(meta_df$viallabel) %in% as.character(curr_outliers),]
  }
  filt_counts = filt_counts[,as.character(meta_df$viallabel)]
  
  # split by sex 
  sex_res = list()
  for(SEX in c('male','female')){
    
    curr_meta = meta_df[meta_df$sex==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    
    full = paste0(c("~ 0", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 0", covariates), collapse=" + ")
    
    # normalize and get voom weights 
    design = model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 0", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 0", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s", label, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = voom(curr_counts, design=design, normalize.method="quantile")
    fit = lmFit(curr_voom,design)
    
    cont.matrix=makeContrasts(
      '1W'='group1w - groupcontrol',
      '2W'='group2w - groupcontrol',
      '4W'='group4w - groupcontrol',
      '8W'='group8w - groupcontrol',
      levels=design)
    
    fit2=contrasts.fit(fit,cont.matrix)
    e=eBayes(fit2)
    
    dea_list = list()
    for(tp in c('1W','2W','4W','8W')){
      res = topTable(e, number=nrow(e), coef=tp, confint=TRUE)
      dt = data.table(assay='epigen-atac-seq',
                      dataset=label, 
                      tissue=gsub(",.*","",label),
                      feature_ID=rownames(res),
                      sex=SEX,
                      comparison_group=tolower(tp),
                      logFC=res$logFC,
                      logFC_se=(res$CI.R - res$CI.L)/3.92,
                      tscore=res$t,
                      p_value = res$P.Value,
                      removed_samples=paste0(curr_outliers, collapse=','),
                      covariates=paste0(curr_cov, collapse=','))
      
      dea_list[[tp]] = dt
    }
    sex_res[[SEX]] = rbindlist(dea_list)
  }
  res = rbindlist(sex_res)
  res[,adj_p_value := p.adjust(p_value, method='BY')]
  return(res)
}

run_gost = function(input, background, 
                    organism='rnorvegicus', 
                    databases=c('REAC','WP','KEGG'),
                    min_pw_set_size=10,
                    max_pw_set_size=200,
                    return_gem=F){
  
  require("gprofiler2")
  require("data.table")
  
  res = gost(input,
             organism=organism,
             significant=FALSE, # return all results
             correction_method='fdr', # it will not return unadjusted p-values. calculate those yourself
             sources=databases,
             custom_bg=background,
             domain_scope='custom',
             evcodes=T)
  
  res_dt = data.table(res$result)
  
  if(nrow(res_dt)==0){
    warning("No results returned by gost.\n")
    return()
  }
  
  # calculate raw p-values
  res_dt[,gost_adj_p_value := p_value]
  res_dt[,p_value := NULL]
  res_dt[,computed_p_value := phyper(intersection_size-1, term_size, effective_domain_size-term_size, query_size,lower.tail= FALSE)]
  res_dt[,BH_adj_p_value := p.adjust(computed_p_value, method='BH')]
  
  # remove tests outside of the parameters
  if(min_pw_set_size>0 | max_pw_set_size<Inf){
    res_dt = res_dt[term_size<=max_pw_set_size & term_size>=min_pw_set_size]
  }
  
  # fix parents column 
  res_dt[,parents := sapply(parents, function(x) paste0(unlist(x), collapse=','))]
  
  res_dt[,BH_adj_p_value := p.adjust(computed_p_value, method='BH')]
  res_dt = res_dt[order(BH_adj_p_value, decreasing=F)]
  
  if(return_gem){
    res_dt[,Phenotype := 1]
    res_dt = res_dt[,.(term_id, term_name, computed_p_value, BH_adj_p_value, Phenotype, intersection)]
    colnames(res_dt) = c('GO.ID','Description','p.Val','FDR','Phenotype','Genes')
  }
  
  return(res_dt)
}
