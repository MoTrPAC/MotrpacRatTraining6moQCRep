#' Implementation of the MOP's section 9 flagging of problematic samples
#' @param pipeline_report_data A matrix of samples as rows and qc metrics as columns
#' @description Assumes that the input is a matrix generated by MoTrPAC's RNA-seq pipeline and flags potentially problematic samples.
rnaseq_pipeline_flagged_samples<-function(pipeline_report_data){
  pct_cols = colnames(pipeline_report_data)[
    grepl("pct_",colnames(pipeline_report_data))]
  pct_mapped_cols = pct_cols[grepl("mapped",pct_cols)]
  pct_unmapped_cols = pct_cols[grepl("unmapped",pct_cols)]
  total_unmapped = rowSums(pipeline_report_data[,pct_unmapped_cols])
  total_exonic = pipeline_report_data$pct_coding + pipeline_report_data$pct_utr
  flagged_sample_comments = rep("",nrow(pipeline_report_data))
  names(flagged_sample_comments) = rownames(pipeline_report_data)
  if("RIN" %in% colnames(pipeline_report_data)){
    currtest = which (pipeline_report_data$RIN < 6)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"RIN<6,")
  }
  if("reads" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$reads < 20000000)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"NumReads<20M,")
  }
  if("pct_GC" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_GC > 80)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"%GC>80,")
    currtest = which(pipeline_report_data$pct_GC < 20)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"%GC<20,")
  }
  if("pct_rRNA" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_rRNA > 20)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"%rRNA>20,")
  }
  if("pct_uniquely_mapped" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_uniquely_mapped < 60)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"%UniqueMap<60,")
  }
  
  currtest = which(total_exonic < 50)
  
  flagged_sample_comments = gsub(",$","",flagged_sample_comments,perl=T)
  flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"total_exonic<50,")
  return(flagged_sample_comments)
}

#' Implementation of the MOP's flagging of problematic samples
#' @param pipeline_report_data A matrix of samples as rows and qc metrics as columns
#' @description Assumes that the input is a matrix generated by MoTrPAC's ATAC-seq pipeline and flags potentially problematic samples.
atacseq_pipeline_flagged_samples<-function(pipeline_report_data){
  flagged_sample_comments = rep("",nrow(pipeline_report_data))
  names(flagged_sample_comments) = rownames(pipeline_report_data)
  if("align.nodup_samstat.total_reads" %in% colnames(pipeline_report_data)){
    currtest = which (pipeline_report_data$align.nodup_samstat.total_reads < 20*10^6)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"NumFiltReads<20M,")
  }else{
    warning("Column 'align.nodup_samstat.total_reads' not found.")
  }
  if("align.samstat.pct_mapped_reads" %in% colnames(pipeline_report_data)){
    currtest = which (pipeline_report_data$align.samstat.pct_mapped_reads < 80)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"%Map<80,")
  }else{
    warning("Column 'align.samstat.pct_mapped_reads' not found.")
  }
  if("peak_enrich.frac_reads_in_peaks.macs2.frip" %in% colnames(pipeline_report_data)){
    currtest = which (pipeline_report_data$peak_enrich.frac_reads_in_peaks.macs2.frip < 0.1)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"MacsFrip<0.1,")
  }else{
    warning("Column 'peak_enrich.frac_reads_in_peaks.macs2.frip' not found.")
  }
  if("align.frag_len_stat.nfr_peak_exists" %in% colnames(pipeline_report_data)){
    currtest = which (pipeline_report_data$align.frag_len_stat.nfr_peak_exists == FALSE)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"NoNfrPeak,")
  }else{
    warning("Column 'align.frag_len_stat.nfr_peak_exists' not found.")
  }
  if(all(c("align.frag_len_stat.mono_nuc_peak_exists","align.frag_len_stat.di_nuc_peak_exists") %in% colnames(pipeline_report_data))){
    currtest = which (pipeline_report_data$align.frag_len_stat.mono_nuc_peak_exists == FALSE & pipeline_report_data$align.frag_len_stat.di_nuc_peak_exists == FALSE)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"NoMonoDiNucPeaks,")
  }else{
    warning("One or both columns not found: 'align.frag_len_stat.mono_nuc_peak_exists', 'align.frag_len_stat.di_nuc_peak_exists'")
  }
  if("align_enrich.tss_enrich.tss_enrich" %in% colnames(pipeline_report_data)){
    currtest = which (pipeline_report_data$align_enrich.tss_enrich.tss_enrich < 4)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"TssEnrich<4,")
  }else{
    warning("Column 'align_enrich.tss_enrich.tss_enrich' not found.")
  }
  flagged_sample_comments = gsub(",$","",flagged_sample_comments,perl=T)
  return(flagged_sample_comments)
}


#' Implementation of the MOP's flagging of problematic samples for RRBS
#' @param pipeline_report_data A matrix of samples as rows and qc metrics as columns
#' @description Assumes that the input is a matrix generated by MoTrPAC's RRBS pipeline and flags potentially problematic samples.
rrbs_pipeline_flagged_samples<-function(pipeline_report_data){
  flagged_sample_comments = rep("",nrow(pipeline_report_data))
  names(flagged_sample_comments) = rownames(pipeline_report_data)
  if("reads" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$reads < 20000000)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"NumReads<20M,")
  }
  if("pct_Uniq" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_Uniq < 50)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"%UniqueMap<50,")
  }
  if("reads" %in% colnames(pipeline_report_data) && "pct_Uniq" %in% colnames(pipeline_report_data)){
    mapped_read_c = pipeline_report_data$reads * pipeline_report_data$pct_Uniq
    for(tissue in pipeline_report_data$Tissue){
      inds = pipeline_report_data$Tissue == tissue
      currtest = which(inds & mapped_read_c < 0.5 * median(mapped_read_c[inds],na.rm=T))
      flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"Num unique reads < median in tissue,")
    }
  }
  if("lambda_pct_Uniq" %in% colnames(pipeline_report_data) &&
     "lambda_pct_CpG" %in% colnames(pipeline_report_data)){
    currtest1 = which(pipeline_report_data$lambda_pct_Uniq > 0.5)
    currtest2 = which(pipeline_report_data$lambda_pct_CpG > 5)
    currtest = intersect(currtest1,currtest2)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],
                                               "lambda_pct_CpG>5 AND lambda_pct_Uniq >0.5,")
  }
  if("pct_OT" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_OT < 30)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"unexpected strands: pct_OT < 30,")
  }
  if("pct_OT" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_OT > 70)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"unexpected strands: pct_OT > 70,")
  }
  if("pct_CTOT" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_CTOT > 10)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"unexpected strands: pct_CTOT > 10,")
  }
  if("pct_CTOB" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_CTOB > 10)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"unexpected strands: pct_CTOB > 10,")
  }
  if("pct_umi_dup" %in% colnames(pipeline_report_data)){
    currtest = which(pipeline_report_data$pct_umi_dup > 20)
    flagged_sample_comments[currtest] = paste0(flagged_sample_comments[currtest],"pct_umi_dup>20,")
  }
  flagged_sample_comments = gsub(",$","",flagged_sample_comments,perl=T)
  return(flagged_sample_comments)
}

pass1a_get_numeric_timepoint<-function(x){
  v = rep(0,length(x))
  tps = c(1,0.5,4,7,24,48)
  for(tp in tps){
    v[grepl(paste0("\\s",tp,"\\s"),x,perl = T)] = tp
  }
  return(v)
}
pass1a_parse_pheno_add_covars<-function(pheno_bucket,pheno_local_dir,gsutil_cmd){
  
  pheno_data = parse_pheno_data(pheno_bucket,local_path = pheno_local_dir,
                                remove_prev_files = T,GSUTIL_PATH=gsutil_cmd)
  # add a tissue variable (for convinience)
  pheno_data$viallabel_data$tissue = 
    pheno_data$viallabel_data$specimen.processing.sampletypedescription
  # add the time to freeze variable ((for convinience))
  pheno_data$viallabel_data$time_to_freeze = 
    pheno_data$viallabel_data$calculated.variables.frozetime_after_acute - 
    pheno_data$viallabel_data$calculated.variables.deathtime_after_acute
  # add a binary is_control variable
  pheno_data$viallabel_data$is_control = as.numeric(grepl("control",
            pheno_data$viallabel_data$key.anirandgroup,ignore.case = T))
  # add the timepoint as a number
  # x - the text description of the group
  pheno_data$viallabel_data$timepoint = pass1a_get_numeric_timepoint(
    pheno_data$viallabel_data$key.anirandgroup
  )
  
  # Update the blood freeze times
  blood_samples = 
    pheno_data$viallabel_data$specimen.processing.sampletypedescription ==
    "EDTA Plasma"
  blood_freeze_times = 
    as.difftime(pheno_data$viallabel_data$specimen.processing.t_freeze,units = "mins") -
    as.difftime(pheno_data$viallabel_data$specimen.processing.t_edtaspin,units="mins")
  blood_freeze_times = as.numeric(blood_freeze_times)
  pheno_data$viallabel_data$time_to_freeze[blood_samples] = 
    blood_freeze_times[blood_samples]
  
  return(pheno_data)
}

