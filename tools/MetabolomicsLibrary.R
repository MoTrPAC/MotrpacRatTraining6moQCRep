#' An ensemble of functions used to assemble PASS1B metabolomics data freezes
#' ###########################################################################

#' AdjustPhenoNames
#' A function to adjust the column names of the phenotypic data in data frame format
#' @param pheno_df the pass1b_6m_viallabel_data.txt data in dataframe format
#' @param pheno_dict_df a data frame of "BICUniqueID" and "FullName" pairs from the merged_dictionary.txt file
#' @return a pheno_df data frame object with column names adjusted from "BICUniqueID" to "FullName"
AdjustPhenoNames <- function(pheno_df,pheno_dict_df){
  # Rename (does not catch all)
  for(i in 1:nrow(pheno_dict_df)){
    bic_name <- pheno_dict_df[i,'BICUniqueID'] %>% as.character()
    full_name <- pheno_dict_df[i,'FullName'] %>% as.character()
    n <- which(names(pheno_df) == bic_name)
    names(pheno_df)[n] <- full_name
  }
  # Correct the names that end in "_1" (just for "NMR.Testing")
  for(i in 1:nrow(pheno_dict_df)){
    bic_name <- pheno_dict_df[i,'BICUniqueID'] %>% as.character()
    full_name <- pheno_dict_df[i,'FullName'] %>% as.character()
    n <- which(names(pheno_df) == paste0(bic_name,'_1'))
    if(!length(names(pheno_df)[n]) == 0){
      if(grepl('_1', names(pheno_df)[n])){
        names(pheno_df)[n] <- paste0(full_name,'_1')
      }
    }
  }
  # Correct the names that end in "_2" (just for "NMR.Testing")
  for(i in 1:nrow(pheno_dict_df)){
    bic_name <- pheno_dict_df[i,'BICUniqueID'] %>% as.character()
    full_name <- pheno_dict_df[i,'FullName'] %>% as.character()
    n <- which(names(pheno_df) == paste0(bic_name,'_2'))
    if(!length(names(pheno_df)[n]) == 0){
      if(grepl('_2', names(pheno_df)[n])){
        names(pheno_df)[n] <- paste0(full_name,'_2')
      }
    }
  }
  # Set a vector for Exercise/Control Levels and Colors
  ec_levels <- c('One-week program',
                 'Two-week program',
                 'Four-week program',
                 'Eight-week program Training Group',
                 'Eight-week program Control Group')
  ec_colors <- c('gold',
                 'darkgoldenrod1',
                 'orange',
                 'darkorange',
                 'darkorange2',
                 'darkorange3',
                 'darkorange4',
                 'steelblue1',
                 'steelblue4')
  # Reorient datatypes
  pheno_df$pid <- as.factor(pheno_df$pid)
  pheno_df$bid <- as.factor(pheno_df$bid)
  pheno_df$labelid <- as.factor(pheno_df$labelid)
  pheno_df$viallabel <- as.factor(pheno_df$viallabel)
  pheno_df$Key.anirandgroup <- factor(pheno_df$Key.anirandgroup,
                                      levels = ec_levels)
  return(pheno_df)
}

#' CountNAs
#' A function to count the NA values in an object, typically a matrix
#' @param x a matrix
#' @return the number of NA values within the matrix
CountNAs <- function(x){
  return(sum(is.na(x)))
}

#' CountZeros
#' A function to count the number of zero values in an object, typically a matrix
#' @param x a matrix
#' @return the number of zero values within the matrix
CountZeros <- function(x){
  return(sum(x==0, na.rm = T))
}

#' CreateCountsMetaTables
#' A function to create a master count (abundance) data dataframe (nested) and a 
#' @param local_data_dir the local directory where the metabolomics data was downloaded
#' @return a list of master count dataframe (1) and meta dataframe (2)
CreateCountsMetaTables <- function(local_data_dir){
  
  #' ExtractValue
  #' Function only used within CreateCountsMetaTables() to extract values from experimentDetails files 
  #' @param string exp_text column name pattern to match 
  #' @return character value from exp_text
  ExtractValue <- function(string){
    val <- as.character(exp_text$VALUE[grepl(string, exp_text$ID)])
    if(length(val)==0){
      return('')
    }
    if(length(val)>1){
      if(!string %in% c('MS:MS_COMMENTS', 'AN:ANALYSIS_DETAILS', 'SP:SAMPLEPREP_SUMMARY', 'CH:CHROMATOGRAPHY_SUMMARY', 'ST:STUDY_SUMMARY')){
        warning(sprintf("More than one value for %s:\n%s", string, paste0(val, collapse='\n')))
      }
      val <- paste0(val, collapse=' ')
    }
    return(val)
  }
  
  # Create empty dataframes
  metadata_df <- countdata_df <- sampledata_df <- metabolitedata_df <- data.frame()
  # Collect targeted metabolites
  OME = 'metabolomics'
  for(DATASET in c('targeted','untargeted')){
    message("\n", toupper(DATASET), "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # Update data directory
    if(DATASET == 'untargeted'){
      tissue_dir <- file.path(local_data_dir,'metabolomics-untargeted')
    }else if(DATASET == 'targeted'){
      tissue_dir <- file.path(local_data_dir,'metabolomics-targeted')
    }
    # List of tissues
    ls <- list.files(tissue_dir) %>% as.vector()
    tissue_dirs <- ls[grepl("t[0-9][0-9]", ls)]
    #t_dir <- tissue_dirs[1]
    # Iterate through the tissues
    for(t_dir in tissue_dirs){
      # Collect Tissue
      TISSUE <- (str_split(t_dir, pattern = "t[0-9]+-") %>% unlist())[2]
      message("\n- ",TISSUE, ": ", appendLF = FALSE)
      # Collect list of technologies
      t_dir <- file.path(tissue_dir, t_dir)
      tech_dirs <- list.files(t_dir) %>% as.vector()
      #tech_dir <- tech_dirs[3]
      for(assay in tech_dirs){
        message(assay,", ", appendLF = FALSE)
        # Collect the METAB_FAMILY
        if(DATASET == 'untargeted'){
          METAB_FAMILY <- (str_split(assay, pattern = "-u-") %>% unlist())[2]
          names <- c('named', 'unnamed')
        }else if(DATASET == 'targeted'){
          METAB_FAMILY <- ((str_split(assay, pattern = "-t-") %>% unlist())[2] %>% str_split(pattern = "-") %>% unlist())[1]
          names <- c('named')
        }
        tech_dir <- file.path(t_dir, assay)
        # Collect variables for named and unnamed compounds
        for(NAMED in names){
          ## Experimental Details-----
          exp_file <- list.files(normalizePath(tech_dir),
                                 pattern=paste0("_", NAMED, "-experimentalDetails.txt"),
                                 ignore.case = TRUE,
                                 full.names = TRUE,
                                 recursive = TRUE)
          # Collect variables
          exp_text <- read.table(file = exp_file, sep = '\t')
          names(exp_text) <- c('ID','VALUE')
          
          #STUDY_TITLE<- filter(exp_text, grepl('ST:STUDY_TITLE', ID)) %>% select(VALUE) %>% str_c() %>% str_squish()
          # above line returns "12" instead of "Molecular Transducers of Physical Activity Consortium (MoTrPAC)"
          STUDY_TITLE<- ExtractValue('ST:STUDY_TITLE')
          STUDY_TYPE<- ExtractValue('ST:STUDY_TYPE')
          STUDY_SUMMARY<- ExtractValue('ST:STUDY_SUMMARY')
          STUDY_INSTITUTE<- ExtractValue('ST:INSTITUTE')
          STUDY_DEPARTMENT<- ExtractValue('ST:DEPARTMENT')
          STUDY_LABORATORY<- ExtractValue('ST:LABORATORY')
          STUDY_LAST_NAME<- ExtractValue('ST:LAST_NAME')
          ST_NUM_GROUPS<- ExtractValue('ST:NUM_GROUPS')
          SUBMIT_DATE<- ExtractValue('ST:SUBMIT_DATE')
          STUDY_COMMENTS<- ExtractValue('ST:STUDY_COMMENTS')
          SUBJECT_TYPE<- ExtractValue('SU:SUBJECT_TYPE')
          SUBJECT_SPECIES<- ExtractValue('SU:SUBJECT_SPECIES')
          SAMPLEPREP_SUMMARY<- ExtractValue('SP:SAMPLEPREP_SUMMARY')
          SAMPLEPREP_PROTOCOL_FILENAME<- ExtractValue('SP:SAMPLEPREP_PROTOCOL_FILENAME')
          CH_CHROMATOGRAPHY_TYPE<- ExtractValue('CH:CHROMATOGRAPHY_TYPE')
          CH_INSTRUMENT_NAME<- ExtractValue('CH:INSTRUMENT_NAME')
          CH_COLUMN_NAME<- ExtractValue('CH:COLUMN_NAME')
          CH_METHODS_FILENAME<- ExtractValue('CH:METHODS_FILENAME')
          CH_SUMMARY<- ExtractValue('CH:CHROMATOGRAPHY_SUMMARY')
          MS_INSTRUMENT_TYPE<- ExtractValue('MS:INSTRUMENT_TYPE')
          MS_INSTRUMENT_NAME<- ExtractValue('MS:INSTRUMENT_NAME')
          MS_TYPE<- ExtractValue('MS:MS_TYPE')
          MS_ION_MODE<- ExtractValue('MS:ION_MODE')
          MS_UNITS<- ExtractValue('MS_METABOLITE_DATA:UNITS')
          MS_COMMENTS<- ExtractValue('MS:MS_COMMENTS')
          MS_RESULTS_FILE<- ExtractValue('MS:MS_RESULTS_FILE')
          ANALYSIS_TYPE<- ExtractValue('AN:ANALYSIS_TYPE')
          ANALYSIS_DETAILS<- ExtractValue('AN:ANALYSIS_DETAILS')
  
          ## Metabolite Metadata-----
          metab_meta_file <- list.files(normalizePath(tech_dir),
                                        pattern=paste0("_", NAMED, "-metadata-metabolites.txt"),
                                        ignore.case = TRUE,
                                        full.names = TRUE,
                                        recursive = TRUE)
          # Load the file
          metab_meta_data <- read.csv(file = metab_meta_file, sep = '\t', 
                                      header = T, fill = T)
          METABOLITE_NAMES <- metab_meta_data %>%
            select(metabolite_name) %>% unlist() %>% 
            unname() %>% paste(collapse = ';')
          METABOLITE_N <- metab_meta_data %>%
            select(metabolite_name) %>% unlist() %>%
            length()
          
          #TODO: Create an annotation column for internal standard N
          ## Sample Metadata-----
          sample_meta_file <- list.files(normalizePath(tech_dir),
                                         pattern=paste0("_", NAMED, "-metadata-samples.txt"),
                                         ignore.case = TRUE,
                                         full.names = TRUE,
                                         recursive = TRUE)
          # Load the file
          sample_meta_data <- read.csv(file = sample_meta_file, 
                                       sep = '\t' , 
                                       header = T) %>% arrange(sample_order)
          sample_meta_data$sample_id <- as.character(sample_meta_data$sample_id)
          sample_meta_data$sample_type <- as.character(sample_meta_data$sample_type)
          SAMPLE_NAMES <- sample_meta_data %>%
            select(sample_id) %>% unlist() %>% 
            unname() %>% paste(collapse = ';')
          SAMPLE_N <- sample_meta_data %>%
            filter(sample_type == 'Sample') %>%
            nrow()
          QC_IS_N <- sample_meta_data %>%
            filter(sample_type == 'QC-InternalStandard') %>%
            nrow()
          QC_PRERUN_N <- sample_meta_data %>%
            filter(sample_type == 'QC-PreRun') %>%
            nrow()
          QC_BLANK_N <- sample_meta_data %>%
            filter(sample_type == 'QC-Blank') %>%
            nrow()
          QC_POOLED_N <- sample_meta_data %>%
            filter(sample_type == 'QC-Pooled') %>%
            nrow()
          QC_REFERENCE_N <- sample_meta_data %>%
            filter(sample_type == 'QC-Reference') %>%
            nrow()
          QC_DRIFT_N <- sample_meta_data %>%
            filter(sample_type == 'QC-DriftCorrection') %>%
            nrow()
          
          ## Collect the counts and append with annotation----
          # Locate the file
          count_file <- list.files(normalizePath(tech_dir),
                                 pattern=paste0("_", NAMED, "-results.txt"),
                                 ignore.case = TRUE,
                                 full.names = TRUE,
                                 recursive = TRUE)
          # Load the file
          count_df <- read.csv(file = count_file, sep = '\t' , 
                               header = T, check.names = F)
          if(sum(duplicated(names(count_df))) > 0){
            print('CAUTION: Duplicate column names')
          }
          names(count_df) <- make.unique(names(count_df))
          
          # Convert the metabolite names to refmet names
          if(NAMED == 'named'){
            metab_meta_names <- metab_meta_data %>% 
              select(metabolite_name, refmet_name)
            count_df <- left_join(count_df, metab_meta_names, by = c('metabolite_name')) %>%
              select(-refmet_name) %>%
              select(metabolite_name, everything())
            names(count_df)[1] <- 'metabolite_name'
          }
          # Gather the dataframe by samples
          gather_names <- names(count_df)[names(count_df)!='metabolite_name']
          
          ## Transform Counts with Different Strategies----
          # Dataframe to matrix
          if(sum(duplicated(count_df$metabolite_name)) > 0){
            dupes <- count_df$metabolite_name[duplicated(count_df$metabolite_name)]
            print(paste0('CAUTION: Duplicate metabolite names in: ',tech_dir))
            print(paste0('duplicate metabolites: ', paste0(dupes, collapse=', ')))
          }
          count_df <- count_df[!duplicated(count_df$metabolite_name),]
          row.names(count_df) <- count_df$metabolite_name %>% as.character()
          count_mat <- count_df %>%
            select(-metabolite_name) %>%
            as.matrix() %>% t()
          
          # Add Sample Annotation to the Count Data
          pheno_join <- pheno_df %>% select(pid,bid,labelid,viallabel)
          pheno_join$pid <- as.character(pheno_join$pid)
          pheno_join$bid <- as.character(pheno_join$bid)
          pheno_join$labelid <- as.character(pheno_join$labelid)
          pheno_join$viallabel <- as.character(pheno_join$viallabel)
          
          ## Original Abundances-----
          count_org <- count_mat %>%
            t() %>% as.data.frame() %>%
            tibble::rownames_to_column(var = "METABOLITE_NAME") %>%
            tidyr::gather(all_of(gather_names), 
                          key = "viallabel", value = "VALUE") %>%
            left_join(y = pheno_join, by = c('viallabel')) %>%
            select(METABOLITE_NAME, labelid, viallabel, pid, bid, VALUE) %>%
            mutate(JOIN_KEY='123') %>%
            group_by(JOIN_KEY) %>%
            nest(COUNT_DATA = -JOIN_KEY) %>%
            ungroup()
          
          ## Check if any viallabels are missing from pheno data----
          original_vl <- sample_meta_data$sample_id[which(sample_meta_data$sample_type == "Sample")]
          tmp_count <- count_org$COUNT_DATA[[1]]
          # current_vl <- unique(tmp_count$viallabel[!is.na(tmp_count$pid)])
          current_vl <- unique(tmp_count$viallabel)
          if(length(current_vl) < length(original_vl)){
            # message("\nLength current: ", length(current_vl), " and length original: ", length(original_vl))
            message(sprintf("\n\tThere are viallabels in %s (%s) that don't exist in the phenotypic or/and results data: %s\n\tThis will not effect the counts tables but will effect the Kruskal-Wallis tests.",
                            assay, 
                            NAMED,
                            paste0(original_vl[!original_vl %in% current_vl], collapse=', ')))
          }
          
          ## Deal with missing variables------
          # If Variables were not collected, make them NA
          vects <- list(OME,DATASET,TISSUE,METAB_FAMILY,NAMED,
                        STUDY_INSTITUTE,STUDY_TITLE,STUDY_TYPE,STUDY_SUMMARY,
                        STUDY_DEPARTMENT,STUDY_LABORATORY,STUDY_LAST_NAME,
                        ST_NUM_GROUPS,SUBMIT_DATE,STUDY_COMMENTS,
                        SUBJECT_TYPE,SUBJECT_SPECIES,SAMPLEPREP_SUMMARY,
                        SAMPLEPREP_PROTOCOL_FILENAME,CH_CHROMATOGRAPHY_TYPE,
                        CH_INSTRUMENT_NAME,CH_COLUMN_NAME,CH_METHODS_FILENAME,
                        CH_SUMMARY,MS_INSTRUMENT_TYPE,MS_INSTRUMENT_NAME,MS_TYPE,
                        MS_ION_MODE,MS_UNITS,MS_COMMENTS,MS_RESULTS_FILE,ANALYSIS_TYPE,
                        ANALYSIS_DETAILS,METABOLITE_NAMES,METABOLITE_N,SAMPLE_NAMES, SAMPLE_N,
                        QC_IS_N,QC_PRERUN_N,QC_BLANK_N,QC_POOLED_N,QC_REFERENCE_N,QC_DRIFT_N)
          
          names(vects) <- c("OME","DATASET","TISSUE","METAB_FAMILY",
                            "NAMED","STUDY_INSTITUTE","STUDY_TITLE","STUDY_TYPE","STUDY_SUMMARY",
                            "STUDY_DEPARTMENT","STUDY_LABORATORY","STUDY_LAST_NAME",
                            "ST_NUM_GROUPS","SUBMIT_DATE","STUDY_COMMENTS","SUBJECT_TYPE",
                            "SUBJECT_SPECIES","SAMPLEPREP_SUMMARY","SAMPLEPREP_PROTOCOL_FILENAME",
                            "CH_CHROMATOGRAPHY_TYPE","CH_INSTRUMENT_NAME","CH_COLUMN_NAME",
                            "CH_METHODS_FILENAME","CH_SUMMARY","MS_INSTRUMENT_TYPE",
                            "MS_INSTRUMENT_NAME","MS_TYPE","MS_ION_MODE","MS_UNITS",
                            "MS_COMMENTS","MS_RESULTS_FILE","ANALYSIS_TYPE",
                            "ANALYSIS_DETAILS","METABOLITE_NAMES","METABOLITE_N",
                            "SAMPLE_NAMES","SAMPLE_N","QC_IS_N","QC_PRERUN_N","QC_BLANK_N",
                            "QC_POOLED_N","QC_REFERENCE_N","QC_DRIFT_N")
          for(i in 1:length(vects)){
            if(length(vects[[i]]) == 0){
              x <- names(vects)[i]
              assign(x, NA)
            }
          }       
          ## Create the Nested Sample Data------
          # # Add the annotation
          sample_meta_data$DATASET <- DATASET
          sample_meta_data$TISSUE <- TISSUE
          sample_meta_data$METAB_FAMILY <- METAB_FAMILY
          sample_meta_data$NAMED <- NAMED
          sample_meta_data$STUDY_INSTITUTE <- STUDY_INSTITUTE
          sample_meta_data$CH_CHROMATOGRAPHY_TYPE <- CH_CHROMATOGRAPHY_TYPE
          sample_meta_data$MS_TYPE <- MS_TYPE
          sample_meta_data$MS_ION_MODE <- MS_ION_MODE
          # Nest sample data in with metadata
          sample_df <- sample_meta_data %>% 
            group_by(DATASET,TISSUE,METAB_FAMILY,NAMED,
                     STUDY_INSTITUTE,CH_CHROMATOGRAPHY_TYPE,
                     MS_TYPE,MS_ION_MODE) %>%
            nest(.key = "SAMPLE_DATA") %>%
            ungroup()
          
          # Potential change
          # sample_df <- sample_meta_data %>% 
          #   group_by(DATASET,TISSUE,METAB_FAMILY,NAMED,
          #            STUDY_INSTITUTE,CH_CHROMATOGRAPHY_TYPE,
          #            MS_TYPE,MS_ION_MODE) %>%
          #   nest(SAMPLE_DATA = c(sample_id, sample_type, sample_order, raw_file)) %>%
          #   ungroup()
          
          sampledata_df <- rbind(sampledata_df, sample_df)
          
          ## Create the Nested Metabolite Data------
          # # Add the annotation
          metab_meta_data$DATASET <- DATASET
          metab_meta_data$TISSUE <- TISSUE
          metab_meta_data$METAB_FAMILY <- METAB_FAMILY
          metab_meta_data$NAMED <- NAMED
          metab_meta_data$STUDY_INSTITUTE <- STUDY_INSTITUTE
          metab_meta_data$CH_CHROMATOGRAPHY_TYPE <- CH_CHROMATOGRAPHY_TYPE
          metab_meta_data$MS_TYPE <- MS_TYPE
          metab_meta_data$MS_ION_MODE <- MS_ION_MODE
          # Nest sample data in with metadata
          metabolite_df <- metab_meta_data %>% 
            group_by(DATASET,TISSUE,METAB_FAMILY,NAMED,
                     STUDY_INSTITUTE,CH_CHROMATOGRAPHY_TYPE,
                     MS_TYPE,MS_ION_MODE) %>%
            nest(.key = "METABOLITE_DATA") %>%
            ungroup()
          # WARNING: to deal with it, this alternative change would require adjustments
          # metabolite_df <- metab_meta_data %>% 
          #   group_by(DATASET,TISSUE,METAB_FAMILY,NAMED,
          #            STUDY_INSTITUTE,CH_CHROMATOGRAPHY_TYPE,
          #            MS_TYPE,MS_ION_MODE) %>%
          #   nest(METABOLITE_DATA = c(metabolite_name, refmet_name, rt, mz, neutral_mass, formula)) %>%
          #   ungroup()
          
          metabolitedata_df <- rbind(metabolitedata_df, metabolite_df)
          
          ## Create the Nested CountData------
          # # Add the annotation
          # Merge the countdata
          # Nest count data in with metadata
          # Join the dataframes and nest
          count_df <- data.frame(DATASET=DATASET, TISSUE=TISSUE,
                                 METAB_FAMILY=METAB_FAMILY, NAMED=NAMED,
                                 STUDY_INSTITUTE=STUDY_INSTITUTE,
                                 CH_CHROMATOGRAPHY_TYPE=CH_CHROMATOGRAPHY_TYPE,
                                 MS_TYPE=MS_TYPE,MS_ION_MODE=MS_ION_MODE,
                                 JOIN_KEY='123') %>%
            left_join(y = count_org, by =c('JOIN_KEY')) %>%
            # left_join(y = count_as, by =c('JOIN_KEY')) %>%
            # left_join(y = count_rs, by =c('JOIN_KEY')) %>%
            # left_join(y = count_ps, by =c('JOIN_KEY')) %>%
            # left_join(y = count_vs, by =c('JOIN_KEY')) %>%
            # left_join(y = count_ls, by =c('JOIN_KEY')) %>%
            # left_join(y = count_lc, by =c('JOIN_KEY')) %>%
            # left_join(y = count_las, by =c('JOIN_KEY')) %>%
            # left_join(y = count_lrs, by =c('JOIN_KEY')) %>%
            # left_join(y = count_lps, by =c('JOIN_KEY')) %>%
            # left_join(y = count_lvs, by =c('JOIN_KEY')) %>%
            # left_join(y = count_lls, by =c('JOIN_KEY')) %>%
          # left_join(y = count_pc, by =c('JOIN_KEY')) %>%
          # left_join(y = count_pas, by =c('JOIN_KEY')) %>%
          # left_join(y = count_prs, by =c('JOIN_KEY')) %>%
          # left_join(y = count_pps, by =c('JOIN_KEY')) %>%
          # left_join(y = count_pvs, by =c('JOIN_KEY')) %>%
          # left_join(y = count_pls, by =c('JOIN_KEY')) %>%
          select(-JOIN_KEY)
          
          countdata_df <- rbind(countdata_df, count_df)
          
          ## Merge the MetaData-----
          int_df <- data.frame(OME,DATASET,TISSUE,METAB_FAMILY,NAMED,
                               STUDY_INSTITUTE,STUDY_TITLE,STUDY_TYPE,STUDY_SUMMARY,
                               STUDY_DEPARTMENT,STUDY_LABORATORY,STUDY_LAST_NAME,
                               ST_NUM_GROUPS,SUBMIT_DATE,STUDY_COMMENTS,
                               SUBJECT_TYPE,SUBJECT_SPECIES,SAMPLEPREP_SUMMARY,
                               SAMPLEPREP_PROTOCOL_FILENAME,CH_CHROMATOGRAPHY_TYPE,
                               CH_INSTRUMENT_NAME,CH_COLUMN_NAME,CH_METHODS_FILENAME,
                               CH_SUMMARY,
                               MS_INSTRUMENT_TYPE,MS_INSTRUMENT_NAME,MS_TYPE,MS_ION_MODE,
                               MS_UNITS,
                               MS_COMMENTS,MS_RESULTS_FILE,ANALYSIS_TYPE,ANALYSIS_DETAILS,
                               METABOLITE_NAMES,METABOLITE_N,SAMPLE_NAMES,SAMPLE_N,
                               QC_IS_N,QC_PRERUN_N,QC_BLANK_N,QC_POOLED_N,QC_REFERENCE_N,QC_DRIFT_N)
          metadata_df <- rbind(metadata_df, int_df)
          
          ## Reassign all variables to NA or empty-----
          resets <- vects[!names(vects) %in% c('OME','DATASET','TISSUE','METAB_FAMILY','NAMED')]
          for(i in 1:length(resets)){
            x <- names(resets)[i]
            assign(x, NA)
          }
          count_df <- sample_df <- metabolite_df <- data.frame()
        }
      }
    }
    message("\n")
  }
  
  # Join the nested count data, sample data, and metabolite data----
  countdata_df <- left_join(countdata_df, sampledata_df, by = c("DATASET","TISSUE","METAB_FAMILY","NAMED",
                                                                "STUDY_INSTITUTE","CH_CHROMATOGRAPHY_TYPE",
                                                                "MS_TYPE","MS_ION_MODE"))
  countdata_df <- left_join(countdata_df, metabolitedata_df, by = c("DATASET","TISSUE","METAB_FAMILY",
                                                                    "NAMED","STUDY_INSTITUTE","CH_CHROMATOGRAPHY_TYPE",
                                                                    "MS_TYPE","MS_ION_MODE"))
  
  # Make custom edits to "Key" columns-----
  # Metadata
  metadata_df <- metadata_df %>%
    mutate(CH_CHROMATOGRAPHY_TYPE = case_when(
      CH_CHROMATOGRAPHY_TYPE == "GC" ~ "GC",
      CH_CHROMATOGRAPHY_TYPE == "gas phase" ~ "GC",
      CH_CHROMATOGRAPHY_TYPE == "Rerverse phase" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "RP" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "Reverse phase" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "Reversed phase" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "flow injection" ~ "FIC",
      CH_CHROMATOGRAPHY_TYPE == "HILIC" ~ "HILIC"))
  metadata_df <- metadata_df %>%
    mutate(MS_ION_MODE = case_when(
      MS_ION_MODE == "negative" ~ "negative",
      MS_ION_MODE == "Negative" ~ "negative",
      MS_ION_MODE == "NEGATIVE" ~ "negative",
      MS_ION_MODE == "positive" ~ "positive",
      MS_ION_MODE == "Positive" ~ "positive",
      MS_ION_MODE == "POSITIVE" ~ "positive"))
  metadata_df <- metadata_df %>%
    mutate(MS_ION_MODE = ifelse(METAB_FAMILY == "rppos", "positive", MS_ION_MODE))
  
  # Count Data----
  countdata_df <- countdata_df %>%
    mutate(CH_CHROMATOGRAPHY_TYPE = case_when(
      CH_CHROMATOGRAPHY_TYPE == "GC" ~ "GC",
      CH_CHROMATOGRAPHY_TYPE == "gas phase" ~ "GC",
      CH_CHROMATOGRAPHY_TYPE == "Rerverse phase" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "RP" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "Reverse phase" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "Reversed phase" ~ "RPC",
      CH_CHROMATOGRAPHY_TYPE == "flow injection" ~ "FIC",
      CH_CHROMATOGRAPHY_TYPE == "HILIC" ~ "HILIC"))
  countdata_df <- countdata_df %>%
    mutate(MS_ION_MODE = case_when(
      MS_ION_MODE == "negative" ~ "negative",
      MS_ION_MODE == "Negative" ~ "negative",
      MS_ION_MODE == "NEGATIVE" ~ "negative",
      MS_ION_MODE == "positive" ~ "positive",
      MS_ION_MODE == "Positive" ~ "positive",
      MS_ION_MODE == "POSITIVE" ~ "positive"))
  countdata_df <- countdata_df %>%
    mutate(MS_ION_MODE = ifelse(METAB_FAMILY == "rppos", "positive", MS_ION_MODE))
  return(list(countdata_df, metadata_df))
}

#' DataGroupMeta
#' A function to extract the data group-specific metadata
#' @param countdata_df the master count (abundance) dataframe
#' @param named named/unnamed
#' @param site where metabolomics measured
#' @param platform metabolomics measured
#' @return a dataframe with sample-specific annotations in a unique annotation (helpful for data visualizations)
DataGroupMeta <- function(countdata_df,named,site,platform){
  meta_df1 <- countdata_df %>%
    ungroup() %>%
    select(-SAMPLE_DATA, -METABOLITE_DATA) %>%
    filter(NAMED == named) %>% filter(STUDY_INSTITUTE == site) %>% 
    filter(METAB_FAMILY == platform) %>%
    unnest(COUNT_DATA) %>%
    select(DATASET,TISSUE,METAB_FAMILY,NAMED,STUDY_INSTITUTE, METABOLITE_NAME,CH_CHROMATOGRAPHY_TYPE,
           MS_TYPE,MS_ION_MODE,METABOLITE_NAME,labelid,viallabel,pid,bid,VALUE) %>%
    unite(TISSUE, viallabel, col = "TISSUE_viallabel", remove = F, sep = ':') %>%
    arrange(TISSUE_viallabel)
  meta_df2 <- countdata_df %>%
    ungroup() %>%
    select(-COUNT_DATA, -METABOLITE_DATA) %>%
    filter(NAMED == named) %>% filter(STUDY_INSTITUTE == site) %>% 
    filter(METAB_FAMILY == platform) %>%
    unnest(SAMPLE_DATA) %>%
    select(TISSUE,sample_id,sample_type,sample_order) %>%
    unite(TISSUE, sample_id, col = "TISSUE_sample_id", remove = T, sep = ':') %>%
    arrange(TISSUE_sample_id)
  if(!all(meta_df1$TISSUE_viallabel %in% meta_df2$TISSUE_sample_id)){
    stop('metadata dataframes are not identical in samples')
  }
  meta_df <- left_join(meta_df1, meta_df2, by = c("TISSUE_viallabel" = "TISSUE_sample_id"))
  # Sample Data
  sample_df <- meta_df %>%
    select(TISSUE_viallabel,TISSUE,viallabel,bid,sample_type,sample_order) %>% unique() %>%
    arrange(TISSUE, sample_order)
  
  # Collect all reference samples and blank samples
  ##########################################################
  #ref_sams <- sample_df %>% filter(sample_type != 'Sample') %>% select(TISSUE_viallabel) %>% unlist() %>% as.character()
  return(sample_df)
}

#' DownloadBucketLocal
#' A function that retrieves data from a google bucket and saves it to a local directory
#' @param bucket the gsutil URI of the google bucket
#' @param local_data_dir the local directory to store the folder/data to be downloaded
#' @param gsutil_cmd the local PATH to the gsutil command
#' Note If you experience problems with multiprocessing on MacOS, they might be related to https://bugs.python.org/issue33725. You can disable multiprocessing by editing your .boto config or by adding the following flag to your command: `-o "GSUtil:parallel_process_count=1"`. Note that multithreading is still available even if you disable multiprocessing.
DownloadBucketLocal <- function(bucket, local_data_dir, gsutil_cmd){
  if(!dir.exists(file.path(local_data_dir))){
    dir.create(file.path(local_data_dir), recursive = TRUE)
  }
  load_cmd <- paste0(gsutil_cmd,' -o  "GSUtil:parallel_process_count=1" -m cp -r ',bucket,' ',local_data_dir)
  system(load_cmd)
  print(paste0('Data saved to : ', local_data_dir))
}

#' KnnImpute
#' A function to apply knn imputation to a numeric matrix
#' @param mat a numeric matrix
#' @param ref_sams2 a vector of reference samples (used to separate reference samples from experimental samples)
#' @param strategy a decision whether to perform knn using rows of matrix or columns of matrix as "neighbors" (currently only coded for 1 decision, second decision can be coded if needed)
#' @param k the number of neighbors used for knn
#' @param impute_bool a boolean whether to impute. If FALSE, the input matrix is returned
#' @param log_bool a boolean whether to log values before applying knn imputation
#' @param plot a boolean wether to plot and visualize before and after imputation (both strategies are visualized).
#' @return the matrix with imputed values
KnnImpute <- function(mat, ref_sams2, strategy = 1, k = 10, impute_bool = T, plot = F, log_bool = T){
  require("impute")
  
  if(impute_bool){
    
    # Only log2-transform before imputation if log_bool
    if(log_bool){
      mat <- log2(mat)
    }
    
    # Calculate features to impute
    fna.mat1 <- apply(mat[!row.names(mat) %in% ref_sams2,], 2, CountNAs)
    (( feat.imp1 <- names(fna.mat1)[fna.mat1 > 0] ))
    print(feat.imp1)
    (( feat.imp2 <- match(feat.imp1, colnames(mat)) ))
    tissue.imp1 <- impute.knn(mat[!rownames(mat) %in% ref_sams2,],k=k)$data # use nearest samples
    tissue.imp2 <- t(impute.knn(t(mat[!rownames(mat) %in% ref_sams2,]),k=k)$data) # use nearest features
    
    # Plot the imputed values
    if(plot){
      redblue100 = colorRampPalette(colors = c("red", "white", "blue"))(100)
      par(mfrow=c(3,1),bg="black")
      image(as.matrix(log(mat[!row.names(mat) %in% ref_sams2,feat.imp2])), col=redblue100,axes=F )
      image(as.matrix(log(tissue.imp1[,feat.imp2])), col=redblue100,axes=F)
      image(as.matrix(log(tissue.imp2[,feat.imp2])), col=redblue100,axes=F)
    }
    # Perform imputation for the reference samples
    if(strategy == 1){
      fna.mat1 <- apply(mat, 2, CountNAs)
      feat.imp1 <- names(fna.mat1)[fna.mat1 > 0]
      feat.imp2 <- match(feat.imp1, colnames(mat))
      tissue.imp3 <- impute.knn(mat,k=k)$data
      n <- dim(mat[!rownames(mat) %in% ref_sams2,])[1]
      if(!all(row.names(tissue.imp3)[1:n] == row.names(tissue.imp1))){
        stop('Cannot combine 2 imputed matrices')
      }else{
        tissue.imp3[1:n,] <-  tissue.imp1
        mat <- tissue.imp3
      }
    }
    
    # Un-log before returning the matrix 
    if(log_bool){
      mat <- 2^mat
    }
  }
  
  return(mat)
}

#' LoadRawMat
#' A function to extract the metabolomics abundance data
#' @param countdata_df the master count (abundance) dataframe
#' @param dataset targeted/untargeted
#' @param named named/unnamed
#' @param site where metabolomics measured
#' @param platform metabolomics measured
#' @param tissue tissue
#' @return a raw abundance matrix with column and row names sorted alphanumerically
LoadRawMat <- function(countdata_df, dataset, named, site, platform, tissue){
  mat <- countdata_df %>% 
    filter(DATASET == dataset) %>%
    filter(NAMED == named) %>%
    filter(STUDY_INSTITUTE == site) %>%
    filter(METAB_FAMILY == platform) %>%
    filter(TISSUE == tissue) %>% 
    ungroup() %>% select(COUNT_DATA) %>% unnest() %>% select(METABOLITE_NAME,viallabel,VALUE) %>%
    pivot_wider(names_from = 'METABOLITE_NAME', values_from = 'VALUE') %>%
    column_to_rownames(var = 'viallabel') %>% as.matrix()
  mat <- mat[sort(row.names(mat)), sort(colnames(mat))]
  return(mat)
}

#' NAFeatureFilter
#' A function to remove features based on the proportion of NA values across features
#' @param mat a numeric matrix
#' @param na_bool a boolean to determine if filtering should be performed
#' @param threshold_freq the threshold frequency of NA values used to filter (remove) features with NA values above the threshold
#' @param tissue_name the name of the tissue
#' @param ref_sams2 a vector of reference samples (used to separate reference samples from experimental samples)
#' @return the matrix with features removed
NAFeatureFilter <- function(mat, na_bool = T, threshold_freq = 0.2, tissue_name, ref_sams2){
  if(na_bool){
    mat2 <- mat[!row.names(mat) %in% ref_sams2,]
    # NA frequency
    fna.mat <- apply(mat2, 2, CountNAs) %>% sort() %>% rev()
    fna.mat_freq <- fna.mat/nrow(mat2)
    # Remove samples
    frm_na <- names(fna.mat_freq[fna.mat_freq > threshold_freq])
    print(paste0(length(frm_na),' Features removed from ',tissue_name,' due to NA values:'))
    print(paste0(frm_na, collapse=', '))
    # Remove the samples from the individual matrices
    mat <- mat[,!colnames(mat) %in% frm_na]
  }
  # Output
  return(mat)
}

#' NormalizeData
#' A function to normalize numeric matrices based on different strategies
#' @param tissue_mat a numeric matrix
#' @param strategy the normalization strategy. Strategies include: none ('raw'), log2 ('log2'), 
#'   autoscaling ('autoscale'), log2 + feature standardization ('log2_fstd'), log2 + feature centering ('log2_fcen'), 
#'   log2 + sample standardization ('log2_sstd'), log2 + sample centering ('log2_scen'), 
#'   log2 + feature standardization + sample centering ('log2_fstd_scen'), 
#'   log2 + feature standardization + sample standardization ('log2_fstd_sstd'), 
#'   log2 + sample standardization + feature standardization ('log2_sstd_fstd'), 
#'   log2 + sample standardization + feature centering ('log2_sstd_fcen'), and 
#'   log2 + sample standardization + feature standardization + sample centering ('log2_sstd_fstd_scen')
#' @param ref_sams2 a vector of reference samples (used to separate reference samples from experimental samples)
#' @param outliers a vector of outlier sample(s) (bid) to remove prior to standardization
#' @param center_f a vector of outlier sample(s) to exclude from feature centering
#' @param scale_f a vector of outlier sample(s) to exclude from feature scaling
#' @return the normalized matrix
NormalizeData <- function(tissue_mat, strategy = 'autoscale', ref_sams2, outliers = NA, center_f = NA, scale_f = NA){
  outliers = outliers[!is.na(outliers)] # remove NA from outlier list 
  outliers = outliers[outliers != '']
  if(strategy == 'raw'){
    freeze3 <- tissue_mat
  }else if(strategy == 'log2'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    # Feature standardization
    # Output data storage
    freeze3 <- log2(tissue_mat1)
  }else if(strategy == 'autoscale'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        exp50 <- exp50[!exp50 %in% outlier]
      }
    }
    # Feature standardization
    if(!is.na(center_f)){
      x <- paste0("^",center_f)
      cen50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.center.tissue <- apply(tissue_mat1[cen50,], 2, median, na.rm = T)
    }else{
      f.center.tissue <- apply(tissue_mat1[exp50,], 2, median, na.rm = T)
    }
    if(!is.na(scale_f)){
      x <- paste0("^",scale_f)
      scl50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.scale.tissue <- apply(tissue_mat1[scl50,], 2, sd, na.rm = T)
    }else{
      f.scale.tissue <- apply(tissue_mat1[exp50,], 2, sd, na.rm = T)
    }
    
    for(i in 1:dim(tissue_mat1)[2]){
      tissue_mat1[,i] <- (tissue_mat1[,i]-f.center.tissue[i])/f.scale.tissue[i]
    }
    # Output data storage
    freeze3 <- tissue_mat1
  }else if(strategy == 'log2_fstd'){
    # log2 + feature-standardization
    ################################################################################
    # Reload the original gastroc mat (to visualize sample medians)
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        exp50 <- exp50[!exp50 %in% outlier]
      }
    }
    tissue_mat3 <- tissue_mat2 <- log2(tissue_mat1)
    # Feature standardization
    if(!is.na(center_f)){
      x <- paste0("^",center_f)
      cen50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.center.tissue <- apply(tissue_mat2[cen50,], 2, median, na.rm = T)
    }else{
      f.center.tissue <- apply(tissue_mat2[exp50,], 2, median, na.rm = T)
    }
    if(!is.na(scale_f)){
      x <- paste0("^",scale_f)
      scl50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.scale.tissue <- apply(tissue_mat2[scl50,], 2, sd, na.rm = T)
    }else{
      f.scale.tissue <- apply(tissue_mat2[exp50,], 2, sd, na.rm = T)
    }
    #f.center.tissue <- apply(tissue_mat2[exp50,], 2, median, na.rm = T)
    #f.scale.tissue <- apply(tissue_mat2[exp50,], 2, sd, na.rm = T)
    for(i in 1:dim(tissue_mat3)[2]){
      tissue_mat3[,i] <- (tissue_mat3[,i]-f.center.tissue[i])/f.scale.tissue[i]
    }
    freeze3 <- tissue_mat3
  }else if(strategy == 'log2_fcen'){
    # log2 + feature-standardization
    ################################################################################
    # Reload the original gastroc mat (to visualize sample medians)
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        exp50 <- exp50[!exp50 %in% outlier]
      }
    }
    tissue_mat3 <- tissue_mat2 <- log2(tissue_mat1)
    # Feature center
    if(!is.na(center_f)){
      x <- paste0("^",center_f)
      cen50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.center.tissue <- apply(tissue_mat2[cen50,], 2, median, na.rm = T)
    }else{
      f.center.tissue <- apply(tissue_mat2[exp50,], 2, median, na.rm = T)
    }
    for(i in 1:dim(tissue_mat3)[2]){
      tissue_mat3[,i] <- (tissue_mat3[,i]-f.center.tissue[i])
    }
    freeze3 <- tissue_mat3
  }else if(strategy == 'log2_scen'){
    # log2 + sample centering
    ################################################################################
    # Reload the original gastroc mat (to visualize sample medians)
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        exp50 <- exp50[!exp50 %in% outlier]
      }
    }
    tissue_mat4 <- tissue_mat2 <- log2(tissue_mat1)
    # Sample standardization (only test samples)
    s.center.tissue <- apply(tissue_mat2, 1, median, na.rm = T)
    for(i in 1:dim(tissue_mat4)[1]){
      tissue_mat4[i,] <- (tissue_mat4[i,]-s.center.tissue[i])
    }
    freeze3 <- tissue_mat4
  }else if(strategy == 'log2_sstd'){
    # log2 + feature-standardization
    ################################################################################
    # Reload the original gastroc mat (to visualize sample medians)
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    tissue_mat2 <- log2(tissue_mat1)
    tissue_mat4 <- tissue_mat2
    s.center.tissue <- apply(tissue_mat4, 1, median, na.rm = T)
    s.scale.tissue <- apply(tissue_mat4, 1, sd, na.rm = T)
    for(i in 1:dim(tissue_mat4)[1]){
      tissue_mat4[i,] <- (tissue_mat4[i,]-s.center.tissue[i])/s.scale.tissue[i]
    }
    freeze3 <- tissue_mat4
  }else if(strategy == 'log2_fstd_scen'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        exp50 <- exp50[!exp50 %in% outlier]
      }
    }
    tissue_mat3 <- tissue_mat2 <- log2(tissue_mat1)
    # Feature standardization
    if(!is.na(center_f)){
      x <- paste0("^",center_f)
      cen50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.center.tissue <- apply(tissue_mat2[cen50,], 2, median, na.rm = T)
    }else{
      f.center.tissue <- apply(tissue_mat2[exp50,], 2, median, na.rm = T)
    }
    if(!is.na(scale_f)){
      x <- paste0("^",scale_f)
      scl50 <- exp50[!grepl(paste(x, collapse = "|"), exp50)]
      f.scale.tissue <- apply(tissue_mat2[scl50,], 2, sd, na.rm = T)
    }else{
      f.scale.tissue <- apply(tissue_mat2[exp50,], 2, sd, na.rm = T)
    }
    #f.center.tissue <- apply(tissue_mat2[exp50,], 2, median, na.rm = T)
    #f.scale.tissue <- apply(tissue_mat2[exp50,], 2, sd, na.rm = T)
    for(i in 1:dim(tissue_mat3)[2]){
      tissue_mat3[,i] <- (tissue_mat3[,i]-f.center.tissue[i])/f.scale.tissue[i]
    }
    # Sample standardization (only test samples)
    tissue_mat4 <- tissue_mat3
    s.center.tissue <- apply(tissue_mat3, 1, median, na.rm = T)
    for(i in 1:dim(tissue_mat4)[1]){
      tissue_mat4[i,] <- (tissue_mat4[i,]-s.center.tissue[i])
    }
    freeze3 <- tissue_mat4
  }else if(strategy == 'log2_fstd_sstd'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        rm_out <- exp50[grepl(paste0("^",as.character(outlier)),exp50)]
        rm_out <- match(rm_out, exp50)
        exp50 <- exp50[-rm_out]
      }
    }
    tissue_mat2 <- log2(tissue_mat1)
    # Feature standardization
    tissue_mat3 <- tissue_mat2
    f.center.tissue <- apply(tissue_mat2[exp50,], 2, median, na.rm = T)
    f.scale.tissue <- apply(tissue_mat2[exp50,], 2, sd, na.rm = T)
    for(i in 1:dim(tissue_mat3)[2]){
      tissue_mat3[,i] <- (tissue_mat3[,i]-f.center.tissue[i])/f.scale.tissue[i]
    }
    # Sample standardization (only test samples)
    tissue_mat4 <- tissue_mat3
    s.center.tissue <- apply(tissue_mat3, 1, median, na.rm = T)
    s.scale.tissue <- apply(tissue_mat3, 1, sd, na.rm = T)
    for(i in 1:dim(tissue_mat4)[1]){
      tissue_mat4[i,] <- (tissue_mat4[i,]-s.center.tissue[i])/s.scale.tissue[i]
    }
    freeze3 <- tissue_mat4
  }else if(strategy == 'log2_sstd_fstd'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    if(length(outliers)>0){
      print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers, collapse = ' '))))
      for(outlier in str_split(outliers, pattern = ';', simplify = T)){
        rm_out <- exp50[grepl(paste0("^",as.character(outlier)),exp50)]
        rm_out <- match(rm_out, exp50)
        exp50 <- exp50[-rm_out]
      }
    }
    tissue_mat2 <- log2(tissue_mat1)
    # Sample standardization
    tissue_mat5 <- tissue_mat2
    s.center.tissue <- apply(tissue_mat5, 1, median, na.rm = T)
    s.scale.tissue <- apply(tissue_mat5, 1, sd, na.rm = T)
    for(i in 1:dim(tissue_mat5)[1]){
      tissue_mat5[i,] <- (tissue_mat5[i,]-s.center.tissue[i])/s.scale.tissue[i]
    }
    
    # Feature standardization (only 40 samples)
    f.center.tissue <- apply(tissue_mat5[exp50,], 2, median, na.rm = T)
    f.scale.tissue <- apply(tissue_mat5[exp50,], 2, sd, na.rm = T)
    for(i in 1:dim(tissue_mat5)[2]){
      tissue_mat5[,i] <- (tissue_mat5[,i]-f.center.tissue[i])/f.scale.tissue[i]
    }
    freeze3 <- tissue_mat5
  }else if(strategy == 'log2_sstd_fcen'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    tissue_mat2 <- log2(tissue_mat1)
    # Sample standardization
    tissue_mat5 <- tissue_mat2
    s.center.tissue <- apply(tissue_mat5, 1, median, na.rm = T)
    s.scale.tissue <- apply(tissue_mat5, 1, sd, na.rm = T)
    for(i in 1:dim(tissue_mat5)[1]){
      tissue_mat5[i,] <- (tissue_mat5[i,]-s.center.tissue[i])/s.scale.tissue[i]
    }
    
    # Feature standardization (only 40 samples)
    f.center.tissue <- apply(tissue_mat5[exp50,], 2, median, na.rm = T)
    for(i in 1:dim(tissue_mat5)[2]){
      tissue_mat5[,i] <- (tissue_mat5[,i]-f.center.tissue[i])
    }
    freeze3 <- tissue_mat5
  }else if(strategy == 'log2_sstd_fstd_scen'){
    tissue_mat1 <- tissue_mat
    exp50 <- row.names(tissue_mat1)[!row.names(tissue_mat1) %in% ref_sams2]
    tissue_mat2 <- log2(tissue_mat1)
    # Sample standardization
    tissue_mat5 <- tissue_mat2
    s.center.tissue <- apply(tissue_mat5, 1, median, na.rm = T)
    s.scale.tissue <- apply(tissue_mat5, 1, sd, na.rm = T)
    for(i in 1:dim(tissue_mat5)[1]){
      tissue_mat5[i,] <- (tissue_mat5[i,]-s.center.tissue[i])/s.scale.tissue[i]
    }
    # Feature standardization (only test samples)
    f.center.tissue <- apply(tissue_mat5[exp50,], 2, median, na.rm = T)
    f.scale.tissue <- apply(tissue_mat5[exp50,], 2, sd, na.rm = T)
    for(i in 1:dim(tissue_mat5)[2]){
      tissue_mat5[,i] <- (tissue_mat5[,i]-f.center.tissue[i])/f.scale.tissue[i]
    }
    s.center.tissue <- apply(tissue_mat5, 1, median, na.rm = T)
    for(i in 1:dim(tissue_mat5)[1]){
      tissue_mat5[i,] <- (tissue_mat5[i,]-s.center.tissue[i])
    }
    freeze3 <- tissue_mat5
  }else{
    accepted_strategies = c('raw', 'log2', 'autoscale', 'log2_fstd', 'log2_fcen', 'log2_sstd', 'log2_scen', 
                            'log2_fstd_scen', 'log2_fstd_sstd', 'log2_sstd_fstd',  'log2_sstd_fcen', 'log2_sstd_fstd_scen')
    error(sprintf("Supplied value '%s' for argument 'strategy' is not one of the accepted values: %s", strategy, paste(accepted_strategies, collapse = ', ')))
  }
  # Output Matrix
  return(freeze3)
}

#' RemoveOutliers
#' A function to remove outlier samples
#' @param mat the count (abundance) matrix
#' @param outliers the outlier sample(s) delimited by ";"
#' @return a count matrix with outlier samples removed
RemoveOutliers <- function(mat, outliers){
  print(paste0('Outlier/flagged sample(s) removed from feature-std: ', as.character(paste(outliers[!is.na(outliers)], collapse = ' '))))
  if(!is.na(outliers)){
    for(outlier in str_split(outliers, pattern = ';', simplify = T)){
      rm_out <- row.names(mat)[grepl(paste0("^",as.character(outlier)),row.names(mat))]
      rm_out <- match(rm_out, row.names(mat))
      mat <- mat[-rm_out,]
    }
  }
  return(mat)
}

#' SaveFreezesOutliers
#' A function that acts as a wrapper around the NormalizeData function and saves the normalized data
#' @param tissue_mat a numeric matrix
#' @param ref_sams2 a vector of reference samples (used to separate reference samples from experimental samples)
#' @param Named named or unnamed
#' @param site_org the site of measurement
#' @param Platform the platform the metabolites were measured on
#' @param Targeted targeted or untargeted
#' @param tissue_name the BIC assigned name of the tissue
#' @param WD the "working directory" that acts as a base directory for analysis (not actual working directory from getwd())
#' @param save_data a boolean whether to save normalized data to disc
#' @return the normalized matrix
SaveFreezesOutliers <-  function(mat, ref_sams2,Named,site_org,Platform,Targeted,tissue, WD, save_data = T){
  source(paste0(WD,'/functions/NormalizeData.R'))
  freeze1 <- NormalizeData(mat, strategy = 'raw', ref_sams2)
  freeze2b <- NormalizeData(mat, strategy = 'log2_sstd', ref_sams2)
  freeze3a <- NormalizeData(mat, strategy = 'autoscale', ref_sams2)
  freeze3b <- NormalizeData(mat, strategy = 'log2_fstd', ref_sams2)
  freeze3c <- NormalizeData(mat, strategy = 'log2_fstd_sstd', ref_sams2)
  freeze3c2 <- NormalizeData(mat, strategy = 'log2_fstd_scen', ref_sams2)
  freeze3d <- NormalizeData(mat, strategy = 'log2_sstd_fstd', ref_sams2)
  freeze3d2 <- NormalizeData(mat, strategy = 'log2_sstd_fstd_scen', ref_sams2)
  if(save_data){
    # Save the freezes
    saveRDS(freeze1, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                   '-',tissue,'-',Targeted,'-freeze1-outliers_steep.RDS'))
    saveRDS(freeze2b, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                    '-',tissue,'-',Targeted,'-freeze2b-outliers_steep.RDS'))
    saveRDS(freeze3a, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                    '-',tissue,'-',Targeted,'-freeze3a-outliers_steep.RDS'))
    saveRDS(freeze3b, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                    '-',tissue,'-',Targeted,'-freeze3b-outliers_steep.RDS'))
    saveRDS(freeze3c, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                    '-',tissue,'-',Targeted,'-freeze3c-outliers_steep.RDS'))
    saveRDS(freeze3c2, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                     '-',tissue,'-',Targeted,'-freeze3c2-outliers_steep.RDS'))
    saveRDS(freeze3d, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                    '-',tissue,'-',Targeted,'-freeze3d-outliers_steep.RDS'))
    saveRDS(freeze3d2, file = paste0(WD,'/freezes/test-ref-',Named,'-',site_org,'-',Platform,
                                     '-',tissue,'-',Targeted,'-freeze3d2-outliers_steep.RDS'))
  }
}

#' Viallabel2Bid
#' A function to convert sample identifiers, viallabel to bid, in data freeze objects
#' @param mat a numeric input matrix
#' @param sample_df a custom dataframe that contains sample identifiers and sample order information
#' @param tissue the BIC assigned name of the tissue
#' @return the matrix with converted sample identifiers in rownames
Viallabel2Bid <- function(mat, sample_df, tissue){
  df <- data.frame(TISSUE_viallabel = paste0(tissue,':',row.names(mat))) %>%
    left_join(y = sample_df)
  if(all(df$viallabel == row.names(mat))){
    idx <- which(!is.na(df$bid))
    row.names(mat)[idx] <- df$bid[idx]
    for(i in 1:length(row.names(mat))){
      if(substr(row.names(mat)[i], 1, 2) == '90'){
        row.names(mat)[i] <- substr(row.names(mat)[i], 1, 5)
      }
    }
    
  }else{
    stop('Conversion of row names is inconsistant')
  }
  return(mat)
}

#' Zero2NA
#' A function to convert all zero values in a numeric matrix to NA values
#' @param mat a numeric matrix
#' @param zero2na a boolean whether to perform conversion
#' @return the numeric matrix with zero values converted to missing values
Zero2NA <- function(mat, zero2na){
  if(zero2na){
    mat[mat == 0] <- NA
  }else{
    mat <- mat
  }
  return(mat)
}

#' Neg2NA
#' A function to convert all negative values in a numeric matrix to NA values
#' @param mat a numeric matrix
#' @param neg2na a boolean whether to perform conversion
#' @return the numeric matrix with negative values converted to missing values
Neg2NA <- function(mat, neg2na){
  if(neg2na){
    mat[mat < 0] <- NA
  }else{
    mat <- mat
  }
  return(mat)
}

#' NA2Min
#' A function to convert all NA values in a numeric matrix to specificied minimum values
#' @param mat a numeric matrix
#' @param na2min a boolean whether to perform conversion
#' @param impute_min the strategy for imputing with a minimum (e.g. half-min)
#' @return the numeric matrix with NA values converted to specificied minimum values
NA2Min <- function(mat, na2min, impute_min){
  if(na2min){
    if(impute_min == 'half-min'){
      for (i in 1:dim(mat)[2]) {
        mat[is.na(mat[,i]),i] <- min(mat[,i], na.rm = T) - log(2,10)
      }
    }else if(impute_min == 'min'){
      for (i in 1:dim(mat)[2]) {
        mat[is.na(mat[,i]),i] <- min(mat[,i], na.rm = T)
      }
    }
  }else{
    mat <- mat
  }
  return(mat)
}

#' SampleStandardizationDecision
#' Perform Kruskal-Wallis tests to decide whether or not a data set should be sample-standardized
#' @param norm_data BID x feature data.frame with 3b normalization
#' @param pheno_df DMAQC phenotypic data
#' @param label dataset label
#' @return a data.table with details of the test results 
SampleStandardizationDecision <- function(norm_data, pheno_df, label){
  require("testit")
  require("data.table")
  
  # add sex and group columns 
  pheno_df$sex = factor(pheno_df$Registration.sex)
  pheno_df$group = factor(paste0(pheno_df$Key.intervention, pheno_df$Key.sacrificetime))
  
  res_list = list()
  
  curr_meta = unique(data.table(pheno_df[as.character(pheno_df$viallabel) %in% rownames(norm_data),c("viallabel","bid","group","sex")]))
  curr_data = data.frame(norm_data, check.names=F)
  curr_data = curr_data[as.character(curr_meta[,viallabel]),]
  assert(all(rownames(curr_data) == curr_meta[,viallabel]))
  
  test_ss = data.table(bid = rownames(norm_data),
                       samp.med = apply(norm_data, 1, median, na.rm=T),
                       quant.75 = apply(norm_data, 1, function(x) quantile(x, 0.75, na.rm=T)))
  
  curr_meta = data.table(cbind(curr_meta, test_ss))
  
  # test sex
  if(length(unique(curr_meta[,sex]))>1){
    sex_pval_1 = kruskal.test(samp.med ~ sex, data=curr_meta)$p.value
    sex_pval_2 = kruskal.test(quant.75 ~ sex, data=curr_meta)$p.value
  }else{
    sex_pval_1 = NA_real_
    sex_pval_2 = NA_real_
  }
  res_list[['sex']] = data.table(dataset=label,
                                 group='sex',
                                 outcome=c('samp.med','quant.75'),
                                 p.value=c(sex_pval_1, sex_pval_2))
  
  # test group
  for(s in unique(curr_meta[,sex])){
    cmet = curr_meta[sex==s]
    group_pval_1 = kruskal.test(samp.med ~ group, data=cmet)$p.value
    group_pval_2 = kruskal.test(quant.75 ~ group, data=cmet)$p.value
    res_list[[s]] = data.table(dataset=label,
                               group=sprintf("group,%s",s),
                               outcome=c('samp.med','quant.75'),
                               p.value=c(group_pval_1, group_pval_2))
  }
  
  return(rbindlist(res_list))
}

#' WriteTData
#' A function to transpose and write a matrix to a file with original column 
#' names as the first column
#' @param mat a numeric matrix
#' @param outfile name of file to write to
#' @param indexname column name for feature column
WriteTData <- function(mat, outfile, indexname="metabolite_name"){
  cn <- colnames(mat)
  tdf <- cbind(data.frame('V1'=cn), as.data.frame(t(mat), check.names=F))
  colnames(tdf)[1] = indexname
  write.table(tdf, file = outfile, quote = F, col.names = T, row.names = F, sep = '\t')
}

#' ZeroFeatureFilter
#' A function to remove features based on the proportion of Zero values across features
#' @param mat a numeric matrix
#' @param threshold_freq the threshold frequency of Zero values used to filter (remove) features with Zero values above the threshold
#' @param tissue_name the name of the tissue
#' @param ref_sams2 a vector of reference samples (used to separate reference samples from experimental samples)
#' @return the matrix with features removed
ZeroFeatureFilter <- function(mat, threshold_freq = 1.0, tissue_name,ref_sams2){
  # CountZeros function
  CountZeros <- function(x){
    sum(x==0, na.rm = T)
  }
  # Remove reference samples
  mat2 <- mat[!row.names(mat) %in% ref_sams2,]
  # NA frequency
  fzero.mat <- apply(mat2, 2, CountZeros) %>% sort() %>% rev()
  fzero.mat_freq <- fzero.mat/nrow(mat2)
  # Remove samples
  frm_zero <- names(fzero.mat_freq[fzero.mat_freq > threshold_freq])
  print(paste0(length(frm_zero),' Features removed from ',tissue_name,' due to zero values:'))
  print(paste0(frm_zero, collapse=', '))
  # Remove the samples from the individual matrices
  mat <- mat[,!colnames(mat) %in% frm_zero]
  # Output
  return(mat)
}

#' TestOnly
#' A function to remove reference samples
#' @param mat a numeric matrix
#' @param 
#' @param ref_sams2 a vector of reference samples (used to separate reference samples from experimental samples)
#' @return the matrix with reference samples removed removed
TestOnly <- function(mat, ref_sams2){
  exp50 <- row.names(mat)[!row.names(mat) %in% ref_sams2]
  return(mat[exp50,])
}

