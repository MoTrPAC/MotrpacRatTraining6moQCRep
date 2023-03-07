#' Matches PTM and proteome by accession number, and tries all possible PTM accession numbers if primary
#' accession number doesn't have proteome match. If there are multiple protein matches, it will pick
#' the one with least missing values using the proteomics quantitative data. 
#' 
#' @param pr_data Proteomics quantitative data in matrix form
#' @param ptm_row_annot PTM row annotations with at least the following columns (protein_id, ptm_id, redundant_id)
#' @param prot_row_annot Protein row annotations with at least the following columns (protein_id, redundant_id)
#' 
#' @description Takes a parsed proteomics dataset object containing both PTM and prot-pr datasets.
#' Assumes that the input is a list object as generated in the proteomics bic_qc_report
#' 
#' @return ptm_row_annotation with an additional column "prot_pr_id_match" indicating the matching identifier in the prot_row_annot table.

ptm_protein_match <- function( pr_data, ptm_row_annot, prot_row_annot){
  
  #Load libraries
  required_libs = c("tidyverse")
  for (lib_name in required_libs){
    tryCatch({library(lib_name,character.only = T)}, error = function(e) {
      print(paste("Cannot load",lib_name,", please install"))
    })
  }

  #convert to data frames for convenience
  ptm_row_annot <- ptm_row_annot %>% as.data.frame() 
  prot_row_annot <- prot_row_annot %>% as.data.frame() 
  
  #Create an additional column with the matching protein_id between prot-pr and PTM
  ptm_row_annot$prot_pr_id_match <- character(length = nrow(ptm_row_annot))
  
  #Expand redundant ids from proteome row annotations
  prot_row_annot.expanded <- prot_row_annot %>% 
    separate_rows(redundant_ids,sep="\\|") %>%
    mutate(redundant_ids = str_extract(redundant_ids,"\\w+")) %>%
    as.data.frame()
  
  #---Match ids across PTMs and proteome -----------------
  #
  for(row in c(1:nrow(ptm_row_annot))){
    
    match <- which(prot_row_annot[,"protein_id"] == ptm_row_annot[,"protein_id"][row])
    
    #If a match is found then save protein id and exit
    if(length(match) > 0){
      ptm_row_annot$prot_pr_id_match[row] <- prot_row_annot[match,"protein_id"]
    } else {
      #First try to match PTM redundant_ids to protein primary ID
      accession.numbers <- unlist(strsplit(ptm_row_annot[,"redundant_ids"][row],'|', fixed = TRUE)) %>%
        str_extract("\\w+")
      
      matches <- prot_row_annot[which(prot_row_annot$protein_id %in% accession.numbers),]
      
      if(nrow(matches) >= 1){
        #If there are multiple matches, find the protein with least missing values
        best <- matches %>% arrange(protein_score) %>% .$protein_id %>% .[1]
        ptm_row_annot$prot_pr_id_match[row] <- best
        
        
      } else {
        #If the above fails, try to match PTM redundant_ids to all Proteome redundant_ids
        matches <- prot_row_annot.expanded[which(prot_row_annot.expanded$redundant_ids %in% accession.numbers),]
        if(nrow(matches) >= 1){
          best <- matches %>% arrange(protein_score) %>% .$protein_id %>% .[1]
          ptm_row_annot$prot_pr_id_match[row] <- best
        } else{
          ptm_row_annot$prot_pr_id_match[row] <- NA
        }
      }
    }
  }
  
  return(as.matrix(ptm_row_annot))
  
}



#' Normalizes PTM to protein level by fitting a global linear model and returning residuals
#' Fits PTM = beta_0 + beta_1*protein to all matched points in a dataset, and returns residuals as
#' protein-corrected PTM values.
#' 
#' @param ptm_data PTM quantitative data in matrix form
#' @param pr_data Proteomics quantitative data in matrix form
#' @param ptm_row_annot PTM row annotations table with at least the following columns (protein_id, ptm_id, redundant_id, prot_pr_id_match)
#' @param prot_row_annot Protein row annotations with at least the following columns (protein_id, redundant_id)
#' @param plot_graphs Logical indicating whether to plot graphs showing PTM-Prot correlation
#' 
#' @description Takes a parsed proteomics dataset object containing both PTM and prot-pr datasets.
#' Assumes that the input is a list object as generated in the proteomics bic_qc_report
#' 
#' @return ptm_data table with protein-corrected quantitative values
ptm_protein_correction <- function(ptm_data, pr_data, ptm_row_annot, prot_row_annot, plot_graphs = T){
  
  #Load libraries
  required_libs = c("tidyverse")
  for (lib_name in required_libs){
    tryCatch({library(lib_name,character.only = T)}, error = function(e) {
      print(paste("Cannot load",lib_name,", please install"))
    })
  }
  
  #convert to data frames for convenience
  ptm_row_annot <- ptm_row_annot %>% as.data.frame() 
  prot_row_annot <- prot_row_annot %>% as.data.frame() 
  
  #warn if not all samples are matched
  matched.samples <- intersect(colnames(pr_data), colnames(ptm_data))
  if(length(matched.samples) != length(colnames(ptm_data))){
    print('WARNING: not all samples in PTM file have matches in proteome.  Unmatched samples removed.')
  }
  
  #create merged data table
  ptm.wide <- full_join(ptm_row_annot,
                        ptm_data %>% as.data.frame() %>% rownames_to_column(var="ptm_id"), by = "ptm_id")
  
  
  ptm.melt <- pivot_longer(ptm.wide,starts_with("9"),names_to= "sample")
  prot.melt <- pivot_longer(pr_data %>% as.data.frame %>% rownames_to_column(var="protein_id_prot_pr"),
                            -protein_id_prot_pr,names_to= "sample")
  
  data <- left_join(ptm.melt, prot.melt, by = c("prot_pr_id_match" = "protein_id_prot_pr","sample" = "sample"),
                    suffix=c(".ptm",".prot"))
  
  #print metrics
  percent <-round(100*sum(!is.na(data$value.prot))/nrow(ptm.melt), digits = 1)
  print(paste(sum(!is.na(data$value.prot)), ' points with proteome match out of ', nrow(ptm.melt),
              ' (', percent, '%).', sep = ''))
  
  #Fit global model
  print("Fitting model...") 
  model <- lm(value.ptm ~ value.prot, data = data,na.action = na.exclude)
  residuals <- residuals(model)
  results <- cbind(data,residuals)
  print("Success.")
  print(summary(model))
  
  #Plot graphs showing PTM values before and after correction
  if(plot_graphs){
    coord_limits <- max(abs(c(data$value.prot,data$value.ptm,results$residuals)),na.rm = T) - 1
    
    g1 <- ggplot(data[sample(1:nrow(data),50000,replace=T),],
                 aes(y= value.ptm, x = value.prot), na.rm=TRUE)+
      geom_point(color="black",alpha=0.05, na.rm = TRUE)+
      geom_smooth(formula = y ~ x, method="lm", na.rm = TRUE)+
      theme_bw()+
      labs(y="PTM Log2 TMT ratio",x = "Prot Log2 TMT ratio",title= paste(tissue,ptm_assay,sep=","),subtitle = "Non-corrected")+
      coord_fixed(xlim = c(-coord_limits,coord_limits),ylim = c(-coord_limits,coord_limits))
    
    g2 <- ggplot(results[sample(1:nrow(results),50000,replace=T),],
                 aes(y= residuals, x = value.prot), na.rm=TRUE)+
      geom_point(color="black",alpha=0.05, na.rm = TRUE)+
      geom_smooth(formula = y ~ x, method="lm", na.rm = TRUE)+
      theme_bw()+
      labs(y="PTM Log2 TMT ratio",x = "Prot Log2 TMT ratio",title= paste(tissue,ptm_assay,sep=","),subtitle = "Protein-corrected PTM")+
      coord_fixed(xlim = c(-coord_limits,coord_limits),ylim = c(-coord_limits,coord_limits))
    
    plot(g1)
    plot(g2)
  }
  
  #Save results
  results.wide <- pivot_wider(results,ptm_id,names_from = sample, values_from = residuals) %>%
    as.data.frame() %>% column_to_rownames("ptm_id") %>% as.matrix()
  return(results.wide)
  
}
