load(file = "rna-seq_meta/rnaseq_meta.RData")
dim(rnaseq_meta)
unique(rnaseq_meta$Tissue)
rnaseq_meta=subset(rnaseq_meta,Tissue %in% c("gastrocnemius","white adipose","liver","Gastrocnemius Powder","PaxGene RNA","Heart Powder"))
dim(rnaseq_meta)
#replace machine error RIN score value with NA
#Update the flag matrix for that sample from TRUE to FALSE
rnaseq_meta[rnaseq_meta$vial_label == "90140017003", "RIN"] <- NA
#rnaseq_meta[rnaseq_meta$vial_label == "90140017003", "IsFlagged"] <- FALSE
save.image(file ="~/work/repo/motrpac-bic-norm-qc/rna-seq_meta/rnaseq_meta_external_release1.RData")
