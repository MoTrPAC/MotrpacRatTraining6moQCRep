
# Two different packages to work with RRBS bismark_results and perform differential methylation analysis, 
# edgeR is very slow , bsseq is multi-threaded Requires R 3.6 or higher

tryCatch(library(edgeR),error=function(){print("cannot run: please install edgeR");q("no")})
tryCatch(library(optparse),error=function(){print("cannot run: please install optparse");q("no")})

option_list = list(
  make_option(c("-p", "--path"), action="store", default=getwd(), type='character',
              help="a path to a directory with the bismark files, defualt: getwd()"),
  make_option(c("-m", "--meta"), action="store", default=NA, type='character',
              help="a path to a file with the RRBS metadata: only bismark files that have vial labels here will be prased, default: NA, which means do not filter bismark files"),
  make_option(c("-o", "--output"), action="store", default="rrbs_data", type='character',
              help="the name of the output RData file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt = parse_args(parser)

print("Parsing bismark files, input parameters are:")
print(paste("directory with bismark files:", opt$p))
print(paste("metadata file path (or NA if not given):",opt$m))
opt$o = paste0(opt$o,".RData")
print(paste("output RData file:",opt$o))


files = list.files(opt$p,full.names = T,recursive = T)
files = files[grepl("_bismark.cov.gz$",perl=T,files)]
print(paste0("scanned for bismark files in input dir, found:",length(files)))
print("first two (for example):")
print(files[1:2])

extract_vial_from_path<-function(p){
  arr = strsplit(p,split='/')[[1]]
  f = arr[length(arr)]
  f = gsub("_bismark.cov.gz","",f)
  return(f)
}

print("Scanning the metadata file for vial labels")
file2sample_name = sapply(files,extract_vial_from_path)
if(!is.na(opt$m)){
  meta=read.csv(opt$m,header=TRUE,sep=",",quote = "\"",dec = ".", fill = TRUE)
  for(vial in meta$vial_label){
    reg = paste0(vial,"_bismark.cov.gz")
    f = files[grepl(reg,files)]
    if(length(f)!=1){next}
    file2sample_name[f] = vial
  }
}
# filter the files object, if meta is NA (no meta file), then this
# line has no effect
files = names(file2sample_name)

print("Done extracting vial names (and filtering by metadata if given)")
print("first two (for example):")
print(file2sample_name[1:2])

print("Parsing the bismark files, this may take a while:")
rrbs_data = readBismark2DGE(files, sample.names= as.character(unname(file2sample_name)))

print("Done!")
print("Saving the RData file")
save(rrbs_data,file = opt$o)

