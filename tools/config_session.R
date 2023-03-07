
# Load packages for visualization and data manipulation
required_libs = c(
  "corrplot","gplots","ggcorrplot","ggplot2","data.table"
)
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}


scatterplot_colors =  scale_colour_manual(
  values=c(liver='#FF00FF',gastroc_stanford='#00cd00',
           gastroc_mssm='#00ee00',kidney='#a020f0',
           'brown adipose'='#ffa500','white adipose'='#1e90ff',
           'paxgene rna'='#ff0000',heart='#40e0d0'))





