load_genExp_data <- function(genExp_path, meta_path){
  
  GE_file_ext <- tools::file_ext(genExp_path)
  if (GE_file_ext %in% c("csv", "txt", "tsv")) {
    # Read data from .csv or .txt file
    genExp <- read.table(genExp_path, header = TRUE, sep = "\t")
  } else if (GE_file_ext == "rds"| GE_file_ext == "RDS") {
    # Read data from .rds file
    genExp <- readRDS(genExp_path)
   # if (class(genExp) %in% c("RangedSummarizedExperiment", "SummarizedExperiment")){
   # genExp <- assay(genExp)
   # }
  } else {
    # Unsupported file type
    stop("Unsupported file type. Only .csv, .txt, and .rds are supported.")
  }
  
  meta_file_ext <- tools::file_ext(meta_path)
  if (meta_file_ext %in% c("csv", "txt", "tsv")) {
    # Read data from .csv or .txt fil
    meta <- read.table(meta_path, header = TRUE, sep = "\t")
  } else if (meta_file_ext %in% c("rds", "RDS")) {
    # Read data from .rds file
    meta <- readRDS(meta_path)
  } else {
    # Unsupported file type
    stop("Unsupported file type. Only .csv, .txt, and .rds are supported.")
  }
  out <- list(raw_data = genExp, raw_meta = meta)
  return(out)
}

#function to fix ID format between different files
id_format_fix <- function(v) {
  v <- gsub("^(\\d)", "X\\1", v)
  out <- gsub("-", ".", v)
  return(out)
}

get_gene_info <- function(genid){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("biomaRt")
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "genes", 
                        dataset = "hsapiens_gene_ensembl")
  bm_gene = getBM(attributes = c('ensembl_gene_id',
                                 'chromosome_name',
                                 'start_position',
                                 'end_position'),
                  filters = c('ensembl_gene_id'),
                  mart = ensembl,
                  values = genid)
  colnames(bm_gene) = c('genid','#Chr','start','end')
  return(bm_gene)
}
