load_data <- function(names = c("data"), ...){
  # Determine the number of arguments
  num_files <- length(list(...))
  num_names <- length(names)
  if (num_files != num_names){
    stop(paste0("You provided ", num_args, " files but only ", num_names, " names. There should be a name matching each file to be loaded"))
  }
  # Initialize an empty list to store the results
  out <- vector("list", length = num_files)
  names(out) <- names
  # Iterate through each argument
  for (i in seq_along(list(...))) {
    #get file type
    file_ext <- tools::file_ext(list(...)[[i]])
    if (file_ext %in% c("csv", "txt", "tsv")) {
      # Read data from .csv or .txt file
      table <- read.table(list(...)[[i]], header = TRUE, sep = "\t")
    } else {
      # Unsupported file type
      stop(paste0(list(...)[[i]]," is an unsupported file type. Only .csv, .txt are supported."))
    }
    # Store the result in the list
    out[[i]] <- table
  }
  return(out)
}
  
  
#  tpm_file_ext <- tools::file_ext(genExp_counts_path)
#   if (tpm_file_ext %in% c("csv", "txt", "tsv")) {
#     # Read data from .csv or .txt file
#     tpm <- read.table(tpm_path, header = TRUE, sep = "\t")
#   } else {
#     # Unsupported file type
#     stop("Unsupported file type. Only .csv, .txt are supported.")
#   }
#   meta_file_ext <- tools::file_ext(meta_path)
#   if (meta_file_ext %in% c("csv", "txt", "tsv")) {
#     # Read data from .csv or .txt fil
#     meta <- read.table(meta_path, header = TRUE, sep = "\t")
#   } else if (meta_file_ext %in% c("rds", "RDS")) {
#     # Read data from .rds file
#     meta <- readRDS(meta_path)
#   } else {
#     # Unsupported file type
#     stop("Unsupported file type. Only .csv, .txt, and .rds are supported.")
#   }
#   out <- list(counts = counts, tpm = tpm)
#   return(out)
# }

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
  colnames(bm_gene) = c('genid','Chr','start','end')
  return(bm_gene)
}

#prep data
prep_datExp <- function(.data, gtf){
  gtf <- gtf %>% filter(!(Chr %in% c("chrX","chrY","chrM")))
  .data <- .data %>% filter(genid %in% gtf$genid)
  gtf <- gtf %>% filter(genid %in% .data$genid)
  rownames(.data) <- .data$genid
  .data <- .data[,!colnames(.data) %in% "genid"]
  .data <- .data[,-1]
  .data <- round(.data)
  return(.data)
}

#function to plot gene density
plot_gen_density <- function(.data, offset = 0.1, title_str = 'Scaled Gene Counts ', xlim = c(-10,30), ylim = c(0,0.5), log = TRUE){
  # View the distribution of expression for each sample.
  # box plot, looking for big differences in read depth (raw counts), symmetry in distribution across samples
  par(mfrow=c(1,2))
  boxplot(.data, range = 0, main = paste(title_str), xlab = 'Samples', xaxt = "n")
  boxplot(log2(offset+.data), range = 0, main = paste('log2(counts+',offset,')'), xlab = 'Samples', xaxt = "n")
  # Histogram/density plot
  # Look for: how well do the distributions line up, outlier samples, zero counts
  par(mfrow=c(1,1))
  i <- 1
  if(log == TRUE){
    plot <- plot(density(log2(offset+.data[,i])), main = paste('Scaled Gene    Counts ',title_str), xlab = paste('log2(counts+',offset,')'), 
                 xlim = xlim, ylim = ylim)
    for(i in 1:ncol(.data)){
      lines(density(log2(.1+.data[,i])), col = i)
    }
  }else{
    plot <- plot(density(.data[,i]), main = paste('Scaled Gene Counts ',title_str),   xlab = 'Counts', 
                 xlim = xlim, ylim = ylim)
    for(i in 1:ncol(.data)){
      lines(density(.data[,i]), col = i)
    }
  }
  plot
  return(plot)
}

#function to remove lowly expressed genes
rm_low <- function(.data,datExpr.tpm, gtf, cutoff = 0.1, percent = 0.25){
  #remove unwanted Chr
  #remove lowly expressed
  # cutoff default 0.1, default 25% subjects
  datExpr.tpm <- datExpr.tpm %>% filter(genid %in% gtf$genid)
  rownames(datExpr.tpm) <- datExpr.tpm$genid
  datExpr.tpm <- datExpr.tpm[,-1]
  # cutoff=0.1, 25% subjects
  keep <- (rowSums(datExpr.tpm > cutoff)) > percent*ncol(datExpr.tpm)
  # # Count filter, use cpm, not raw or scaled from TPM counts!
  # keep <- filterByExpr(datExpr, min.count = (cutoff*100), min.prop=percent)
  datExpr <- .data[keep,]
  return(.data)
}
