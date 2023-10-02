## File of useful functions to load at top of script

#Function to create custom size & step batches from df, sorted via given column
create_batch <- function(df, batch_size, method = split, batch_step = NULL, save = FALSE, folder_name = NULL, loc = NULL){
  methods <- c("split", "sliding_window")
  method <- match.arg(method, methods)
  # Load or install necessary libraries
  required_pkg <- c("dplyr", "purr")
  load_install_pkg(required_pkg)
  # Split the dataframe into regular batches
  if (method == "split"){
    batches <- df %>%
    group_by(batch = cumsum(row_number() %% batch_size == 1)) %>%
    ungroup() %>%
    nest()
    }
    #Split into sliding window batches
    else if (method == "sliding_window"){
      #check if batch step argument is given
      if (!is.null(batch_step)) {
        for (i in seq(1, ncol(data), batch_step)) {
            end <- min(i + batch_size - 1, ncol(data))
            batch <- data[, i:end]
            batches[[length(batches) + 1]] <- batch
            }
        } else {stop("Batch Step size is missing. Cannot perform sliding window analysis without step size.")}
      
      }
    
  # Save all batchs if save == TRUE
    if (save == TRUE){
    setwd(loc)
    # Create a directory to store the batches
    dir.create(folder_name, showWarnings = FALSE)
    walk2(batches$data, seq_along(batches$data), batch_dir, save_batch)
    } else{return(batches$data)}
 }


# Function to save a batch to given directory or as a TEMP file
save_batch <- function(batch_df, batch_index, batch_dir, file_extension = "csv", temp = FALSE) {
  #library(VariantAnnotation)
  
  # Get directory from arguments
  # Generate the file name with the specified file extension
  batch_file_name <- paste0("batch.", batch_index, file_extension)
  batch_file_path <- file.path(batch_dir, batch_file_name)
  
  #if fct specifies creating temp files
  if (temp == TRUE){
    batch_file_path <- tempfile(batch_file_name, batch_dir, fileext = file_extension)
  }
  #if saving non-temp files
  else{
    path <- file.path(batch_dir_name, batch_file_name) 
  }
  
  # Save the batch dataframe with the specified file extension
  if (file_extension == "csv") {
    write.csv(batch_df$data, file = batch_file_path, row.names = FALSE)
  } else if (file_extension == "vcf") {
    writeVcf(batch_df$data, batch_file_path)
  } else {
    stop("Unsupported file extension. Use 'csv' or 'vcf'.")
  }
}

# Function to load or install required packages
load_install_pkg <- function(required_packages) {
  # Check if packages are installed, and install missing ones
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
  # Load all required packages
  lapply(required_packages, require, character.only = TRUE)
}

`%notin%` <- function(x, y) {
  !(x %in% y)
}
