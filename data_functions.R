## File of useful functions to load at top of script

#Function to create custom size & step batches from df, sorted via given column
create_batch <- function(df, sort_col, batch_size, batch_step, save = FALSE){
  # Load necessary libraries
  library(dplyr)
  library(purrr)
  # Specify the batch size
  batch_size <- batch_size  #input variable
  
  # Split the dataframe into batches
  batches <- df %>%
    arrange(sort_col) %>%
    group_by(batch = cumsum(row_number() %% batch_size == 1)) %>%
    ungroup() %>%
    nest()
  
  return(batches$data)
  
    # Save all batchs if save == TRUE
    if (save == TRUE){
    # Create a directory to store the batches
    batch_dir <- "batches"  # Change this to your desired directory name
    dir.create(batch_dir, showWarnings = FALSE)
    walk2(batches$data, seq_along(batches$data), batch_dir, save_batch)
    }
 }


# Function to save a batch to given directory or as a TEMP file
save_batch <- function(batch_df, batch_index, batch_dir, file_extension = "csv", temp = FALSE) {
  library(VariantAnnotation)
  
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