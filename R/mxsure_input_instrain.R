#' Process inStrain compare input data
#'
#' Efficiently reads or processes the output from the inStrain compare command to extract required columns
#' for MxSure It can calculate date differences between samples it by joining a provided 'dates'.
#'
#' @param input A data.frame or data.table. If NULL, file_path is used.
#' @param file_path A string path to a file to be read by fread().
#' @param dates A data.frame or data.table with sample IDs and their collection dates.
#' @param dates_sample_col A string, the name of the sample ID column in 'dates'.
#' @param dates_date_col A string, the name of the date column in 'dates'.
#'
#' @return A data.table with columns: snp_dist, time_dist (if found or if sampling dates are provided), and sites (if found).
#'
#' @import data.table
#'
#' @export
mxsure_input_instrain <- function(input = NULL,
                               file_path = NULL,
                               dates = NULL,
                               dates_sample_col = "sample_id",
                               dates_date_col = "date") {

  library(data.table)


  #  Input Handling
  if (is.null(input) && is.null(file_path)) {
    stop("Either 'input' (a data.frame/data.table) or 'file_path' (a string) must be provided.")
  }

  local_input <- NULL

  if (!is.null(file_path)) {
    if (!is.null(input)) {
      warning("Both 'input' and 'file_path' provided. 'file_path' will override 'input'.")
    }
    local_input <- fread(file=file_path)
  } else {
    if (!is.data.table(input)) {
      local_input <- as.data.table(input)
    } else {
      # Use a copy to avoid modifying the original data.table by reference
      local_input <- copy(input)
    }
  }

  # Column Extraction
  output_list <- list()
  input_names <- names(local_input)

  #SNP distance
  if ("population_SNPs" %in% input_names) {
    output_list$snp_dist <- local_input$`population_SNPs`
  } else {
    stop("Required column not found. No 'population_SNPs' column.")
  }

  # Date difference
  if (!is.null(dates)) {
    # Option A: Calculate from dates

    # Check if required columns exist in the input
    if (!all(c("name1", "name2") %in% input_names)) {
      stop("To use 'dates', the 'input' data must contain 'name1' and 'name2' columns.")
    }

    # Ensure dates is a data.table and prepare it for joining
    if (!is.data.table(dates)) {
      dates <- as.data.table(dates)
    }

    cols_to_keep <- c(dates_sample_col, dates_date_col)
    dates_lookup <- dates[, cols_to_keep, with=F]

    # Ensure date column is Date type
    tryCatch({
      dates_lookup[, (dates_date_col) := as.Date(get(dates_date_col))]
    }, error = function(e) {
      stop(paste("Failed to convert date column:", dates_date_col, ". Ensure it's in a standard format (e.g., YYYY-MM-DD)."))
    })

    # Set standard names for joining
    setnames(dates_lookup, old = c(dates_sample_col, dates_date_col), new = c("sample_key", "date_val"))
    setkey(dates_lookup, sample_key)

    # Perform two update-joins (very fast, no copies made)
    local_input[dates_lookup, on = .(name1 = sample_key), dateA := i.date_val]
    local_input[dates_lookup, on = .(name2 = sample_key), dateB := i.date_val]

    # Calculate time difference in days
    local_input[, time_dist_calc := abs(as.numeric(difftime(dateA, dateB, units = "days")))]

    output_list$time_dist <- local_input$time_dist_calc

  }

  # --- Sites considered (Optional) ---
  if ("compared_bases_count" %in% input_names) {
    output_list$sites <- local_input$`compared_bases_count`
  } else {
    warning("Column 'compared_bases_count' not found. Skipping.")
  }

  # 3. --- Return ---
  return(as.data.frame(output_list))
}
