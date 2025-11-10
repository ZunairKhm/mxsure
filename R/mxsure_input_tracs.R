#' Process TRACS distance input data
#'
#' Efficiently reads or processes output from the TRACS distance command to extract required columns
#' for MxSure. It can either extract an existing "date difference"
#' column or calculate it by joining a provided 'dates'.
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
mxsure_input_tracs <- function(input = NULL,
                               file_path = NULL,
                               dates = NULL,
                               dates_sample_col = "sample_id",
                               dates_date_col = "date") {


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
    if (!(is.data.table(input))) {
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
  if ("filtered SNP distance" %in% input_names) {
    output_list$snp_dist <- local_input$`filtered SNP distance`
  } else if ("SNP distance" %in% input_names) {
    output_list$snp_dist <- local_input$`SNP distance`
  } else {
    stop("Required column not found. No 'filtered SNP distance' or 'SNP distance' column.")
  }

  # Date difference
  if (!is.null(dates)) {
    # Option A: Calculate from dates

    # Check if required columns exist in the input
    if (!all(c("sampleA", "sampleB") %in% input_names)) {
      stop("To use 'dates', the 'input' data must contain 'sampleA' and 'sampleB' columns.")
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
    local_input[dates_lookup, on = .(sampleA = sample_key), dateA := i.date_val]
    local_input[dates_lookup, on = .(sampleB = sample_key), dateB := i.date_val]

    # Calculate time difference in days
    local_input[, time_dist_calc := abs(as.numeric(difftime(dateA, dateB, units = "days")))]

    output_list$time_dist <- local_input$time_dist_calc

  } else if ("date difference" %in% input_names) {
    # Option B: Fallback to existing column
    output_list$time_dist <- local_input$`date difference`
  } else {
    # Option C: No date info provided
    warning("Optional column 'date difference' not found and 'dates' not provided. Skipping time_dist.")
  }

  # --- Sites considered (Optional) ---
  if ("sites considered" %in% input_names) {
    output_list$sites <- local_input$`sites considered`
  } else {
    warning("Optional column 'sites considered' not found. Skipping.")
  }

  # 3. --- Return ---
  return(as.data.frame(output_list))
}
