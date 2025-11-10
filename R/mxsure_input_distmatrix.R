#' Process a TRACS distance matrix for MxSure
#'
#' Efficiently reads a distance matrix, "melts" it into a long-form pairwise
#' table, and optionally joins a 'dates' table to calculate time differences.
#'
#' @param input A data.frame or data.table (distance matrix). If NULL, file_path is used.
#' @param file_path A string path to a file to be read by fread() (distance matrix).
#' @param rownames If input has row names (TRUE) or if column one are sample names (FALSE)
#' @param dates A data.frame or data.table with sample IDs and their collection dates.
#' @param dates_sample_col A string, the name of the sample ID column in 'dates'.
#' @param dates_date_col A string, the name of the date column in 'dates'.
#'
#' @return A data.table with columns: sampleA, sampleB, snp_dist,
#'   and time_dist (if 'dates' was provided).
#'
#' @import data.table
#'
#' @export
mxsure_input_distmatrix <- function(input = NULL,
                                    file_path = NULL,
                                    dates = NULL,
                                    dates_sample_col = "sample_id",
                                    dates_date_col = "date") {
  library(data.table)


  # Input Handling
  if (is.null(input) && is.null(file_path)) {
    stop("Either 'input' or 'file_path' must be provided.")
  }

  local_input <- NULL

  if (!is.null(file_path)) {
    if (!is.null(input)) {
      warning("Both 'input' and 'file_path' provided. 'file_path' will override 'input'.")
    }
    local_input <- fread(file=file_path, header = T)
    setnames(local_input, 1, "sampleA")
  } else {
    if(sum(row.names(input)!=1:(nrow(input)))==0){ #detects if there are row names
      local_input <- as.data.table(input)
      setnames(local_input, 1, "sampleA")
    }else{
      local_input <- as.data.table(input, keep.rownames = "sampleA")
      }
  }

  # Melt the matrix from wide to long
  output_dt <- melt(local_input,
                                id.vars = "sampleA",
                                variable.name = "sampleB",
                                value.name = "snp_dist",
                                variable.factor = FALSE)

  # 3. --- Clean Melted Data ---
  # Remove self-comparisons (e.g., sample1 vs sample1)
  output_dt <- output_dt[sampleA != sampleB]

  # Remove duplicate pairs (e.g., keep A vs B, drop B vs A)
  # Create a sorted key for each pair
  output_dt[, c("pairA", "pairB") := .(pmin(sampleA, sampleB), pmax(sampleA, sampleB))]
  # Keep only the first unique row for each sorted pair
  output_dt <- unique(output_dt, by = c("pairA", "pairB"))

  # Tidy up - remove the helper columns
  output_dt[, c("pairA", "pairB") := NULL]

  # 4. --- Date difference ---
  if (!is.null(dates)) {
    # Option A: Calculate from dates_df

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

    # Perform two update-joins on our new 'output_dt'
    output_dt[dates_lookup, on = .(sampleA = sample_key), dateA := i.date_val]
    output_dt[dates_lookup, on = .(sampleB = sample_key), dateB := i.date_val]

    # Calculate time difference in days
    output_dt[, time_dist := abs(as.numeric(difftime(dateA, dateB, units = "days")))]

    # Remove temporary date columns
    output_dt[, c("dateA", "dateB") := NULL]

  } else {
    # Option B: No date info provided
    warning("No 'dates' data.frame provided. 'time_dist' column will not be created.")
  }


  # 6. --- Return ---
  # Select final columns (in case time_dist wasn't created)
  final_cols <- c("sampleA", "sampleB", "snp_dist")
  if ("time_dist" %in% names(output_dt)) {
    final_cols <- c(final_cols, "time_dist")
  }

  return(output_dt[, ..final_cols])
}
