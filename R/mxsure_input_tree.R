#' Process a Phylogenetic Tree for MxSure
#'
#' Efficiently reads a phylogenetic tree, calculates pairwise tip distances,
#' checks for valid SNP-scale branch lengths, rounds them to integers,
#' and optionally maps a vector of dates to calculate time differences.
#'
#' @param tree A 'phylo' object (from the 'ape' package). If NULL, file_path is used.
#' @param file_path A string path to a tree file (e.g., Newick format) to be read by ape::read.tree().
#' @param tip_dates A vector of dates corresponding to the order of tips in the tree.
#' @param tip_person_ids A vector of person/patient IDs corresponding to the order of tips in the tree.
#'
#' @return A data.table with columns: sampleA, sampleB, snp_dist, time_dist (if 'tip_dates' was provided), and same_person (if 'tip_person_ids' provided).
#'
#' @import data.table
#' @importFrom ape read.tree cophenetic.phylo
#'
#' @export
mxsure_input_tree <- function(tree = NULL,
                              file_path = NULL,
                              tip_dates = NULL,
                              tip_person_ids = NULL) {

  library(data.table)

  # Input Handling
  if (is.null(tree) && is.null(file_path)) {
    stop("Either 'tree' or 'file_path' must be provided.")
  }

  local_tree <- NULL

  if (!is.null(file_path)) {
    if (!is.null(tree)) {
      warning("Both 'tree' and 'file_path' provided. 'file_path' will override 'tree'.")
    }
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("The 'ape' package is required to read trees. Please install it.")
    }
    local_tree <- ape::read.tree(file_path)
  } else {
    local_tree <- tree
  }

  if (!inherits(local_tree, "phylo")) {
    stop("Input must be a phylogenetic tree of class 'phylo'.")
  }

  if (is.null(local_tree$edge.length)) {
    stop("The provided tree does not have branch lengths.")
  }

  # Extract Pairwise Distances
  dist_mat <- ape::cophenetic.phylo(local_tree)
  off_diag_dists <- dist_mat[upper.tri(dist_mat)]

  # Branch Length Checks & Rounding
  max_dist <- max(off_diag_dists, na.rm = TRUE)
  mean_dist <- mean(off_diag_dists, na.rm = TRUE)

  if (max_dist < 2) {
    stop("All pairwise tip distances are below 2. Do not proceed: ensure your tree's branch lengths represent total SNPs, not substitutions per site.")
  }

  if (mean_dist < 5) {
    warning("The mean pairwise tip distance is very low (< 5). Double-check that your tree represents total SNPs.")
  }

  # Force distances to integers (whole SNPs)
  dist_mat <- round(dist_mat)

  # 4. --- Melt Matrix from Wide to Long ---
  local_input <- as.data.table(dist_mat, keep.rownames = "sampleA")

  output_dt <- melt(local_input,
                    id.vars = "sampleA",
                    variable.name = "sampleB",
                    value.name = "snp_dist",
                    variable.factor = FALSE)

  # Clean Melted Data
  output_dt <- output_dt[sampleA != sampleB]
  output_dt[, c("pairA", "pairB") := .(pmin(sampleA, sampleB), pmax(sampleA, sampleB))]
  output_dt <- unique(output_dt, by = c("pairA", "pairB"))
  output_dt[, c("pairA", "pairB") := NULL]

  # Date Difference
  if (!is.null(tip_dates)) {

    if (is.data.frame(tip_dates)) {
      if (ncol(tip_dates) == 1) {
        tip_dates <- tip_dates[[1]]
      } else {
        stop("'tip_dates' is a data.frame with multiple columns. Please provide a single vector.")
      }
    }

    # Check length match
    if (length(tip_dates) != length(local_tree$tip.label)) {
      stop("The 'tip_dates' vector must have the exact same length as the number of tips in the tree.")
    }

    # Safely convert to dates and bind to tip labels
    tryCatch({
      dates_lookup <- data.table(
        sample_key = local_tree$tip.label,
        date_val = as.Date(as.vector(tip_dates))
      )
    }, error = function(e) {
      stop("Failed to convert 'tip_dates' to Date format. Ensure it contains valid dates (e.g., 'YYYY-MM-DD').")
    })

    setkey(dates_lookup, sample_key)

    # Perform update-joins
    output_dt[dates_lookup, on = .(sampleA = sample_key), dateA := i.date_val]
    output_dt[dates_lookup, on = .(sampleB = sample_key), dateB := i.date_val]

    # Calculate time difference in days
    output_dt[, time_dist := abs(as.numeric(difftime(dateA, dateB, units = "days")))]
    output_dt[, c("dateA", "dateB") := NULL]

  } else {
    warning("No 'tip_dates' vector provided. 'time_dist' column will not be created.")
  }

  # Person ID Tagging (Within-Host vs Between-Host)
  if (!is.null(tip_person_ids)) {
    if (is.data.frame(tip_person_ids)) {
      if (ncol(tip_person_ids) == 1) {
        tip_person_ids <- tip_person_ids[[1]]
      } else {
        stop("'tip_person_ids' is a data.frame with multiple columns. Please provide a single vector.")
      }
    }


    if (length(tip_person_ids) != length(local_tree$tip.label)) {
      stop("The 'tip_person_ids' vector must have the exact same length as the number of tips in the tree.")
    }

    person_lookup <- data.table(
      sample_key = local_tree$tip.label,
      person_val = tip_person_ids
    )
    setkey(person_lookup, sample_key)

    # Map person IDs to sampleA and sampleB
    output_dt[person_lookup, on = .(sampleA = sample_key), personA := i.person_val]
    output_dt[person_lookup, on = .(sampleB = sample_key), personB := i.person_val]

    # Evaluate if they are the same person
    output_dt[, same_person := (personA == personB)]

    # Clean up temporary columns
    output_dt[, c("personA", "personB") := NULL]
  }

  # 8. --- Return ---
  final_cols <- c("sampleA", "sampleB", "snp_dist")
  if ("time_dist" %in% names(output_dt)) {
    final_cols <- c(final_cols, "time_dist")
  }
  if ("same_person" %in% names(output_dt)) {
    final_cols <- c(final_cols, "same_person")
  }

  return(output_dt[, ..final_cols])

}
