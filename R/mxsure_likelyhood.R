#' Reports Likelyhoods for related and unrelated distributions for each datapoint provided
#'
#'
#'
#' @param mixed_snp_dist vector of SNP distances from mixed dataset
#' @param unrelated_snp_dist vector of SNP distances from unrelated dataset
#' @param mixed_time_dist vector of time differences for each SNP distacne in the mixed dataset
#' @param mixed_sites vector of sites considered for each SNP distance in the mixed dataset
#' @param truncation_point maximum limit of SNP distances to consider
#'
#' @importFrom dplyr distinct mutate
#' @importFrom tibble tibble
#' @importFrom stats dpois dnbinom
#'
#' @return a dataframe with SNP distances, time differences, sites considered and the likeyhoods of related and unrelated models fitting for each datapoint
#' @export
mxsure_likelyhood <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point=2000,original_result=NA,start_params=NA,
                              tree=NA, sampleA=NA, sampleB=NA,  branch_offset=NA){

  snp_dist<- time_dist<- rel_loglh<- unrel_loglh<- logLHR <-  NULL


  if(is.na(truncation_point)){
    truncation_point <- Inf
  }
  unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]

  if(anyNA(original_result)){
  mix_res <- suppressWarnings(mxsure_estimate(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point = 2000,start_params=start_params,
                                              tree=tree, sampleA=sampleA, sampleB=sampleB, branch_offset=branch_offset))
  }else{
  mix_res <- original_result
}


  if(is.na(mean(mixed_sites, na.rm=TRUE))){
    mixed_sites <- 1
  }

  #### tree likelyhood ####
  if(!anyNA(tree)|!anyNA(sampleA)|!anyNA(sampleB)){
    # Ensure tree is a list (even if single tree provided)
    tree_list <- if (inherits(tree, "phylo")) list(tree) else tree

    # Internal function to compute corrected time for a single pair
    compute_shared_snps <- function(tip1, tip2, snp_dist, time_diff) {
      for (tree_i in tree_list) {
        tips <- tree_i$tip.label
        if (tip1 %in% tips && tip2 %in% tips) {
          mrca_node <- ape::getMRCA(tree_i, c(tip1, tip2))
          if (is.null(mrca_node)) return(snp_dist)

          tip1_node <- which(tree_i$tip.label == tip1)
          tip2_node <- which(tree_i$tip.label == tip2)

          # Get distances from MRCA to each tip
          edge_dist <- ape::dist.nodes(tree_i)

          mrca_to_tip1 <- edge_dist[mrca_node, tip1_node]
          mrca_to_tip2 <- edge_dist[mrca_node, tip2_node]
          root_to_mrca <- phytools::fastHeight(tree_i, tip1, tip2)
          #shared_snps <- min(mrca_to_tip1, mrca_to_tip2)

          #if (distance_since_sample<0  ) return(NA)
          #if( full_distance / distance_since_sample > max_correction_factor) return(NA)

          return(tibble(mrca_to_tip1=mrca_to_tip1, mrca_to_tip2=mrca_to_tip2, root_to_mrca=root_to_mrca))  #shared_snps) #time_diff * full_distance / distance_since_sample)
        }
      }
      # If no matching tree or MRCA, return original time
      return(tibble(mrca_to_tip1=NA, mrca_to_tip2=NA, root_to_mrca=NA))
    }

    # Vectorised using pmap
    branch_lengths <- purrr::pmap_dfr(
      list(sampleA, sampleB, mixed_snp_dist, mixed_time_dist),
      compute_shared_snps
    )


    # if (is.na(max_time)) {
    #   max_time <- max(mixed_time_dist)
    # }

    mixed_snp_dist <- mixed_snp_dist[!anyNA(branch_lengths)]
    mixed_time_dist <- mixed_time_dist[!anyNA(branch_lengths)]
    mixed_sites <- mixed_sites[!anyNA(branch_lengths)]
    sampleA <- sampleA[!anyNA(branch_lengths)]
    sampleB <- sampleB[!anyNA(branch_lengths)]
    if(!is.na(branch_offset)){
      branch_lengths$root_to_mrca <- branch_offset
    }

    LH <-  distinct(
      tibble(
        snp_dist = mixed_snp_dist,
        time_dist = mixed_time_dist,
        sites = mixed_sites,
        mrca_to_tip1 = branch_lengths$mrca_to_tip1,
        mrca_to_tip2 = branch_lengths$mrca_to_tip2,
        root_to_mrca = branch_lengths$root_to_mrca
      )
    )

    LH <- LH |>
      rowwise()|>
      mutate(rel_loglh = (
        skellam::dskellam(mrca_to_tip1 - mrca_to_tip2,
          (mix_res$lambda * (time_dist / 365.25) * sites )+ mix_res$intercept + mrca_to_tip2 + root_to_mrca,
          mrca_to_tip2+ root_to_mrca,
          log = TRUE) +
          dpois(x = mrca_to_tip1,
                lambda = mix_res$single_branch_lambda,
                log=TRUE)+
          dpois(x = mrca_to_tip2,
                lambda = mix_res$single_branch_lambda,
                log=TRUE)
      ),
      unrel_loglh = (
        dnbinom(
          snp_dist,
          mu = mix_res$nb_mu,
          size = mix_res$nb_size,
          log = TRUE
        ) - pnbinom(
          truncation_point,
          mu = mix_res$nb_mu,
          size = mix_res$nb_size,
          log.p = TRUE
        )
      )) |>
      mutate(logLHR = (rel_loglh - unrel_loglh)) |>
      mutate(
        rel_lh = exp(rel_loglh),
        unrel_lh = exp(unrel_loglh),
        LHR = exp(logLHR)
      )


    return(LH)

#### non tree likelyhood #####
}else{
  LH <-
    tibble(
      snp_dist = mixed_snp_dist,
      time_dist = mixed_time_dist,
      sites = mixed_sites
    )


  LH <- LH |>
    mutate(rel_loglh = (
      dpois(snp_dist, (
        mix_res$lambda * (time_dist / 365.25) * sites + mix_res$intercept
      ), log = TRUE) #/ ppois(truncation_point, (mix_res$lambda*(time_dist/365.25)*mean(mixed_sites)+mix_res$intercept))
    ),
    unrel_loglh = (
      dnbinom(
        snp_dist,
        mu = mix_res$nb_mu,
        size = mix_res$nb_size,
        log = TRUE
      ) - pnbinom(
        truncation_point,
        mu = mix_res$nb_mu,
        size = mix_res$nb_size,
        log.p = TRUE
      )
    )) |>
    mutate(logLHR = (rel_loglh - unrel_loglh)) |>
    mutate(
      rel_lh = exp(rel_loglh),
      unrel_lh = exp(unrel_loglh),
      LHR = exp(logLHR)
    )


  return(LH)

  }

}
