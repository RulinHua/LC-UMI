#' @title Integrate Multiple LC-UMI Matrices Across Samples
#'
#' @description
#' This function integrates multiple LOESS-corrected UMI count matrices (LC-UMI) across samples
#' by correcting for sample-level variability in sequencing depth and cell number.
#' Specifically, a normalization factor is computed for each sample based on the average total UMI count per cell
#' in the corresponding raw UMI matrix. These factors are then used to rescale the LC-UMI matrices,
#' ensuring comparability of UMI counts across samples prior to integration.
#' This approach ensures that variation in total RNA abundance across samples due to differences in sequencing depth and cell recovery is minimized.
#'
#' @param list_rawUMI A named list of raw UMI count matrices (genes × cells), one matrix per sample.
#' The names of the list must correspond exactly to those in \code{list_LcUMI}.
#' For each sample, the raw and LC-UMI matrices must contain the same set of genes (rows) and cells (columns),
#' with consistent row and column names.
#'
#' @param list_LcUMI A named list of LOESS-corrected UMI count matrices (genes × cells), one for each sample.
#' The names of this list must exactly match those in \code{list_rawUMI}.
#' For each sample, the corresponding raw and corrected matrices must have identical gene sets (row names)
#' and cell barcodes (column names).
#' \strong{Note:} Cell barcodes (i.e., column names) must be unique across all samples to avoid duplication during matrix integration.
#'
#' @return
#' A single gene expression matrix obtained by column-wise concatenation of all corrected LC-UMI matrices,
#' in which each sample has been scaled by a sample-specific normalization factor derived from its raw UMI matrix.
#' This adjustment compensates for differences in sequencing depth and cell number across samples,
#' enabling fair comparison of gene expression values at the single-cell level.
#'
#' @details
#' For each sample, a normalization factor is computed as the average total UMI count per cell
#' (i.e., \code{sum(UMIs) / number of cells}). These factors are normalized to their median across all samples,
#' and each LC-UMI matrix is scaled accordingly. This ensures that technical variation due to sequencing depth
#' and cell recovery is minimized when integrating multiple samples.
#' It is required that cell barcodes (i.e., column names) are unique across all samples.
#' If different samples contain overlapping cell names, they should be prefixed or modified beforehand
#' to ensure uniqueness, as duplicate column names will cause errors or data loss during integration.
#'
#' @examples
#' \dontrun{
#' # Assume you have 3 samples with matching raw and LC-UMI matrices
#' # Make sure that all cell barcodes are unique across samples
#' raw_list <- list(S1 = cnt_mtx1, S2 = cnt_mtx2, S3 = cnt_mtx3)
#' lcumi_list <- list(S1 = LcUMI_matrix1, S2 = LcUMI_matrix2, S3 = LcUMI_matrix3)
#' combined <- integrate_samples(raw_list, lcumi_list)
#' }
#'
#' @export
integrate_samples <- function(list_rawUMI,list_LcUMI) {
  if (any(duplicated(unlist(lapply(list_LcUMI, colnames))))) {
    stop("Cell names (column names) in list_LcUMI must be unique across all samples.")
  }
  norm_factor <- unlist(lapply(list_rawUMI, function(x){sum(x)/ncol(x)}))
  names(norm_factor) <- names(list_rawUMI)
  norm_factor <- norm_factor / median(norm_factor)
  lst <- list()
  for (i in names(list_LcUMI)) {
    lst[[i]] <- list_LcUMI[[i]] / norm_factor[i]
  }
  combined_matrix <- do.call(cbind,lst)
  return(combined_matrix)
}
