#' @title LOESS-based Correction of UMI Amplification Bias (LC-UMI)
#'
#' @description
#' Performs LOESS regression-based correction of PCR amplification bias in UMI-based single-cell RNA-seq data.
#' The method models the nonlinear relationship between each cell’s average reads per UMI and its total UMI count,
#' and adjusts total UMI counts accordingly to reduce cross-cell variability introduced by non-uniform amplification.
#'
#' @param cnt_mtx A sparse matrix or dense matrix of raw UMI counts with genes in rows and cells in columns (e.g., a \code{dgCMatrix} or \code{matrix}).
#' The column names of \code{cnt_mtx} must correspond to cell barcodes present in \code{mol_info$cell}.
#'
#' @param mol_info A data.frame with three columns: \code{cell} (barcode), \code{umi} (UMI identifier), and \code{read} (read count per UMI).
#' Only cells that match the column names of \code{cnt_mtx} will be used for correction.
#'
#' @param span A numeric value controlling the LOESS smoothing span (default is \code{0.5}). Smaller values capture finer details but may overfit.
#' @param plot_path Optional file path (character string). If provided, a diagnostic LOESS fit plot will be saved to this location (e.g., "fit_plot.png").
#'
#' @return
#' A UMI count matrix of the same dimensions as \code{cnt_mtx}, where each cell’s gene expression values have been scaled to correct for amplification bias.
#'
#' @details
#' This function computes the average read count per UMI for each cell, fits a LOESS curve between this value and the log10-transformed total UMI counts,
#' and uses the residuals to estimate corrected UMI totals. The adjusted UMI matrix preserves the original count structure while reducing technical variability due to PCR amplification.
#' The corrected matrix can subsequently be log-transformed (e.g., \code{log1p}) for use in downstream normalization and analysis.
#'
#' @examples
#' # Define dimensions
#' set.seed(123)
#' n_cells <- 100
#' n_genes <- 500
#' cell_names <- paste0("Cell", seq_len(n_cells))
#' gene_names <- paste0("Gene", seq_len(n_genes))
#'
#' # Simulate total UMI count per cell (log-normal distribution)
#' total_umi_per_cell <- rlnorm(n_cells, meanlog = 4, sdlog = 0.4)
#' total_umi_per_cell <- round(total_umi_per_cell)
#'
#' # Simulate raw count matrix (sparse, randomly distributed UMIs)
#' cnt_mtx <- matrix(0, nrow = n_genes, ncol = n_cells,
#'                   dimnames = list(gene_names, cell_names))
#'
#' for (i in seq_len(n_cells)) {
#'   umi_count <- total_umi_per_cell[i]
#'   umi_genes <- sample(n_genes, umi_count, replace = TRUE)
#'   gene_table <- table(umi_genes)
#'   cnt_mtx[as.integer(names(gene_table)), i] <- as.integer(gene_table)
#' }
#'
#' # Simulate mol_info: each row is one UMI, with read count positively correlated to cell UMI total
#' mol_info <- data.frame()
#' for (i in seq_len(n_cells)) {
#'   cell_id <- cell_names[i]
#'   n_umis <- sum(cnt_mtx[, i])
#'   umi_ids <- paste0("UMI", seq_len(n_umis))
#'
#'   # Simulate read count per UMI: Poisson with lambda increasing with total_UMI
#'   base_reads <- rpois(n_umis, lambda = 3 + 0.006 * total_umi_per_cell[i])
#'
#'   cell_dt <- data.frame(
#'     cell = rep(cell_id, n_umis),
#'     umi = umi_ids,
#'     read = base_reads + 1  # ensure minimum read = 1
#'   )
#'   mol_info <- rbind(mol_info, cell_dt)
#' }
#'
#' # Run LC-UMI correction and save diagnostic plot
#' plot_path <- tempfile(fileext = ".png")
#' corrected <- caculate_LcUMI(cnt_mtx, mol_info, span = 0.7, plot_path = plot_path)
#' print(paste("Plot saved to:", plot_path))
#'
#' @import data.table
#' @importFrom Matrix colSums
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw ggsave
#' @export
caculate_LcUMI <- function(cnt_mtx,mol_info,span=0.5,plot_path=NULL) {
  colnames(mol_info) <- c("cell","umi","read")
  mol_info <- as.data.table(mol_info)
  mol_info <- mol_info[cell %in% colnames(cnt_mtx)]
  mol_info <- mol_info[,.(mean_read=mean(read)),by=cell]
  dt <- data.table(cell=colnames(cnt_mtx),total_umi=colSums(cnt_mtx))
  dt <- merge(dt,mol_info,by="cell")
  dt <- dt[match(colnames(cnt_mtx),dt$cell)]
  dt$log10_total_umi <- log10(dt$total_umi)
  model <- loess(log10_total_umi ~ mean_read, data = dt,span = span)
  dt$fitted_log10_umi <- predict(model)
  dt$corrected_total_umi <- 10 ^ (dt$log10_total_umi - dt$fitted_log10_umi)
  dt$corrected_total_umi <- dt$corrected_total_umi * median(dt$total_umi) / median(dt$corrected_total_umi)
  LcUMI_matrix <- t(t(cnt_mtx) * (dt$corrected_total_umi / dt$total_umi))

  if (!is.null(plot_path)) {
    p <- ggplot(dt, aes(x = mean_read, y = log10_total_umi)) +
      geom_point(alpha = 0.3, size = 0.5, color = "grey") +
      geom_smooth(method = "loess", span = span, se = FALSE, color = "blue") +
      theme_bw() +
      labs(x = "Mean UMI Amplification (Reads per UMI)",
           y = "log10(Total UMI Count)",
           title = "LOESS Fit: Mean Read vs. log10 Total UMI")
    ggsave(filename = plot_path, plot = p, width = 5, height = 5, dpi = 300)
  }
  return(LcUMI_matrix)
}
