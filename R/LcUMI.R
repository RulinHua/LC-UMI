
caculate_LcUMI <- function(cnt_mtx,mol_info,span=0.5) {
    library(data.table)
    library(Matrix)

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
    return(LcUMI_matrix)
  }
#' @export

integrate_samples <- function(list_rawUMI,list_LcUMI) {
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
#' @export

