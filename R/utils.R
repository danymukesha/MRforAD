# all necessary functions

#' Loading GWAS summary statistics
#'
#' This f(x) loads GWAS summary statistics from a file and standardizes
#' the column names.
#'
#' @param file_path `a character string` specifying the path to the GWAS
#'      summary statistics file (tab-limited).
#' @param col_map (Optional) `named list` mapping required column names
#'      (`SNP`, `beta`, `se`, `pval`, `effect_allele`, `other_allele`,
#'      `eaf`) to actual column names in the file.
#'      If NULL, it assumes that the columns are already correctly named.
#'
#' @return a dataframe with standardized column names: `SNP`, `beta`,
#' `se`, `pval`, `effect_allele`, `other_allele`, `eaf`.
#'
#' @examples
#' \dontrun{
#' #e.g. with custom column mapping
#' col_map <- list(SNP = "rsid", beta "b", se = "standard_error",
#'                 pval = "p", effect_allele = "a", other_allele = "a2",
#'                 eaf = "freq")
#' gwas_data <- load_gwas_data("path/to/gwas.txt", col_map)
#' }
#'
#' @importFrom readr read_delim
#' @importFrom dplyr rename
#'
#' @export
load_gwas_data <- function(file_path, col_map = NULL) {
    data <- readr::read_delim(file_path, delim = "\t", progress = FALSE)

    if (!is.null(col_map)) {
        data <- dplyr::rename(data, !!!col_map)
    }

    required_cols <- c("SNP", "beta", "se", "pval",
                       "effect_allele", "other_allele", "eaf")

    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        stop("We are missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }

    data
}
