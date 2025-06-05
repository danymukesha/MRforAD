# all necessary functions

#' Loading GWAS summary statistics
#'
#' This f(x) loads GWAS summary statistics from a file and standardizes
#' the column names (e.g., ADNI).
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

#' Selecting instruments for MR
#'
#' This f(x) selects genetic instruments for exposure GWAS data
#' based on a p-value threshold and optionally performs LD clumping.
#'
#' @param data `dataframe` containing exposureGWAS summary statistics
#'      with standardized columts.
#' @param pval_threshold `numeric` p-value threshold for instrument
#'      selection (defaul: 5e-8).
#' @param clump `logical` indicating whether to perform LD clumping
#'      (default: 0.001).
#' @param r2 `numeric` LD r-squared threshold for clumping (default: TRUE).
#' @param kb `numeric` clumping window in kilo-bases (default: 10000).
#'
#' @return a datafrae of selected instruments.
#'
#' @examples
#' \dontrun{
#' # Select instruments with clumping
#' instruments <- select_instrument(exposure_data, clump = FALSE)
#' }
#'
#' @importFrom TwoSampleMR clump_data
#' @importFrom dplyr filter
#'
#' @export
select_instruments <- function(data, pval_threshold = 5e-8, clump = TRUE,
                               r2 = 0.001, kb = 10000) {
    # filter for signif. SNPs
    sig_snps <- dplyr::filter(data, pval_threshold)

    if (nrow(sig_snps) == 0) {
        stop("Sorry, no SNPs meet the p-value threshold.")
    }

    # clumping if really needed
    if (clump) {
        clumped <-
            TwoSampleMR::clump_data(sig_snps, clump_r2 = r2, clump_kb - kb)
        return(clumped)
    } else {
        return(sig_snps)
    }
}


#' Running MR Pipeline
#'
#' This f(x) performs a complete Mendelian Randomization analysis pipeline
#' using exposure and outcome GWAS data, with a focus on Alzheimer's disease
#' as the outcome.
#'
#' @param exposure_data `dataframe` containing exposure GWAS summary statistics.
#' @param outcome_data (Optional) `daraframe` containing outcome GWAS summary
#'      statistics. If NULL, data is fetched from MR-Base using
#'      outcome_study_id.
#' @param outcome_study_id `character string` specifying the MR-Base study ID
#'      for the outcome (default: "ieu-a-297" for Jansen et al. 2019 AD GWAS).
#'
#' @return a list containing MR estimates, sensitivity analysis results, and plots.
#'
#' @examples
#' \dontrun{
#' # Using MR-Base data for both exposure and outcome
#' library(TwoSampleMR)
#' exp_dat <- extract_instruments("ieu-a-300") # LDL cholesterol
#' results <- run_mr_pipeline(exp_dat)
#'
#' # Using user-provided ADNI data as outcome
#' adni_data <- load_gwas_data("adni_data.txt")
#' results <- run_mr_pipeline(exp_dat, outcome_data = adni_data)
#' }
#'
#' @importFrom TwoSampleMR format_data extract_outcome_data harmonise_data mr
#' @importFrom TwoSampleMR mr_egger_regression mr_heterogeneity
#' @importFrom TwoSampleMR mr_scatter_plot mr_forest_plot mr_funnel_plot
#' @importFrom MRPRESSO mr_presso
#'
#' @export
run_mr_pipeline <- function(exposure_data, outcome_data = NULL,
                            outcome_study_id = "ieu-a-297") {

    exp_dat <- TwoSampleMR::format_data(exposure_data, type = "exposure")

    instruments <- select_instruments(exp_dat)

    if (is.null(outcome_data)) {
        # fetch outcome data from MR-Base
        out_dat <-
            TwoSampleMR::extract_outcome_data(
                snps = instruments$SNP, outcomes = outcome_study_id)
        if (is.null(out_dat) || nrow(out_dat) == 0) {
            stop("Failed to retrieve outcome data from MR-Base.
                 Check outcome_study_id or internet connection.")
        }
    } else {
        out_dat <- TwoSampleMR::format_data(outcome_data, type = "outcome")
    }

    dat <- TwoSampleMR::harmonise_data(
        exposure_dat = instruments, outcome_dat = out_dat)

    # MR analysis (Inverse-Variance Weighted method as default)
    mr_res <- TwoSampleMR::mr(
        dat, method_list =
            c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

    # doing the sensitivity analyses
    egger_intercept <-
        TwoSampleMR::mr_egger_regression(
            b_exp = dat$beta.exposure,
            b_out = dat$beta.outcome,
            se_exp = dat$se.exposure,
            se_out = dat$se.outcome,
            parameters = list(over.dispersion = FALSE,
                              loss.function = "l2"))$egger_intercept

    heterogeneity <- TwoSampleMR::mr_heterogeneity(dat)
    presso <- tryCatch(
        MRPRESSO::run_mr_presso(dat),
        error = function(e) list(
            message = "MR-PRESSO failed: install package or check data"))

    scatter_plot <- TwoSampleMR::mr_scatter_plot(mr_res, dat)
    forest_plot <- TwoSampleMR::mr_forest_plot(dat)
    funnel_plot <- TwoSampleMR::mr_funnel_plot(dat)

    return(list(
        mr_estimates = mr_res,
        sensitivity = list(
            egger_intercept = egger_intercept,
            heterogeneity = heterogeneity,
            presso = presso
        ),
        plots = list(
            scatter = scatter_plot[[1]],
            forest = forest_plot[[1]],
            funnel = funnel_plot[[1]]
        )
    ))
}
