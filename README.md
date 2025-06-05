---
editor_options: 
  markdown: 
    wrap: 72
---

# MRforAD

## Package overview

The `MRforAD` is an R package designed for Mendelian Randomization (MR)
analysis, focusing on Alzheimer's Disease (AD). Its main purpose is to
be used to study causal relationships in AD using genetics data. It uses
publicly available GWAS summary statistics, such as those from Jansen et
al. (2019), and supports integration with user-provided ADNI data from
those with approved access.

## Introduction

Mendelian Randomization (MR) is a statistical method in genetic
epidemiology used to infer causal relationships between casual
relationships between an exposure and an outcome by using genetic
variants as instrumental variables. This approach is particularly
interesting for studying complex disease like Alzheimer's disease (AD),
where observational studies may be confounded by environmental factors.

The Alzheimer's Disease Neuroimaging Initiative (ADNI) provides a rich
dataset including genetic, clinical, and imaging data; making it a
primary candidate for MR analyses. However, ADNI data access is
restricted, requiring approval through the LONI Image and Data Archive.
With that being said, the development of an R package focuses on using
publicly available summary statistics for AD, while providing
functionality for any user with ADNI access.

## Bioinformatic analyses background

MR relies on three key assumptions: the genetic instrument is associated
with the exposure (relevance), is independent of confounders
(independence), and affects the outcome only through the exposure
(exclusion restriction). For AD, this method can help to identify causal
risk factors, such as genetic variants influencing AD risk through
pathways like immune response or lipid metabolism.

The ADNI cohort, with is longitudinal data on cognitively normal, mild
cognitive impairment and AD/dementia participants, is ideal for such
analyses. But its restricted access poses a challenge. To address this,
the package uses publicly available AD GWAS summary statistics, ensuring
broad accessibility while maintaining relevance for ADNI users.

## Data and functionality

Since ADNI data requires approved access through the LONI Image and Data
Archive (<http://adni.loni.usc.edu>), this tool cannot include it
directly. However, i used directly summary statistics from Jansen et al.
(2019), available at <https://ctg.cncr.nl/software/summary_statistics>.
I included functionns to download these datasets, perform MR analyses
(including instrument selection, harmonization, and sensitivity checks),
and visualize results by making is suitably for the users to present the
results.

**Core functions:**

-   `download_gwas_data_url(url)`: this f(x) downloads GWAS summary
    statistics from specified URLs.
-   `get_ad_gwas_summary(study_id)`: this fx() retrieves data from the
    GWAS/NIAGADS API, which provides programmatic access to publicly AD
    GWAS summary statistic.
-   `run_mr_pipeline(exposure_data, outcome_data)`: this f(x) performs
    the full MR pipeline, including instrument selection, harmonization,
    MR analysis (e.g., IVW, MR-Egger, weighted median), sensitivity
    analyses (e.g., MR-egger intercept, heterogeneity tests), and
    visualization.
-   `plot_mr_results(mr_results)`: this f(x) generates standard MR plots
    like scatter, forest and funnel plots.
-   `prepare_mr_data(exposure_file, outcome_file)`: this f(x) loads and
    prepares user-provided data for analysis.
