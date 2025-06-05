# MRforAD (Mendelian Randomization for Alzheimer's Disease)

## Package overview
The `MRforAD` is an R package designed for Mendelian Randomizaion (MR) analysis, focusing on Alzheimer's Disease (AD). 
Its main purpose is to be used to study causal relationships in AD using genetics data. 
It uses publicly available GWAS summary statistics, such as those from Jansen et al. (2019), 
and supports integration with user-pprovided ADNI data from those with approved access.

## Introduction
Mendelian Randomization (MR) is a statistical method in genetic epigemiology used to infer causal relationships 
between casual relationships between an exposure and an outcome by using genetic variants as intrumental variables. 
This approach is particularly interesting for studying complex disease like Alzheimer's disease (AD), 
where observational studies mey be confounded by environmental factors. 

The Alzheimer's Disease Neuroimaging Initiative (ADNI) provides a rich dataset including genetic, clinical,
and imaging data; making it a primary candidate for MR analyses. However, ADNI data access is restricted, 
requiring approval through the LONI Image and Data Archive. 
With that being said, the development of an R package docuses on using publicaly available summary statisics for AD,
while providing functionality for any user with ADNI access. 

## Bioinformatic analyses background
MR relies on three key assumptions: the genetic instrument is associated with the exposure (relevance), is independent 
of confounders (independence), and affects the outcome only through the exposure (exlusion restriction). 
For AD, this method can hellp to identify causal risk factors, such as genetic variants influencing AD risk
through pathways like immune response or lipid metabolism. 

The ADNI cohort, with is longitudinal data on cognitively normal, mild cognitive impairment and AD/dementia participants, 
is ideal for such analyses. But its restricted access poses a challengee. 
To addresss this, the package uses publically available AD GWAS summary statistics, ensuring broad accessibility 
while maintaining relevance for ADNI users.
