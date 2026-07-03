# Negative-controls-in-MR

This repository contains R code used to implement the real-data analysis of the paper "Using Negative Control Outcomes to Detect Selection Bias in Mendelian Randomization Studies". The paper is available as a pre-print [here](https://www.medrxiv.org/content/10.64898/2026.01.30.26345215v1).

The analysis uses summary-level genetic data from [Schoeler et al. (2023)](https://www.nature.com/articles/s41562-023-01579-9). The datasets are freely available for download [here](https://drive.google.com/drive/folders/1QMq9Y_BDK-9o0TlaJL-fwaqQl4BbR6Ek).

The R code used to run the analysis is included in the file "Negative_Controls_Application.R". The main analysis requires the user to have downloaded summary-level data from the above repository. For the additional analysis using natural hair colour as a negative control outcome, the data for the binary "hair colour" variables can be accessed directly from [OpenGWAS](https://opengwas.io/), as described in our R code. The summary-level data for the categorical "natural hair colour" variable come from (Sanderson et al. (2021))[https://academic.oup.com/ije/article/50/4/1350/6133127]. These data are not publicly available; access was granted to us by the authors of that study, in accordance with ongoing UK Biobank applications.

The results of our analysis are contained in the file "Negative_Controls_Application.RData".
