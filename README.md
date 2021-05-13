[![DOI](https://zenodo.org/badge/360270100.svg)](https://zenodo.org/badge/latestdoi/360270100)

A repository for: 
================

Sánchez CA, Venkatachalam-Vaz J, Drake JM. 2021. Spillover of zoonotic pathogens: A review of reviews. Zoonoses & Public Health. doi: 10.1111/zph.12846


**This repository has 2 R scripts used for this analysis:** The first script 
(searchConcatenation.R) assembles, cleans, and removes duplicates from the results 
of a systematic search across multiple databases to identify review papers about 
zoonotic disease spillover. Authors CAS and JV-V screened the abstract and titles
of these papers according to exclusion criteria, followed by a full-text screening 
of papers retained after the initial screening. These processes are described in
detail in the manuscript. Qualitative codes were assigned to text snippets from 
review papers, then hierarchical cluster analysis was performed in R on these code 
data. The second script (analysis.R) performs all analyses and generates accompanying 
figures. Figures were then post-processed by Éric Marty to generate the published 
versions.
