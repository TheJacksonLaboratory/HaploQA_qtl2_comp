# HaploQA_qtl2_comp
Input data received from HaploQA - https://haploqa.jax.org/, software platform that provides interpretation of Mouse Universal Genotyping array platforms (MUGA) genotype calls available from Neogen Genomics (Morgan et al., 2015).

Computation implementations based on R/qtl2 - Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen Ś, Yandell BS, Churchill GA (2018) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multi-parent populations. Genetics 211:495-502

# Disclaimer
This pipeline is still in draft/dev state as of EOD 7/18/2022, especially the test_script branch...

# Overview

HaploQA_qtl2_output.R - Launch script
web_scrape_functions.R - script with functions

haploqa_collab_cross - input data of individuals, output of 'sample_individual_scrape' function
cc_qtl2_test - input data of qtl2, output of 'get_qtl2_input' function

annotations_config.csv, founder_strains_table.csv - dictionaries used for mapping

# Usage
The complete pipeline from HaploQA data retrieval to qtl2 computations is in HaploQA_qtl2_output.R, change the toggles on top of the script to select whether to regenerate files and run the entire script.
