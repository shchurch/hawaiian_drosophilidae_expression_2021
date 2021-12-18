# hawaiian_drosophilidae_expression_2021
This repostiory contains the data and code for a future manuscript on tissue-specific expression evolution in Hawaiian _Drosophila_ and _Scaptomyza_ flies.

There is an accompanying data visualization that can be found here:
https://shchurch.shinyapps.io/hawaiian_fly_dataviz_2021

and here:
https://github.com/shchurch/hawaiian_fly_dataviz_2021

## Contents

### Manuscript files

The main text and supplementary methods can be reproduced from the Rmarkdown files

```
Hawaiian_transcriptome_manuscript.Rmd
Hawaiian_transcriptome_supplement.Rmd
```

In addition, the supplementary methods Rmd file contains all the script calls required to reproduce all analyses.

### Analyses

All analyses script files are including in the directory `analysis`.
This includes input data files downloaded from https://flybase.org/ and other online repositories, under `analysis/data`.
It also includes scripts for each stage of the analysis in the following directories:

The code used to map and annotate RNA reads is located in the directories

```
analysis/agalma/
analysis/BLAST/
```

There are also several files containing collection and sequencing data in the directory ```expression_methods_data/```
Code to import the results of these analysis into R is located in the files

```
analysis/read_in_expression.R
```

Code to perform all comparative analyses are located in the directories

```
analysis/differential_expression/
analysis/ANOVA/
analysis/phylogenetic_expression/
analysis/head_analyses/
```

### Results

Figures are printed to the directory ```figures_and_panels/```

### Environment

Information on the R environment used to perform the comprative analyses is in the file

```
r_session_info.txt
```
