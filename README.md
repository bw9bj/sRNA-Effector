# sRNA-Effector: a tool to expedite discovery of novel small RNA regulators

sRNA-Effector is a random forest based algorithm trained on enhanced Crosslinking and Immunoprecipitation sequencing and knockdown RNA sequencing data. sRNA-Effector can accurately identify known miRNA biogenesis and effector proteins and has been used to identify novel miRNA effector proteins. sRNA-Effector takes as input a gene of interest and log2 fold change data from RNA sequencing or microarray data. Optionally, a miRNA of interest can also be used as input. Output includes a text file of predictions, a stacked bar plot showing the proportion of miRNAs classified as either "BinderEffectors", "BinderNotEffectors", "NotBinderEffector", and "NotBinderNotEffector", and a ranked plot for predicted "BinderEffectors". If the optional miRNA flag is passed, sRNA-Effector will also output a cumulative distribution function plot for the predicted targets of that miRNA.

sRNA-Effector requires the following R packages:
tidyverse,
caret,
ranger,
PepTools,
Peptides,
fgsea,
optparse,
biomaRt,
data.table,
dplyr,
readr,
ggrepel

Example use:

        Rscript sRNAeffector_v1.0.R -g BMI1 -f shBMI1_GSE7578.top.table.tsv

Options:

        -g CHARACTER, --gene=CHARACTER
                Gene of interest

        -f CHARACTER, --file=CHARACTER
                Log2FC input file name

        -m CHARACTER, --miRNA=CHARACTER
                MiRNA of interest. Example format: miR-21-5p

        -h, --help
                Show this help message and exit
