<p align="center">
<img src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/0482c749-7692-4f0b-9233-dbeea0dbe991" width="400" />
</p>


# sRNA-Effector: a tool to expedite discovery of novel small RNA regulators

sRNA-Effector is a random forest based algorithm trained on enhanced Crosslinking and Immunoprecipitation sequencing and knockdown RNA sequencing data. sRNA-Effector can accurately identify known miRNA biogenesis and effector proteins and has been used to identify novel miRNA effector proteins. sRNA-Effector takes as input a gene of interest and log2 fold change data from RNA sequencing or microarray data. Optionally, a miRNA of interest can also be used as input. 

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

The following files are also required to run sRNA-Effector: <br>
c3.mir.mirdb.v7.4.symbols.gmt <br>
edited-nkim_precursor_miRNAseqs.csv <br>
high_conf_mature_miRs.txt <br>
miRmetProts.txt <br>
multiclass_nnetfit_mirprtints10042022.rds <br>

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


Output includes a text file of predictions ("predictions.txt"), a stacked bar plot showing the proportion of miRNAs classified as either "BinderEffectors", "BinderNotEffectors", "NotBinderEffector", and "NotBinderNotEffector" ("inputGene_stackedbarplot.pdf"), and a ranked plot for predicted "BinderEffectors" ("BErankplot.pdf"). If the optional miRNA flag is passed, sRNA-Effector will also output a cumulative distribution function plot for the predicted targets of that miRNA.

Output from example use: <br>
<p align="center">

<img src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/08b8475b-5107-4801-96e1-4237bed63188" width="500" /> <img src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/e9b42434-1265-467d-94d5-775780796457" width="400" /> <br>
</p>







