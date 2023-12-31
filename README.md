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
c3.mir.mirdb.v7.4.symbols.gmt (file with miRNA targets) <br>
edited-nkim_precursor_miRNAseqs.csv (file with primary and precursor miRNA sequences) <br>
high_conf_mature_miRs.txt (list of high confidence miRNAs) <br>
miRmetProts.txt (list of proteins that interact with known miRNA biogenesis proteins) <br>
multiclass_nnetfit_mirprtints10042022.rds (final trained sRNA-Effector model) <br>

To run the example, the file shBMI1_GSE7578.top.table.tsv is required. This file can be replaced with any other file with the same format. Specifically, the input file should be a text file with gene names as the first column and log2 fold changes with the knockdown condition as the numerator and the control as the denominator as the second column. Here is a screenshot of a portion of the shBMI1_GSE7578.top.table.tsv file:

<img width="126" alt="image" src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/4cf16623-91e8-4d13-a884-105c7215ab2f">


Install with git (or download zipped file):
        
        git clone https://github.com/bw9bj/sRNA-Effector.git

        cd sRNA-Effector 


Example use:

        Rscript sRNAeffector_v1.0.R -g BMI1 -f shBMI1_GSE7578.top.table.tsv -d miRBase

Options:

        -g CHARACTER, --gene=CHARACTER
                Gene of interest

        -f CHARACTER, --file=CHARACTER
                Log2FC input file name

        -m CHARACTER, --miRNA=CHARACTER
                MiRNA of interest. Example format: miR-21-5p. Optional.

        -d CHARACTER, --miRdb=CHARACTER
                Database of high confidence miRNAs. If using option either input 'miRBase' or 'MirGeneDB'.

        -h, --help
                Show this help message and exit

sRNA-Effector has run successfully if the following output appears:
<img width="567" alt="image" src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/0695bc21-2ad9-4bff-b174-7f7a15995166">



Output includes a text file of predictions ("predictions.txt"), a stacked bar plot showing the proportion of miRNAs classified as either "BinderEffectors", "BinderNotEffectors", "NotBinderEffector", and "NotBinderNotEffector" ("inputGene_stackedbarplot.pdf"), and a ranked plot for predicted "BinderEffectors" ("BErankplot.pdf"). If the optional miRNA flag is passed, sRNA-Effector will also output a cumulative distribution function plot for the predicted targets of that miRNA.

Output from example use: <br>
<p align="center">
<img width="500" alt="image" src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/bd9f9385-faef-413b-8099-5b385429f0ff"> <img width="400" alt="image" src="https://github.com/bw9bj/sRNA-Effector/assets/42174020/e7a3f001-b936-40ba-8453-412db07b1d5e"> <br>
</p>







