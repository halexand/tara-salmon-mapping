# tara-salmon-mapping

Snakemake framework for mapping all salmon metaG and metaT reads against a set of transcripts. 

Pulls and creates salmon indicies for all transcriptome files located in `transcripts/` that end with `.fasta`. It then will map all metaT and metaG data from Tara located in directory indidcated in the `config.yaml` file. 

To use: create a conda environment and install `snakemake`. 

```
conda create -n snakemake -c bioconda -c conda-forge snakemake
conda activate snakemake
```
