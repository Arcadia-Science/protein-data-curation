# Downloading and processing proteomes and transcriptomes
### This repository contains a Snakemake workflow to complete the following:
#### - Download a proteome or transcriptome file using a URL link
#### - Combine proteomes/transcriptomes if multiple are present for the same species
#### - For assembled transcriptomes, run [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki) to predict coding regions within transcripts
#### - Format files to be compatible with [NovelTree](https://github.com/Arcadia-Science/noveltree)
#### - Annotate proteins with [Eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper), [KofamScan](https://github.com/takaram/kofam_scan)
#### - Predict signal peptides with [DeepSig](https://github.com/BolognaBiocomp/deepsig)

* Note: this pipeline was successfully run on Ubuntu 22.04. It currently does not work on MacOS because of incompatibility with DeepSig requirements 

### Set-up
#### First create a conda environment using the ```envs/environment.yml``` file and activate your new environment
This repository uses snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html). After installing conda, run the following commands to create and activate an environment for running this pipeline:
```
conda env create -f envs/environment.yml -n curate
conda activate curate
```

#### Choose what you want to download
The ```inputs/samples.tsv``` file contains all of the necessary information for the samples that will go through this workflow. If you want to run this workflow on a different set of files, you will need to edit this table. Each transcriptome/proteome should have a single row in the table. These four columns are required for running the pipeline:

```link```: Link to the transcriptome or proteome fasta file that will be downloaded and used as input in the pipeline. For example, https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000001610/UP000001610_983644.fasta.gz <br>
```sample```: The name of the sample. There should be a unique entry for each link/data file pulled. <br>
```genus_species```: The genus and species of the organism associated with the sample, separated by '_'. The final protein fasta will contain this information in the headers and this will be used to name output files.<br>
```transdecoder```: Should be 'yes' or 'no'. If yes, transdecoder will be run on the input file for this genus_species.<br>
```isoform```: Should be 'yes' or 'no'. If yes, there is isoform information available in the header of the fasta file and isoform filtering will be performed to only keep the longest isoform.<br>
```uniprot```: Whether or not the proteome is sourced from Uniprot. Should be 'yes' or 'no'. This will change how proteins are clustered in the NovelTree processing step. <br>

The pipeline is set up so that if a single species has multiple input transcriptomes/proteomes, these will be combined prior to any other processing. Downloaded files that have the same information in the genus_species column of the input table will be processed as a single sample.


#### Run Snakemake from your activated environment
```
snakemake -p -j 1 --use-conda
```

```-p``` tells Snakemake to print the shell commands that it runs for each rule, which can be useful for understanding what is happening/debugging.<br>
```-j``` tells Snakemake how many cores to use for executing jobs in parallel.<br>
```--use-conda``` tells Snakemake to utilize Conda environments to manage dependencies for individual rules.<br>
The ```--rerun-incomplete``` flag can be appended to the above command if the pipeline terminates prematurely for some reason and needs to be restarted.<br>

