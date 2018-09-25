# SNPs2ChIP
## Latent Factors of ChIP-seq to infer functions of non-coding SNPs

### Publication
Please check our conference proceedings at PubMed. If you use SNPs2ChIP in your publication, please cite our proceedings.

### Installation

You can install our source code by cloning our Git repository and download the reference data.

```
$ git clone git@github.com:lkalesinskas/bigtranscription.git
$ cd bigtranscription
$ wget http://bejerano.stanford.edu/SNPs2ChIP/SNPs2ChIP_data.tar
$ tar xvf SNPs2ChIP_data.tar
```


### Getting started

<<This part needs to be updated once we have a wrapper>>

`src/validation_template_notebook.ipynb` has some examples.


### Description of the codebase

All codes are in relevant sub-directories

### `src`

This directory contains Python codes of the pipeline. If you are interested in using the software without updating the reference dataset, this directory should contain everything you need.

#### ``snp2chip.py``

This wrapper is used to query SNPS.

The following gets the top 10 (k) principal components and their respective scores for SNPS 1105,1124,1245,1233 on chr 12, 15,12 on chr 8, and SNP 12 on chr 10.
```
snps2chip.py -i 12:1105,1124,1245,1233 8:15,12 10:12 -k 10 -d "PATH/TO/DATA/"
```


#### `enrichment`

This sub-directory contains the code to apply genomic region enrichment analysis tool (GREAT) for latent components.

#### `validation_for_class_project_write-up`

This directory will be deleted once we prepare the camera-ready version of PSB paper.
