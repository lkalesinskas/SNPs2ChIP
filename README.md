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

### ``snp2chip.py``

This wrapper is used to query SNPS.

The following is an example used to query principal components for the Vitamin D Receptor (VDR).

```
python snp2chip.py -i 10:6390000 13:42952000 2:33701000 5:1315000 -k 3 --verbose -M -o vdr_snps_out.txt -d PATH/TO/DATA/DIR
```

Each input (``-i``) takes in a chromosome:SNP-location pair. One can query multiple SNPs from the same chromosome as such: ``2:132000,142000``: this queries at the SNP locations 132000 and 142000 at chromosome 2. This queries for 3 (``-k``) components and runs on multi-run mode, toggled by the ``-M``, which runs each SNP independently. If multiple SNPs are given to ``snp1chip.py``, it will compute summary values over all SNPs. One may also provide weights (``-w``) if running on summary mode (default). The output file (``-o``) is the path to a text file where results will be appended to or the file will be created if it does not exist. The data directory (``-d``) is also required but will default to ``../private_data``.

See ``vdr_snps_out.txt`` for the output from multi-run mode. See ``vdr_snps_summary.txt`` for output from summary (default) mode.

#### `enrichment`

This sub-directory contains the code to apply genomic region enrichment analysis tool (GREAT) for latent components.

#### `validation_for_class_project_write-up`

This directory will be deleted once we prepare the camera-ready version of PSB paper.
