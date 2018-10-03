# SNPs2ChIP
## Latent Factors of ChIP-seq to infer functions of non-coding SNPs

### Publication
Please check our conference proceedings at PubMed. If you use SNPs2ChIP in your publication, please cite our proceedings.

### Installation

You can install our source code by cloning our Git repository and download the reference data.

```
$ git clone git@github.com:lkalesinskas/bigtranscription.git
$ cd bigtranscription
$ wget http://cs.stanford.edu/people/ytanigaw/SNPs2ChIP_GREAT_data.tar
$ tar xvf SNPs2ChIP_GREAT_data.tar
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

##### Input
Each input (``-i``) takes in a chromosome:SNP-location pair. One can query multiple SNPs from the same chromosome as such: ``2:132000,142000``: this queries at the SNP locations 132000 and 142000 at chromosome 2. This queries for 3 (``-k``) components and runs on multi-run mode, toggled by the ``-M``, which runs each SNP independently and outputs results for each. The data directory (``-d``) is also required but will default to ``../private_data``.

##### Multi vs. Summary Mode
By default, if multiple SNPs are given to ``snp1chip.py``, it will compute summary values over all SNPs running in __summary mode__. One may also provide weights (``-w``) if running on default summary mode, but it must be a list of the same length as the number of SNPs or you'll receive an error. Toggle __multi-mode__ by the ``-M`` flag and each SNP will be run independently and PCs/scores will be given as such.

##### Output
The output file (``-o``) is the path to a text file where results will be appended to or the file will be created if it does not exist.

* Summary mode format: ``chromosomes, SNPs, k-components, PCs, scores, weights``
  * ex. ``./src/vdr_snps_summary.txt``
* Multi-mode output format: ``chr:SNP, k-components, PCs, scores``
  * ex. ``./src/vdr_snps_out.txt``

#### `enrichment`

This sub-directory contains the code to apply genomic region enrichment analysis tool (GREAT) for latent components.

#### `validation_for_class_project_write-up`

This directory will be deleted once we prepare the camera-ready version of PSB paper.
