# Enrichment analysis

## objective
- Given latent components, we would like to perform enrichment analysis to investigate the functional relevance of each of the component.

## directory structure
- `src`: source code
- `figs`: figure output
- `notebook`: Jupyter notebook to browse the data
- `private_data`: data directory. It is ignored by git version control (`.gitignore`).

## method description
1. Given the loci score files (PCA0.csv, etc), we will generate query BED files for GREAT enrichment analysis [cite]. This is done by `src/generate_GREAT_query.py`. The usage of the program will be shown if you type `$python src/generate_GREAT_query.py`.
1. Perform GREAT enrichment analysis. For each of the query BED files, we apply GREAT. Yosuke has an access to private commandline tool that perform the same functionality as the web version of GREAT (http://great.stanford.edu/). `src/GREAT_wrapper.sh` is the wrapper to run this analysis. We apply enrichment analysis for the following ontologies:
  - GOBiologicalProcess
  - GOCellularComponent
  - GOMolecularFunction
  - MGIPhenotype
  - MGIPhenoSingleKO
  - HumanPhenotypeOntology
1. Manual inspection of enrichment analysis.
  - `notebook/GREAT_enrichment_res.ipynb`
  
