# Enrichment analysis

Yosuke Tanigawa (ytanigaw@stanford.edu)

## objective
- Given latent components, we would like to perform enrichment analysis to investigate the functional relevance of each of the component.

## directory structure
- `src`: source code
- `figs`: figure output
- `notebook`: Jupyter notebook to browse the data
- `private_data`: data directory. It is ignored by git version control (`.gitignore`).

## method description
1. Given the results of decomposition (`loadingSquared.csv`), we will generate query BED files for GREAT enrichment analysis [cite]. This is done by `notebook/GREAT_query_generation.ipynb`. 
1. Perform GREAT enrichment analysis. For each of the query BED files, we apply GREAT. Yosuke has an access to private commandline tool that perform the same functionality as the web version of GREAT (http://great.stanford.edu/). `src/GREAT_wrapper.sh` is the wrapper to run this analysis through the tools in Bejerano lab (it will produce the same results as the ones through the web version). We apply enrichment analysis for the following ontologies:

  - GOBiologicalProcess
  - GOCellularComponent
  - GOMolecularFunction
  - MGIPhenotype
  - MGIPhenoSingleKO
  - HumanPhenotypeOntology

1. Manual inspection of enrichment analysis.
  - The results of the enrichment analysis is uploaded to Google Drive:
    - https://drive.google.com/file/d/1j2F0w_r9Y-cRdAh0UrDy-QZl4qWndK8W/view?usp=sharing
  - `notebook/GREAT_enrichment_res.ipynb`
  

