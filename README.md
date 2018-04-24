# bigtranscription
Let's make transcription BIG!

# Code Files
**createFakeBED.py:** Creates fake BED file with bins for a certain length

**intersectRunner.py:** Runs intersect script on a folder of files to find intersect between filtered files and a chromosome/bin file

**downloadCHiPAtlas.py:** Download script to grab .BED files from the CHiP-Atlas database

**filterTopHits.py:** Filtration Script to filter .BED files to make relevant files. 
Options: 
1) chromosomeFilter: Gets all peaks corresponding to a certain chromosome
2) scoreFilter: Gets top k peaks from .bed file

**matrixMaker.py:** Makes the matrix for PCA/SVD/whatever. TODO: Get rid of blank stuff...


# Example Files
**ERX329611.bed:** Example .bed file from CHiP-Atlas

**filtered_chr_ERX329611.bed:** Chromosome filtered .bed file

**filtered_score_ERX329611.bed:** Score filtered .bed file

**final_filtered_chr_ERX329611.bed:** Filtered and intersected BED file

**chr19.bed:** Chromosome 19 file for BED intersection

# CHiP-Atlas Files
**chip_atlas_celltype_list:** List of IDs assosciated with different cell lines. Currently using Lymphoblastoid cell line
