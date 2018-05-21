_README_ = '''
-------------------------------------------------------------------------
generate_GREAT_query.py

Given a collection of csv file with importance score of genomic region and
the genomic loci definition file, generate bed files to perform 
GREAT enrichment analysis

Usage: python generate_GREAT_query.py -l loci.bed -i pca{}.csv -o pca{}.bed -n 50 -k 5000 

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2018/5/21
-------------------------------------------------------------------------
'''


import argparse, os, logging
import numpy as np
import pandas as pd
from functools import reduce

from logging.config import dictConfig
from logging import getLogger

dictConfig(dict(
    version = 1,
    formatters = {'f': {'format': '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}},
    handlers = {
        'h': {'class': 'logging.StreamHandler','formatter': 'f',
              'level': logging.DEBUG}},
    root = {'handlers': ['h'], 'level': logging.DEBUG,},
))

logger = getLogger('generate_GREAT_query')



def nlargest_with_ties(df, n, by='score'):
    if(n >= len(df)):
        return df
    else:
        thr=list(df.sort_values(by=by, ascending=False)[by])[n]
        return df[df[by] >= thr]
    
def generate_GREAT_query(loci_df, score_df, topk):
    '''
    Given the definition of loci and scores from decomposition, 
    generate a BED file (data frame) for GREAT genomic region enrichment analysis.
    [Input]
      - loci_df  : data frame of the loci definition bed file
      - score_df : data frame from decomposition
    [Output]
      - df: data frame
    '''
    df = nlargest_with_ties(score_df, topk).merge(
        loci_df,
        on='name',
    )
    df['name']=['{}:{:.3e}'.format(x[0], x[1]) for x in zip(df['name'], df['score'])]
    return df[['chrom', 'chromStart', 'chromEnd', 'name']]
    
    
def generate_GREAT_query_main(loci_bed, score_csv_template, nPCs, topk, out_bed_template):
    loci_df = pd.read_table(
        loci_bed,
        names=['chrom', 'chromStart', 'chromEnd', 'name']
    )
    for i in range(nPCs):
        if(not os.path.exists(os.path.dirname(out_bed_template.format(i)))):
            os.makedirs(os.path.dirname(out_bed_template.format(i)))        
        
        in_df=pd.read_csv(
            score_csv_template.format(i),
            names=['score']
        )
        in_df['name'] = np.arange(len(in_df))
        
        generate_GREAT_query(
            loci_df, 
            in_df,
            topk
        ).to_csv(
            out_bed_template.format(i), 
            sep='\t', index=False, header=False
        )    

def main():        
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )

    parser.add_argument('-i', metavar='i', required=True,
                        help='input csv file name template (use {} as a placeholder for PC index)')    
    parser.add_argument('-o', metavar='o', required=True,
                        help='output bed file name template (use {} as a placeholder for PC index)')
    parser.add_argument('-l', metavar='l', required=True,
                        help='the loci definition file')
    parser.add_argument('-k', metavar='k', default=5000,
                        help='the number of loci to keep (default: 5000)')
    parser.add_argument('-n', metavar='n', default=50,
                        help='the number of components (default: 50)')

    args = parser.parse_args()

    loci_bed = args.l
    score_csv_template = args.i 
    nPCs = int(args.n)
    topk = int(args.k)
    out_bed_template = args.o
    
    for key, value in zip(
        ['loci', 'in file', 'out_file', 'nPCs', 'topk'],
        [loci_bed, score_csv_template, out_bed_template, nPCs, topk],
    ):        
        logger_main.info('{}:\t {}'.format(key, value))
    
    generate_GREAT_query_main(loci_bed, score_csv_template, nPCs, topk, out_bed_template)
    
    
if __name__ == "__main__":
    main()     
    
    