import pandas as pd
import os


def read_great_res_all(filename, nrows=None, BFold=None, BPval=None, topk=None, sort_by='BPval', ascending=False, compression=None):
    '''
    This function reads GREAT (Genomic Region Enrichment Analysis Tool) 
    results file (tsv or tsv.gz) file. 
    
    - sort_by, ascending: By default, the results are sorted by 
                          binomial p-value (BPval) column with increasing order.
    - BFold, BPvalu: One may impose filter based on binomial fold and/or 
                     binomial p-value to mimic the default setting in GREAT web app [1].

    - nrows: Read the GREAT results file up to the given number of rows.
             The binomial fold and/or p-value filter will be applied after 
             this truncation. This is useful when you want to speed up the data load.
             
    - topk:  Return the top k hits after imposing the filters.
    
    For more info, please check GREAT web app and GREAT paper.
    [1]: http://great.stanford.edu
    [2]: https://www.ncbi.nlm.nih.gov/pubmed/20436461    
    '''
    
    if((compression is None) and (len(filename) > 3) and (filename[-3:] == '.gz')):
        compression='gzip'
    df = pd.read_table(
        filename, sep='\t', skiprows=7, nrows=nrows, compression=compression
    )    
    
    trues = pd.Series([True] * len(df))
    
    if (BFold is not None):
        f_BFold = df['BFold'] >= BFold
    else:
        f_BFold = trues
    
    if (BPval is not None):
        f_BPval = df['BPval'] <= BPval
    else:
        f_BPval = trues
        
    df_filtered = df[f_BFold & f_BPval]
    
    if (topk is None):
        topk = len(df_filtered)
    
    if (sort_by == 'BFold'):
        ascending=False
    else:
        ascending=True
        
    return df_filtered.sort_values(by=sort_by, ascending=ascending).head(topk)

def read_great_res_wrapper(data_dir, pc, ontology, BFold=2.0, **argv):
    '''
    Wrapper function to read GREAT results file for SNPs2ChIP.
    By default, we impose binomial fold >= 2.0 filter.
    This is the default behavior of the GREAT web app.
    '''
    
    def read_great_res(
        filename, nrows=None, BFold=None, BPval=None, topk=None, sort_by='BPval'
    ):
        '''
        We only care about binomial fold, p-value, and FDR in our application.
        '''
        return read_great_res_all(
            filename, nrows, BFold, BPval, topk, sort_by
        )[
            ['# ID', 'Desc', 'BFold', 'BPval', 'BFDR']
        ]
    
    return read_great_res(
        os.path.join(
            data_dir, 'GREAT', '{}'.format(pc), '{}.tsv.gz'.format(ontology)
        ), 
        BFold=BFold, 
        **argv
    )
