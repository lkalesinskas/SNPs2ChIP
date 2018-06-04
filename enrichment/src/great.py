import os, pandas

def read_great_res_all(filename, nrows=None, BFold=None, BPval=None, topk=None, sort_by='BPval', compression=None):
    '''Read GREAT results file (tsv) file.
    '''
    if((compression is None) and (len(filename) > 3) and (filename[-3:] == '.gz')):
        compression='gzip'
    df = pandas.read_table(
        filename, sep='\t', skiprows=7, nrows=nrows, compression=compression
    )    
    
    trues = pandas.Series([True] * len(df))
    
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


def read_great_res(filename, nrows=None, BFold=None, BPval=None, topk=None, sort_by='BPval'):
    return read_great_res_all(filename, nrows, BFold, BPval, topk, sort_by)[['# ID', 'Desc', 'BPval', 'BFold']]


def read_great_res_wrapper(data_dir, pc, ontology):
    return read_great_res(
        os.path.join(
            data_dir, 'GREAT', '{}'.format(pc), '{}.tsv.gz'.format(ontology)
        ), 
        BFold=2.0
    )
