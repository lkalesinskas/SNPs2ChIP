import pandas as pd

class Genomic_bins:
    def read_genomic_bin_bed_file(self, genomic_bin_bed_file):
        '''
        Read bed file that defines genomic bins
        '''
        genomic_bin_df=pd.read_csv(
            genomic_bin_bed_file,
            names=['chr', 'chromStart', 'chromEnd', 'name'],
            sep='\t'
        )
        return genomic_bin_df
    
    def __init__(self, genomic_bin_bed_file, bin_size = 1000):
        '''
        Read the definition of the genomic loci BED file and load to memory
        '''
        self.genomic_bin = self.read_genomic_bin_bed_file(genomic_bin_bed_file)
        self.name_to_loci = dict(zip(self.genomic_bin['name'], self.genomic_bin.index))
        self.bin_size = bin_size
    
    def find_locus_given_snp(self, chrom, position):
        '''
        Given a genomic coordinate a pair of (chrom, position),
        return the index of the corresponding genomic locus.
        '''
        name='chr{}_{}'.format(chrom, int(position / self.bin_size))
        return self.name_to_loci[name]

    def find_loci_given_snps(self, chroms, positions):
        '''
        Given a list of genomic coordinates 
        (one need to pass as two lists of chromosome and positions),
        return the list of genomic loci
        '''
        return [
            self.find_locus_given_snp(x[0], x[1]) for x in
            zip(chroms, positions)
        ]
