import numpy as np
import pandas as pd
import argparse
import os
import json
import sys

from Decomposed_matrices import Decomposed_matrices
from Genomic_bins import Genomic_bins
from great import read_great_res_wrapper

# TODO: implement something for find X given Y
# TODO: figure out best output files
#-----------------------------------------
#   Parser
#-----------------------------------------
class StoreDictKeyPair(argparse.Action):
     def __init__(self, option_strings, dest, nargs=None, **kwargs):
         self._nargs = nargs
         super(StoreDictKeyPair, self).__init__(option_strings, dest, nargs=nargs, **kwargs)
     def __call__(self, parser, namespace, values, option_string=None):
         my_dict = {}
         try:
             for value in values:
                 chr,snps = value.split(":")
                 if int(chr) < 0 or int(chr) > 26:
                     raise ValueError("Please provide valid chormosome number.")
                 my_dict[int(chr)] = [int(snp.strip()) for snp in snps.split(",")]
         except:
            raise ValueError("Please provide input as chr1:snp1,snp2,snp3... etc.")
         setattr(namespace, self.dest, my_dict)

parser = argparse.ArgumentParser(prog='SNPs2ChIPs', description='A tool that uses latent factors of ChIP-seq to infer functions of non-coding SNPs.')

parser.add_argument('-i', '--input_snps', action=StoreDictKeyPair, nargs="+", metavar="KEY:VAL", default=None,
                    help='input dictionary of chr:SNPs of interest; ex. {chr1: snp1, snp2, snp3 .., chr6: snp4, ...}')

parser.add_argument('-a', '--input_assays', nargs="+", type=int, default=None,
                    help='assay indices')

parser.add_argument('-o', '--output_dir', type=str,
                    default='.', help='output directory')

parser.add_argument('-w', '--weights', type=float, nargs="+",
                    default=None, help='weights for each intput snp/assay')

parser.add_argument('-k', '--top_k_compoments', type=int,
                    default=5, help='top K principal components to find')

parser.add_argument('-d', '--data_dir', type=str,
                    default='../private_data',
                    help='data directory with genomic bins and decomposed matrices extracted form ChIP peaks')

parser.add_argument('-v', '--verbose', type=bool,
                    default=True, help='verbose run')

args = parser.parse_args()

#-----------------------------------------
#   Process/Check Inputs
#-----------------------------------------
if args.input_snps is None and args.input_assays is None:
    raise ValueError("Please provide valid assay IDs or SNP indices.")
if args.input_snps is not None and args.input_assays is not None:
    raise ValueError("Please provide either assays or SNPs, not both.")

def loc_to_list(input_dict):
    '''
    Take input dict of chr:snps and convert to lists.
    '''
    chrs = []
    snps = []
    for key in input_dict.keys():
        snps += input_dict[key]
        chrs += len(input_dict[key])*[key]
    return chrs,snps

def check_weights(inputs,weights):
    '''
    Ensure correct no. weights input.
    '''
    if len(weights)!=len(inputs):
        raise ValueError("Please ensure no. inputs ({}) matches no. weights ({})".format(len(inputs),len(weights)))

#-----------------------------------------
#   Build Genomic Bins and Decomposed Matrices
#-----------------------------------------
genomic_bins = Genomic_bins(os.path.join(args.data_dir, 'loci_def.bed'))
decomposed_mats = Decomposed_matrices(os.path.join(args.data_dir, 'diagonalScore.csv.gz'), os.path.join(args.data_dir, 'uScore.csv.gz'), os.path.join(args.data_dir, 'vScore.csv.gz'))

#-----------------------------------------
#   Output
#-----------------------------------------
if args.input_snps is not None:
    snps, chps = loc_to_list(args.input_snps)

    if args.weights is not None:
        check_weights(snps, args.weights)

    top_pcs, scores = decomposed_mats.find_pcs_given_loci_list(genomic_bins.find_loci_given_snps(snps,chps), args.weights, topk=args.top_k_compoments)
    print(top_pcs, scores)

else:
    # TODO: Needs to be tested
    if args.weights is not None:
        check_weights(args.input_assays, args.weights)

    top_pcs, scores = decomposed_mats.find_pcs_given_loci_list(genomic_bins.find_loci_given_snps(snps,chps), args.weights, topk=args.top_k_compoments)
    print(top_pcs, scores)
