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
# TODO: test assay indices
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

parser.add_argument('-o', '--output_file', type=str,
                    default='./snp2chip_out.txt', help='output directory')

parser.add_argument('-w', '--weights', type=float, nargs="+",
                    default=None, help='weights for each intput snp/assay')

parser.add_argument('-k', '--top_k_compoments', type=int,
                    default=5, help='top K principal components to find')

parser.add_argument('-d', '--data_dir', type=str,
                    default='../private_data',
                    help='data directory with genomic bins and decomposed matrices extracted form ChIP peaks')

parser.add_argument('-M', '--multi', action='store_true',
                    help='multi-run: multiple input SNPs are treated separately and not queried for summary statistics.')

parser.add_argument('--verbose', action='store_true',
                    help='verbose run')

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
if args.verbose:
    print('Loading data...')
genomic_bins = Genomic_bins(os.path.join(args.data_dir, 'loci_def.bed'))
decomposed_mats = Decomposed_matrices(os.path.join(args.data_dir, 'diagonalScore.csv.gz'), os.path.join(args.data_dir, 'uScore.csv.gz'), os.path.join(args.data_dir, 'vScore.csv.gz'))

#-----------------------------------------
#   Output
#-----------------------------------------
if args.verbose:
    print('Computing results...')

if args.input_snps is not None:
    # if user inputs SNPs
    chrs,snps = loc_to_list(args.input_snps)

    if args.weights is not None:
        check_weights(snps, args.weights)

    if args.multi:
        if args.weights is not None and args.verbose:
            print("Ignoring weights, multi-run mode in use.")

        if args.verbose:
            print("-------------------")
        for idx, snp in enumerate(snps):
            if args.verbose:
                print("Outputting PCs for chr{}:{}".format(chrs[idx], snp))
            top_pcs, scores = decomposed_mats.find_pcs_given_locus(genomic_bins.find_locus_given_snp(chrs[idx], snp), topk=args.top_k_compoments)
            with open(args.output_file, 'a+') as f:
                f.write('{}:{},{},{},{}\n'.format(chrs[idx], snp, args.top_k_compoments, top_pcs, scores))
        if args.verbose:
            print("-------------------")
    else:
        top_pcs, scores = decomposed_mats.find_pcs_given_loci_list(genomic_bins.find_loci_given_snps(chrs,snps), args.weights, topk=args.top_k_compoments)
        with open(args.output_file, 'a+') as f:
            f.write('{},{},{},{},{}\n'.format(chrs, snps, args.top_k_compoments, top_pcs, scores, args.weights))
else:
    # if user inputs assay indices
    if args.weights is not None:
        check_weights(args.input_assays, args.weights)
    top_pcs, scores = decomposed_mats.find_pcs_given_loci_list(genomic_bins.find_loci_given_snps(snps,chps), args.weights, topk=args.top_k_compoments)
    print(top_pcs, scores)

if args.verbose:
    print("____COMPLETE____")
