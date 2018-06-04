#!/bin/bash
set -beEuo pipefail

inBed=$1

ontologies="GOBiologicalProcess,GOCellularComponent,GOMolecularFunction,MGIPhenotype,MGIPhenoSingleKO,HumanPhenotypeOntology"

out_dir_name="GREAT"
out=$(dirname $(dirname $inBed))/${out_dir_name}/$(basename $(basename ${inBed} .gz) .bed)

if [ ! -d $(dirname $out) ] ; then mkdir -p $(dirname $out) ; fi


cat $inBed \
| grep -v "chr17_81195" \
| GREATER --date=20171028 --minAnnotCount=5 --maxAnnotCount=1000 --requiredTests=neither --ontologies=${ontologies} hg19 /dev/stdin $out

