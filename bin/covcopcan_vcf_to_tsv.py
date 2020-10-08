#!/usr/bin/env python3
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import re
import io
import pandas as pd
import numpy as np
from cnv_functions import open_exon_annotations, exons_in_range
# Read VCF as dataframe:
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def return_normalized_depth(read_depth, sample, chrom, start, end):
    depth = []
    for idx, row in read_depth[(read_depth.Chromosome == chrom)].iterrows():
        if start <= row.End and end >= row.Start:
            depth.append(row[sample])
    return depth


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reads all VCF for CovCopCan and output a tab-sep file with all calls.")
    # Inputs:
    parser.add_argument("vcf_files", nargs="+", type=str, help="VCF files with CovCopCan CNV calls.")
    parser.add_argument("-e", "--exons_file",  type=str, help="Bed file with all exon annotations to format the output.")
    parser.add_argument("-m", "--matrix_file",  type=str, help="Normalized matrix depth from CovCopCan.")
    # Outputs:
    parser.add_argument("-o", "--output", type=str, help="Tab-separated file with all CovCopCan calls, formatted with common columns between callers.")
    # Optional parameters
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    # Read results
    data = []
    for vcf_file in param.vcf_files:
        df = read_vcf(vcf_file)[["ALT","CHROM","INFO","POS"]]
        df["Run"] = vcf_file.split("/")[-3]
        df["Sample"] = vcf_file.split("/")[-1].replace(".vcf","")
        data.append(df)
    covcop = pd.concat(data)
    # add end coordinate:
    covcop["END"] = covcop.apply(lambda r: int(re.search("END=(\d+);", r.INFO).group(1)), axis=1) # Possibly change pattern depending on how samples are named;
    # exon annotations:
    exons = open_exon_annotations(param.exons_file)
    covcop["Exons"] = covcop.apply(lambda r: exons_in_range(exons, r.CHROM, r.POS, r.END), axis=1)
    covcop["ExonN"] = covcop.Exons.str.len()
    covcop["Gene"] = covcop.Exons.map(lambda exons: tuple(set(map(lambda s:re.sub("\_.+", "", s), exons))))

    # Split multiple gene calls:
    dup = covcop[covcop.Gene.str.len() > 1]
    single = covcop[covcop.Gene.str.len() == 1]
    rows = []
    # For each row, break it into the different genes:
    for idx, row in dup.iterrows():
        for gene in row.Gene:
            new_row = row.copy()
            new_row.Exons = [ex for ex in new_row.Exons if ex.startswith(gene)]
            new_row.Gene = (gene,) # has to be a tuple, for consistency
            rows.append(new_row)
    # merge all:
    covcop = pd.concat([single, pd.DataFrame(rows)])
    # tuples to single gene:
    covcop.Gene = covcop.Gene.map(lambda x: x[0])

    # Add estimated copy number from the matrix file:
    read_depth = pd.read_csv(param.matrix_file, sep="\t")
    covcop["Depth"] = 0 # Creates the column if dataframe is empty
    covcop["Depth"] = covcop.apply(lambda r: return_normalized_depth(read_depth, r.Sample, r.CHROM, r.POS, r.END), axis=1)
    covcop["MedianCN"] = covcop.Depth.map(lambda d: 2 * np.median(d))

    ## Formatted output:
    cols = ["Sample","CNVtype","CN", "Chr","Start","End","Exons","ExonN","Gene"]
    cc_fmt = covcop[["Sample","ALT","MedianCN","CHROM","POS","END","Exons","ExonN","Gene"]].sort_values(["Sample"])
    cc_fmt.columns = cols
    cc_fmt["CNVtype"].replace({"<DUP>": "duplication", "<DEL>": "deletion"}, inplace=True)
    cc_fmt["Caller"] = "CovCopCan"
    cc_fmt.to_csv(param.output, sep="\t", index=False, header=True)
