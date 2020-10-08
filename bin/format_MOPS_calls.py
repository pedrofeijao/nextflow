#!/usr/bin/env python
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import re
import pandas as pd
from cnv_functions import open_exon_annotations, exons_in_range
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Format panelcn.MOPS calls.")
    # Inputs:
    parser.add_argument("mops", type=str, help="panelcn.MOPS call file.")
    parser.add_argument("-e", "--exons_file",  type=str, help="Bed file with all exon annotations to format the output.")
    # Outputs:
    parser.add_argument("-o", "--output", type=str, help="Formatted output.")
    # Optional parameters
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    # Format sample name
    mops = pd.read_csv(param.mops, sep="\t")
    mops.Sample = mops.apply(lambda r: r.Sample.replace(".bam",""), axis=1)
    mops = mops[~(mops.CN.isin(["CN2"]))]
    mops["CNVtype"] = mops.apply(lambda r:"deletion" if r.CN=="CN1" else "duplication", axis=1)
    # add chr:
    mops.Chr = 'chr' + mops.Chr.astype(str)
    # exon annotations:
    exons = open_exon_annotations(param.exons_file)
    mops["Exons"] = mops.apply(lambda r: exons_in_range(exons, r.Chr, r.Start, r.End), axis=1)

    # Simple Sample N:
    mops["SampleN"] = mops.apply(lambda r: int(re.search("_S(\d+)", r.Sample).group(1)), axis=1)
    # At least 100 reads?
    mops = mops[(mops.RC > 100)].copy()

    # Agreggate calls on the same gene as one:
    mops_agg = mops.groupby(['Sample','Gene','Chr','CNVtype'], as_index = False).agg({'Exons': list, 'Start': min, 'End': max, 'RC':list, 'medRC':list, 'RC.norm':list,
        'medRC.norm':list, 'CN':list})
    mops_agg["ExonN"] = mops_agg.apply(lambda r: len(r.Exons), axis=1)

    # Filter ExonN > 2
    mops_f = mops_agg[(mops_agg.ExonN > 2)].copy()

    mops_f["CNList"] = mops_f.CN.map(lambda cn: [int(c[2:]) for c in cn])
    mops_f["MedianCN"] = mops_f.CNList.map(lambda cn: np.median(cn))

    cols = ["Sample","CNVtype","CN", "Chr","Start","End","Exons","ExonN","Gene"]
    mopsmt = mops_f[["Sample","CNVtype","MedianCN","Chr","Start","End","Exons","ExonN","Gene"]].copy()
    mopsmt.columns = cols
    mopsmt["Caller"] = "MOPS"

    mopsmt.to_csv(param.output, sep="\t", index=False, header=True)
