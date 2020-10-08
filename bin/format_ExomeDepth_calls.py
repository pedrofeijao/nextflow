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

    parser = argparse.ArgumentParser(description="Format ExomeDepth calls.")
    # Inputs:
    parser.add_argument("exomedepth", type=str, help="ExomeDepth call file.")
    # Outputs:
    parser.add_argument("-o", "--output", type=str, help="Formatted output.")
    # Optional parameters
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    exd = pd.read_csv(param.exomedepth, sep="\t")

    # formatting and filtering:
    exd["CN"] = 2 * exd["Reads.ratio"]
    exd["Exons"] = exd.Gene
    exd["Gene"] = exd.apply(lambda r: "/".join(set(map(lambda s:re.sub("\_.+", "", s), r.Exons.split(", ")))), axis=1)
    exd["Sample"] = exd.Sample.map(lambda sample: sample.replace(".bam","").replace(".","-"))

    cols = ["Sample","CNVtype","CN", "Chr","Start","End","Exons","ExonN","Gene"]
    exd_fmt = exd[(exd["N.exons"]>2)][["Sample","CNV.type","CN","Chromosome","Start","End","Exons","N.exons","Gene"]].copy()
    exd_fmt.columns = cols
    exd_fmt["Caller"] = "ExomeDepth"
    exd_fmt.to_csv(param.output, sep="\t", index=False, header=True)
