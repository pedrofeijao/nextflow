#!/usr/bin/env python
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import re
import pandas as pd
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Takes a file with read depth for all targets outputs a matrix file for CovCopCan.")
    # Inputs:
    parser.add_argument("depth_files", nargs="+", type=str, help="Depth file, obtained from the mapped bam file and a bed file, with bedtools.")
    # Outputs:
    parser.add_argument("-m", "--covcopcan_matrix", type=str, help="Matrix file output in CovCopCan format.")
    # Optional parameters
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    depth_data = []
    for depth_file in param.depth_files:
        depth = pd.read_csv(depth_file, sep="\t", header=None, names=["Chr","Start", "End", "Gene", "AmpliconID", "Depth", "L1", "L2", "L3"])
        depth["Sample"] = os.path.basename(re.sub("\.bam.+", "", depth_file))
        depth_data.append(depth)
    all_depth = pd.concat(depth_data)
    # Pivot for the matrix format:
    cov = all_depth.pivot_table(index=["Gene","AmpliconID"],columns="Sample",values="Depth").reset_index().rename(columns = {'AmpliconID':'Target'})
    cov.to_csv(param.covcopcan_matrix, header=True, sep="\t", index=False)
