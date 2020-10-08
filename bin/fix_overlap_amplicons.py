#!/usr/bin/env python3
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import re
import pandas as pd
from functools import reduce

def agg_list(series):
    return reduce(lambda a,b: "{},{}".format(a,b), series)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Corrects coordinates of all overlapping amplicons, changing the start position so no overlap occurs.")
    # Inputs:
    parser.add_argument("manifest_input", type=str, help="Input amplicon manifest.")
    # Outputs:
    parser.add_argument("manifest_output", type=str, help="Output manifest, with corrected amplicon coordinates.")
    # Optional parameters
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    manifest = pd.read_csv(param.manifest_input, sep="\t").sort_values(['Chr','Start'])
    manifest["Exon"] = manifest.Amplicon_name.str.replace(r'(Exon[0-9]+).+',r'\1')
    # Group if they overlap:
    manifest['Group']=((manifest.End.rolling(window=2,min_periods=1).min()-manifest.Start.rolling(window=2,min_periods=1).max())<0).cumsum()

    previous_row = None
    row_list = []
    for idx, current_row in manifest.iterrows():
        if previous_row is None:
            previous_row = current_row
        elif previous_row.Chr == current_row.Chr and previous_row.End > current_row.Start:
            start = current_row.Start
            current_row.Start = previous_row.End + 1
            previous_row.End = start - 1
        row_list.append(current_row)
        previous_row = current_row

    corrected_manifest = pd.DataFrame(row_list)
    corrected_manifest = corrected_manifest[(corrected_manifest.Start < corrected_manifest.End)]
    print("Removed %d rows." % (len(manifest) - len(corrected_manifest)))
    corrected_manifest.to_csv(param.manifest_output, sep="\t", index=False)
