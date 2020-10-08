#!/usr/bin/env python
import logging
logger = logging.getLogger()

import argparse
import sys
import os
import re
import pandas as pd
import numpy as np
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reads calls from all callers and creates an aggregated call file.")
    # Inputs:
    parser.add_argument("exomedepth", type=str, help="ExomeDepth call file.")
    parser.add_argument("mops", type=str, help="panelcn.MOPS call file.")
    parser.add_argument("covcopcan", type=str, help="CovCopCan call file.")
    # Outputs:
    parser.add_argument("-o", "--output", type=str, help="Aggreagated calls.")
    # Optional parameters
    parser.add_argument('-r','--remove_suffix',type=str, help='Regex for a suffix to be removed from the sample name, if given.')
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()
    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')


    # Read formatted CovCopCan:
    cc_fmt = pd.read_csv(param.covcopcan, sep="\t")

    # Read ExomeDepth:
    exd_fmt = pd.read_csv(param.exomedepth, sep="\t")

    # Format cnpanel.MOPS:
    mops_fmt = pd.read_csv(param.mops, sep="\t")

    # Merge all calls
    calls = pd.concat([cc_fmt,exd_fmt,mops_fmt]).sort_values(["Chr","Gene"])
    # Format sample names:
    if param.remove_suffix is not None:
        suffix_re = re.compile(param.remove_suffix)
        calls.Sample = calls.Sample.map(lambda sample: suffix_re.sub("", sample))


    # Agreggate calls:
    calls_agg = calls.groupby(['Sample','CNVtype','Chr','Gene'], as_index = False).agg({'Start': list, 'End': list, 'Exons':list, 'ExonN':tuple,
        'Caller':list, 'CN':list })
    calls_agg["Caller"] = calls_agg.apply(lambda r: set(r.Caller), axis=1)
    calls_agg["CallerN"] = calls_agg.Caller.str.len()
    calls_agg["ExonAvg"] = calls_agg.ExonN.map(lambda r: np.mean(r))
    calls_agg["MedianCN"] = calls_agg.CN.map(lambda r: np.median(r))

    calls_agg.to_csv(param.output, sep="\t", index=False, header=True)
