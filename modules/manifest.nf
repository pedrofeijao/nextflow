/****************************
** Manifest pre-processing **
*****************************/
process merge_overlapping_amplicons {
    /* merge overlapping amplicons into a single amplicon. */
    container "quay.io/biocontainers/pandas:0.24.1"
    input:
        path manifest
    output:
        path "${manifest}.merged.txt"
    script:
    """
    merge_overlap_amplicons.py $manifest ${manifest}.merged.txt
    """
}

process fix_overlapping_amplicons {
    /* fix overlapping amplicons, adjusting coordinates to they don't overlap. */
    container "quay.io/biocontainers/pandas:0.24.1"
    input:
        path manifest
    output:
        path "${manifest}.fixed.txt"
    script:
    """
    fix_overlap_amplicons.py $manifest ${manifest}.fixed.txt
    """
}

process manifest_to_bed {
    /*
    Reads the amplicon manifest and outputs a bed file with the coordinates, gene and amplicon ID for each amplicon.
    This file be used in some of the processes.
     */
     container "quay.io/biocontainers/pandas:0.24.1"
    input:
        file manifest
    output:
        file "${manifest}.bed"
    """
    #!/usr/bin/env python
    import pandas as pd
    manifest = pd.read_csv("$manifest", sep="\t")
    bed = manifest.sort_values(["Chr","Start"]).reset_index()
    bed.to_csv("${manifest}.bed",  columns=["Chr","Start","End","HGNC_gene","Amplicon_ID"],  header=False, sep="\t", index=False)
    """
}

process manifest_to_covcopcan_design {
    // Reads the amplicon manifest and outputs a "CovCopCan design" file, that is required for CovCopCan to run.
    container "quay.io/biocontainers/pandas:0.24.1"
    input:
        file manifest
    output:
        file "covcopcan_design.txt"
    """
    #!/usr/bin/env python
    import pandas as pd
    manifest = pd.read_csv("$manifest", sep="\t")
    copcon = manifest.sort_values(["Chr","Start"]).reset_index()
    copcon["Filter"] = copcon.apply(lambda r:"EXON" if "Exon" in r.Amplicon_name else "NOT", axis=1)
    copcon = copcon[["Pool","Chr","Start","End","Amplicon_ID","HGNC_gene","Filter","Amplicon"]]
    copcon.columns = ["#pool","chromosome","start","end","amplicon","gene","filter","sequence"]
    copcon.to_csv("covcopcan_design.txt", header=True, sep="\t", index=False)
    """
}
