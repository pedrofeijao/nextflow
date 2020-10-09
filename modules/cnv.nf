// ExomeDepth
params.trans_prob = 0.15
params.exomedepth_output = "exome_depth_cnv.txt"
params.exomedepth_fmt_output = "exome_depth_cnv_fmt.txt"
params.exomedepth_count_output = "exome_depth_counts.txt"

// panelcn.MOPS
params.panelcn_MOPS_output = "panelcn_MOPS_cnv.csv"
params.panelcn_MOPS_fmt_output = "panelcn_MOPS_cnv_fmt.csv"

// CovCopCan
params.covcopcan_matrix_file = "covcopcan_matrix.txt"
params.referenceAmpliconNumber = 45


/**************
** CovCopCan **
***************/
process depth_to_covcopcan_matrix {
    // CovCopCan does not use bam files, it needs a "matrix file" with all the read counts for all samples and all amplicons.
    // Here we get the depth from the "get_read_depth_per_amplicon" process and create the matrix file for CovCopCan.
    input:
        path bam_depths
    output:
        path params.covcopcan_matrix_file
    """
    depth_to_covcopcan_matrix.py $bam_depths -m $params.covcopcan_matrix_file
    """
}

process covcopcan_cnv {
    /*
    Main CovCopCan CNV caller step. Uses a design file (similar to the amplicon manifest) and a matrix file (with read depth per amplicon per sample)
    to make CNV calls.
    */
    input:
        path covcopcan_jar
        path covcopcan_design
        path covcopcan_matrix
    output:
        path "covcopcan_output"
    """
    xvfb-run --auto-servernum --server-num=1 java -Dprism.order=sw -jar $covcopcan_jar -g \\
    -d $covcopcan_design -m $covcopcan_matrix -o covcopcan_output --gcCorrection false --ampLenCorrection false \\
    --minCNVLength 2  --exportRawData false --referenceAmpliconNumber $params.referenceAmpliconNumber \\
    --deviationFromAverage 2 --zScoreDetection true
    """
}
process covcopcan_format {
    // Format output VCFs by CovCopCan into a .txt (tab-separated) file. This is the final CNV output for CovCopCan.
    input:
        path exome_bed
        path covcopcan_output_folder
    output:
        path "covcopcan_cnv.txt"
    """
    covcopcan_vcf_to_tsv.py $covcopcan_output_folder/VCF/*vcf -e $exome_bed -m $covcopcan_output_folder/covcopcan_matrix_normalized_data.tsv -o covcopcan_cnv.txt
    """

}
/***************
** ExomeDepth **
****************/
process exome_depth_cnv {
    /*
    Main ExomeDepth step. An R script takes all bam files and calls CNVs, outputting a .txt (tab-separated) file with all calls.
    Also outputs a file with read counts per exon (optional file, might be good for troubleshooting)
    */
    input:
        path  exomedepth_path
        path  reference_genome
        path  bam_files
        path  bai_files
        path  bed_file
    output:
        path params.exomedepth_output, emit: exomedepth_output
        path params.exomedepth_count_output, emit: exomedepth_count_output
    """
    # Rscript $exomedepth_path --bamdir . --ref_genome $reference_genome --trans_prob $params.trans_prob --out $params.exomedepth_output --counts_out $params.exomedepth_count_output --bed_file $bed_file
    Rscript $exomedepth_path --bamdir . --ref_genome $reference_genome --trans_prob $params.trans_prob --out $params.exomedepth_output --counts_out $params.exomedepth_count_output
    """
}

process format_exome_depth_output {
    input:
        path exomedepth_output
    output:
        path params.exomedepth_fmt_output
    """
    format_ExomeDepth_calls.py $exomedepth_output -o $params.exomedepth_fmt_output
    """

}
/*****************
** panelcn.MOPS **
******************/
process panelcn_MOPS_cnv {
    /*
    R script with the panelcn.MOPS caller. Uses all .bam files and outputs a .txt file (tab-separated) with the CNV calls.
    */
    input:
        path  panelcn_MOPS_path
        path  all_bam_files
        path  all_bai_files
        path  bed_file
    output:
        path params.panelcn_MOPS_output
    """
    Rscript $panelcn_MOPS_path --bamdir . --ref_genome $params.reference_genome --out $params.panelcn_MOPS_output --bedfile $bed_file
    """
}

process format_MOPS_calls {
    /*
    panelcn.MOPS calls need some formatting. This is done here with a python script.
    */
    input:
        path exome_bed
        path panelcn_MOPS_output
    output:
        path params.panelcn_MOPS_fmt_output
    """
        format_MOPS_calls.py $panelcn_MOPS_output -e $exome_bed -o $params.panelcn_MOPS_fmt_output
    """
}
/*****************
** Final Output **
******************/
process aggregate_cnv_calls {
    // Takes the results from all CNV callers and produces a file with combined calls.
    input:
        path exomedepth_fmt_output
        path panelcn_MOPS_fmt_output
        path covcopcan_cnv_formatted_output
    output:
        path "aggregated_cnv_calls.txt"
    """
    aggregate_calls.py $exomedepth_fmt_output $panelcn_MOPS_fmt_output $covcopcan_cnv_formatted_output -o aggregated_cnv_calls.txt -r ".primers.hardclipped"
    """
}
