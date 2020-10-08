#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
- To run:

nextflow run cnv_pipeline.nf --run_id [RUN_ID]

- where run_id is the folder with the sequencer output (f.i., 190129_M03829_0188_000000000-C6M8H)
*/

/* Include required modules: */
include { merge_overlapping_amplicons; manifest_to_bed; manifest_to_covcopcan_design} from './modules/manifest'
include { get_read_depth_per_amplicon } from  './modules/mapping'
include {depth_to_covcopcan_matrix; covcopcan_cnv; covcopcan_format; exome_depth_cnv;
    format_exome_depth_output; panelcn_MOPS_cnv; format_MOPS_calls; aggregate_cnv_calls } from './modules/cnv'

// Set up the main input channel: BAM files
bam_files = Channel.fromPath(params.bam_files)
bai_files = Channel.fromPath(params.bai_files)

/*****************
** CNV Workflow **
*****************/
workflow  {
    // PreProcessing
    merged_manifest = merge_overlapping_amplicons(params.manifest_file)
    bed_file = manifest_to_bed(merged_manifest)

    read_depth = get_read_depth_per_amplicon(bam_files, bed_file)

    // CovCopCan
    covcopcan_design = manifest_to_covcopcan_design(merged_manifest)
    covcopcan_matrix = depth_to_covcopcan_matrix(read_depth.collect())
    covcopcan_output_folder = covcopcan_cnv(covcopcan_design, covcopcan_matrix)
    covcopcan_cnv_calls = covcopcan_format(covcopcan_output_folder)

    // ExomeDepth
    exome_depth_cnv(bam_files.collect(), bai_files.collect(), bed_file)
    exome_depth_cnv_calls = format_exome_depth_output(exome_depth_cnv.out.exomedepth_output)

    // MOPS
    panelcn_MOPS_cnv_calls = format_MOPS_calls(panelcn_MOPS_cnv(bam_files.collect(), bai_files.collect(), bed_file))

    // Aggregate calls:
    aggregate_cnv_calls(exome_depth_cnv_calls, panelcn_MOPS_cnv_calls, covcopcan_cnv_calls)

}
