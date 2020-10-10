#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
- To run:
nextflow run cnv_pipeline.nf --run_id [RUN_ID]
- where run_id is the folder with the sequencer output (f.i., 190129_M03829_0188_000000000-C6M8H)
*/

/* Include required modules: */
include { fix_overlapping_amplicons; manifest_to_bed; manifest_to_covcopcan_design} from './modules/manifest'
include { get_read_depth_per_amplicon } from  './modules/mapping'
include {depth_to_covcopcan_matrix; covcopcan_cnv; covcopcan_format; exome_depth_cnv;
    format_exome_depth_output; panelcn_MOPS_cnv; format_MOPS_calls; aggregate_cnv_calls } from './modules/cnv'


/*****************
** CNV Workflow **
*****************/
workflow  {
    // Set up the main input channel: BAM files
    bam_files = Channel.fromPath(params.bam_files).filter( ~/.+(DNA|NF).+/)
    bai_files = Channel.fromPath(params.bai_files).filter( ~/.+(DNA|NF).+/)
    reference_genome = Channel.fromPath(params.reference_genome)
    reference_genome_idx = Channel.fromPath(params.reference_genome_idx)
    // bam_files.view()


    // // PreProcessing
    manifest_file = "$workflow.projectDir/$params.manifest_file"
    exome_bed = "$workflow.projectDir/$params.exome_bed"
    //
    fixed_manifest = fix_overlapping_amplicons(manifest_file)
    bed_file = manifest_to_bed(fixed_manifest)

    // Read depths from BAM files
    read_depth = get_read_depth_per_amplicon(bam_files, bed_file)

    // CovCopCan
    // covcopcan_jar = "$workflow.projectDir/$params.covcopcan_jar"
    // covcopcan_design = manifest_to_covcopcan_design(fixed_manifest)
    // covcopcan_matrix = depth_to_covcopcan_matrix(read_depth.collect())
    // covcopcan_output_folder = covcopcan_cnv(covcopcan_jar, covcopcan_design, covcopcan_matrix)
    // covcopcan_cnv_calls = covcopcan_format(exome_bed, covcopcan_output_folder)

    // ExomeDepth
    exomedepth_path = "$workflow.projectDir/$params.exomedepth_path"
    exome_depth_cnv(exomedepth_path, reference_genome, reference_genome_idx, bam_files.collect(), bai_files.collect(), bed_file)
    exome_depth_cnv_calls = format_exome_depth_output(exome_depth_cnv.out.exomedepth_output)

    // MOPS
    panelcn_MOPS_path = "$workflow.projectDir/$params.panelcn_MOPS_path"
    panelcn_MOPS_cnv_calls = format_MOPS_calls(exome_bed, panelcn_MOPS_cnv(panelcn_MOPS_path, bam_files.collect(), bai_files.collect(), bed_file))

    // Aggregate calls:
    // aggregate_cnv_calls(exome_depth_cnv_calls, panelcn_MOPS_cnv_calls, covcopcan_cnv_calls)

}
