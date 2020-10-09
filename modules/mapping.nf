/******************
** Read Mapping ***
******************/
process map_reads {
    // Maps the FastQ files to the reference genome, then filters for quality > 30 and indexes the .bam files.
    cpus 5
    input:
        tuple val(sampleId), path(fastq)
    output:
        path "${sampleId}.bam", emit: bam_file
        path "${sampleId}.bam.bai", emit: bai_file
    """
    bwa mem -T 50 -t 4 $params.reference_genome ${fastq[0]} ${fastq[1]} -o bwa.sam
    samtools view -F 256 -q 30 -h -bS bwa.sam | samtools sort - -o ${sampleId}.bam
    samtools index ${sampleId}.bam ${sampleId}.bam.bai
    """
}


process get_read_depth_per_amplicon {
    /*
    Here we use bedtools to count the number of reads per amplicon.
    */
    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'
    input:
        path bam_file
        path bed_file
    output:
        path "${bam_file}.depth.txt"
    """
    bedtools coverage -a $bed_file -b $bam_file > ${bam_file}.depth.txt
    """
}
