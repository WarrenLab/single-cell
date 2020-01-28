#!/usr/bin/env nextflow

/*
 * --------------- Parameters to edit start here. ------------------
 */

// Replace this with the path to a directory containing raw fastq files
params.fastqs_dir = '/path/to/fastqs'

// Replace this with the path to the directory containing the CellRanger
// reference. If you haven't made this yet, see the script
// 'make_cellranger_reference.sh'
params.ref_dir = '/path/to/cellranger/ref'

// Replace this with a list of the ID's of the samples to analyze
Channel.of('61C', '62C', '63C', '72C', '74C', '75C').set { ids }
/*
 * ---------------- Parameters to edit end here. -------------------
 *
 * Once you are done setting the parameters, you can run this pipeline
 * with 'nextflow run cellranger.nf'. A nextflow.config file with your
 * preferred settings should be in the directory where you run it.
 */


process cellranger_count {
    publishDir 'molecule_info'
    cpus 24
    memory '240 GB'

    input:
    val id from ids

    output:
    tuple val(id), file("molecule_info.${id}.h5") into molecule_info

    """
    cellranger count \
        --id=${id} \
        --fastqs=${params.fastqs_dir} \
        --sample=${id} \
        --transcriptome=${params.ref_dir} \
        --localcores=${task.cpus} \
        --localmem=240 \
        --disable-ui
    ln -s \$PWD/${id}/outs/molecule_info.h5 molecule_info.${id}.h5
    """
}

molecule_info.map { "${it[0]},${it[1]}" }
    .collectFile(name: "samples.csv",
                 newLine: true,
                 seed: "library_id,molecule_h5")
    .set { molecule_info_csv }

process cellranger_aggregate {
    publishDir 'aggregate'
    cpus 16

    input:
    file "libraries.csv" from molecule_info_csv

    output:
    file "aggregated/outs"

    """
    cellranger aggr --id=aggregated --csv=libraries.csv
    """
}
