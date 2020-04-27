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

// Replace this with a table of the ID's of the samples to analyze, plus any
// additional relevant information about each sample. The sample ID must be
// the first field, and the header should call this field 'library_id'. For
// example, the first few lines might look like this:
//
// library_id,treated
// C3,control
// C4,treatment
// C5,control
// C6,treatment
params.sample_sheet = 'samples.csv'

// path to the run_seurat.R script. If it's in your PATH, you can leave it as
// is; otherwise, set this to explicitly point to its location.
params.path_to_run_seurat = 'run_seurat.R'
/*
 * ---------------- Parameters to edit end here. -------------------
 *
 * Once you are done setting the parameters, you can run this pipeline
 * with 'nextflow run cellranger.nf'. A nextflow.config file with your
 * preferred settings should be in the directory where you run it.
 */

// read the sample sheet and make three channels from it:
Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header:true)
    .into { sampleSheet1; sampleSheet2; sampleSheet3 }
// one to extract the header from...
keys = sampleSheet1.first().keySet().value
// one to get a list of library ids from for the count process...
sampleSheet2.map { it.library_id }.set { ids }
// and one to use in the aggregate process
sampleSheet3.map { tuple(it.library_id, it) }.set { sampleSheetRows }

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

// use the sample sheet and the output of the count process to make
// a new sample sheet for the aggregate process
molecule_info.join(sampleSheetRows).map {
    it[2].remove('library_id')
    values = it[2].values().join(',')
    return [it[0], it[1], values].join(',')
}.collectFile(
    name: 'molecule_info.csv',
    newLine: true,
    seed: "library_id,molecule_h5," + keys.drop(1).join(',')
).set { molecule_info_csv }

process cellranger_aggregate {
    publishDir 'aggregate'
    cpus 16

    input:
    file "molecule_info.csv" from molecule_info_csv

    output:
    file "aggregated/outs" into aggregated

    """
    cellranger aggr --id=aggregated --csv=molecule_info.csv
    """
}

process seurat {
    publishDir 'seurat_out'

    input:
    file "aggregated/outs" from aggregated

    output:
    file 'seurat_out' into seurat_out

    """
    ${params.path_to_run_seurat} \
        --output-dir seurat_out \
        --aggregation aggregated/outs/aggregation.csv \
        aggregated/outs/filtered_feature_bc_matrix
        
    """
}

