# single-cell
Pipelines and scripts for single-cell best practices

## Installation
You'll need the following software:
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
* [Seurat](https://satijalab.org/seurat/install.html)

If you want to use nextflow to run the Cell Ranger pipeline (recommended),
you'll also need to install nextflow with the following command if you haven't
already:
```
curl -s https://get.nextflow.io | bash
```
This creates a binary named 'nextflow' in your current directory; you probably
want to move it to somewhere in your path.

### Nextflow on lewis
Nextflow needs to be run from `htc` due to file-locking issues, but is much
faster when it's doing file operations on `hpc`. To get around this, you can set
it to always use a work directory on `hpc` by adding something like this to your
`.bashrc` (changing the location, of course):
```bash
export NXF_WORK=/storage/hpc/group/warrenlab/users/esrbhb/nx_work
```
A good nextflow configuration for lewis is in `nextflow.config`. Copy this file
to the directory you are running nextflow from.

## Cell Ranger
### Creating a reference
Cell Ranger needs reference indices in its own special format. If you are
running it on a reference genome you haven't run it on before, you'll need to
create those indices. Check whether there are already is already a reference
available, though. For example, there is a Cell Ranger chicken reference at
```
/storage/hpc/group/warrenlab/reference_genomes/GalGal5/GalGal5_cellranger
```
The commands I used to generate this reference can be found in this repository
in `make_cellranger_reference.sh`. Take a look at the
[https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references](official documentation)
on making a Cell Ranger reference for more information. As recommended in the
documentation, this example script keeps only protein-coding genes (including
those in mtDNA) and lncRNAs, discarding other types of features like rRNAs that
could confound the analysis.

### Creating a sample sheet for comparison experiments
(in progress)

### Running with nextflow
The whole Cell Ranger pipeline can be run for all your samples with a single
nextflow command. Just copy the files `cellranger.nf` and `nextflow.config` from
this repository into your project directory (must be on htc), edit the first few
lines of `cellranger.nf` to point it to the right reference and fastq files,
and start the pipeline:
```
nextflow run cellranger.nf
```
This pipeline runs the [count]() command for each sample separately on its own
node, and then runs the [aggregate]() command to normalize across samples once
all the count commands are done. Cell Ranger's batch capabilities do not play
well with SLURM in my experience, so this pipeline runs each command in local
mode, but parallelizes the running of the commands across different nodes.

Once the pipeline is finished running, there should be a directory called
`aggregate` with various outputs from cellranger.
