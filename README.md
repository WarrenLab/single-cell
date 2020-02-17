# single-cell
Pipelines and scripts for single-cell best practices

## Installation
You'll need the following software:
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
* [R](https://www.r-project.org/) with the following packages, all available
  through CRAN unless otherwise noted:
    - [Seurat](https://satijalab.org/seurat/install.html)
    - dplyr
    - argparser
    - ggplot2
    - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
      (get through Bioconductor rather than CRAN)
    - [MAST](https://github.com/RGLab/MAST) (get through BioConductor
      rather than CRAN)

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
[official documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references)
on making a Cell Ranger reference for more information. As recommended in the
documentation, this example script keeps only protein-coding genes (including
those in mtDNA) and lncRNAs, discarding other types of features like rRNAs that
could confound the analysis.

### Creating a sample sheet
Cell Ranger needs a table of information about your samples as input. This is a
good place to note information about, e.g., which samples are control and which
are treatment, so that these data are associated with the samples in the output
tables. It takes this as a comma-separated value (csv) file with the first line
as a header. The first column must be called `library_id` and contain the id
of the sample as encoded in the fastq filename (e.g., the library ID of a file
called `61C_S7_L001_R1_001.fastq.gz` is just `61C`). All other columns are
optional, and can contain whatever data you'd like. Here's an example:
```
library_id,disease_resistance,challenge_status
61C,resistant,control
62C,resistant,control
63C,resistant,control
72C,susceptible,control
74C,susceptible,control
75C,susceptible,control
62M,resistant,challenged
63M,resistant,challenged
64M,resistant,challenged
71M,susceptible,challenged
72M,susceptible,challenged
73M,susceptible,challenged
```

### Running with nextflow
The whole Cell Ranger pipeline can be run for all your samples with a single
nextflow command. Just copy the files `cellranger.nf` and `nextflow.config` from
this repository into your project directory (must be on htc), edit the first few
lines of `cellranger.nf` to point it to the right reference and fastq files,
and start the pipeline:
```
nextflow run cellranger.nf
```
This pipeline runs the `count` command for each sample separately on its own
node, and then runs the `aggregate` command to normalize across samples once
all the count commands are done. Cell Ranger's batch capabilities do not play
well with SLURM in my experience, so this pipeline runs each command in local
mode, but parallelizes the running of the commands across different nodes.

Once the pipeline is finished running, there should be a directory called
`aggregate` with various outputs from cellranger.
