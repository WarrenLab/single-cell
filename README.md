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

You can install the R package contained in this repository if you've got R
installed by downloading the source from github and then telling R to install
it:
```bash
git clone https://github.com/WarrenLab/single-cell.git
cd single-cell
R CMD INSTALL .
```
This command should in theory be able to install all the dependencies so you
don't have to do it yourself, but usually installing Seurat requires some manual
intervention to get some of _its_ dependencies.

If you want to use nextflow to run the Cell Ranger pipeline (recommended),
you'll also need to install nextflow with the following command if you haven't
already:
```
curl -s https://get.nextflow.io | bash
```
This creates a binary named 'nextflow' in your current directory; you probably
want to move it to somewhere in your path.

### Nextflow on lewis
Nextflow needs to be run from `htc` due to file-locking issues*, but is much
faster when it's doing file operations on `hpc`. To get around this, you can set
it to always use a work directory on `hpc` by adding something like this to your
`.bashrc` (changing the location, of course):
```bash
export NXF_WORK=/storage/hpc/group/warrenlab/users/esrbhb/nx_work

# Change work directory to group ownership and make it "stick" 
chgrp -R warrenlab-gropu
chmod -R ug+rwX $NXF_WORK
chmod g+x $NXF_WORK
```

Since files in this directory count towards the `hpc` quota associated with the 
group ownership of that directory, this example sets the group ownership of this
directory and its subdirectories to "warrenlab-group" (assuming that you are 
performing work for said group).

A good nextflow configuration for lewis is in `nextflow.config`. Copy this file
to the directory you are running nextflow from.

* WARNING:  A subdirectory called `.nextflow` is created (or updated) inside the
directory from which you run Nextflow. This automatically created `.nextflow` 
directory is what needs to exist on `htc`. We have seen cases where running 
from a directory on `hpc` worked *some of the time*, but would later fail 
with hard-to-diagnose errors, sometimes including Nextflow apparently 
"forgetting" previously successful hours-long steps and rerunning them 
needlessly.

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
nextflow run pipeline.nf
```
This pipeline runs the `count` command for each sample separately on its own
node, and then runs the `aggregate` command to normalize across samples once
all the count commands are done. Cell Ranger's batch capabilities do not play
well with SLURM in my experience, so this pipeline runs each command in local
mode, but parallelizes the running of the commands across different nodes.

Once all the Cell Ranger stuff is done, the pipeline runs some basic analyses
in Seurat with very permissive cutoffs and (hopefully) sensible defaults. This
should not be the only Seurat run you do. See the section below for more
defaults.

Once the pipeline is finished running, there should be a directory called
`aggregate` with various outputs from cellranger, and `seurat_out` for some
plots from Seurat.

## Seurat
Cell Ranger aligns all the scRNA-seq reads to the reference, and then tabulates
a matrix of how many times each gene (or rather, "feature," as you can also
count non-gene things like lncRNAs) shows up in each cell. You probably want to
actually use this big dataset to learn things, though, which is where Seurat
comes in. Seurat makes it easy to do some common data analysis tasks with this
big matrix, such as:
* filtering out things that aren't useful, like cells that are undergoing
  apoptosis, or genes that don't get expressed often enough to tell you anything
* normalizing and scaling the counts to make them comparable to each other
  across cells and features
* integrating data from different experiments or conditions
* reducing the dimensionality of the data so you can look at it in two
  dimensions instead of two thousand dimensions
* clustering the cells into groups with similar expression profiles
* finding marker genes for each cluster so that you can label them with cell
  identities
* visualizing the data in a bunch of different ways

This package contains a wrapper script in R that uses Seurat to do all of these
things, `run_seurat.R`. The nextflow pipeline runs this script on the Cell
Ranger output with permissive default parameters. That is, very few things are
filtered out. It also makes some diagnostic plots that will hopefully be helpful
for choosing better parameters than the defaults, such as `feature_plot.pdf`,
`violin_plot.pdf`, and `elbow.pdf`.

How to choose good parameters is beyond the scope of this README, but the
[Seurat vignettes](https://satijalab.org/seurat/vignettes.html) and
[this paper](https://www.embopress.org/doi/10.15252/msb.20188746) are good
places to learn about that. We **strongly** recommend running `run_seurat.R`
with a couple different sets of parameters. You can see the options with
```
run_seurat.R -h
```
You can also use this script for comparative analysis between different sets
of conditions (e.g., treated and untreated) with the `--integrate` and
`--group-var` options. For example, with the example sample sheet from this
README, you could set `--group-var=challenge_status` to compare the challenged
and control cells.

Finally, although we've included in this script all of the analyses we always
run on a new single-cell data set, there are likely some other things you want
to do based on the questions you generated your data to answer, so you can do
this by starting an interactive R session and loading the seurat object, which
the script saves as an RDS file, like this:
```R
library(Seurat)
seurat <- readRDS('seurat_out/seurat.rds')
```
and then plot away!
