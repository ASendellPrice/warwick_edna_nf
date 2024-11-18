# warwick_edna_nf
A nextflow pipeline for the Warwick eDNA project.

## Dependencies
The following dependencies are required:
- Java v11 or later (required by nextflow)
- Nextflow (see installation instructions [here](https://nextflow.io/docs/latest/install.html))
- A conda installation e.g. [miniconda](https://docs.anaconda.com/miniconda/)
- PIA (see installation instructions [here](https://github.com/Allaby-lab/PIA/))
- A local blast database

## Installation & environment setup
```
git clone https://github.com/ASendellPrice/warwick_edna_nf.git
cd warwick_edna
conda env create -f env/warwick_edna.yaml
```

## Updating the configuration file
Before running the pipeline the configuration file "config/warwick_edna.config" will need to be updated to include the correct paths to the blast database and PIA installation directory.

## Running the pipeline
The pipeline (which runs FastQC, multiQC, trim_galore, seqtk, mmseq2, blast, & pia) can be initiated as follows:
```
# Activate conda environment
conda activate warwick_edna

# Launch nextflow pipeline
nextflow run warwick_edna.nf -config config/warwick_edna.config \
--input_csv test.csv --outdir results_test -process.cpus 16
```
**Note:** There is an optional parameter "--taxa" which can be parsed to restrict blast search to a specific taxanomic group. For example, if you want to search against metazoa only add the following: --taxa 33208

## To do
- update blast database (perl update_blastdb.pl --decompress nt)
- add flag to change percent identity
- taxaranks auto install - currently installed using pip
- update how resources are allocated
- add step which generates fastas for each taxa