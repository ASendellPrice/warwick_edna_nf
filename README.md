# warwick_edna_nf
A nextflow pipeline for the Warwick eDNA project.

## Dependencies
The following dependencies are required:
- Java v11 or later
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
--input_csv invertebrate_samples.csv \
--outdir results_invertebrate_primers_full_blast \
-process.cpus 16 -resume
```







conda activate warwick_edna
nextflow run warwick_edna.nf \
--input_csv original_V_samples.csv \
--blast_db /home/u2271009/eDNA/bin/blast_db/nt \
--outdir results_original_V_metazoa_blast \
--taxa 33208 -process.cpus 16 -resume



## To do

- pia will need to be installed locally
- update blast database (perl update_blastdb.pl --decompress nt)
- add flag to change percent identity
- taxaranks auto install - currently installed using pip



