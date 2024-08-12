# warwick_edna_nf
nextflow pipeline for Warwick eDNA project

## Conda environment should be created as follows
```
conda env create -f env/warwick_edna.yaml
```

## Running the pipeline
```
conda activate warwick_edna

nextflow run warwick_edna.nf \
--R1 test_data/M10-211123_S10_L001_R1_001.fastq.gz \
--R2 test_data/M10-211123_S10_L001_R2_001.fastq.gz
--forward_primer GGGTTGGTAAATTTCGTGCCAGC \
--reverse_prime CATAGTGGGGTATCTAATCCCAGTTTG \
--prefix M10-211123_S10
```


## To do

- pia will need to be installed locally
- update blast database (perl update_blastdb.pl --decompress nt)




