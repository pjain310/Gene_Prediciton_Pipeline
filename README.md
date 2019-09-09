# Team3-GenePrediction

## Bacterial ab-initio Gene Prediction Pipeline

Gene prediction is the process of finding which regions of genomic DNA encodes genes and non-coding RNA.This pipeline is meant to predict coding and non-coding regions in assembled contigs of a bacterial genome. Tool and paramter selection is carried out to ensure best performance for de-novo assembled Listeria monocytogenes genomes.

## Installation and Setup 

This pipeline uses as conda based environment to ensure you have the appropriate dependencies. We recommend that you download and install Miniconda from https://conda.io/en/latest/miniconda.html

Example installation for Miniconda on Linux:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
rm  Miniconda3-latest-Linux-x86_64.sh
```

Next, clone the repository into your local system:

```
git clone  https://github.gatech.edu/compgenomics2019/Team3-GenePrediciton.git
```

Create and activate a conda environment using the yml file provided in our lib folder:

```
#Create environment after downloading yml file
conda-env create -f lib/gp_env.yml -n myenv
source activate myenv
```

From within the hmmer-2.2 folder in lib, compile binaries for hmmer-2.2 (dependency for rnammer which is part of our pipeline):
```
cd hmmer-2.2
./configure
make install
```

Export path to 'lib' to path variable (lib contains precompiled binaries for GenemarkS-2, Glimmer, RNAmmer which are part of the pipeline)
```
export PATH=$PATH:<path to lib>
```

## Running the pipeline

To run our pipeline with sample data provided in our repository (check sample_input folder)

```
./gp_pipeline.sh -i sample_input -o sample_output
```

For each input genome, the list of generated outputs is as follows:
1. gff file containing the coordinates for the coding sequences
2. fna file for coding nucleotide sequences
3. faa file for coded protein sequences
3. fna file for RNA predictions
