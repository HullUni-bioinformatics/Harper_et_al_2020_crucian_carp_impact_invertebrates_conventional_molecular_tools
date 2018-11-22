# Harper *et al.* (2018) Invertebrate metabarcoding

Data processing workflow and supplementary data for:

Harper *et al.* Assessing the impact of the threatened crucian carp (*Carassius carassius*) on pond invertebrate diversity - a comparison of conventional and molecular tools.
Permanently archived at: **TBC**


## Contents

Curated reference databases used in analyses (GenBank/fasta format) [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Reference_database)

Notebooks to run metaBEAT pipeline [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Jupyter_notebooks)

NCBI Sequence Read Archive (SRA) accession numbers for raw Illumina data [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Data/Sample_accessions.tsv)

Taxonomic assignment results [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Data/Taxonomic_Assignment_Results)

R scripts used to analyse metaBEAT output and produce figures [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/R_scripts)

Sample metadata needed to run analyses in R [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Data/Sample_Metadata)


## Instructions to set up dependencies for data processing and analyses

To facilitate full reproducibility of our analyses, we provide Jupyter notebooks illustrating our workflow and all necessary associated data in this repository.

Illumina data was processed (from raw reads to taxonomic assignment) using the [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT) pipeline. The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self-contained docker image that includes all necessary dependencies [here](https://hub.docker.com/r/chrishah/metabeat/).


## Setting up the environment

In order to retrieve scripts and associated data (reference sequences, sample metadata etc.), start by cloning this repository to your current directory:

```
git clone --recursive https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_netabarcoding.git
```

In order to make use of our self contained analysis environment, you will have to install Docker on your computer. Docker is compatible with all major operating systems, but see the Docker documentation for details. On Ubuntu, installing Docker should be as easy as:

```
sudo apt-get install docker.io
```

Once Docker is installed, you can enter the environment by typing:

```
sudo docker run -i -t --net=host --name metaBEAT -v $(pwd):/home/working chrishah/metabeat /bin/bash
```

This will download the metaBEAT image (if not yet present on your computer) and enter the 'container' i.e. the self contained environment (**NB:** ```sudo``` may be necessary in some cases). With the above command, the container's directory ```/home/working``` will be mounted to your current working directory (as instructed by ```$(pwd)```). In other words, anything you do in the container's ```/home/working``` directory will be synced with your current working directory on your local machine.


## Data processing workflow as Jupyter notebooks

Raw illumina data has been deposited on the NCBI SRA:
- Study: SRP163672
- BioProject: PRJNA494857
- BioSample accessions: SAMN10181701 - SAMN10182084 (bulk tissue DNA) and SAMN10187732 - SAMN10188115 (eDNA)
- SRA accessions: SRR7969394 - SRR796977 (bulk tissue DNA) and SRR7985814 - SRR7986197 (eDNA)

The sample specific accessions can be found [here](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Data/Sample_accessions.tsv). Before following the workflow for data processing, you'll need to download the raw reads from the SRA. To download the raw read data, you can follow the steps in this [Jupyter notebook](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/raw_reads/How_to_download_from_SRA.ipynb).

With the data in place, you should be able to fully reproduce our analyses by following the steps outlined in the [Jupyter notebooks](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Jupyter_notebooks).

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory ```raw_reads``` at the base of the repository structure and that the files are named according to the following convention: 'sampleID-marker', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. SampleID must correspond to the first column in the file ```Sample_accessions.tsv``` [here](https://github.com/HullUni-bioinformatics/Harper_et_al_2018_invertebrate_metabarcoding/tree/master/Data/Sample_accessions.tsv).
