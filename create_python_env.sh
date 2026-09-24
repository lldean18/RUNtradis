#!/bin/bash
# 21/9/26

# script to create python virtual env for installing software for the RUNtradis script

# make a hidden directory in the shared folder to install the python environment
mkdir -p /share/bryant_lab/.python_3.12.3_runtradis

# activate the python module
module load python-uoneasy/3.12.3-GCCcore-13.3.0

# create the virtual environment
python -m venv /share/bryant_lab/.python_3.12.3_runtradis

# activate the virtual environment
source /share/bryant_lab/.python_3.12.3_runtradis/bin/activate

# install matplotlib
pip install matplotlib

#install biopython
pip install biopython

# to use the environment interactively in TMUX:
conda activate tmux
tmux attach
srun --partition defq --cpus-per-task 4 --mem 20g --time 06:00:00 --pty bash
source $HOME/.bash_profile
module load python-uoneasy/3.12.3-GCCcore-13.3.0
source /share/bryant_lab/.python_3.12.3_runtradis/bin/activate

# to run the script on the ecoli dataset
python /gpfs01/home/mbzlld/github/RUNtradis/insertion_sites_per_gene.py \
    /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.fna \
    /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.embl \
    /share/bryant_lab/essential_genes.txt \
    /share/bryant_lab/non-essential_genes.txt \
    BWtacXpress_merge123_1.fq.gz.NZ_CP009273.1.insert_site_plot.gz

# to run the script on the photorhabdus dataset
python /gpfs01/home/mbzlld/github/RUNtradis/insertion_sites_per_gene.py \
    /share/bryant_lab/reference_genomes/323630L_Photorhabduskhanii.fna \
    /share/bryant_lab/reference_genomes/323630L_Photorhabduskhanii.embl \
    /share/bryant_lab/essential_genes.txt \
    /share/bryant_lab/non-essential_genes.txt \
    

# deactivate the python virtual env
deactivate

