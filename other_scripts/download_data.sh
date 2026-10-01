#!/bin/bash

# script to download the short read data for Photorhabdus from
# the onedrive zip folder shared with me after I copied it to my sharepoint

#SBATCH --time=24:00:00
#SBATCH --job-name=rclone
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16g


# setup env
module load rclone-uon/1.65.2

# copy the directory with rclone
rclone --transfers 4 --checkers 4 --bwlimit 100M --onedrive-chunk-size 5M \
--checksum copy Laura2:HPC_data_dirs_backup/tradis/rawdata/Photorhabdus_khanii /gpfs01/home/mbzlld/data/tradis/rawdata/Photorhabdus_khanii

# Check the directory has copied successfully
rclone check --one-way Laura2:HPC_data_dirs_backup/tradis/rawdata/Photorhabdus_khanii /gpfs01/home/mbzlld/data/tradis/rawdata/Photorhabdus_khanii

# unload module
module unload rclone-uon/1.65.2
