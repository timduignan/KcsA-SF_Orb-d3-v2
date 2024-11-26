#!/bin/bash
#PBS -l ncpus=16
#PBS -l ngpus=1
#PBS -l mem=96GB
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/tn51+gdata/tn51
#PBS -l wd
#PBS -q dgxa100
#PBS -P tn51
#PBS -N LiBr-run1
#PBS -l jobfs=10GB

module load python3/3.12.1
source /g/data/tn51/shared/orb-env/bin/activate
#wget https://storage.googleapis.com/orbitalmaterials-public-models/forcefields/orb-d3-v1-20240902.ckpt -O orb-d3-v1-20240902.ckpt
python3 KcsA.py

