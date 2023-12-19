#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=13
#SBATCH --partition=ventana
#SBATCH --time=24:00:00
#SBATCH --mem=250G
#SBATCH --output=CM1stitch-%j.txt

# Optional lines. Yours will be different. Change as necessary or remove/comment out.
source /etc/profile
source ~/.bashrc
conda activate ventana

# Execute the code.
time python stitchCM1.py
