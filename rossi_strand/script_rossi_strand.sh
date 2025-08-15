#!/bin/bash
##SBATCH --gres=gpu:a5000:1
#SBATCH --mem=64G
#SBATCH --output=rossi_strand.out
#SBATCH --error=rossi_strand.err
##SBATCH --exclusive
##SBATCH -w compsci-cluster-fitz-14



conda activate pyranges_env2

python3 rossi_strand.py






