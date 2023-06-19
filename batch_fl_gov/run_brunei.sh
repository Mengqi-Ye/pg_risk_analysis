#!/bin/bash
#
#SBATCH --job-name=run_brunei
#SBATCH --partition=ivm
#SBATCH --nodelist node240
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=4
#SBATCH --output=out_brunei
#SBATCH --error=err_brunei
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/risk.py BRN fl gov