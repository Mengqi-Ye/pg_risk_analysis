#!/bin/bash
#
#SBATCH --job-name=run_singapore
#SBATCH --partition=ivm
#SBATCH --nodelist node242
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=2
#SBATCH --output=out_singapore
#SBATCH --error=err_singapore
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/risk.py SGP tc gov