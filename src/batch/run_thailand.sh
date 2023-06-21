#!/bin/bash
#
#SBATCH --job-name=run_thailand
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=4
#SBATCH --output=out_thailand
#SBATCH --error=err_thailand
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/exposure.py THA tc
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/exposure.py THA fl