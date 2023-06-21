#!/bin/bash
#
#SBATCH --job-name=run_laos
#SBATCH --partition=ivm
#SBATCH --nodelist node242
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=2
#SBATCH --output=out_laos
#SBATCH --error=err_laos
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/exposure.py LAO tc
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/exposure.py LAO fl