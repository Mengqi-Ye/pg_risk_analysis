#!/bin/bash
#
#SBATCH --job-name=run_malaysia
#SBATCH --partition=ivm
#SBATCH --nodelist node240
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=2
#SBATCH --output=out_malaysia
#SBATCH --error=err_malaysia
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/exposure.py MYS tc
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/exposure.py MYS fl