#!/bin/bash
#
#SBATCH --job-name=run_south_korea
#SBATCH --partition=ivm
#SBATCH --nodelist node001
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=2
#SBATCH --output=out_south_korea
#SBATCH --error=err_south_korea
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/risk.py KOR fl gov