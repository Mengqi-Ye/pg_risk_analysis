#!/bin/bash
#
#SBATCH --job-name=run_myanmar
#SBATCH --partition=ivm
#SBATCH --nodelist node241
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=2
#SBATCH --output=out_myanmar
#SBATCH --error=err_myanmar
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/risk.py MMR tc osm