#!/bin/bash
#
#SBATCH --job-name=run_china
#SBATCH --partition=ivm
#SBATCH --nodelist node241
#SBATCH --ntasks=1
#SBATCH --time=10-00
#SBATCH --cpus-per-task=4
#SBATCH --output=out_china
#SBATCH --error=err_china
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengqi.ye@vu.nl

source activate py310
python /scistor/ivm/mye500/projects/pg_risk_analysis/src/risk.py CHN tc osm