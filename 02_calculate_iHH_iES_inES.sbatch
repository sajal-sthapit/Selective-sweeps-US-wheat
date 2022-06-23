#!/bin/bash -l
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END --mail-user=sthapit@landinstitute.org
#SBATCH --time=0-04:00:00
#SBATCH --output="02_calculate_iHH_iES_inES_console.output"
#SBATCH --error="02_calculate_iHH_iES_inES_console.error"

module purge
module load R
R --no-save -q < 02_calculate_iHH_iES_inES.R