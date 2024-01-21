#!/bin/bash

#SBATCH --job-name=scIntegrate
#SBATCH --mail-user=makpabey@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=10g
#SBATCH --time=24:00:00
#SBATCH --account= XXXXXX
#SBATCH --partition=standard
#SBATCH --output=~/Logs/sc_covid_integrate.log

ml python3.10-anaconda/2023.03

python integrate.py

