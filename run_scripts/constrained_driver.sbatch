#!/bin/bash

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

###########
### Run EvoGeneX to test for for constrained evolution on all genes within the input file.
### Params: 
### Arg 1: File of tree in Newick format with branch lengths
### Arg 2: File of the single regime
### Arg 3: Path to the data file
### Arg 4: Path to save the output
### See EvoGeneX documentation for requirements for these files
##########

module load R/4.1.2/3.14

TREE_FILE=$1
REGIME_FILE=$2
DATA_PATH=$3
OUTPUT_FILE=$4

Rscript run_scripts/constrained_evogenex.R \
  $TREE_FILE \
  $REGIME_FILE \
  $DATA_PATH \
  $OUTPUT_FILE
