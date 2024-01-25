#!/bin/bash

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=05:00:00
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

###########
### Run EvoGeneX to test for for adaptivity on all genes within the input file.
### Params: 
### Arg 1: File of tree in Newick format with branch lengths
### Arg 2: File of the single regime
### Arg 3: File for the two-regime
### Arg 4: Path to the data file
### Arg 5: Path to save the output
### See EvoGeneX documentation for requirements for these files
##########

module load R/4.1.2/3.14

TREE_FILE=$1
SINGLE_REGIME_FILE=$2
TWO_REGIME_FILE=$3
DATA_PATH=$4
OUTPUT_FILE=$5

Rscript run_scripts/adaptive_evogenex.R \
  $TREE_FILE \
  $SINGLE_REGIME_FILE \
  $TWO_REGIME_FILE \
  $DATA_PATH \
  $OUTPUT_FILE

