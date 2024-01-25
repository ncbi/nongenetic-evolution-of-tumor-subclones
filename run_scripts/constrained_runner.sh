#!/bin/bash

########
### Loop over all files in given directory and run EvoGeneX to test for 
### constrained evolution over all genes within each file.
### Params: 
### -t File of tree in Newick format with branch lengths
### -r File of the single regime
### -d Path to the data files
### -o Path to save the output
### -p Prefix for the output files (optional)
### -s Flag to use SLURM
### See EvoGeneX documentation for requirements for the files.
### Saves the output for each data file separately under files with the same name.
##########

USE_SLURM="false"
OUTPUT_PREFIX=""

while getopts t:r:d:o:p:s flag; do
  case "$flag" in
    t) TREE_FILE=${OPTARG};;
    r) REGIME_FILE=${OPTARG};;
    d) DATA_PATH=${OPTARG};;
    o) OUTPUT_PATH=${OPTARG};;
    p) OUTPUT_PREFIX=${OPTARG};;
    s) USE_SLURM="true"
  esac
done

if [ ! -f $TREE_FILE ]; then
  echo "Tree file ${TREE_FILE} doesn't exist."
  exit
fi  
  
if [ ! -f $REGIME_FILE ]; then
  echo "Regime file ${REGIME_FILE} doesn't exist."
  exit
fi

if [ ! -d $DATA_PATH ]; then
  echo "Data path ${DATA_PATH} doesn't exist."
  exit
fi  

if [ ! -d $OUTPUT_PATH ]; then
  echo "Output path ${OUTPUT_PATH} doesn't exist."
  exit
fi

if [ "${OUTPUT_PREFIX}" != "" ]; then
  if [ ${OUTPUT_PREFIX: -1} != "_" ]; then
    $OUTPUT_PREFIX="${OUTPUT_PREFIX}_"
  fi
fi

for data_file in $DATA_PATH*; do
  data_name="${data_file%.csv}"
  data_name="${data_name##*/}"
  
  output_file="${OUTPUT_PATH}${OUTPUT_PREFIX}${data_name}.csv"

  if [ "$USE_SLURM" == "true" ]; then
    job_name="${data_name}"
    echo "Using SLURM with job name $job_name"
    sbatch \
      --job-name=$job_name \
      ./run_scripts/constrained_driver.sbatch $TREE_FILE $REGIME_FILE $data_file $output_file  
  else
    echo "Using CPU"
    Rscript run_scripts/constrained_evogenex.R \
      $TREE_FILE \
      $REGIME_FILE \
      $data_file \
      $output_file
  fi
done
