#!/bin/sh

#Usage:
#./runJob jobScript inputDir outputDir

source $HOME/setup.sh

if [[ "$*" == "" || "$#" -lt "1" ]]; then
    echo "Usage: ./runJob.sh jobScript inputDir (outputDir)"
    exit 1
fi

JOB_SCRIPT=$1

for iFolder in $(seq 0 23); do
    INPUT_DIR=$HOME/timeChunck/data/$iFolder
    OUTPUT_DIR=$HOME/timeChunck/data/$iFolder/analysis

    mkdir -p $OUTPUT_DIR

    for file in $(ls $INPUT_DIR/*.root); do
	CMD="qsub -v 'OUTPUT_ID=$OUTPUT_DIR' -F '$file $OUTPUT_DIR' $JOB_SCRIPT"
	echo "Submitting ... $CMD"

	eval $CMD

	sleep 1
    done
done
