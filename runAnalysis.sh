#!/bin/bash

function readTxtFile() {
	counter=0
	FILELIST=$1

    # TODO: Specify where the PBS script is
    PBS_SCRIPT=$HOME/BoostedDM/runAnalysis.pbs

	while read CMD; do
        echo "Submiting job $counter - Analysis of $CMD ..."

		qsub -F "$CMD $2" $PBS_SCRIPT
		counter=$((counter + 1))

		if (( $counter % 11 == 0)); then
			sleep 60
			counter=0
		fi
	done < "$FILELIST"
}


function checkLiveTime(){
# Create script to launch ROOT command
SCRIPT="/tmp/readRunTime.C"

cat > ${SCRIPT}<< EOF
//------------FILE-CONTENT------------//

void readRunTime(const char* filename ){
    TFile *_file = new TFile(filename);
	
	auto runtime = (TVectorT<double>*) _file->Get("runtime");

	if (runtime->GetNrows() != 1) {
		printf("Multiple runtime? Seriously?");
		exit(1);
	}

	std::cout << std::setprecision(20) << runtime->Max() << std::endl;

    _file->Close();
}

//------------FILE-CONTENT------------//
EOF

	touch $2/runtimeCheck.txt

    FILELIST=$1

	while read FILE; do
        echo "Reading $FILE ..."

		CMD="root -l -q '/tmp/readRunTime.C(\"$FILE\")'"	
		echo $(eval $CMD | tail -1) >> $2/runtimeCheck.txt
	done < "$FILELIST"

	rm -rf /tmp/readRunTime.C
}

################################
#         Main Program         #
################################
# Check the right number of args
if [ "$#" -ne 2 ]; then
    echo "Usage: ./runAnalysis.sh listOfFilesToAnalyze.txt path/to/output_dir"
    exit 1
fi

# Prompt to get user's choice of analysis
CONF_="n"
ANALYSIS_=10

LISTOFFILES=$1

readTxtFile "${LISTOFFILES}" $2
# Calling function
#if [ $ANALYSIS_ -eq 1 ]; then
#elif [ $ANALYSIS_ -eq 2 ]; then
    # Overwrite outputfile
#    rm -f $2/runtimeCheck.txt

#    checkLiveTime "${LISTOFFILES}" $2
#fi

