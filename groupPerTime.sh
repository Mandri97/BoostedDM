#!/bin/bash


# Make directories
for i in $(seq 0 23); do
    mkdir -p $HOME/timeChunck/data/$i
done

FILE=listOfRootFiles.txt

# Create script to launch ROOT command
SCRIPT="/tmp/readAbsTime.C"

cat > ${SCRIPT}<< EOF
//------------FILE-CONTENT------------//

void readAbsTime(const char* filename ){
    TFile *_file = new TFile(filename);
	    
    std::cout << std::setprecision(20) <<
	((TVectorT<double>*) _file->Get("abstime"))->Max() << std::endl;

    _file->Close();
}

//------------FILE-CONTENT------------//
EOF

while read CMD; do
    COMMAND="root -l -q '/tmp/readAbsTime.C(\"$CMD\")'"

    EPOCH=$(eval $COMMAND | tail -1)
    TIME=$(date -d @${EPOCH} +"%T")

	echo $TIME


    HOUR=$(echo $TIME | cut -d ':' -f 1)

    MIN_HOUR=$(echo $TIME | cut -d ':' -f 2)

	# Convert time zone to EST
    ROUNDED_HOUR=$(awk 'BEGIN {printf "%.0f", '$HOUR'+1 + '$MIN_HOUR'/60}')

	echo "$TIME -> $HOUR:$MIN_HOUR  = $ROUNDED_HOUR"

    if [ $ROUNDED_HOUR -eq 24 ]; then ROUNDED_HOUR=0; fi

    echo "Moving $CMD to $ROUNDED_HOUR dataset"
    LINKING_FILE="ln -s $CMD $HOME/timeChunck/data/$ROUNDED_HOUR"

    eval $LINKING_FILE

done < "$FILE"

rm ${SCRIPT}
