#!/bin/bash

usage="$(basename "$0") [-h] [-g gain] -n name
Program to take samples at a fixed frequency from the direction antenna.
where: 
    -h|--help              show this help text
    -g|--gain     (=62)    set the gain of the receiver
    -r|--rate     (=7.68MHz)  set the sampling rate of the receiver
    -d|--duration (=8)   set the duration (seconds) during which to take samples
    -i|--nfiles   (=230)   number of files of samples to save
    -f|--freq     (=320)   set the frequency (MHz) at which to take samples
    -n|--name              set the name of the folder and archive in which samples will be saved"

while [[ $# -gt 0 ]] 
do
key="$1"
case $key in
    -n|--name)
    NAME="$2"
    shift
    ;;
    -g|--gain)
    GAIN="$2"
    shift
    ;;
    -r|--rate)
    RATE="$2"
    shift
    ;;
    -d|--duration)
    DUR="$2"
    shift
    ;;
    -i|--nfiles)
    NFILES="$2"
    shift
    ;;
    -f|--freq)
    FREQ="$2"
    shift
    ;;
    -h|--help)
    echo "$usage"
    exit
    ;;
  esac
  shift
done

NAME=${NAME:-"NA"}
if [ "$NAME" == "NA" ]; then
    echo "Please input name (-n)"
    exit 0
fi

GAIN=${GAIN:-62}
RATE=${RATE:-15360000}
DUR=${DUR:-8}
NFILES=${NFILES:-99999999}
FREQ=${FREQ:-320}
FILETYPE=".dat"

M=1000000       # MHz
freq=$[FREQ*M]
FILE="Thwaites2019_E312_wGPS_freq"$FREQ"_gain"$GAIN"_BW"$RATE"_ChirpLen"$DUR"_SN310725C_"
mkdir /media/usb0/$NAME
cd /tmp
rm /tmp/*.dat
COUNTER=0
LASTCOUNTER=0
while [ $COUNTER -lt $NFILES ]; do
 	echo Counter $COUNTER
	#sleep 3
	#cd /tmp
	if [ $COUNTER == 0 ]
	then	
		$HOME/DeploymentCode/rx_samples_to_file_stampGPS_for_E312_ver5/build/rx_samples_to_file_stampGPS --file /tmp/$FILE$COUNTER$FILETYPE --gain $GAIN --rate $RATE --bw $RATE --duration $DUR --freq $freq --ant RX2 --subdev A:A --args="master_clock_rate=61440000" --stats --progress
	else
		let LASTCOUNTER=$COUNTER-1
		#echo Last File: $FILE$LASTCOUNTER$FILETYPE
		#ls
		mv $FILE$LASTCOUNTER$FILETYPE /media/usb0/$NAME/$FILE$LASTCOUNTER$FILETYPE &
		$HOME/DeploymentCode/rx_samples_to_file_stampGPS_for_E312_ver5/build/rx_samples_to_file_stampGPS --file /tmp/$FILE$COUNTER$FILETYPE --gain $GAIN --rate $RATE --bw $RATE --duration $DUR --freq $freq --ant RX2 --subdev A:A --args="master_clock_rate=61440000" --stats --progress &
		wait
		#ls
		#rm $FILE$LASTCOUNTER$FILETYPE
	fi
	#mv $FILE$COUNTER /media/usb0/$NAME
	#rm $FILE$COUNTER
	
	let COUNTER+=1
done

let COUNTER-=1
mv $FILE$COUNTER$FILETYPE /media/usb0/$NAME/$FILE$COUNTER$FILETYPE

echo bash script done
#~/grc/applications/take_samples_fixed/build/take_samples_fixed --file ~/grc/samples/$NAME/ --gain $GAIN --nfiles $NFILES --freq $FREQ --rate $RATE --duration $DUR


# #make archive and clean up
# echo "Compressing..."
# ARCHIVENAME="dir_fixed_"
# ARCHIVENAME+=$NAME
# ARCHIVENAME+=".tar.xz"
#


#
# #tar -cf ~/$ARCHIVENAME $NAME/*
# tar -cf /media/usb0/$ARCHIVENAME $NAME/*
# rm -r ~/grc/samples/$NAME
#
# cd ~/grc/scripts
