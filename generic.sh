#!/bin/bash
#A script for running Balestasis

START=$(date +%s.%N)

if [ ! -f ./balestasis ]

	then

	gcc ./main.c -fopenmp -o balestasis -lm

	echo "Balestasis is ready." 

fi

#set the mode. Remember to include the quotation marks.
#'d': raw data in, net out
#'c': caches (parent sets) in, net out 
#'h': raw data in, caches out
mode='d'

#if mode='d'/'h'
#put the data file in the data folder
#data file should be in the format shown as follows: 
#	the first line is the number of variables, 
#	the second line is the number of data points/instances, 
#	the third line is variable names, 
#	the fourth line is variables' #state, 
#	then instances line by line.

#if mode='c'
#put the cache (scored parent sets) file in the caches folder, the file extension should be .dat
#the file should be in the format as follows: 
#	the first line is the number of variables, 
#	the second line is variable names
#	the third line is variables' #state, 
#	then parent sets in the following format
#	<node> <# of parent sets>
#	<score of the set> <# of parents in the set> <parent> <parent> ......
#	......

#put here the name of the data file, the output net file will have the same name plus the .dsc extension.
datafile=<filename>

#set the name of the log file
logfile=<filename>

#set the name of the net file, files appear only by setting output=1
netfile=<filename>

#set the time (per node) for cache construction (searching local strutures/parent sets). 
#This time limit is not absolute. For ordinarily complicated networks with less than 1000 nodes, 5-10 seconds is usually enough. 
#not needed in mode 'c'.
psiTime=5

#set the time for building global structures. 
#Only the networks found within the limited will be recorded or saved. 
#node needed in mode 'h'
soTime=2

#set the number of ASOBS restarts. The total ASOBS running time will be restart*soTime.
restart=1

#If you want to save the learnt networks to a file, set it to 1, otherwise 0.
#Saved networks will be placed in the network folder.
output=0


printf "Input: ./data/$datafile\n\n"

echo "started at $(date)" >> ./log/time.txt

#--------------1------2-------3------------4---------------5----=----6------7-----8-----9-------------10----------11
./balestasis $mode $psiTime $soTime ./log/$logfile ./data/$datafile 'l' $restart None $output ./network/$netfile None

echo "completed at $(date)" >> ./log/time.txt


END=$(date +%s.%N)

TIME=$(echo "$END - $START" | bc)

echo "time passed: $TIME" >> ./log/time.txt

exit 0


