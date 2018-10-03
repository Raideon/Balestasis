#!/bin/bash
#A script for running Balestasis

START=$(date +%s.%N)

if [ ! -f ./balestasis ]

	then

	gcc ./main.c -fopenmp -o balestasis -lm

	echo "Balestasis is ready." 

fi


#put the data file in the data folder, the extension should be .dat
#data file should be in the format shown as follows: 
#	the first line is the number of variables, 
#	the second line is the number of data points/instances, 
#	the third line is variable names, 
#	the fourth line is variables' #state, 
#	then instances line by line.

#put here the name of the data file without the extension and change the extension of the file to .dat, for example, data from data.dat
datafile=<filename>

#set the time (per node) for cache construction (searching local strutures/parent sets). 
#This time limit is not absolute. For ordinarily complex networks with less than 1000 nodes, 5-10 seconds is usually enough. 
psiTime=5

#set the time for building global structures. 
#Only the network found within the limited will be recorded or saved.
soTime=2

#set the number of ASOBS restarts. The total ASOBS running time will be restart*soTime. (Sometimes, for example, it is better to run five 2-second runs than a single 10-second run.)
restart=1

#if you want to save the learnt network to a file, set it to 1, otherwise 0. (The log will always say that the networks are saved.)
#Networks will be saved in the network folder.
output=0


printf "Input: ./data/$datafile\n\n"

echo "started at $(date)" >> ./log/time.txt

#-------------1------2-------3------------4-------------------5-------------6------7------8-----9-------------10--------------11
./balestasis 'd' $psiTime $soTime ./log/$datafile.txt ./data/$datafile.dat 'l' $restart None $output ./network/$datafile.dsc None

echo "completed at $(date)" >> ./log/time.txt


END=$(date +%s.%N)

TIME=$(echo "$END - $START" | bc)

echo "time passed: $TIME" >> ./log/time.txt

exit 0


