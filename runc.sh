#!/bin/bash
#A script for running Balestasis

START=$(date +%s.%N)

if [ ! -f ./balestasis ]

	then

	gcc ./main.c -fopenmp -o balestasis -lm

	echo "Balestasis is ready." 

fi


#put the data file in the data folder, the extension should be .dat
#data file should be in the format as follows: 
#	the first line is the number of variables, 
#	the second line is the number of data points/instances, 
#	the third line is variable names, 
#	the fourth line is variables' #state, 
#	then instances line by line.

#set the name of the data file without the extension, for example, data from data.dat
datafile=alarmCache

#set the time (per node) for cache construction (searching local strutures). 
#This time limit isn't strict. For ordinarily complex networks with less than 1000 nodes, 5-10 seconds is usually enough. 
psiTime=5

#set the time for building global structures. 
#Only the network found within the limited will be recorded or saved, which might not be necessary outside experimental usages. 
soTime=2

#set the number of ASOBS restarts. The total ASOBS running time will be restart*soTime. Sometimes, for example, it is better to run five 2-second runs than a single 10-second run.
restart=1

#if you want to save the learnt network to a file, set it as 1, otherwise 0. (The log will always say the network has been saved.)
#Saved network will be placed in the network folder. This won't work with preloaed cache.
saveNet=0


printf "Data: $datafile\n\n"

echo "started at $(date)" >> ./log/time.txt

#			  1	   2	   3			4					  5			    6	 7	     8	   9----------------10------------11
./balestasis 'c' $psiTime $soTime ./log/$datafile.txt ./data/$datafile.dat 'l' $restart None $saveNet ./network/$datafile.dsc None

echo "completed at $(date)" >> ./log/time.txt


END=$(date +%s.%N)

TIME=$(echo "$END - $START" | bc)

echo "time passed: $TIME" >> ./log/time.txt

exit 0


