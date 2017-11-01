#!/bin/bash
#A script for running Balestasis

START=$(date +%s.%N)

if [ ! -f ./balestasis ]

	then

	gcc ./main.c -fopenmp -o balestasis -lm

	echo "Balestasis is ready." 

fi


#put the cache (scored parent sets) file in the caches folder
#the file should be in the format as follows: 
#	the first line is the number of variables, 
#	the second line is variable names
#	the third line is variables' #state, 
#	then parent sets in the following format
#	<node> <# of parent sets>
#	<score of the set> <# of parents in the set> <parent> <parent> ......
#	......

#set the name of the data file
datafile=<filename>

#set the time for building global structures. 
#Only the network found within the limited will be recorded or saved.
#This strict limit might not be necessary outside experimental usages. 
soTime=2

#set the number of ASOBS restarts. The total ASOBS running time will be restart*soTime. Sometimes, it's better to do multiple short runs than a single long run.
restart=1

#If you want to save the learnt structure to a file, set it as 1, otherwise 0. (The log will always say the network has been saved.)
#The node order follows the ordering used to form the structure. 
#In 'c' mode, where there aren't raw data, the output is structure only in a file of selected parent sets.
output=0


printf "Input: ./caches/$datafile\n\n"

echo "started at $(date)" >> ./log/time.txt

#-------------1--2----3------------4-------------------5-----------6-----7------8------9-----------10----------------11
./balestasis 'c' 0 $soTime ./log/$datafile.txt ./caches/$datafile 'l' $restart None $output ./network/$datafile.txt None

echo "completed at $(date)" >> ./log/time.txt


END=$(date +%s.%N)

TIME=$(echo "$END - $START" | bc)

echo "time passed: $TIME" >> ./log/time.txt

exit 0


