#!/bin/bash
#A script for running Balestasis

START=$(date +%s.%N)

if [ ! -f ./balestasis ]

	then

	gcc ./main.c -fopenmp -o balestasis -lm

	echo "Balestasis is ready." 

fi


#put the cache (scored parent sets) file in the caches folder
#the file should be in the format shown as follows: 
#	the first line is the number of variables, 
#	the second line is variable names
#	the third line is variables' #state, 
#	then parent sets in the following format
#	<node> <# of parent sets>
#	<the score of the set> <# of parents in the set> <parent> <parent> ......
#	......

#put here the name of the data file
datafile=<filename>

#set the time for building global structures. 
#Only the networks found within the time limit will be recorded or saved.
#This limit is absolute, which might not be necessary in practical uses. 
soTime=2

#set the number of ASOBS restarts. The total ASOBS running time will be restart*soTime. (Sometimes, it is better to do multiple short runs than a single long run.)
restart=1

#If you want to save learnt network structures to a file, set it to 1, otherwise 0. (But the log will always say that the networks has been saved.)
#The node order follows the ordering used in finding the structure. 
#In 'c' mode, where there is no raw datum, the output is a structure specified by parent sets.
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


