#!/bin/bash
cat $1 | while read line
do
	tax=$(awk -v OFS='\t' '{print $NF}' $line)
	if [ "$tax" == "NA" ]
	then
		awk 'NF{NF-=1};1' $line
	else 
		echo $line
	fi
done > $2