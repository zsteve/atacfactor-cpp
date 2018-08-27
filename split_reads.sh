#!/bin/bash

while read x; do
	read_id=$(echo $x | cut -d' ' -f 1);
	echo $x >> "test_data/split/$read_id"
done < $1 
