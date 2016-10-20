#!/bin/bash

while read CASE RANK
do
	echo $CASE
	cd pyclone/$CASE
	cd trace/
	cp ../../../pyclone-post-processing.R pyclone-post-processing.R
	$HOME/R/bin/R CMD BATCH --no-save --no-restore pyclone-post-processing.R
	cd ../../..
done < "caselist.txt"

