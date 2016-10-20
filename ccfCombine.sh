#!/bin/sh

for file in *_ccf.txt
do
	echo $file
	python ccfcombine.py $file
done