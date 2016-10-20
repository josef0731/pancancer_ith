#!/bin/bash

for file in *_ccfplot_only_nonsynonymous.pdf
do
	trimmed=${file%.pdf}
	convert -verbose -density 150 $file -quality 100 -sharpen 0x1.0 `echo $trimmed`.jpg
done