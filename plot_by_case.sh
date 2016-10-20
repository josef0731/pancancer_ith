#!/bin/sh

while read CASE FILENAME
do
	echo $CASE
	Rscript plot_ccf_by_case.R $CASE "gene_col" ccf_by_case/$FILENAME COAD.driver.txt hypermutated.txt --no-save --no-restore
done < "filename.txt"