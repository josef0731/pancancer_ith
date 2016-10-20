#!/bin/bash

while read CASE
do
	echo $CASE
	maf=`echo $CASE`parsed.maf
	append=APPENDED_`echo $CASE`parsedTRIMMED.vcf
	segdat=`echo $CASE`.txt
	mkdir pycloneFINAL/$CASE
	cd maf
	perl $HOME/vcf2maf/maf2vcf.pl --input-maf `echo $CASE`parsedTRIMMED.maf --output-dir ../vcfFINAL --ref-fasta $HOME/vep/homo_sapiens/82_GRCh37/human_g1k_v37.fasta.gz
	cd ../vcfFINAL
	vcf="`echo $CASE`parsedTRIMMED.vcf"
	Rscript append_read_depth.R $vcf --no-save --no-restore
	cd ..
	cp vcfFINAL/$append pycloneFINAL/$CASE/$append
	cp pyclone-parser.py pycloneFINAL/$CASE/pyclone-parser.py
	cp pyclone-parserNO-CNV.py pycloneFINAL/$CASE/pyclone-parserNO-CNV.py
	cp pyclone-yaml-parser.py pycloneFINAL/$CASE/pyclone-yaml-parser.py
	cp pyclone-config.yaml pycloneFINAL/$CASE/pyclone-config.yaml
	cp cnv/`echo $CASE`.txt pycloneFINAL/$CASE/`echo $CASE`.txt
	cd pycloneFINAL/$CASE
#	python pyclone-parser.py $CASE
	python pyclone-parserNO-CNV.py $CASE
	python pyclone-yaml-parser.py $CASE
	python ~/pyclone-0.12.9/PyClone build_mutations_file --ref_prior normal --var_prior parental_copy_number pyclone.mutation.tsv `echo $CASE`.yaml
	python ~/pyclone-0.12.9/PyClone analyse pyclone-config.yaml
	cp ../../pyclone-cluster.R trace/
	cd trace/
	bzip2 -d labels.tsv.bz2
	$HOME/R/bin/R CMD BATCH --no-save --no-restore pyclone-cluster.R
	mv `echo $CASE`.cellular_frequencies.tsv.bz2 cellular.frequencies.tsv.bz2
	bzip2 -d cellular.frequencies.tsv.bz2
	cp ../../../pyclone-post-processing.R pyclone-post-processing.R
	$HOME/R/bin/R CMD BATCH --no-save --no-restore pyclone-post-processing.R
	cd ../../..
done < "msi`echo $1`.txt"

