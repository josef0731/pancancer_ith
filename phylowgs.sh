while read case rank
do
	mkdir phylowgsFINAL/$case
	echo $case
	append="APPENDED_`echo $case`T.vcf"
	cp vcf/$append phylowgsFINAL/$case/$append
	cnv="`echo $case`_oncosnp.txt"
	cp oncosnp/cnv/$cnv phylowgsFINAL/$case/$cnv
	cd phylowgsFINAL/$case
	python ~/phylowgs/parser/phylowgs/create_phylowgs_inputs.py $append -v vardict --cnvs $cnv
	python ~/phylowgs/evolve.py ssm_data.txt cnv_data.txt -s 500
	python ~/phylowgs/posterior_trees.py ssm_data.txt cnv_data.txt
	cd ../..
	cp clonal.py phylowgsFINAL/$case/top_trees/clonal.py
	cp tree_collapse.py phylowgsFINAL/$case/posterior_trees/tree_collapse.py
	cp phylowgs_summary.R phylowgsFINAL/$case/posterior_trees/phylowgs_summary.R
	cp cnv.py phylowgsFINAL/$case/cnv.py
	cd phylowgsFINAL/$case/top_trees
	python clonal.py
	cd ../posterior_trees
	python tree_collapse.py
	cd ..
	python ~/phylowgs/parser/create_phylowgs_inputs.py $append -v vardict --cnvs $cnv
	python cnv.py
	cd posterior_trees
	$HOME/R/bin/Rscript phylowgs_summary.R WES $case --no-save --no-restore
	cd ../../..
done < "caselist`echo $1`.txt"