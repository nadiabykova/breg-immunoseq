cd /mnt/mapr/group/ms-solid/immuno/
path=$1
pathout=$2
script_path_1="scripts/util.filter_chain_reads.R"
mixcr_path=/mnt/mapr/group/ms-solid/immuno/mixcr-3.0.10/mixcr_n

declare -a chains=("IGH" "IGK" "IGL")
for f in $(ls $path*_R1_001.fastq.gz)
do
	x=$(basename $f)
	name="${x%_R*_001.fastq.gz}"
	echo $name
	f1=$path$name"_R1_001.fastq.gz"
	f2=$path$name"_R2_001.fastq.gz"
	f1_fq=$path$name"_R1_001.fastq"
	f2_fq=$path$name"_R2_001.fastq"
	y=$pathout$name
	echo $y
	#rm $y.align.report
	#$mixcr_path align -b imgt.201822-5.sv4 -s hs -f -r $y.align.report -OsaveOriginalReads=true $f1 $f2 $y.vdjca
	#gunzip -c $f1  > $f1_fq
	#gunzip -c $f2  > $f2_fq
	for ch in "${chains[@]}"
	do
		#$mixcr_path  exportAlignments -c $ch -f -readIds $y.vdjca $y.$ch.readIds.txt
		#Rscript $script_path_1 $y.$ch.readIds.txt $f1_fq $y.$ch.R1.fastq
		#Rscript $script_path_1 $y.$ch.readIds.txt $f2_fq $y.$ch.R2.fastq
		#$mixcr_path align -b imgt.201822-5.sv4 -s hs -f -r $y.$ch.align.report -OsaveOriginalReads=true $y.$ch.R1.fastq $y.$ch.R2.fastq $y.$ch.vdjca
		#$mixcr_path assemble -f -OaddReadsCountOnClustering=true -r $y.$ch.assemble.report -a $y.$ch.vdjca $y.$ch.clna
		#$mixcr_path exportReadsForClones -f $y.$ch.clna $y.$ch.reads.fastq
		#$mixcr_path exportClones -c $ch -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nFeature CDR3 $y.$ch.clna $y.$ch.clones.txt
		$mixcr_path exportAlignments -f -cloneId -readIds -targets -targetSequences -defaultAnchorPoints -vGene -vHit -vHitScore -vAlignment -jGene -aaFeature CDR3 -nFeature CDR3 -nFeature CDR1 -nMutationsRelative CDR1 VRegion -nFeature FR2 -nMutationsRelative FR2 VRegion -nFeature CDR2 -nMutationsRelative CDR2 VRegion -nFeature FR3 -nMutationsRelative FR3 VRegion $y.$ch.clna $y.$ch.aln.txt
	done
done 
