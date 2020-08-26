#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l hostname=hpc03


#Java 8 is requires for mixcr. It's installed (only?) on the hpc03 node, so we limit our request by "-l hostname=hpc03" to this node

seed=$1
path=$2
path2=$3
mixcr_path=/mnt/mapr/group/ms-solid/immuno/mixcr-3.0.10/mixcr_n

declare -a chains=("IGH" "IGK" "IGL")

cd /mnt/mapr/group/ms-solid/immuno/

for ch in "${chains[@]}"
do
	for f in $path/*$ch".sampled_seed"$seed"_R1.fastq"
	do
		x=$(basename $f)
		echo $x
		name="${x%.sampled_seed*}"
		f1=$path$name".sampled_seed"$seed"_R1.fastq"
		f2=$path$name".sampled_seed"$seed"_R2.fastq"
		y=$path2$name".seed"$seed
		echo $y
		rm $y.align.report
		$mixcr_path align -b imgt.201822-5.sv4 -s hs -f -r $y.align.report $f1 $f2 $y.vdjca
		rm $y.assemble.report
		$mixcr_path assemble -f -a -OaddReadsCountOnClustering=true -r $y.assemble.report $y.vdjca $y.clna
		$mixcr_path exportClones -c $ch -f $y.clna $y.clones.txt
		$mixcr_path exportAlignments -f -cloneId -readIds -targets -targetSequences -defaultAnchorPoints -vGene -vHit -vHitScore -vAlignment -jGene -aaFeature CDR3 -nFeature CDR3 -nFeature CDR1 -nMutationsRelative CDR1 VRegion -nFeature FR2 -nMutationsRelative FR2 VRegion -nFeature CDR2 -nMutationsRelative CDR2 VRegion -nFeature FR3 -nMutationsRelative FR3 VRegion $y.clna $y.aln.txt
	done
done
