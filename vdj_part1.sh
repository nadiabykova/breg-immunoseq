source ./vdj_config.sh 

#### bash scripts ##############
vdj_conv_bash="./vdj_conv.sh"
vdj_analyse_bash="./vdj_analyse.sh"

#0
Rscript $cdr13_script $clones_in $fixer_script
Rscript $decont_script $clones_in $antibodies

#1 make meta data
Rscript $metadata_script $clones_in $sample_info $m1

#2 convert to vdi + correct, filter nc - decontaminate!!
rm -r $vdj_out
mkdir $vdj_out
$vdj_conv_bash $m1 $vdj_out

#3 separate into folders -orig,naive and memory
Rscript $sep_script $vdj_out $clones_in $o $of $n $mem $n0 $mem0

#4 vdjtools downsample
declare -a sep=($o $n $mem $n0 $mem0)
declare -a chains=("IGH" "IGK" "IGL")

for i in "${!sep[@]}"
	do
		fold=${sep[$i]}
		echo $fold
		cp $vdj_out"metadata.txt" $vdj_out$fold"metadata.txt"
		cd $vdj_out$fold
		java -jar $vdjpath SplitMetadata -c chain metadata.txt ./
		for ch in "${chains[@]}"
			do
				echo $ch
				xmin=$(cat min.count.$ch.txt)
				java -jar $vdjpath Downsample -m metadata.$ch.txt -x $xmin vdj.corr.func.csd.down
				mv metadata.txt metadata.$ch.txt
			done
		ls metadata.IG* | xargs -n 1 tail -n +2 > metadata-1.txt
		head -n 1 metadata.IGH.txt > header.txt
		cat header.txt metadata-1.txt > metadata.txt
		rm metadata-1.txt
		rm header.txt
		cd ../../
		$vdj_analyse_bash $vdj_out$fold
	done

