#### paths #####################
vdjpath=/Users/koshka/Documents/favorov2/final/vdjtools-1.2.1/vdjtools-1.2.1.jar

#### init files and folders ####
clones_in="fastq_final_results/"
sample_info="sample_info_final.txt"

#### metadata file #############
m1="metadata_final.convert.txt"
antibodies="antibodies.aln.txt"

echo $antibodies

#### out folders ###############
#level 1
vdj_out=fastq_final_vdj/
results=IBH_29112019/
#level 2
o="orig/"
n="naive3/"
mem="mem3/"
n0="naive0/"
mem0="mem0/"
rc_stat="read_statistics/"

# level 3
naive_p_folder=naive_p/
naive_p_v_folder=naive_p_v/
clonality_folder=clonality/
aa_f_folder=aa_features/
vgenes_folder=vgenes/


#### R scripts #################
metadata_script=scripts_final/make_metadata.R
sep_script=scripts_final/sep_naive.R
stat_source=scripts_final/stat_test.util.R
naive_p_script=scripts_final/naive_p.R
naive_p_v_script=scripts_final/naive_p_v.R
aa_f_script=scripts_final/aaFeatures.R
clonality_script=scripts_final/clonality.R
vgenes_script=scripts_final/vgenes.R
combine_script1=scripts_final/combine1.R
cdr13_script=scripts_final/filter_fullTarget.R
fixer_script=scripts_final/aln_fixer.R
decont_script=scripts_final/decontamination.R
count_clones_script=scripts_final/countClones1.R

declare -a chains=("IGH" "IGK" "IGL")
declare -a sep=($o $n $mem $n0 $mem0)
