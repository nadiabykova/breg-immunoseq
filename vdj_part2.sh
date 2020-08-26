source ./vdj_config.sh

declare -a sep=($o $n)

rm -r $results
mkdir $results

#4.1 naive_p analysis
mkdir $results$o
mkdir $results$o$naive_p_folder
Rscript $naive_p_script $stat_source $vdj_out$o $clones_in $sample_info $results$o$naive_p_folder
mkdir $results$o$naive_p_v_folder
Rscript $naive_p_v_script $stat_source $vdj_out$o $clones_in $sample_info $results$o$naive_p_v_folder

#4.2 in this folders do analysis: aa_features, clonality, genes
for i in "${!sep[@]}"
        do
                fold=${sep[$i]}
                echo $fold
                mkdir $results$fold

                mkdir $results$fold$clonality_folder
                Rscript $clonality_script $stat_source $vdj_out$fold $results$fold$clonality_folder

                mkdir $results$fold$aa_f_folder
                Rscript $aa_f_script $stat_source $vdj_out$fold $results$fold$aa_f_folder

                mkdir $results$fold$vgenes_folder
                Rscript $vgenes_script $stat_source $vdj_out$fold $sample_info $results$fold$vgenes_folder
                
                Rscript $combine_script1 $results$fold
        done
Rscript $combine_script1 $results

mkdir $results$rc_stat
Rscript $count_clones_script $sample_info $clones_in $vdj_out $results$rc_stat
