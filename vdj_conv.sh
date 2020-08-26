vdjpath=/Users/koshka/Documents/favorov2/final/vdjtools-1.2.1/vdjtools-1.2.1.jar
m1=$1
dir=$2
java -jar $vdjpath Convert -S mixcr -m $m1 $dir/vdj
cd $dir
java -jar $vdjpath Correct -m metadata.txt vdj.corr
java -jar $vdjpath FilterNonFunctional -m metadata.txt vdj.corr.func
java -jar $vdjpath Decontaminate -m metadata.txt vdj.corr.func.csd
