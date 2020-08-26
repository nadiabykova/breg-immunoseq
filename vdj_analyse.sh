vdjpath=/Users/koshka/Documents/favorov2/final/vdjtools-1.2.1/vdjtools-1.2.1.jar
cd $1
java -jar $vdjpath CalcDiversityStats -m metadata.txt stat
java -jar $vdjpath CalcBasicStats -m metadata.txt basic
java -jar $vdjpath CalcCdrAaStats -r CDR3-full,cdr3-center-5,cdr3-center-3,VJ-junc,V-germ,J-germ -a strength,kf10,turn,cdr3contact,rim,alpha,beta,polarity,charge,surface,hydropathy,count,mjenergy,volume,core,disorder,kf2,kf1,kf4,kf3,kf6,kf5,kf8,kf7,kf9 -n -w -m metadata.txt aas
java -jar $vdjpath CalcCdrAaStats -r CDR3-full,cdr3-center-5,cdr3-center-3,VJ-junc,V-germ,J-germ -a strength,kf10,turn,cdr3contact,rim,alpha,beta,polarity,charge,surface,hydropathy,count,mjenergy,volume,core,disorder,kf2,kf1,kf4,kf3,kf6,kf5,kf8,kf7,kf9 -w -m metadata.txt aas
