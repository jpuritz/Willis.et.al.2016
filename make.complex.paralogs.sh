#!/bin/bash

#Function to create fixed differences replacing the 5th G in each sequence with an A in even sequences and T in odd sequences
fixdiff(){
j=0
while read NAME;
do
RL=$(echo $NAME | mawk '/[AAT,CCG]/' )
if [ -n "$RL" ]; 
then
	if [ $(($j%2)) -eq 0 ];
        	then
		echo $NAME | sed -e 's/G/A/5'
        	j=$(($j + 1))
	else
		echo $NAME | sed -e 's/G/T/4'
		j=$(($j + 1))
	fi
else 
echo $NAME
fi
done < $1
}

#Loop through 20 individuals and create paralogs by copying sequences from "paralog" populations for the first 50 loci
for i in {0..19};
do
rm p.fq q.fq p1.fq q1.fq 2>/dev/null
	for j in {1..50};
	do
	zcat PopC_${i}.F.fq.gz | grep -A3 "data${j}_"  >> p.fq &
	zcat PopC_${i}.R.fq.gz | grep -A3 "data${j}_"  >> q.fq &
	zcat PopD_${i}.F.fq.gz | grep -A3 "data${j}_"  >> p1.fq &
	zcat PopD_${i}.R.fq.gz | grep -A3 "data${j}_"  >> q1.fq
	wait
	done

#Create fixed differences in the paralog loci
fixdiff p.fq | sed 's/fake/faaake/g' > pp.fq
fixdiff q.fq | sed 's/fake/faaake/g' > qq.fq
fixdiff p1.fq | sed 's/fake/faaake/g' > pp1.fq
fixdiff q1.fq | sed 's/fake/faaake/g' > qq1.fq

wait
#Change name of sequences to differentiate paralogs
cat p.fq | sed 's/fake/faaake/g' > pp.fq &
cat q.fq | sed 's/fake/faaake/g' > qq.fq &
cat p1.fq | sed 's/fake/faaake/g' > pp1.fq &
cat q1.fq | sed 's/fake/faaake/g' > qq1.fq

wait 
#Copy sequences to original individuals
gzip -c pp.fq >> PopA_${i}.F.fq.gz 	
gzip -c qq.fq >> PopA_${i}.R.fq.gz
gzip -c pp1.fq >> PopB_${i}.F.fq.gz 	
gzip -c qq1.fq >> PopB_${i}.R.fq.gz
wait
done

rm *.fq
