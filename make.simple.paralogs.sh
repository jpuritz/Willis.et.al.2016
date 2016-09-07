#!/bin/bash

#Loop through 20 individuals and create paralogs by copying sequences from "paralog" populations for the first 50 loci
for i in {0..19};
do
rm p.fq q.fq p1.fq q1.fq 2>/dev/null
	for j in {1..50};
	do
	zcat PopC_${i}.F.fq.gz | grep -A3 "data${j}_" | head -100 >> p.fq &
	zcat PopC_${i}.R.fq.gz | grep -A3 "data${j}_" | head -100 >> q.fq &
	zcat PopD_${i}.F.fq.gz | grep -A3 "data${j}_" | head -100 >> p1.fq &
	zcat PopD_${i}.R.fq.gz | grep -A3 "data${j}_" | head -100 >> q1.fq
	wait
	done


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
