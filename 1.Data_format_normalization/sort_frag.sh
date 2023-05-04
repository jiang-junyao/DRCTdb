for i in *.gz
do
bgzip -@ 10 -d $i
done
for i in *.tsv
do
sort -k 1,1 -k 2,2n $i > $i.sorted.tsv
bgzip -@ 10 $i.sorted.tsv
tabix -p bed $i.sorted.tsv.gz
done


