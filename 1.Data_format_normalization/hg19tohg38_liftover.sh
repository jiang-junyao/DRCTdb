for i in *.bed
do
~/software/liftOver $i ~/software/hg19ToHg38.over.chain.gz ../bg38/$i ../bg38/$i_umapped.bed
gzip ../bg38/$i
done