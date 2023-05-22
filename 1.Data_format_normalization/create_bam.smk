import glob
import os

directory_path = "bam_pe"

pattern = "*.bedpe.gz"
bedpe_files = []
for file_path in glob.glob(os.path.join(directory_path, pattern)):
    bedpe_files.append(os.path.splitext(os.path.basename(file_path))[0])

genome = '~/sample20_test/hg38.genome'
rule all:
  input:
    expand("bam_pe/{sample}.gz",sample=bedpe_files),
    expand("sortedbam/{sample}.bam",sample=bedpe_files),
    expand("sortedbam/{sample}.bam.bai",sample=bedpe_files)


rule bedpe2bam:
  input:
    'bam_pe/{sample}.gz'
  output:
    temp('bam/{sample}.bam')
  threads: 4
  shell:
    'bedtools bedpetobam -i {input} -g {genome} > {output}'


rule samtools_sort:
  input:
    temp('bam/{sample}.bam')
  output:
    'sortedbam/{sample}.bam'
  threads: 4
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule samtools_index:
  input:
    'sortedbam/{sample}.bam'
  output:
    'sortedbam/{sample}.bam.bai'
  shell:
    "samtools index {input}"
