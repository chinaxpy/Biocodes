#!/bin/bash
##  ## build index

##  ln -s  ../data/genome.fasta  ./genome.fasta
##  samtools faidx genome.fasta 
##  samtools index -@ 10 -c pa1K.sort.mkdup.bam
##  bwa index  genome.fasta  
##  gatk CreateSequenceDictionary R=iwgsc_refseqv2.1_assembly.fa O=iwgsc_refseqv2.1_assembly.dic 
##  hisat2-build -p 4 genome.fa genome
## 1. trim reads ###  fastp
## 2. bwa mem align to genome
## 3. add group and mark duplicates
########################       align to genome
hisat2_db=
bwa_db=
sampledir=/usrdata/users/xfli/work/wheat_Rye/01_cleandate/rye/
outdir=/usrdata/users/xfli/work/wheat_Rye/03_hisat2_align/
# mkdir -p $outdir

cat $1 | while read id
do
  arr=($id)
  fq1=${arr[0]}
  fq2=${arr[1]}
  sample=${arr[2]}
  echo "***********" Now analysis on $sample  "***********"
	date
	echo hisat2 mapping
hisat2 -p 12 --dta \
    -x $hisat2_db \
    -1 ${sampledir}$fq1 -2 ${sampledir}$fq2  \
    -S $outdir${sample}.sam  &> ${sampledir}.hisat2.log
	samtools sort -@ 12 -O bam  -o $outdir${sample}.sort.bam  $outdir${sample}.sam  
	echo hisat2 done 
    rm -f $outdir${sample}.sam  
	date
done


cat $1 | while read id
do
  arr=($id)
  fq1=${arr[0]}
  fq2=${arr[1]}
  sample=${arr[2]}
  echo "***********" Now analysis on $sample  "***********"
	date
	echo bwa mapping
bwa mem -t 5 -R '@RG\tID:baisehunchi\tSM:baisehunchi\tPL:illumina'  \
  $bwa_db \
  ${sampledir}$fq1   ${sampledir}$fq2 \
 |  samtools sort -@ 10 -m 1G -o $outdir${sample}.sort.bam  - 
	echo bwa done 
	date
done
#############################     align to genome done

## rnaseq测序应该加上标签
gatk  AddOrReplaceReadGroups -I ./pa1K.sort.bam  -O ./pa1K.sort.RG.bam  -SO coordinate -ID  pa1K  -LB  pa1K  -PL illumina -PU machine -SM pa1K
gatk  AddOrReplaceReadGroups -I ./pa2G.sort.bam  -O ./pa2G.sort.RG.bam  -SO coordinate -ID  pa2G  -LB  pa2G  -PL illumina -PU machine -SM pa2G
gatk  AddOrReplaceReadGroups -I ./po6K.sort.bam  -O ./po6K.sort.RG.bam  -SO coordinate -ID  po6K  -LB  po6K  -PL illumina -PU machine -SM po6K
gatk  AddOrReplaceReadGroups -I ./po6G.sort.bam  -O ./po6G.sort.RG.bam  -SO coordinate -ID  po6G  -LB  po6G  -PL illumina -PU machine -SM po6G

gatk --java-options "-Djava.io.tmpdir=./tmp -Xmx20G -Xms20G"  MarkDuplicates -I  ./pa1K.sort.RG.bam    -O  ./pa1K.sort.RG.mkrp.bam   --VALIDATION_STRINGENCY SILENT   -M ./pa1K._marked_dup_metrics.txt > pa1K.mkdup.log 
gatk --java-options "-Djava.io.tmpdir=./tmp -Xmx20G -Xms20G"  MarkDuplicates -I  ./pa2G.sort.RG.bam    -O  ./pa2G.sort.RG.mkrp.bam   --VALIDATION_STRINGENCY SILENT   -M ./pa2G._marked_dup_metrics.txt > pa2G.mkdup.log 
gatk --java-options "-Djava.io.tmpdir=./tmp -Xmx20G -Xms20G"  MarkDuplicates -I  ./po6K.sort.RG.bam    -O  ./po6K.sort.RG.mkrp.bam   --VALIDATION_STRINGENCY SILENT   -M ./po6K._marked_dup_metrics.txt > po6K.mkdup.log 
gatk --java-options "-Djava.io.tmpdir=./tmp -Xmx20G -Xms20G"  MarkDuplicates -I  ./po6G.sort.RG.bam    -O  ./po6G.sort.RG.mkrp.bam   --VALIDATION_STRINGENCY SILENT   -M ./po6G._marked_dup_metrics.txt > po6G.mkdup.log 

samtools index -c -@ 10  ./pa1K.sort.RG.mkrp.bam 
samtools index -c -@ 10  ./pa2G.sort.RG.mkrp.bam 
samtools index -c -@ 10  ./po6K.sort.RG.mkrp.bam 
samtools index -c -@ 10  ./po6G.sort.RG.mkrp.bam 

## 4. use HaplotypeCaller to generate gvcf
python3 gatk4.py --process 22 --input  ./pa1K.sort.RG.mkrp.bam   --ref ~/data/genome/cs2.1/iwgsc_refseqv2.1_assembly.fa --output ./pa1K 
python3 gatk4.py --process 22 --input  ./pa2G.sort.RG.mkrp.bam   --ref ~/data/genome/cs2.1/iwgsc_refseqv2.1_assembly.fa --output ./pa2G 
python3 gatk4.py --process 22 --input  ./po6K.sort.RG.mkrp.bam   --ref ~/data/genome/cs2.1/iwgsc_refseqv2.1_assembly.fa --output ./po6K 
python3 gatk4.py --process 22 --input  ./po6G.sort.RG.mkrp.bam   --ref ~/data/genome/cs2.1/iwgsc_refseqv2.1_assembly.fa --output ./po6G 


#####
## gatk4.py
#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor
def gatk4(chrom):
   subprocess.run(['/share/software/gatk-4.1.8.1/gatk','--java-options','-Xmx8G','HaplotypeCaller','-I',args.input,
                    '-O',args.output + '.' + chrom + '.g.vcf.gz','-R',args.ref,'--emit-ref-confidence','GVCF','-OVI','False', '-L',chrom],shell=False)
parser = argparse.ArgumentParser(description='输入参数如下:')
parser.add_argument('--process', '-P', type=int, help='进程数量，必要参数', required=True)
parser.add_argument('--input', '-I', help='bam file，必要参数', required=True)
parser.add_argument('--ref', '-R', help='reference genome，必要参数', required=True)
parser.add_argument('--output', '-O', help='sample name，必要参数', required=True)
args = parser.parse_args()
if __name__ == '__main__':
    commom_wheat_chrom = ['Chr1A','Chr2A','Chr3A','Chr4A','Chr5A','Chr6A','Chr7A',
                         'Chr1B','Chr2B','Chr3B','Chr4B','Chr5B','Chr6B','Chr7B',
                         'Chr1D','Chr2D','Chr3D','Chr4D','Chr5D','Chr6D','Chr7D','ChrUnknown']
    chrom_gvcf = []
    for chrom in commom_wheat_chrom:
        chrom_gvcf.append( args.output + '.' + chrom + '.g.vcf.gz')
    try:
        with ProcessPoolExecutor(max_workers=args.process) as pool:
            future1 = pool.map(gatk4,commom_wheat_chrom)
            print(list(future1))
    except Exception as e: 
        print(e)

    subprocess.run('/share/software/gatk-4.1.8.1/gatk GatherVcfs -RI TRUE -I ' + ' -I '.join(chrom_gvcf) + ' -O ' + args.output + '.g.vcf.gz',shell=True)
    # subprocess.run('rm -fr ' + ' '.join(chrom_gvcf),shell=True)
    subprocess.run('rm -fr ' + '*Chr*g.vcf.gz',shell=True)
#####


## 4. use HaplotypeCaller to generate gvcf
gatk --java-options "-Xmx10g" HaplotypeCaller -R ref/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -I bam/B97Y55-3.dedupped.bam    -O vcf/B97Y55-3.g.vcf.gz   -ERC GVCF
gatk --java-options "-Xmx10g" HaplotypeCaller -R ref/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -I bam/B97Y55-R3.dedupped.bam   -O vcf/B97Y55-R3.g.vcf.gz  -ERC GVCF
gatk --java-options "-Xmx10g" HaplotypeCaller -R ref/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -I bam/GK53-3.dedupped.bam      -O vcf/GK53-3.g.vcf.gz     -ERC GVCF
gatk --java-options "-Xmx10g" HaplotypeCaller -R ref/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -I bam/GK53-R5.dedupped.bam     -O vcf/GK53-R5.g.vcf.gz    -ERC GVCF
gatk --java-options "-Xmx10g" HaplotypeCaller -R ref/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -I bam/W4713.dedupped.bam       -O vcf/W4713.g.vcf.gz      -ERC GVCF

## 5. CombineGVCFs
gunzip *gz
gatk IndexFeatureFile -I pa2G.g.vcf
gatk IndexFeatureFile -I pa2G.g.vcf
ls *g.vcf > gvcf.list
#  gatk --java-options "-Xmx10g" CombineGVCFs -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa --variant pa1K.g.vcf --variant pa2G.g.vcf --variant po6G.g.vcf --variant  po6K.g.vcf -O ../03-gvcf/merged.g.vcf
gatk --java-options "-Xmx10g" CombineGVCFs -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa  -V gvcf.list -O ../03-gvcf/merged.g.vcf
## 6.1 GenotypeGVCFs on the combined gvcf
gatk --java-options "-Xmx10g" GenotypeGVCFs -R  /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa -V  ./merged.g.vcf     -O ./merged.raw.vcf
## 6.2 filter vcf file
#SNP

gatk --java-options "-Xmx10g" SelectVariants -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa -V  ./merged.raw.vcf --select-type SNP  -O  ./raw.snp.vcf
gatk --java-options "-Xmx10g" VariantFiltration -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa  -V  ./raw.snp.vcf -O ./filter.snp.vcf   --filter-name "SNP_Filter" --filter-expression "DP < 5.0 || FS > 60.0 || MQ < 40.0 || QD < 2.0 " --cluster-size 3 --cluster-window-size 100 
gatk --java-options "-Xmx10g" SelectVariants -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa  -V ./filter.snp.vcf --exclude-filtered  -O ./filtered.snp.vcf
#INDEL

gatk --java-options "-Xmx10g" SelectVariants -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa -V  ./merged.raw.vcf --select-type INDEL  -O  ./raw.indel.vcf
gatk --java-options "-Xmx10g" VariantFiltration -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa  -V  ./raw.indel.vcf -O ./filter.indel.vcf   --filter-name "INDEL_Filter" --filter-expression "DP < 5.0 || FS > 60.0 || MQ < 40.0 || QD < 2.0 " --cluster-size 3 --cluster-window-size 100 
gatk --java-options "-Xmx10g" SelectVariants -R /usrdata/users/xfli/data/genome/cs/iwgsc_refseqv2.1_assembly.fa  -V ./filter.indel.vcf --exclude-filtered  -O ./filtered.indel.vcf

######## 7_filt_variants.pl

 
###USAGE: perl $0 output.converted.vcf >output.final.vcf
####Perl script START
#!/usr/bin/perl -w
while(<>){
    if (/\A#CHROM.*/){
        print;
        last;
    }
}

while(<>){
    my @line = split;
    #next if (($line[6] ne 'PASS') or ($line[8] ne 'GT:AD:DP:GQ:PL'));
    next if ($line[6] ne 'PASS');
    my $mis = 0;
	for (9..$#line){
		my $j = $_;
        if ($line[$j] =~ /\.\/\./){
            $line[$j] = "./.";
            $mis ++;
            next;
        }   
		my @sample = split (/:/, $line[$j]);
        next if $sample[2] == 0;
		if ($sample[2] =~ /\d+/){
            if ($sample[2] < 5){
                $line[$j] = "./.";
                $mis ++;
            }
        }else{
            $line[$j] = "./.";
            $mis ++;
        }
	}
    next if $mis > 2;
	my $print = join("\t", @line);
	print "$print\n";
}
####Perl script END



# 8_修改vcf文件表头，##fileformat=VCFv4.2
# 使vcf文件材料名对应，去除下面的./ 以及 .sort.bam	
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	./pa1K.sort.bam	./pa2G.sort.bam	./po6G.sort.bam	./po6K.sort.bam
# Chr1A	1167	.	T	C	72.91	PASS	AC=2;AF=0.333;AN=6;DP=9;ExcessHet=0.4576;FS=0.000;MLEAC=2;MLEAF=0.333;MQ=60.00;QD=25.36;SOR=2.303	GT:AD:DP:GQ:PL	./.:0,0:0:.:0,0,0	0/0:1,0:1:3:0,3,15	1/1:0,2:2:6:84,6,0	0/0:6,0:6:3:0,3,193
