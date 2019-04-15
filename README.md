# ttgeno
transmissible tumor genotyper

This tool is designed to genotype each genomic site of a transmissible tumor unbiasedly.


Assuming WGS of one transmissible tumor and the corresponding host were performed, and reads were aligned to reference genome respectively. Two bam files were then manipulated by sorting, duplication mark, indel realign and BQSR in GATK. Two bam files' final names are T.srt.rmdup.realign.bqsr.bam (Tumor) and H.srt.rmdup.realign.bqsr.bam (Host).

All input files for ttgeno should be prepared in the following steps.

1) GATK HaplotypeCaller for host bam

      Call SNPs
   ```shell
   java -jar GenomeAnalysisTK.jar -R ref.fa -T HaplotypeCaller -I H.srt.rmdup.realign.bqsr.bam -o H.vcf.gz
   java -jar GenomeAnalysisTK.jar -R ref.fa -T SelectVariants -V H.vcf.gz -selectType SNP -o H.snp.vcf.gz
   java -jar GenomeAnalysisTK.jar -R ref.fa -T VariantFiltration -V H.snp.vcf.gz -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "GATK experience hard filtration" -o H.snp.filter.vcf.gz
   bcftools view -f PASS -O z -o H.snp.filter.PASS.vcf.gz H.snp.filter.vcf.gz
   bcftools index -t H.snp.filter.PASS.vcf.gz
    ```   
      Call INDELs
   ```shell
   java -jar GenomeAnalysisTK.jar -R ref.fa -T SelectVariants -V H.vcf.gz -selectType INDEL -o H.indel.vcf.gz
   java -jar GenomeAnalysisTK.jar -R ref.fa -T VariantFiltration -V H.indel.vcf.gz -filter "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0 " --filterName "GATK experience hard filtration" -o H.indel.filter.vcf.gz
   bcftools view -f PASS -O z -o H.indel.filter.PASS.vcf.gz H.indel.filter.vcf.gz
   bcftools index -t H.indel.filter.PASS.vcf.gz
   ```
      Extract complex loci around INDELs in 3 bps
   ```shell
   bcftools concat -O z -o H.snpNindel.filter.PASS.vcf.gz -a H.snp.filter.PASS.vcf.gz H.indel.filter.PASS.vcf.gz
   bcftools index -t H.snpNindel.filter.PASS.vcf.gz
   bcftools filter -O z -o H.snpNindel.filter.PASS.SnpGap3.vcf.gz -g 3 H.snpNindel.filter.PASS.vcf.gz
   bcftools index -t H.snpNindel.filter.PASS.SnpGap3.vcf.gz
   bcftools view -O z -o H.snp.filter.PASS.SnpGap3.vcf.gz -v snps H.snpNindel.filter.PASS.SnpGap3.vcf.gz
   bcftools index -t H.snp.filter.PASS.SnpGap3.vcf.gz
   bcftools isec -o H.snp.filter.PASS.AroundIndel3 -C H.snp.filter.PASS.vcf.gz H.snp.filter.PASS.SnpGap3.vcf.gz
   bgzip H.snp.filter.PASS.AroundIndel3
   tabix -f -s 1 -b 2 -e 2 H.snp.filter.PASS.AroundIndel3.gz
   ```

2) GATK MuTect2 for tumor bam

      Call SNV and INDEL
   ```shell
   java -jar GenomeAnalysisTK.jar -R ref.fa -T MuTect2 -I:tumor T.srt.rmdup.realign.bqsr.bam -I:normal H.srt.rmdup.realign.bqsr.bam -o T.m2.vcf.gz
   ```

3) bcftools mpileup for tumor bam

      Call INDEL
   ```shell
   bcftools mpileup -Q 20 -d 500 -L 500 -a AD -Ou -f ref.fa T.srt.rmdup.realign.bqsr.bam | bcftools call -mv -Oz -V snps -o T.indel.vcf.gz
   bcftools index -t T.indel.vcf.gz
   ```

4) cnvnator for host bam

   ```shell
   cnvnator -root H.400.root -chrom 1 2 3 4 ... X -unique -tree H.srt.rmdup.realign.bqsr.bam
   cnvnator -root H.400.root -chrom 1 2 3 4 ... X -d ref.fa.split.dir -his 400
   cnvnator -root H.400.root -chrom 1 2 3 4 ... X -stat 400
   cnvnator -root H.400.root -chrom 1 2 3 4 ... X -partition 400
   cnvnator -root H.400.root -chrom 1 2 3 4 ... X -call 400 | gzip - > H.400.cnv.gz
   ```
      Filter by p and q value
   ```shell
   gunzip H.400.cnv.gz | awk '{if($5<0.01 && $9<0.5) print $0}'| gzip - > H.400.cnv.p001q05.gz
   ```
      Filter by overlap with gaps in reference
      Gap anotation can be downloaded from UCSC ftp, and changed to bed format by yourself.
   ```shell
   gunzip H.400.cnv.p001q05.gz | awk '{print $2":"$1}'|awk -F '[:-]' '{if($3-$2>999) print $1"\t"$2-1"\t"$3"\t"$4}'| bedops -n 50% - ref.gap.bed > H.400.cnv.p001q05.gap50.bed
   ```
      Filter by overlap with repeatmask regions in reference
      Repeatmask anotation can be downloaded from UCSC ftp, as well as nestedrepeat anotation. I merged them and changed to bed format.
   ```shell
   bedops -n 50% H.400.cnv.p001q05.gap50.bed ref.repeatmask.bed > H.400.cnv.p001q05.gap50.rpmk50.bed
   ```
      Regenotype the filtered region
   ```shell
   awk '{print $1":"$2"-"$3}END{print "exit"}' H.400.cnv.p001q05.gap50.rpmk50.bed | cnvnator -root H.400.root -genotype 400|cut -d " " -f 2,4|sed 's/[-: ]/\t/g' | bgzip > H.400.cnv.p001q05.gap50.rpmk50.cn.gz
   tabix -f -0 -s 1 -b 2 -e 3 -c male,X H.400.cnv.p001q05.gap50.rpmk50.cn.gz
   ```

5) Sequenza and sequenza-utils
   
   stat gc content wiggle
   ```shell
   sequenza-utils gc_wiggle -f ref.fa -w 50 -o ref.fa.gc50base.txt.gz
   ```
   generate seqz file
   ```shell
   sequenza-utils bam2seqz -n H.srt.rmdup.realign.bqsr.bam -T T.srt.rmdup.realign.bqsr.bam -gc ref.fa.gc50base.txt.gz -F ref.fa |bgzip > T.seqz.gz
   tabix -f -s 1 -b 2 -e 2 -S 1 T.seq.gz
   sequenza-utils seqz_binning -w 50 -s T.seqz.gz | bgzip > T.bin50.seqz.gz
   ```
   Estimate ploidy, contamination and local cnv in R.
   Set the parameter of sex `female=T or F`
   ```R
   library("sequenza")
   test<- sequenza.extract("T.bin50.seqz.gz",normalization.method = "median",breaks.method= "fast",max.mut.types=3)
   CP.example <- sequenza.fit(test,mc.cores=getOption("mc.cores",8L),female=T,chromosome.list=test$chromosomes)
   sequenza.results(sequenza.extract=test,cp.table=CP.example,sample.id="T",out.dir="T.f.mutmax3",chromosome.list=test$chromosomes,female=T)
   ```
   Check results whether the best solution is correspond to prior knowledge or result from immunohistochemistry. Then you need compress the cnv file.
   ```shell
   sed 's/"//g' T_segments.txt| bgzip > T_segments.txt.gz
   tabix -f -S 1 -s 1 -b 2 -e 3 T_segments.txt.gz
   ```
   Generate the acgt file for tumor
   ```shell
   samtools mpileup -Q 20 -f ref.fa T.srt.rmdup.realign.bqsr.bam |sequenza-utils pileup2acgt -p -|bgzip > T.acgt.gz
   tabix -S 1 -s 1 -b 2 -e 2 T.acgt.gz
   ```
   
When you finish all steps above, you can use the ttgeno.pl to generate the per-site genotype of the tumor now.
All parameters of ttgeno must be as the same order as the description below.
```shell
perl ttgeno.pl H.snp.filter.PASS.AroundIndel3.gz H.snp.filter.PASS.SnpGap3.vcf.gz H.indel.filter.PASS.vcf.gz H.400.cnv.p001q05.gap50.rpmk50.cn.gz T.seqz.gz T_segments.txt.gz T.m2.vcf.gz T.indel.vcf.gz T.acgt.gz cellularity_estimated_in_sequenza chromosome_name
```
Results are two files. The T.sampleid.chr.baseCN.gz contains per-site copy number of tumor. The T.sampleid.chr.amb.gz contains amphibolous sites, which can be abandoned if sites are not too much.
