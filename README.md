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
   tabix -s 1 -b 2 -e 2 H.snp.filter.PASS.AroundIndel3.gz
   ```

2) GATK MuTect2 for tumor bam

      Call SNV and INDEL
   ```shell
   java -jar GenomeAnalysisTK.jar -R ref.fa -T MuTect2 -I:tumor T.srt.rmdup.realign.bqsr.bam -I:normal H.srt.rmdup.realign.bqsr.bam -o T.m2.vcf.gz
   ```

3) bcftools
