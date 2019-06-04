#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 
use Math::Round;

my $chrnum=$ARGV[-1];
print $chrnum,"\n";
open my $ai3,"tabix $ARGV[0] $chrnum|" or die ("Can't open $ARGV[0]\n");
my %ai3;
print "Host ai3\n";
while (<$ai3>) {
    my @tmp=split /\t/;
    $ai3{$tmp[1]}++;
}
close $ai3;
open my $gatkhom,"bcftools view -r $chrnum -g hom -H $ARGV[1]|" or die ("Can't open $ARGV[1]\n");
print "Host hom\n";
my %gatkhom;
while (<$gatkhom>) {
    my @tmp=split /\t/;
    $gatkhom{$tmp[1]}=$tmp[4];
}
close $gatkhom;
open my $gatkhet,"bcftools view -r $chrnum -g het -H $ARGV[1]|" or die ("Can't open $ARGV[1]\n");
print "Host het\n";
my %gatkhet;
while (<$gatkhet>) {
    my @tmp=split /\t/;
    if ($tmp[4]=~/(\w),(\w)/) {
        $gatkhet{$tmp[1]}{$1}++;
        $gatkhet{$tmp[1]}{$2}++;
    }else {
        $gatkhet{$tmp[1]}{$tmp[3]}++;
        $gatkhet{$tmp[1]}{$tmp[4]}++;
    }
}
close $gatkhet;
open my $gatkdel,"bcftools view -H -v indels -r $chrnum $ARGV[2]|vcf2bed -n -d|cut -d \";\" -f 1|";
print "Host del\n";
my %gatkdel;
while (<$gatkdel>) {
    chomp;
    my @tmp=split /\t/;
    if ($tmp[-1] eq "AC=2"){
        for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
            $gatkdel{$i}{c}+=2;
            $gatkdel{$i}{g}=1;
        }
    }else {
        for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
            $gatkdel{$i}{c}++;
            $gatkdel{$i}{g}=1;
        }
    }
}
close $gatkdel;
open my $tdel,"bcftools view -H -r $chrnum -v indels $ARGV[6]|vcf2bed -n -d|" or die ("Can't open $ARGV[6]\n");
my %tumordel;
print "Tumor m2.del\n";
while (<$tdel>) {
    my @tmp=split /\t/;
    my @format=split /:/,$tmp[10];
    $format[1]=~/,(\d+)/;
    my $tumor_num=$1;
    @format=split /:/,$tmp[11];
    $format[1]=~/(\d+),(\d+)/;
    my $host_c=0;
    $host_c=round($2*2/($1+$2)) unless ($1+$2)==0;
    for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
        $tumordel{$i}{c}+=$tumor_num;
        $tumordel{$i}{g}=1;
        if ($host_c != 0) {
            $gatkdel{$tmp[1]}{c}+=$host_c unless exists $gatkdel{$tmp[1]}{g};
        }
    }
}
close $tdel;
open my $tbdel2,"bcftools view -H -M 2 -r $chrnum $ARGV[7]|vcf2bed -n -d|" or die ("Can't open $ARGV[7]\n");
print "Tumor bcftools.del2\n";
while (<$tbdel2>) {
    my @tmp=split /\t/;
    my @format=split /:/,$tmp[10];
    $format[2]=~/,(\d+)/;
    my $tumor_num=$1;
    for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
        unless (exists $tumordel{$i}) {
            $tumordel{$i}{c}+=$tumor_num;
            $tumordel{$i}{g}=1;
        }
    }
}
close $tbdel2;
open my $tbdelm,"bcftools view -H -m 3 -r $chrnum $ARGV[7]|" or die ("Can't open $ARGV[7]\n");
print "Tumor bcftools.delm\n";
while (<$tbdelm>) {
    my @tmp=split /\t/;
    my $reflength=length($tmp[3]);
    my @alt=split /,/,$tmp[4];
    my %index;
    for (my $i=0;$i<$#alt;$i++) {
        my $altlength=length($alt[$i]);
        $index{$i}=$altlength if ($altlength<$reflength);
    }
    my @format=split /:/,$tmp[9];
    my @AD=split /,/,$format[2];
    for my $index (sort keys %index) {
        for (my $i=($tmp[1]+$index{$index});$i<($tmp[1]+$reflength);$i++) {
            unless (exists $tumordel{$i}{g}) {
	$tumordel{$i}{c}+=$AD[$index+1];
            }
        }
    }
}
close $tbdelm;
open my $hcnv,"tabix $ARGV[3] $chrnum|" or die("Can't open $ARGV[3]\n");
print "Host cnv\n";
my %hostcnv;
while (<$hcnv>) {
    my @tmp=split /\t/;
    chomp $tmp[3];
    for (my $i=($tmp[1]+1);$i<=$tmp[2];$i++) {
        $hostcnv{$i}=$tmp[3];
    }
}
close $hcnv;
my %final;
my %amb;
open my $seqz,"tabix $ARGV[4] $chrnum|" or die("Can't open $ARGV[4]\n");
print "seqz\n";
while (<$seqz>) {
    my @tmp=split /\t/;
    if (exists $ai3{$tmp[1]}) {
        if (exists $gatkdel{$tmp[1]}) {
            next unless $gatkdel{$tmp[1]}{c}==1;
        }else {
            next;
        }
    }
    if ($tmp[8] eq "hom") {
        if (exists $gatkhom{$tmp[1]}) {
            if ($tmp[11] eq $gatkhom{$tmp[1]}) {
	$final{$tmp[1]}{h}=$tmp[11]; # alt allele hom
            }else {
	$amb{$tmp[1]}{0}++; # amb0
            }
        }elsif (exists $gatkhet{$tmp[1]}) {
            $amb{$tmp[1]}{2}++; # amb2
        }elsif (exists $gatkdel{$tmp[1]}) {
            $gatkdel{$tmp[1]}{b}=$tmp[11]; # het del remeaning base
        }else {
            $final{$tmp[1]}{h}=$tmp[11]; # ref allele hom
        }
    }else {		# final het
        if (exists $gatkhet{$tmp[1]}) {
            $tmp[11]=~/(\w)(\w)/;
            if ((exists $gatkhet{$tmp[1]}{$1}) and (exists $gatkhet{$tmp[1]}{$2})) {
	$final{$tmp[1]}{e}{$1}++; # het
	$final{$tmp[1]}{e}{$2}++;
            }else {
	$amb{$tmp[1]}{1}++; # amb1
            }
        }elsif (exists $gatkdel{$tmp[1]}) {
            $amb{$tmp[1]}{5}++;	# amb5
        }else {
            $amb{$tmp[1]}{3}++; # amb3
        }
    }
}
close $seqz;
for my $pos (sort keys %gatkhom) {
    $amb{$pos}{4}++ unless (exists $final{$pos} or exists $amb{$pos}); # amb4
}
for my $pos (sort keys %gatkhet) {
    $amb{$pos}{4}++ unless (exists $final{$pos} or exists $amb{$pos}); # amb4
}
undef %gatkhet;
undef %gatkhom;
open my $cnv,"tabix $ARGV[5] $chrnum|" or die ("Can't open $ARGV[5]\n");
print "Tumor cnv\n";
my %cnv;
while (<$cnv>) {
    my @tmp=split /\t/;
    for (my $i=$tmp[1];$i<$tmp[2];$i++) {
        $cnv{$i}=$tmp[9];
    }
}
close $cnv;
open my $atcg,"tabix $ARGV[8] $chrnum|" or die ("Can't open $ARGV[8]\n");
my $sample=`basename $ARGV[4] .seqz.gz|cut -d "." -f 2`;
chomp $sample;
open my $out,"|bgzip >T.$sample.$chrnum.baseCN.gz";
print "Tumor atcg\n";
my $cellularity=$ARGV[9];
while (<$atcg>) {
    my @tmp=split /\t/;
    next if exists $amb{$tmp[1]};
    next unless (exists $final{$tmp[1]} or exists $gatkdel{$tmp[1]});
    if (exists $cnv{$tmp[1]}) {
        my $hostRN;
        my $hostcnv;
        if (exists $hostcnv{$tmp[1]}) {
            $hostcnv=$hostcnv{$tmp[1]};
        }else {
            $hostcnv=2;
        }
        my $total=$tmp[4]+$tmp[5]+$tmp[6]+$tmp[7];
        if (exists $tumordel{$tmp[1]}) {
            $total+=$tumordel{$tmp[1]}{c};
        }
        if (exists $gatkdel{$tmp[1]}) {
            if ($gatkdel{$tmp[1]}{c}==1) {
	$hostRN=$total*($hostcnv/2)*(1-$cellularity)/(($hostcnv/2)*(1-$cellularity)+$cnv{$tmp[1]}*$cellularity);
            }
        }else {
            $hostRN=$total*$hostcnv*(1-$cellularity)/($hostcnv*(1-$cellularity)+$cnv{$tmp[1]}*$cellularity);
        }
        my @hostgeno;
        if (exists $gatkdel{$tmp[1]}) {
            if ($gatkdel{$tmp[1]}{c}==1) {
            	if (exists $gatkdel{$tmp[1]}{b}) {
	    @hostgeno="d";
	}else {
	    $amb{$tmp[1]}{6}++;# amb6
	    next;
	}
            }else {
	@hostgeno="d";
            }
        }else {
            @hostgeno=keys $final{$tmp[1]};
        }
        my %tumor_rawRN;
        $tumor_rawRN{A}=$tmp[4];
        $tumor_rawRN{C}=$tmp[5];
        $tumor_rawRN{G}=$tmp[6];
        $tumor_rawRN{T}=$tmp[7];
        if ($hostgeno[0] eq 'h') {
            $tumor_rawRN{$final{$tmp[1]}{h}}=$tumor_rawRN{$final{$tmp[1]}{h}}-$hostRN;
        }elsif ($hostgeno[0] eq 'e') {
            for my $allele (sort keys $final{$tmp[1]}{e}) {
	$tumor_rawRN{$allele}=$tumor_rawRN{$allele}-($hostRN/2);
            }
        }elsif ($hostgeno[0] eq 'd') {
            if ($gatkdel{$tmp[1]}{c}==1) {
	$tumor_rawRN{$gatkdel{$tmp[1]}{b}}=$tumor_rawRN{$gatkdel{$tmp[1]}{b}}-$hostRN;
            }
        }
        my $sum;
        for my $allele (sort keys %tumor_rawRN) {
            $sum+=$tumor_rawRN{$allele};
        }
        my $realcnv=$cnv{$tmp[1]};
        if ($sum==0) {
            $realcnv=0;
        }
        if (exists $tumordel{$tmp[1]}) {
            $realcnv=round($realcnv*$sum/($sum+$tumordel{$tmp[1]}{c}));
        }
        print $out $chrnum,"\t",$tmp[1],"\t",$tmp[2],"\t";
        if ($sum != 0) {
            for my $allele (sort keys %tumor_rawRN) {
	$tumor_rawRN{$allele}=round($tumor_rawRN{$allele}/$sum*$realcnv);
	print $out $tumor_rawRN{$allele},"\t";
            }
        }else {
            print $out "0\t0\t0\t0\t";
        }
        print $out $realcnv,"\n";
    }else {
        print $out $chrnum,"\t",$tmp[1],"\t",$tmp[2],"\t0\t0\t0\t0\t0\n";
    }
}
close $atcg;
close $out;
open my $amb,"|bgzip >T.$sample.$chrnum.amb.gz";
print "amb\n";
for my $pos (sort keys %amb) {
    for my $type (sort keys $amb{$pos}) {
        print $amb $type,"\t",$chrnum,"\t",$pos,"\n";
    }
}
close $amb;
