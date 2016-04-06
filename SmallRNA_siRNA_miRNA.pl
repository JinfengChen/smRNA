#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"ref:s","repeat:s","rfam:s","mirbase:s","step:s","read:s","cpu:s","linker3:s","linker5:s","help");


my $help=<<USAGE;
perl $0 --read ERRR000077.fastq.gz
Prepare sequence (0) and trim quality and adaptor for smRNA reads (1), filter r/t/sn/snoRNA (2), filter miRNA (3) and map reads to genome (4) and non-redundency repeat (4).
1.The resulting fastq should be siRNA in plant. We map this onto genome and TE to analyze siRNA and TE interaction.
2.miRNA expression level can be estimate by other software, like miRNAexpress, miRNAkey or mirDeep2.
3.novel miRNA can be predicted by mirDeep2.

OUTPUT:
Trimed fastq files of each trims
SRR771499_1.trim3.fastq.gz
SRR771499_1.trim3_5.fastq.gz
SRR771499_1.trim3_5.rfam.fastq.gz
SRR771499_1.trim3_5.rfam.mirbase.fastq.gz

Bam files map SRR771499_1.trim3_5.rfam.mirbase.fastq.gz to genome and/or TE
SRR771499_1.trim3_5.rfam.mirbase.genome.bam
SRR771499_1.trim3_5.rfam.mirbase.repeat.bam

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{ref}   ||="/rhome/cjinfeng/BigData/03.SmallRNA/reference/MSU_r7.fa";
$opt{repeat}||="/rhome/cjinfeng/BigData/03.SmallRNA/reference/Ping_Pong.fa";
$opt{rfam}  ||="$Bin/database/Rfam/Rfam.fasta";
$opt{mirbase} ||="$Bin/database/mirbase/hairpin_osa.fa";
$opt{cpu}   ||=1;
$opt{step}  ||="0";
=adaptor
yizhou 3' adaptor TGGAATTCTCGGGTGCCAAGGC
Zimberman 3' adaptor CTGTAGGCACCATCAAT  5' adaptor ACACTCTTTCCCTACACGACGCTGTTCCATCT
Illumina 3'(Illumina_Small_RNA_3p_Adapter_1) ATCTCGTATGCCGTCTTCTGCTTG 5' (Illumina_Small_RNA_Adapter_1) GTTCAGAGTTCTACAGTCCGACGATC
Rice pollen 3' TCGTATGCCGTCTTCTGCTTG 5' GTTCAGAGTTCTACAGTCCGACGATC
YijunQi 3' adaptor CTGTAGGCACCATCAAT 5' adaptor GTTCAGAGTTCTACAGTCCGACGATC
=cut
$opt{linker3} ||="CTGTAGGCACCATCAAT"; ###Illumina_Small_RNA_3p_Adapter_1
$opt{linker5} ||="GTTCAGAGTTCTACAGTCCGACGATC"; ###Illumina_Small_RNA_Adapter_1
my $bwa="~/BigData/00.RD/RelocaTE2/tools/bwa-0.6.2/bwa";
my $SAMtool="/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools";
my $bam="/rhome/cjinfeng/BigData/software/bamUtil_1.0.9/bin/bam";
my $fastx_clipper="/opt/linux/centos/7.x/x86_64/pkgs/fastx_toolkit/0.0.13/bin/fastx_clipper";
my $prefix="";

if ($opt{read}=~/gz$/){
   $prefix=basename($opt{read},".fastq.gz");
}else{
   $prefix=basename($opt{read},".fastq");
}

#######Step 0###############################################
#Index reference genome, Rfam and mirbase hairpin
#clean Rfam, only keep rRNA/tRNA/snoRNA
###########################################################
if($opt{step}=~/0/){
  print "Step0: index database\n";
  if (-e $opt{ref} and !(-e "$opt{ref}.sa")){
     print "index reference genome\n";
     `$bwa index $opt{ref}`;
  }
  if (-e $opt{repeat} and !(-e "$opt{repeat}.sa")){
     print "index reference repeat\n";
     `$bwa index $opt{repeat}`;
  }
  if (-e $opt{rfam}){
     my $headrfam = $1 if ($opt{rfam}=~/(.*)\.fa*/);
     unless (-e "$headrfam.clean.fasta.sa"){
        print "index Rfam genome\n";
        cleanRfam($opt{rfam},"$headrfam.clean.fasta");
        `$bwa index $headrfam.clean.fasta`;
     }
     $opt{rfam}="$headrfam.clean.fasta";
  }
  if (-e $opt{mirbase} and !(-e "$opt{mirbase}.sa")){
     print "index mirbase genome\n";
     `$bwa index $opt{mirbase}`;
  }
}

#######Step 1###############################################
#Quality and adaptor trim
############################################################
if($opt{step}=~/1/ and !(-e "$prefix.trim3_5.fastq.gz")){
  print "Step1: Trim adaptor: $opt{read} ......\n";
  if ($opt{read}=~/gz$/){
      $prefix=basename($opt{read},".fastq.gz");
      `zcat $opt{read} | $fastx_clipper -a $opt{linker3} -l 16 -Q 33 -v -c -z -o $prefix.trim3.fastq.gz > $prefix.trim3.log 2> $prefix.trim3.log2` unless (-e "$prefix.trim3.fastq.gz");
      `zcat $prefix.trim3.fastq.gz | $fastx_clipper -a $opt{linker5} -l 16 -Q 33 -v -z -o $prefix.trim3_5.fastq.gz > $prefix.trim3_5.log 2> $prefix.trim3_5.log2` unless (-e "$prefix.trim3_5.fastq.gz");
      print "Done\n";
  }else{
      $prefix=basename($opt{read},".fastq");
      `cat $opt{read} | $fastx_clipper -a $opt{linker3} -l 16 -Q 33 -v -c -z -o $prefix.trim3.fastq.gz > $prefix.trim3.log 2> $prefix.trim3.log2` unless (-e "$prefix.trim3.fastq.gz");
      `zcat $prefix.trim3.fastq.gz | $fastx_clipper -a $opt{linker5} -l 16 -Q 33 -v -z -o $prefix.trim3_5.fastq.gz > $prefix.trim3_5.log 2> $prefix.trim3_5.log2` unless (-e "$prefix.trim3_5.fastq.gz");
  }
}

#######Step 2##############################################
#Rfam trim
###########################################################
if($opt{step}=~/2/ and !(-e "$prefix.trim3_5.rfam.fastq.gz")){
  print "Step2: Filter r/t/sn/snoRNA by mapping to Rfam: $opt{read} ......\n";
  `$bwa aln -t $opt{cpu} $opt{rfam} $prefix.trim3_5.fastq.gz > $prefix.trim3_5.rfam.sai` unless (-e "$prefix.trim3_5.rfam.sai");
  `$bwa samse $opt{rfam} $prefix.trim3_5.rfam.sai $prefix.trim3_5.fastq.gz > $prefix.trim3_5.rfam.sam` unless (-e "$prefix.trim3_5.rfam.sam");
  `awk '\$2==4' $prefix.trim3_5.rfam.sam | $bam bam2FastQ --in - --outBase $prefix.trim3_5.rfam --readname` unless (-e "$prefix.trim3_5.rfam.fastq.gz");
  `gzip $prefix.trim3_5.rfam.fastq` unless (-e "$prefix.trim3_5.rfam.fastq.gz");
  #`rm $prefix.trim3_5.rfam_1.fastq $prefix.trim3_5.rfam_2.fastq $prefix.trim3_5.rfam.sa*`;
  print "Done\n";
}

#######Step 3##############################################
#mirbase trim
###########################################################
if($opt{step}=~/3/ and !(-e "$prefix.trim3_5.rfam.mirbase.fastq.gz")){
  print "Step3: Filter miRNA by mapping to mirbase hairpin: $opt{read} ......\n";
  `$bwa aln -t $opt{cpu} $opt{mirbase} $prefix.trim3_5.rfam.fastq.gz > $prefix.trim3_5.rfam.mirbase.sai` unless (-e "$prefix.trim3_5.rfam.mirbase.sai");
  `$bwa samse $opt{mirbase} $prefix.trim3_5.rfam.mirbase.sai $prefix.trim3_5.rfam.fastq.gz > $prefix.trim3_5.rfam.mirbase.sam` unless (-e "$prefix.trim3_5.rfam.mirbase.sam");
  `awk '\$2==4' $prefix.trim3_5.rfam.mirbase.sam | $bam bam2FastQ --in - --outBase $prefix.trim3_5.rfam.mirbase --readname` unless (-e "$prefix.trim3_5.rfam.mirbase.fastq.gz");
  `gzip $prefix.trim3_5.rfam.mirbase.fastq` unless (-e "$prefix.trim3_5.rfam.mirbase.fastq.gz");
  #`rm $prefix.trim3_5.rfam.mirbase_1.fastq $prefix.trim3_5.rfam.mirbase_2.fastq $prefix.trim3_5.rfam.mirbase.sa*`;
  print "Done\n";
}

#######Step 4##############################################
#map the clean fastq to reference and/or TE
###########################################################
if ($opt{step}=~/4/){
  ###genome mapping
  if (-e $opt{ref} and !(-e "$prefix.trim3_5.rfam.genome.bam")){
    print "Step4: Mapping reads to genome ......\n";
    `$bwa aln -t $opt{cpu} $opt{ref} $prefix.trim3_5.rfam.fastq.gz > $prefix.trim3_5.rfam.genome.sai` unless (-e "$prefix.trim3_5.rfam.genome.sai");
    `$bwa samse $opt{ref} $prefix.trim3_5.rfam.genome.sai $prefix.trim3_5.rfam.fastq.gz > $prefix.trim3_5.rfam.genome.sam` unless (-e "$prefix.trim3_5.rfam.genome.sam");
    `$SAMtool view -bS -o $prefix.trim3_5.rfam.genome.new.bam $prefix.trim3_5.rfam.genome.sam` unless (-e "$prefix.trim3_5.rfam.genome.new.bam");
    `$SAMtool sort $prefix.trim3_5.rfam.genome.new.bam $prefix.trim3_5.rfam.genome` unless (-e "$prefix.trim3_5.rfam.genome.bam");
    `$SAMtool index $prefix.trim3_5.rfam.genome.bam`;
    `rm $prefix.trim3_5.rfam.genome.sa* $prefix.trim3_5.rfam.genome.new.bam $prefix.trim3_5.rfam_*.fastq $prefix.trim3_5.rfam.sa*`;
    print "Done\n";
  }
  ###transposon mapping
  if (-e $opt{repeat} and !(-e "$prefix.trim3_5.rfam.repeat.bam")){
    print "Step4: Mapping reads to repeat ......\n";
    `$bwa aln -n 0 -t $opt{cpu} $opt{repeat} $prefix.trim3_5.rfam.fastq.gz > $prefix.trim3_5.rfam.repeat.sai` unless (-e "$prefix.trim3_5.rfam.repeat.sai");
    `$bwa samse $opt{repeat} $prefix.trim3_5.rfam.repeat.sai $prefix.trim3_5.rfam.fastq.gz > $prefix.trim3_5.rfam.repeat.sam` unless (-e "$prefix.trim3_5.rfam.repeat.sam");
    `$SAMtool view -bS -o $prefix.trim3_5.rfam.repeat.new.bam $prefix.trim3_5.rfam.repeat.sam` unless (-e "$prefix.trim3_5.rfam.repeat.new.bam");
    `$SAMtool sort $prefix.trim3_5.rfam.repeat.new.bam $prefix.trim3_5.rfam.repeat` unless (-e "$prefix.trim3_5.rfam.repeat.bam");
    `$SAMtool index $prefix.trim3_5.rfam.repeat.bam`;
    `rm $prefix.trim3_5.rfam.repeat.sa* $prefix.trim3_5.rfam.repeat.new.bam $prefix.trim3_5.rfam_*.fastq $prefix.trim3_5.rfam.sa*`;
    print "Done\n";
  }
}



###########################################################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}


sub cleanRfam
{
$/=">";
my %hash;
my ($file,$out)=@_;
print "in cleanRfam: $file\t$out\n";
open OUT, ">$out" or die "$!";
open IN,"$file" or die "$!";
while (<IN>){
    print "$_\n";
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\s//g;
    print "$head\n";
    if ($head=~/rRNA/ or $head =~/tRNA/ or $head =~/sno/i or $head=~/;U\d{1,2}.*?\;/){
       print OUT ">$head\n$seq\n";
    }
}
close IN;
close OUT;
$/="\n";
}
 
