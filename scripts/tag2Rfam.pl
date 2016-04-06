#!/usr/bin/perl
# for small RNA solexa sequencing
use FindBin qw ($Bin);

if (@ARGV != 3) {
	print "usage: $0 <tag.fas> <Rfam.fas> <result_dir>\n";
	print "* Please formatdb Rfam first.\n";
	print "* Request match percent >= 90%\n";
	print "* only consider rRNAs, tRNAs ,snRNAs and snoRNAs\n";
	exit;
}

my $tag_file=shift;   # fa
my $Rfam_file=shift;  # fa
my $work_dir=shift;   # result dir

# work dir
$work_dir.="/" unless ($work_dir=~/\/$/);
mkdir $work_dir unless (-e $work_dir);

my $formatdb="formatdb";
my $blastall="blastall";
my $blast2table="$Bin/blast2table.pl";

my $blast_file=$work_dir."match_Rfam.blast";
my $table_file=$work_dir."match_Rfam.blast.tab";
my $result_file=$work_dir."match_Rfam.txt";
my $stat_file=$work_dir."match_Rfam.stat";

system("$formatdb -i $Rfam_file -p F -o") unless (-e "$Rfam_file.nhr");
system("$blastall -p blastn -i $tag_file -d $Rfam_file -F F -e 0.01 -o $blast_file");
system("$blast2table -D $blast_file > $table_file");

my %tagId_tag;
my %tagId_number;
&read_tag_file($tag_file,\%tagId_tag,\%tagId_number);

# filter random and redundancy match
my %tagId_type;
open IN, $table_file || die $!;
open OUT, ">$result_file" || die $!;
while (<IN>) {
	chomp;
	my @d=split;
	next unless ($d[4]/$d[5] >= 0.9); # identity < 0.9, skip
#	next if ($d[9] > $d[10]); # match minus strand, skip
	
	my $tagId=$d[0];
	if (not defined $tagId_type{$tagId}) {
		my $type;
		if ($d[14]=~/mir|let-7|lin-4|bantam|lsy|iab/i) {
			$type="miRNA";
		}
		elsif ($d[14]=~/tRNA/) {
			$type="tRNA";
		}
		elsif ($d[14]=~/rRNA/) {
			$type="rRNA";
		}
		elsif ($d[14]=~/sno/i) {
			$type="snoRNA";
		}
		elsif ($d[14]=~/U\d{1,2}/) {
			$type="snRNA"; # ???
		}
		else {
			$type="undef";
		}
		$tagId_type{$d[0]}=$type;
		my $length=length $tagId_tag{$d[0]};

		# only output non-miRNA
		next if ($type eq "miRNA" || $type eq "undef");
		print OUT join("\t",$d[0],$length,$tagId_number{$d[0]},$tagId_tag{$d[0]},$type,$d[1]),"\n";
	}
}
close IN;
close OUT;

# stat
my $match_Rfam_c=0;
my $match_Rfam_n=0;
my %type_c;
my %type_n;
foreach my $tagId (keys %tagId_type) {
	my $type=$tagId_type{$tagId};
	next if ($type eq "miRNA" || $type eq "undef");
	++$match_Rfam_c;
	$match_Rfam_n+=$tagId_number{$tagId};
	++$type_c{$type};
	$type_n{$type}+=$tagId_number{$tagId};
}
open OUT, ">$stat_file" || die $!;
foreach my $type (sort keys %type_c) {
	next if ($type eq "miRNA" || $type eq "undef");
	print OUT join("\t",$type,$type_c{$type},$type_n{$type}),"\n";
}
#print OUT "only consider rRNA,tRNA,snRNA,snoRNA\n";
print OUT "total\t$match_Rfam_c\t$match_Rfam_n\n";
close OUT;


sub read_tag_file{
	my $infile=shift;
	my $tagId_tag=shift;
	my $tagId_number=shift;
	
	open IN, $infile || die $!;
	while (<IN>) {
		chomp;
		if(/^>(\S+)\s+(\d+)/){
			my $id=$1;
			my $number=$2;
			my $tag = <IN>;
			chomp $tag;
			$tagId_tag->{$id}=$tag;
			$tagId_number->{$id}=$number;
		}
	}
	close IN;
}

