#!/usr/bin/perl
# blast2table.pl
# version 4.0
# liqb@genomics.org.cn

use strict;

##### options
use Getopt::Std;
use vars qw($opt_p $opt_s $opt_e $opt_a $opt_D $opt_h);
getopts('p:s:e:aDh');
my $pi_cut = $opt_p ? $opt_p : 0;
my $bits_cut = $opt_s ? $opt_s : 0;
my $expect_cut = $opt_e ? $opt_e : 1e30;
my $ALIGN = $opt_a ? $opt_a : 0;
my $DESC = $opt_D ? $opt_D : 0;
my $help     = $opt_h ? $opt_h : 0;
#####


##### usage
my $usage= << "USAGE";

blast2table: convert standard blast result into tabular format
Usage: blast2table.pl options *.blast
Version: 4.0, 2007-9-7
Options:
	-p <int>   indentity percent cutoff, default: 0
	-s <int>   bit score cutoff, default: 0
	-e <float> evalue cutoff value, default: 1e30
	-a         print align
	-D         print query and target description
	-h         output help information
Author: liqb\@genomics.org.cn

USAGE

if ($help || @ARGV==0) {
	print $usage;
	exit;
}
#####

##### parse
my $qry_id; # query id
my $qry_d; # query description
my $qry_l; # query length
my $sbj_id;
my $sbj_d;
my $sbj_l;
my $HSP; 

my $report_c=0;
while (<>) {
	if (/^T?BLAST[NPX]/) { # start of a new report
		++$report_c;
		print_HSP($HSP) if (defined $HSP);
		undef $HSP;
		undef $qry_id;
		undef $qry_d;
		undef $qry_l;
		undef $HSP;
	}
	elsif (/^Query=\s+(\S+)\s*(.*)/){ # query info
		$qry_id=$1;
		#### query description
		$qry_d=$2;
		$qry_d.=" ";
		while (<>) {
			if (/([\d,]+) letters/) {
				$qry_l=$1;
				$qry_l=~s/\D//g;
				last;
			}
			$qry_d.=$_;
		}
		$qry_d =~ s/\s+/ /g;
		$qry_d=~s/\s*$//;
		####

#		print "qry_id: $qry_id\n";
#		print "qry_d: $qry_d\n";
#		print "qry_l: $qry_l\n";
	}
	elsif(/^>(\S+)\s*(.*)/){ # start of a new subject
		print_HSP($HSP) if (defined $HSP);
		undef $HSP;
		undef $sbj_id;
		undef $sbj_d;
		undef $sbj_l;
		
		$sbj_id=$1;
		#### subject description
		$sbj_d=$2;
		$sbj_d.=" ";
		while (<>) {
			if (/Length = ([\d,]+)/) {
				$sbj_l=$1;
				$sbj_l=~s/\D//g;
				last;
			}
			$sbj_d.=$_;
		}
		$sbj_d =~ s/\s+/ /g;
		$sbj_d=~s/\s*$//;
		####

#		print "sbj_id: $sbj_id\n";
#		print "sbj_d: $sbj_d\n";
#		print "sbj_l: $sbj_l\n";
	}
	elsif (/^ Score = /) { # start of a new HSP
		print_HSP($HSP) if (defined $HSP);
		undef $HSP;
		
		#### HSP stat
		my $stat=$_;
		while (<>) {
			last unless (/\S/);
			$stat.=$_;
		}
        my ($bits)=$stat=~/(\d\S+) bits/;
        my ($expect)=$stat=~/Expect\S* = ([\d\-\.e]+)/;
        my ($match, $total, $percent)=$stat=~/Identities = (\d+)\/(\d+) \((\d+)%\)/;
		#my $mismatch = $total - $match;
#		print $stat;
#		print "bits:$bits, expect:$expect, match:$match, total:$total, percent:$percent\n\n";
		####

		# strand or frame info
		# Strand = Plus / Plus  [blastn]
		# Frame = -2 / -3       [tblastx]
		#						[blastp]
		# Frame = -2			[blastx]
		# Frame = +2            [tblastn]
		#my($strand)=$stat=~/Strand = (\w+ \/ \w+)/;		

		#### create a new HSP
		$HSP = {qid => $qry_id, ql => $qry_l, qd => $qry_d,
			sid => $sbj_id, sl => $sbj_l, sd => $sbj_d,
			bits => $bits, expect => $expect, match => $match,
            percent => $percent, total => $total, 
			qb => 0, qe => 0, qa => "",
            sb => 0, se => 0, sa => ""};
	}
	elsif (/^Query:\s+(\d+)\s+(\S+)\s+(\d+)/) {
		$HSP->{qb}  = $1 unless ($HSP->{qb});
		$HSP->{qe}  = $3;
		$HSP->{qa} .= $2;
	}
    elsif (/^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/) {
		$HSP->{sb}  = $1 unless ($HSP->{se});
		$HSP->{se}  = $3;
		$HSP->{sa} .= $2;
	}
	elsif (/^\s+Database:/){ # end of each report or end of entire blast
		print_HSP($HSP) if (defined $HSP); undef $HSP;
	}
}

#### print_HSP
sub print_HSP{
	my $hsp=shift;
	return if ($hsp->{percent} < $pi_cut);
	return if ($hsp->{bits} < $bits_cut);
	return if ($hsp->{expect} > $expect_cut);
	if ($DESC) {
		print join("\t", $hsp->{qid},$hsp->{sid},
			$hsp->{percent},$hsp->{total},$hsp->{match},
			$hsp->{ql},$hsp->{qb},$hsp->{qe},
			$hsp->{sl},$hsp->{sb},$hsp->{se},
			$hsp->{expect},$hsp->{bits},
			$hsp->{qd},$hsp->{sd}), "\n";
	}
	else {
		print join("\t", $hsp->{qid},$hsp->{sid},
			$hsp->{percent},$hsp->{total},$hsp->{match},
			$hsp->{ql},$hsp->{qb},$hsp->{qe},
			$hsp->{sl},$hsp->{sb},$hsp->{se},
			$hsp->{expect},$hsp->{bits}), "\n";
	}
	if ($ALIGN) {
		printf "Query: %10d $hsp->{qa} %-10d\n", $hsp->{qb},$hsp->{qe};
		printf "Sbjct: %10d $hsp->{sa} %-10d\n", $hsp->{sb},$hsp->{se};
		print "//\n";
	}
	undef $HSP;
}
