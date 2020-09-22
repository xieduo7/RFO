#!/usr/bin/perl
use strict;
use warnings;

die "Usage:$0 <muscle.out> \n" if @ARGV < 1;

my $out = shift;

#print STDERR "$out\n";
my $Pep1;
my $Pep2;

open IN,$out or die "$!";
$/=">";<IN>;$/="\n";

my $id1 = <IN>;
$id1 =~ s/^(\S+)\s*.*\n/$1/;
$/=">";
$Pep1 = <IN>;
chomp $Pep1;
$Pep1 =~ s/\s//g;
$/="\n";

my $id2 = <IN>;
$id2 =~ s/^(\S+)\s*.*\n/$1/;
$/=">";
$Pep2 = <IN>;
chomp $Pep2;
$Pep2 =~ s/\s//g;
$/="\n";

close IN;

my $Len = length($Pep1);

my $headGap;
my $tailGap;

$headGap=$1 if($Pep1 =~ /^([-]+)/);
$tailGap=$1 if($Pep1 =~ /([-]+)$/);

#($headGap,$tailGap)=(length $headGap,length $tailGap);
$headGap=$headGap?length $headGap:0;
$tailGap=$tailGap?length $tailGap:0;

#my $sum_pep=len
my $Pep_len1 = 0;
my $Pep_len2 = 0;
my $match = 0;
my $not_gap = 0;
for(my $i=0;$i<$Len;$i++){
	my $c1 = uc(substr($Pep1,$i,1));
	my $c2 = uc(substr($Pep2,$i,1));
	if($c1 ne '-'){
		$Pep_len1++;
	}
	if($c2 ne '-'){
		$Pep_len2++;
	}
	if ( $c1 ne '-' && $c2 ne '-'){

		if ($c1 eq $c2){
			$match++;
		}
		$not_gap++;
	}
}
unless($id1 && $id2 && $Pep_len1 && $Pep_len2){
    print STDERR "$out may be incompleted!\n";
}
#print "pep1\tpep2\tpep1_len\tpep2_len\tnot_gap\tnot_gap/pep2_len(%)\tmatch\tmatch/pep1_len(%)\tmatch/pep2_len(%)\theadGap\ttailGap\n";
print join("\t",$id1,$id2,$Pep_len1,$Pep_len2,$not_gap,$not_gap/$Pep_len2*100,$match,$match/$Len*100,$match/$Pep_len1*100,$match/$Pep_len2*100)."\t".$headGap."\t".$tailGap."\n";
