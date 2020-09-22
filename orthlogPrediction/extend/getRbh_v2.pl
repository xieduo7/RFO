#!/usr/bin/perl -w
use strict;
unless (@ARGV == 4) {
	die <<End;

Description:
    This script is used to get ortholog pairs.

Usage:
    perl $0 <idAdd> <ref_species> <hit_species> > ref_hit.idAdd.ort

End
}
my $tab_file = shift;
my $species1 = shift;
my $species2 = shift;
my $dest=shift;
my @speciesOrder = ($species1, $species2);
##################################################################################################################

## get the gene pairs which are best hit to each other.
my %hit; ## mark the hits of a query
my %bestPair; ## mark the RBH pairs
&getRBH($tab_file, $species1, $species2, \%hit, \%bestPair);

my (@result,@copy_result); ## rank the hit pairs
my %genePair; ## mark the gene pairs determine by mcscan.
## output result and mark the gene pairs added.
foreach my $g1 (keys %bestPair) {
	foreach my $g2 (keys %{$bestPair{$g1}}) {
		next if $genePair{$g1} || $genePair{$g2};
		($g1, $g2) = &getGeneOrder($g1, $g2);
		my ($align1, $align2, $id, $align, $score) = @{$bestPair{$g1}{$g2}};
		push @result, "$g1\t$g2\t$align1\t$align2\t$id\t$score\tL1\n";
		push @copy_result, "$g2\t$g1\t$align2\t$align1\t$id\t$score\tL1\n";
		$genePair{$g1} ++;
		$genePair{$g2} ++;
	}
}

## add the gene pairs which past the cut off but not best hit the each otheri

## output result
my $out1=$dest.$species1."-".$species2.".gff.intersect.ortholog.overlap.tab.best";
my $out2=$dest.$species2."-".$species1.".gff.intersect.ortholog.overlap.tab.best";
open OUT,">$out1" or die"Cannot find the file!\n";
print OUT "#Referrence \tSubject \tR_ratio \tS_ratio \tIdentity \tSocre \tLevel \n";  
print OUT "$_" foreach @result;
close OUT;

open OUT1,">$out2" or die"Cannot find the file!\n";
print OUT1 "#Referrence \tSubject \tR_ratio \tS_ratio \tIdentity \tSocre \tLevel \n";  
print OUT1 "$_" foreach @copy_result;
close OUT1;
##subroutine
##################################################################
sub getRBH{
my ($tab_file, $species1, $species2, $hit_ref, $bestPair_ref) = @_;
my %species = ($species1=>1, $species2=>1);
## restore the hits for each query
my %hit;
open IN, $tab_file;
while (<IN>) {
    chomp;
	next if /^#/;
	my ($g1, $align1, $g2, $align2, $score, $id) = (split /\s+/)[0,5,1,6,3,4];
    next if $g1 eq $g2;


#    my $sp1 = (split /_/, $g1)[0];
    my $sp1 = (split /-/, $g1)[0];
#    my $sp2 = (split /_/, $g2)[0];
    my $sp2 = (split /-/, $g2)[0];
    next unless $species{$sp1} && $species{$sp2};
	my $align = ($align1 + $align2) / 2; ## take the average alignRate for sort
	push @{$hit{$g1}}, [$g2, $align1, $align2, $id, $align, $score];
    push @{$hit{$g2}}, [$g1, $align2, $align1, $id, $align, $score];
	push @{$hit_ref->{$g1}}, [$g2, $align1, $align2, $id, $align, $score];
	push @{$hit_ref->{$g2}}, [$g1, $align2, $align1, $id, $align, $score];
}
close IN;
## mark the best hit for each query
my %pair;
foreach my $g1 (keys %hit) {
	@{$hit{$g1}} = sort {$b->[-1] <=> $a->[-1] or $b->[-3] <=> $a->[-3] or $b->[-4] <=> $a->[-4]} @{$hit{$g1}};
	foreach my $p (@{$hit{$g1}}) {
		my ($g2, $align1, $align2, $id, $align, $score) = @$p;
		$pair{$g1}{$g2} = [$align1, $align2, $id, $align, $score];
		last;
	}
}
## restore RBH pairs
my %printMark;
foreach my $g1 (sort keys %pair) {
	foreach my $g2 (keys %{$pair{$g1}}) {
		next unless $pair{$g2}{$g1};
		($g1, $g2) = sort ($g1, $g2);

		my ($align1_1, $align1_2, $id11, $align11, $score11) = @{$pair{$g1}{$g2}};
		my ($align2_1, $align2_2, $id22, $align22, $score22) = @{$pair{$g2}{$g1}};
		my ($align1, $align2, $id, $align, $score);
		if ($score22 > $score11) {
			($align1, $align2, $id, $align, $score) = ($align2_2, $align2_1, $id22, $align22, $score22);	
		} else {
			($align1, $align2, $id, $align, $score) = ($align1_1, $align1_2, $id11, $align11, $score11);
		}

		
		#my ($align1, $align2, $id, $align, $score) = @{$pair{$g1}{$g2}};
		unless ($printMark{$g1}{$g2}) {
			$bestPair_ref->{$g1}{$g2} = [$align1, $align2, $id, $align, $score];
			$bestPair_ref->{$g2}{$g1} = [$align2, $align1, $id, $align, $score];
		}
		$printMark{$g1}{$g2} ++;
	}
}
}

###############################
sub getGeneOrder {
my ($g1, $g2) = @_;
#my $sp1 = (split /_/, $g1)[0];
my $sp1 = (split /-/, $g1)[0];
#my $sp2 = (split /_/, $g2)[0];
my $sp2 = (split /-/, $g2)[0];
my %sp_to_gene = ($sp1=>$g1, $sp2=>$g2);
($g1, $g2) = ($sp_to_gene{$speciesOrder[0]}, $sp_to_gene{$speciesOrder[1]});
return ($g1, $g2);
}
