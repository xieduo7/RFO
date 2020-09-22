#!/usr/bin/perl -w
use strict;
use Data::Dumper;
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
#my $dest=shift;
my @speciesOrder = ($species1, $species2);
##################################################################################################################

## get the gene pairs which are best hit to each other.
my %hit; ## mark the hits of a query
my %bestPair; ## mark the RBH pairs
&getRBH($tab_file, $species1, $species2, \%hit, \%bestPair);

#my (@result,@copy_result); ## rank the hit pairs
my (@result); ## rank the hit pairs
my %genePair; ## mark the gene pairs determine by mcscan.
## output result and mark the gene pairs added.
foreach my $g1 (keys %bestPair) {
	foreach my $g2 (keys %{$bestPair{$g1}}) {
		next if ($genePair{$g1} || $genePair{$g2});
		($g1, $g2) = &getGeneOrder($g1, $g2);
		my ($align1, $align2, $id, $align, $score,$gene_syn) = @{$bestPair{$g1}{$g2}};
#        next if($gene_syn ne "Perfect");
		push @result, "$g1\t$g2\t$align1\t$align2\t$id\t$score\t$gene_syn\tL1\n";
#		push @copy_result, "$g2\t$g1\t$align2\t$align1\t$id\t$score\t$gene_syn\tL1\n";
# i       print "$g1\t$g2\n";
        delete $hit{$g1}; 
        delete $hit{$g2}; 
		$genePair{$g1} ++;
		$genePair{$g2} ++;
	}
}
my (%bestPairL2);
&getRBH2(\%hit,\%bestPairL2);
## add the gene pairs which past the cut off but not best hit the each other
foreach my $g21 (keys %bestPairL2){
    foreach my $g22 (keys %{$bestPairL2{$g21}}){
        next if ($genePair{$g21} || $genePair{$g22});
        ($g21, $g22) = &getGeneOrder($g21, $g22);
        my ($align1, $align2, $id, $align, $score,$gene_syn) = @{$bestPairL2{$g21}{$g22}};
        push @result, "$g21\t$g22\t$align1\t$align2\t$id\t$score\t$gene_syn\tL2\n";
 #       push @copy_result, "$g22\t$g21\t$align2\t$align1\t$id\t$score\t$gene_syn\tL2\n";               
        delete $hit{$g21};
        delete $hit{$g22};
        $genePair{$g21} ++;
        $genePair{$g22} ++;
    }
}
my (%bestPairL3);
&getRBH3(\%hit,\%bestPairL3);
foreach my $g31 (keys %bestPairL3){
    foreach my $g32 (keys %{$bestPairL3{$g31}}){
        next if ($genePair{$g31} || $genePair{$g32});
        ($g31, $g32) = &getGeneOrder($g31, $g32);
        my ($align1, $align2, $id, $align, $score,$gene_syn) = @{$bestPairL3{$g31}{$g32}};
#        print "$g31\t$g32\n";
        push @result, "$g31\t$g32\t$align1\t$align2\t$id\t$score\t$gene_syn\tL3\n";
 #       push @copy_result, "$g32\t$g31\t$align2\t$align1\t$id\t$score\t$gene_syn\tL3\n";
        delete $hit{$g31};
        delete $hit{$g32};
        $genePair{$g31} ++;
        $genePair{$g32} ++;        
    }
}
## add the gene pairs which past the cut off but not best hit the each otheri

## output result
#my $out1=$dest.$species1."-".$species2.".";
#my $out2=$dest.$species2."-".$species1.".gff.intersect.ortholog.overlap.tab.best";
#open OUT,">$out1" or die"Cannot find the file!\n";
#print OUT "#Referrence \tSubject \tR_ratio \tS_ratio \tIdentity \tSocre \tLevel \n";  
#print "#Referrence \tSubject \tR_ratio \tS_ratio \tIdentity \tSocre \tLevel \n";  
print "$_" foreach @result;
#close OUT;

#open OUT1,">$out2" or die"Cannot find the file!\n";
#print OUT1 "#Referrence \tSubject \tR_ratio \tS_ratio \tIdentity \tSocre \tLevel \n";  
#print OUT1 "$_" foreach @copy_result;
#close OUT1;
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
	my ($g1, $align1, $g2, $align2, $score, $id,$gene_syn,$state) = (split /\s+/)[0,5,1,6,3,4,8,9];
	next if ($state ne "duplication" || $gene_syn eq "NO_gene_synteny" || $g1 eq $g2);
#    next if $g1 eq $g2;


#    my $sp1 = (split /_/, $g1)[0];
    my $sp1 = (split /-/, $g1)[0];
#    my $sp2 = (split /_/, $g2)[0];
    my $sp2 = (split /-/, $g2)[0];
    next unless $species{$sp1} && $species{$sp2};
	my $align = ($align1 + $align2) / 2; ## take the average alignRate for sort
	push @{$hit{$g1}}, [$g2, $align1, $align2, $id, $align, $score, $gene_syn];
    push @{$hit{$g2}}, [$g1, $align2, $align1, $id, $align, $score, $gene_syn];
	push @{$hit_ref->{$g1}}, [$g2, $align1, $align2, $id, $align, $score, $gene_syn];
	push @{$hit_ref->{$g2}}, [$g1, $align2, $align1, $id, $align, $score, $gene_syn];
}
close IN;
## mark the best hit for each query
my %pair;
foreach my $g1 (keys %hit) {
#	@{$hit{$g1}} = sort {$b->[-4] <=> $a->[-4] or $b->[-2] <=> $a->[-2] or $b->[-5] <=> $a->[-5]} @{$hit{$g1}};
	@{$hit{$g1}} = sort {$b->[-2] <=> $a->[-2] or $b->[-4] <=> $a->[-4] or $b->[-5] <=> $a->[-5]} @{$hit{$g1}};
	foreach my $p (@{$hit{$g1}}) {
		my ($g2, $align1, $align2, $id, $align, $score,$gene_syn) = @$p;
        next if($gene_syn ne "Perfect" || $id < 30);
		$pair{$g1}{$g2} = [$align1, $align2, $id, $align, $score, $gene_syn];
		last;
	}
}
## restore RBH pairs
my %printMark;
foreach my $g1 (sort keys %pair) {
	foreach my $g2 (keys %{$pair{$g1}}) {
		next unless $pair{$g2}{$g1};
		($g1, $g2) = sort ($g1, $g2);

		my ($align1_1, $align1_2, $id11, $align11, $score11,$gene_syn11) = @{$pair{$g1}{$g2}};
		my ($align2_1, $align2_2, $id22, $align22, $score22,$gene_syn22) = @{$pair{$g2}{$g1}};
		my ($align1, $align2, $id, $align, $score,$gene_syn);
		if ($score22 > $score11) {
			($align1, $align2, $id, $align, $score,$gene_syn) = ($align2_2, $align2_1, $id22, $align22, $score22,$gene_syn22);	
		} else {
			($align1, $align2, $id, $align, $score,$gene_syn) = ($align1_1, $align1_2, $id11, $align11, $score11,$gene_syn11);
		}

		
		#my ($align1, $align2, $id, $align, $score) = @{$pair{$g1}{$g2}};
		unless ($printMark{$g1}{$g2}) {
			$bestPair_ref->{$g1}{$g2} = [$align1, $align2, $id, $align, $score,$gene_syn];
			$bestPair_ref->{$g2}{$g1} = [$align2, $align1, $id, $align, $score,$gene_syn];
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
####################
sub getRBH2{
my($hit,$bestPair_ref)=@_;
my %pair;
foreach my $g1 (keys %$hit) {
#       @{$hit{$g1}} = sort {$b->[-4] <=> $a->[-4] or $b->[-2] <=> $a->[-2] or $b->[-5] <=> $a->[-5]} @{$hit{$g1}};
        @{$hit->{$g1}} = sort {$b->[-2] <=> $a->[-2] or $b->[-4] <=> $a->[-4] or $b->[-5] <=> $a->[-5]} @{$hit->{$g1}};
        foreach my $p (@{$hit->{$g1}}) {
                my ($g2, $align1, $align2, $id, $align, $score,$gene_syn) = @$p;
        next if($gene_syn eq "Orphan");
                $pair{$g1}{$g2} = [$align1, $align2, $id, $align, $score, $gene_syn];
                last;
        }
}
## restore RBH pairs
my %printMark;
foreach my $g1 (sort keys %pair) {
        foreach my $g2 (keys %{$pair{$g1}}) {
                next unless $pair{$g2}{$g1};
                ($g1, $g2) = sort ($g1, $g2);

                my ($align1_1, $align1_2, $id11, $align11, $score11,$gene_syn11) = @{$pair{$g1}{$g2}};
                my ($align2_1, $align2_2, $id22, $align22, $score22,$gene_syn22) = @{$pair{$g2}{$g1}};
                my ($align1, $align2, $id, $align, $score,$gene_syn);
                if ($score22 > $score11) {
                        ($align1, $align2, $id, $align, $score,$gene_syn) = ($align2_2, $align2_1, $id22, $align22, $score22,$gene_syn22);
                } else {
                        ($align1, $align2, $id, $align, $score,$gene_syn) = ($align1_1, $align1_2, $id11, $align11, $score11,$gene_syn11);
                }


 #       next if($gene_syn ne "Perfect");
                #my ($align1, $align2, $id, $align, $score) = @{$pair{$g1}{$g2}};
                unless ($printMark{$g1}{$g2}) {
                        $bestPair_ref->{$g1}{$g2} = [$align1, $align2, $id, $align, $score,$gene_syn];
                        $bestPair_ref->{$g2}{$g1} = [$align2, $align1, $id, $align, $score,$gene_syn];
                }
                $printMark{$g1}{$g2} ++;
        }
}
}
######
sub getRBH3{
my($hit,$bestPair_ref)=@_;
my %pair;
foreach my $g1 (keys %$hit) {
#       @{$hit{$g1}} = sort {$b->[-4] <=> $a->[-4] or $b->[-2] <=> $a->[-2] or $b->[-5] <=> $a->[-5]} @{$hit{$g1}};
        @{$hit->{$g1}} = sort {$b->[-4] <=> $a->[-4] or $b->[-2] <=> $a->[-2] or $b->[-5] <=> $a->[-5]} @{$hit->{$g1}};
        foreach my $p (@{$hit->{$g1}}) {
                my ($g2, $align1, $align2, $id, $align, $score,$gene_syn) = @$p;
                $pair{$g1}{$g2} = [$align1, $align2, $id, $align, $score, $gene_syn];
        next if($gene_syn eq "DownNone" || $gene_syn eq "NO_gene_synteny" || $gene_syn eq "UpNone" ||$gene_syn eq "Orphan");
                last;
        }
}
## restore RBH pairs
my %printMark;
foreach my $g1 (sort keys %pair) {
        foreach my $g2 (keys %{$pair{$g1}}) {
                next unless $pair{$g2}{$g1};
                ($g1, $g2) = sort ($g1, $g2);

                my ($align1_1, $align1_2, $id11, $align11, $score11,$gene_syn11) = @{$pair{$g1}{$g2}};
                my ($align2_1, $align2_2, $id22, $align22, $score22,$gene_syn22) = @{$pair{$g2}{$g1}};
                my ($align1, $align2, $id, $align, $score,$gene_syn);
                if ($score22 > $score11) {
                        ($align1, $align2, $id, $align, $score,$gene_syn) = ($align2_2, $align2_1, $id22, $align22, $score22,$gene_syn22);
                } else {
                        ($align1, $align2, $id, $align, $score,$gene_syn) = ($align1_1, $align1_2, $id11, $align11, $score11,$gene_syn11);
                }


        #                print "$g1\t$g2\n";
     #   foreach my $p (@{$hit{$g1}}) {
                #my ($align1, $align2, $id, $align, $score) = @{$pair{$g1}{$g2}};
                unless ($printMark{$g1}{$g2}) {
                        $bestPair_ref->{$g1}{$g2} = [$align1, $align2, $id, $align, $score,$gene_syn];
                        $bestPair_ref->{$g2}{$g1} = [$align2, $align1, $id, $align, $score,$gene_syn];
                }
                $printMark{$g1}{$g2} ++;
        }
}


}

