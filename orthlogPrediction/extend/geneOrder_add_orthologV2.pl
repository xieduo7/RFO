#!/usr/bin/perl

if (@ARGV != 3)
{
	print "perl $0 <.gff1> <.gff2> <.ortholog>\n";
	exit;
}

use strict;

my $gff1 = shift;
my $gff2 = shift;
my $orth = shift;
#my $cande = shift;
my (%chrGene1, %geneIndex,%can_gene);
&read_orth($orth,\%can_gene);
&read_gff_ChrGene($gff1, \%chrGene1, \%geneIndex,\%can_gene); # sub 1

my %chrGene2;
&read_gff_ChrGene($gff2, \%chrGene2, \%geneIndex,\%can_gene); #sub 1

my %ortholog_geneSyn;
&gene_synteny_add($orth,\%ortholog_geneSyn,\%geneIndex,\%chrGene1); # sub 2

foreach my $chr(sort keys %ortholog_geneSyn)
{
	my @tmp = sort {$a->[2] <=> $b->[2]} @{$ortholog_geneSyn{$chr}};
	for (my $i = 0; $i < scalar @tmp; $i++)
	{
		my $ln;
		for (my $j = 0; $j < scalar @{$tmp[$i]};$j++)
		{
			$ln .= "$tmp[$i][$j]\t" 
		}
		$ln .= "\n";
		print $ln;
	}
}




## sub 1
##
sub read_orth{
    my($table,$gene_h)=@_;
    open IN,"<$table" or die"!";
    while(<IN>){
        chomp;
        my @t=split;
        #print "$t[0]\t$t[1]\n";
        my $gene_id1 = (split(/-/,$t[0]))[1];
        my $gene_id2 = (split(/-/,$t[1]))[1];
        $gene_h->{$gene_id1} = 1;
        $gene_h->{$gene_id2} = 1;
    }
    close IN;

}
sub read_gff_ChrGene
{

	my $gff_file = shift;
	my $chr_hshp = shift;
	my $ind_hshp = shift;
    my $gene_hs =shift;
	open IN, $gff_file or die "Fail to open $gff_file\n";
	while (<IN>) 
	{
        	chomp;
        	my @info = split /\t/;
        	next unless ($info[2] eq 'mRNA');
        	die unless $info[8] =~ /ID=(\S+?);/;
        	my $gene = $1;
 #           my $gene_id = (split(/-/,$gene))[1];
            next unless(exists $gene_hs->{$gene});
        	push @{$chr_hshp->{$info[0]}}, [$gene, $info[3], $info[4], $info[6]];
	}
	close IN;

	foreach my $chr (sort keys %{$chr_hshp}) 
	{
        	@{$chr_hshp->{$chr}} = sort {$a->[1] <=> $b->[1]} @{$chr_hshp->{$chr}};
        	for (my $i = 0; $i < @{$chr_hshp->{$chr}}; $i ++) 
		{
                	my ($gene, $bg, $ed, $strand) = @{$chr_hshp->{$chr}->[$i]};
                	my $index = $i + 1;
                	$ind_hshp->{$gene} = [$chr, $bg, $ed, $strand, $index];
        	}
	}

}

## sub 2
## $orth,\%ortholog_geneSyn,\%geneIndex,\%chrGene1
sub gene_synteny_add
{
	my $file = shift;
	my $hshp = shift;
	my $hsp2 = shift;
	my $hsp3 = shift;
	
	my %ind;
	open IN, $file or die "Fail to open $file\n";
    my($spe1,$spe2);
	while (<IN>)
	{
		next if ($_ =~ m/#/);
#		my ($gn1, $gn2, $level) = (split /\s+/)[0,1,-1];
		my ($gn1, $gn2) = (split /\s+/)[0,1];
		my ($g1,$g2) = ($gn1,$gn2);
        $g1 =~ /(\w+)-(\S+)/;
       # $g1 =~ /(\w+)-(\w+)/;
        my $g1_id=$2;
       # print "$g1_id\n";
        $spe1=$1;
        $g2 =~ /(\w+)-(\w+)/;
        my $g2_id=$2;
        #print "$g1_id\t$g2_id\n";
        $spe2=$1;
		#$g1 =~ s/_[A-Z]{5}//;
		#$g2 =~ s/_[A-Z]{5}//;
		my ($chr1, $bg1, $ed1, $strand1, $index1) = @{$hsp2->{$g1_id}};
#		my ($chr1, $bg1, $ed1, $strand1, $index1) = @{$hsp2->{$g1}};
#		die "$g2" unless ($hsp2->{$g2});
		die "$g2_id" unless ($hsp2->{$g2_id});
#		my ($chr2, $bg2, $ed2, $strand2, $index2) = @{$hsp2->{$g2}};
		my ($chr2, $bg2, $ed2, $strand2, $index2) = @{$hsp2->{$g2_id}};
		$ind{$chr1}{$index1} = 1;
#		push @{$hshp->{$chr1}}, [$gn1,$chr1,$index1,$strand1,$bg1,$ed1,$gn2,$chr2,$index2,$strand2,$bg2,$ed2,$level];
		push @{$hshp->{$chr1}}, [$gn1,$chr1,$index1,$strand1,$bg1,$ed1,$gn2,$chr2,$index2,$strand2,$bg2,$ed2];
	}
	close IN;

	foreach my $chr1 (keys %{$hsp3}) 
	{
	        foreach my $p (@{$hsp3->{$chr1}}) {
        	        my ($g1, $bg1, $ed1, $strand1) = @$p;
                	my $index1 = $hsp2->{$g1}->[-1];
	                next if ($ind{$chr1}{$index1});
        	        push @{$hshp->{$chr1}}, ["$spe1-$g1",$chr1,$index1,$strand1,$bg1,$ed1,"NA","NA","NA","NA","NA","NA","NA"];
        	}
	}
}


