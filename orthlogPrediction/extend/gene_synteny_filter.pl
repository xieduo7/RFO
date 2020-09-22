#!/usr/bin/env  perl
#####################################
## Duo Xie                         ##
## xieduo1907@gmail.com            ##
############################################################################################################
## Identifying orthologs and lineage-specific genes using the reference-free whole-genome alignment       ##
## Duo Xie                                                                                                ##
## Version:2.0 2020-04-05                                                                                 ##
############################################################################################################
=head1 NAME

gene_synteny_filter  --  filter the gene pair with gene synteny

=head1 Usage

    perl gene_synteny_filter.pl ortholog genesynteny

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
die `pod2text $0` if (@ARGV == 0 );
my ($ortholog,$genesynteny)=@ARGV;
###########################
####### MAIN SCRIPT #######
###########################
my %pair;
read_synteny($genesynteny,\%pair);
open INN,"<$ortholog" or die"!";
while(<INN>){
    chomp;
    my @arr=split;
    if($arr[-1] eq "ortholog_one2one"){
        print "$_\tortholog_one2one\n";
    }else{
        if(exists $pair{$arr[0]}){
            if(exists $pair{$arr[0]}{$arr[1]}){
                print"$_\t$pair{$arr[0]}{$arr[1]}\n";
            }else{
                print "$_\tNO_gene_synteny\n";
            }
        }else{
            print "$_\tNO_gene_synteny\n";
        }
    }
}
close INN;
#########################
####### FUNCTIONS #######
#########################
###########################################################################
sub read_synteny{
    my($syn,$pair_h)=@_;
    open IN,"<$syn" or die"!";
    while(<IN>){
    chomp;
    my @t=split;
    $pair_h->{$t[0]}{$t[6]}=$t[-1];
#    print "$t[0]\t$t[6]\n";
    }
}
