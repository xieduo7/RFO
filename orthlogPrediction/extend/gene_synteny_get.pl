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

gene_synteny_get.pl  --  classfiy the orthologs into different categories

=head1 Usage

    perl gene_synteny_get.pl inputTable

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
die `pod2text $0` if (@ARGV == 0 );
my $table=shift;
my (%left,%right);
count_input($table,\%left,\%right);
open INN,"<$table" or die"!";
while(<INN>){
    chomp;
    my @arr=split;
    if($arr[-1] eq "ortholog_one2one"){
        print "$_\tortholog_one2one\n";
    }elsif($arr[-1] ne "NO_gene_synteny"){
        if($left{$arr[0]}==1 && $right{$arr[1]}==1){
            print "$_\tanc_copy\n";
        }else{
            print "$_\tduplication\n";
        }
    }else{
            print "$_\tduplication\n";
    }
}
close INN;



#########################
####### FUNCTIONS #######
#########################
###########################################################################
sub count_input{
    my ($ortholog,$left_h,$right_h)=@_;
    open IN,"<$ortholog" or die"!";
    while(<IN>){
        chomp;
        my @t=split;
        if($t[-1] ne "NO_gene_synteny"){
        $left_h->{$t[0]}+=1;
        $right_h->{$t[1]}+=1;
        }

    }
    close IN;
}

