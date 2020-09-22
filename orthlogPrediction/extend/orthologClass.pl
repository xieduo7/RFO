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

orthologClass.pl  --  classfiy the orthologs into different categories

=head1 Usage

    perl COP.pl inputTable

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
die `pod2text $0` if (@ARGV == 0 );
###########################
####### MAIN SCRIPT #######
###########################
my $table=shift;
my (%left,%right);
count_input($table,\%left,\%right);
open INN,"<$table" or die"!";
while(<INN>){
    chomp;
    my @arr=split;
    if($left{$arr[0]}==1 && $right{$arr[1]}==1){
        print "$_\tortholog_one2one\n";
    }elsif(($left{$arr[0]}==1 && $right{$arr[1]}>1)){
        print "$_\tortholog_many2one\n";
    }elsif($left{$arr[0]}>1 && $right{$arr[1]}==1){
        print "$_\tortholog_one2many\n";
    }elsif($left{$arr[0]}>1 && $right{$arr[1]}>1){
        print "$_\tortholog_many2many\n";
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
        $left_h->{$t[0]}+=1;
        $right_h->{$t[1]}+=1;

    }
    close IN;
}

