#!user/bin/perl -w
use strict;

my $class=shift;
my $dist1=shift;
my $dist2=shift;
my $rbh=shift;
my(%overlap,%best);
&readdist($dist1,\%overlap);
&readdist($dist2,\%overlap);
&readdist_1($rbh,\%best);
my $ref;
my $query;
open IN,"<$class" or die"!";
while(<IN>){
    chomp;
    my @arr = split;
    $ref = (split("-",$arr[0]))[1];
    $query = (split("-",$arr[1]))[1];
    if((exists $overlap{$ref} || exists $overlap{$query}) && $arr[-1] ne "ortholog_one2one"){
        if(exists $best{$ref}){
            if($query eq $best{$ref}){
    #            $arr[-1]="ortholog_one2one";
#                print "$_i\n";
                print join("\t", @arr), "\n";
            }
        }
    }else{
        print "$_\n";
    }
}
close IN;

sub readdist{
    my $file=shift;
    my $hash=shift;
    open IN,"<$file" or die"!";
    while(<IN>){
        chomp;
        my @t=split;
        $hash->{$t[0]}=$t[1];
    }
    close IN;
}
sub readdist_1{
    my $file=shift;
    my $hash=shift;
    open IN,"<$file" or die"!";
    while(<IN>){
        next if(/^#/);
        chomp;
        my @t=split;
        my $ref=(split("-",$t[0]))[1];
        my $query=(split("-",$t[1]))[1];
        $hash->{$ref}=$query;
    }
    close IN;
}

