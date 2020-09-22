#!user/bin/perl -w
use strict;
use File::Basename qw(basename dirname);
if (@ARGV != 3)
{
        print "perl $0 <ortholog_files.list> <prefix for output: ortholog> <ref>\n";
        exit;
}

my $list = shift;
my $out  = shift;
my $ref = shift;
my %ortholog;
my @taxa;
push @taxa,$ref;
&read_files($list,\%ortholog,\@taxa); # sub 1
my @taxa2 = @taxa;
shift @taxa2;
my $i = 1;
my ($out1,$out2);
foreach my $ref_gene(sort keys %ortholog)
{
#        $out1 .= "$i\t$ref\t";
        $out1 .= "$ref\t$ref_gene\t";
 #       $out2 .= "$i\t1\t";
#        $out2 .= "$i\t1\t";
        my $num = 1;
#        @taxa2 = sort @taxa2;
        foreach (@taxa2)
        {
                my $abb_name=$_;
                if (exists $ortholog{$ref_gene}{$abb_name})
                {
 #                       $out2 .= "1\t";
                        $out1 .= "$ortholog{$ref_gene}{$abb_name}\t";
                        $num ++;
                }
                else
                {
  #                      $out2 .= "0\t";
                        $out1 .= "-\t"
                }
        }
        $out1 .= "\n";
  #      $out2 .= "$num\n";
  #      $i++;
}


open OUT, ">${out}_ortholog_group" or die;
print OUT "reference\t".join ("\t", @taxa)."\n";
print OUT $out1;
close OUT;

#open OUT, ">${out}_ortholog_group.stat" or die;
#print OUT "ID\t".join ("\t", @taxa)."\tGeneNum\n";
#print OUT $out2;
#close OUT;



sub read_files
{
        my $file = shift;
        my $hshp = shift;
        my $hsha = shift; 
        my $specie="";
#        push @$hsha,$ref;
        open IN,$file or die "Fail to open $file\n";
        while (my $ln = <IN>)
        {
                chomp $ln;
                my $taxa = basename($ln);
                #print $ln;
                my $tem=(split(/\./,$taxa))[0];
              #  print $tem;
                $specie=(split(/-/,$tem))[1];
                push @$hsha,$specie;
                open ORTH,$ln or die "Fail to open $ln\n";
                while (<ORTH>)
                {
                        next if(/^#/);
                        my ($id1,$id2) = (split)[0,1];
                #        print "$specie\n";
 #                       $id2 =~ /^([\w,-]+)_/;
                        $hshp->{$id1}{$specie} = $id2;
                }
                close ORTH;

        }
        close IN;
}

