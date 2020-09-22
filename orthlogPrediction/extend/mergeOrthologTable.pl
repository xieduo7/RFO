#!/usr/bin/env  perl
use strict;
use warnings;
use Data::Dumper;
#####################################
## Duo Xie                         ##
## xieduo1907     AT gmail DOT com ##
############################################################################################################
## Predicting orthologs using the reference-free whole-genome alignment                                   ##
## Duo Xie                                                                                                ##
############################################################################################################





###########################
####### MAIN SCRIPT #######
###########################
### command line arguments
# 1: infile
my $ortholog_list = "";
if(scalar @ARGV == 0){
        die "Please provide a list file\n";
} else {
        $ortholog_list = shift @ARGV;
}

# 2: outfile
my $outdir = "";
if(scalar @ARGV == 0){
        die "Please provide a output directory\n";
} else {
        $outdir = shift @ARGV;
}

#3: read ortholog file list
my %hash = readOrthologList($ortholog_list);

#4: merge the ortholog table one by one
my @container;
my $count = 0;
my @intial_table;
foreach my $key(sort{$a <=> $b} keys %hash){
        if($key == 1){
           push(@intial_table,$hash{$key}[1]);
           push(@container,$hash{$key}[0])
        }
        else{
           mergeOrthologTable(\@intial_table,$hash{$key}[0],$hash{$key}[1],$outdir,\@container); 
        }
}



#########################
####### FUNCTIONS #######
#########################
###########################################################################
###########################################################################
sub readOrthologList{
        my $file = shift @_;
        eval{ checkIfExists($file) };
        die if($@);
        open IN,"<$file";
        my %result;
        while(<IN>){
            chomp;
            my @arr = split;
            $result{$arr[0]} = [$arr[1],$arr[2]];
        }
        close IN;
        return %result;
}
###########################################################################
sub mergeOrthologTable{
        my ($left,$taxon2,$ortholog_table2,$dest,$priority) = @_;
        testmkdir($dest);
      #  my $taxon1 = $left->[0];
        my $ortholog_table1 = $left->[0];
        eval{ checkIfExists($ortholog_table1) };
        die if($@);
        eval{ checkIfExists($ortholog_table2) };
        die if($@);
        my (@left,%ortholog1,%all_genes_left,%ortholog2,@order);
        fetchLeftTable($ortholog_table1,\%ortholog1,\%all_genes_left,\@order);

        fetchRightTable($ortholog_table2,\%ortholog2);
#=cut
        foreach my $reference(keys %ortholog2){
#            print Dumper( keys %{$ortholog2{$reference}});

            foreach my $gene(keys %{$ortholog2{$reference}}){
            #    print "$gene\n";
                my $switch=1;
                foreach my $rank(@$priority){
                    my @get1 = keys %{$ortholog1{$rank}};
                    #################
                    foreach my $left_gene(@get1){
#                    print "$left_gene\n";
                    if($ortholog2{$reference}{$gene}{$rank} eq $left_gene && $left_gene ne "-"){
                        foreach my $species_check(sort { $ortholog1{$rank}{$left_gene}{$b} cmp $ortholog1{$rank}{$left_gene}{$a} } keys %{$ortholog1{$rank}{$left_gene}}){

                            if($ortholog1{$rank}{$left_gene}{$species_check} eq $ortholog2{$reference}{$gene}{$species_check}){
                                next;
                        } 
                        else{
                            if($ortholog1{$rank}{$left_gene}{$species_check} eq "-"){
                                if(exists $all_genes_left{$ortholog2{$reference}{$gene}{$species_check}}){
                                    $switch=0;
                                    last;
                                }
                                else{
                     $ortholog1{$rank}{$left_gene}{$species_check}=$ortholog2{$reference}{$gene}{$species_check};                                                     $switch=0;
                                }
                            }
                            else{
                                $switch=0;
                                last;   
                            }
                        }
                        
                    }
                 }
                 }
                 #############################################
                
            }
#            print "$switch\n";
            if($switch==1){
#              print "OK\n"; 
               my @add=values %{$ortholog2{$reference}{$gene}};
                    my $flag=1;
                    foreach my $add_gene(@add){
                        if(exists $all_genes_left{$add_gene}){
                            $flag=0;
                            next;
         #                   last;
                        }
                        else{
                            next;
                        }
                    }
                    if($flag==0){
                    next;
        #                last;
                    }
                    else{
                        foreach my $key_add(keys %{$ortholog2{$reference}{$gene}}){
                            $ortholog1{$reference}{$gene}{$key_add}=$ortholog2{$reference}{$gene}{$key_add};
                            }
                    }
 
            }
                    else{
                        next;
                    }
            ####################################
        }
}
    push(@$priority,$taxon2);
    my $file_name = join("-",@$priority);
    my $first_line = join("\t",@order);
#    print "$dest/$file_name.ortholog.table\n";
#    my $tem_dest="$dest/$file_name.ortholog.table";
    open OUT,"> $dest/$file_name.ortholog.table" or die "Cannot open the file!";
    #open OUT, "> $tem_dest" or die $!;
    print OUT "reference\t$first_line\n";
        foreach my $ref_print(keys %ortholog1){
            foreach my $gene_print(keys%{$ortholog1{$ref_print}}){
            print OUT "$ref_print\t";
                foreach my $species(@order){
    #              print OUT"$species\t";
                  print OUT"$ortholog1{$ref_print}{$gene_print}{$species}\t";
                }
                print OUT"\n";
            }
        }

    close OUT;
    $left->[0]="$dest/$file_name.ortholog.table";
#        push(@$priority,$taxon2);
#=cut
}
###########################################################################
sub fetchRightTable{
        my $table =shift;
        my $ortholog_hash = shift;
        open TABLE,"<$table" or die"Cannot open the ortholog table!";
        my $first_line = <TABLE>;
        chomp($first_line);
        my @arr = split(/\s/,$first_line);
        splice @arr,0,1;
        my $count = 0;
        my %column;
        my %column_rev;
        foreach my $species(@arr){
            $column{$species} = $count;
            $column_rev{$count} = $species;
            $count +=1;
            }
        while(<TABLE>){
            chomp;
            my @t =split;
            my $index=0;
            my $ref = shift @t;
            foreach my $gene(@t){
 #           if($index == $column{$ref}){
   #                     $index +=1;
  #                      next;
    #                }
     #               else{
                        $ortholog_hash->{$ref}{$t[int($column{$ref})]}{$column_rev{$index}} = $gene;
                        #print "$column_rev{$index}\n";
                      #  print "$gene\n";
                        $index +=1;
       #             }
            }

        }
        close TABLE;
            


}

###########################################################################
sub fetchLeftTable{
        my($table,$ortholog_hash,$all_gene,$order) = @_;
#        my (%ortholog);
        open TABLE,"<$table" or die"Cannot open the ortholog table!";
        my $first_line = <TABLE>;
        chomp($first_line);
        @$order = split(/\s/,$first_line);
        splice @$order,0,1;
        my $count = 0;
#        my $ref_index = 0;
        my %column;
        my %column_rev;
#####get the index of different species
        foreach my $species(@$order){
                $column{$species} = $count;
                $column_rev{$count} = $species;
                $count +=1;
            }
        while(<TABLE>){
            chomp;
            my @t =split;
            my $index=0;
            my $ref = shift @t;
        #    print("@t\n");
            foreach my $gene(@t){
                if($gene ne "-"){
                    $all_gene->{$gene} =1;
                    }
               
#                    $all_genes{$gene} = 1;
               #     if($index == $column{$ref}){
                #        $index +=1;
                 #       next;
              #      }
            #        else{
                        $ortholog_hash->{$ref}{$t[int($column{$ref})]}{$column_rev{$index}} = $gene;           
                        $index +=1;
             #       }
                }
            }
        close TABLE;
#        return (%ortholog);
        }



###########################################################################
sub checkIfExists{
        my $file = shift @_;
        if(! -f $file){
                die "File $file does not exist";
        }
}


##########################################################################
sub testmkdir{
        my $dir=shift;
        if (-e $dir){
#20170830               warn "Warning: Folder ($dir) exists! all files in it will be deleted!\n";
#20170830               `rm -r $dir`;
        }
        else{
        `mkdir -p $dir`;
        }

}

