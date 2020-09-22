#!/usr/bin/env  perl
#####################################
## Duo Xie                         ##
## xieduo1907@gmail.com            ##
############################################################################################################
## Identifying orthologs and lineage-specific genes using the reference-free whole-genome alignment       ##
## Duo Xie                                                                                                ## 
## Version:1.0 2019-09-05                                                                                 ##
############################################################################################################
=head1 NAME

RFO_pipeline.pl  --  the pipeline of identifying ortholog using reference-free whole-genome alignment

=head1 DESCRIPTION

This pipeline prepares data,identifies ortholog,and constructs gene family using annotation,HAL file 
and genomes.

=head1 Usage

    perl COP.pl [options] inputTable
    --step <str>         set the start step for data prepartion,ortholog prediction or gene family construction
    --hal  <hal>         the path of HAL alignment
    --outdir <str>       set the result directory, default="./"

=head1 Example

    perl COP.pl --step datapre  --hal 200m-v1.hal inputTable
    perl COP.pl --step ortholog --hal 200m-v1.hal inputTable
    perl COP.pl --step genefam  --hal 200m-v1.hal inputTable

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my($step,$hal,$outdir,$outgroup,$Help);
GetOptions(
    "step:s"=>\$step,
    "hal:s" =>\$hal,
    "outdir:s" =>\$outdir,
#    "outgroup:s"=>\$outgroup,
    "help!"=>\$Help
);
$outdir ||= ".";
die `pod2text $0` if (@ARGV == 0 || $Help || (not defined $step));


###########################
####### MAIN SCRIPT #######
###########################
my $table = shift;
my (%info,@order);
parse_input($table,\%info,\@order);
my (@temp,@combine);
push(@temp,[$_->[0],$_->[1]]) for combine(\@order, 2);
foreach my $temp_pair(@temp){
    my $name="$temp_pair->[0] $temp_pair->[1]";    
    $name =~ s/\s+/-/;
    push(@combine,$name);
}
############################################
##############data preparation##############
############################################
if ($step eq "datapre"){
    mkdir "$outdir/data",0755;
    mkdir "$outdir/data/cds",0755;
    mkdir "$outdir/data/pep",0755;
    mkdir "$outdir/data/gff",0755;
    mkdir "$outdir/data/genome",0755;
    open OUT,">$outdir/data/category" or die"!";
    foreach my $species(@order){
         print OUT "$species\t1\n";
        `ln -s $info{$species}{"annotation"} $outdir/data/gff/$species.gff`;
        `ln -s $info{$species}{"genome"} $outdir/data/genome/$species.fa`;
        `perl $Bin/dataPre/getGene.pl $outdir/data/gff/$species.gff $outdir/data/genome/$species.fa|perl $Bin/dataPre/simplify_pep_fasta.pl - >$outdir/data/cds/$species.cds`;
        `perl $Bin/dataPre/cds2aa_sim.pl -species $species $outdir/data/cds/$species.cds >$outdir/data/pep/$species.pep`;
    }
    close OUT;
}
############################################
##############ortholog prediction###########
############################################
if ($step eq "ortholog"){
    mkdir "$outdir/ortholog",0755;
    my($cmdStep1,$cmdStep2,$cmdStep3);
    foreach my $pair(@combine){
        my @species = split(/-/,$pair); 
        mkdir "$outdir/ortholog/$pair",0755;
        mkdir "$outdir/ortholog/$pair/genome_files",0755;
        mkdir "$outdir/ortholog/$pair/annotation",0755;
        `ln -s $outdir/data/genome/$species[0].fa $outdir/ortholog/$pair/genome_files/$species[0].fa`;
        `ln -s  $outdir/data/genome/$species[1].fa $outdir/ortholog/$pair/genome_files/$species[1].fa`;
        `ln -s  $outdir/data/gff/$species[0].gff $outdir/ortholog/$pair/annotation/$species[0].gff`;
        `ln -s  $outdir/data/gff/$species[1].gff $outdir/ortholog/$pair/annotation/$species[1].gff`;
        `ln -s $outdir/data/cds/$species[0].cds $outdir/ortholog/$pair/annotation/$species[0].fa`;
        `ln -s $outdir/data/cds/$species[1].cds $outdir/ortholog/$pair/annotation/$species[1].fa`;
        $cmdStep1 .= "cd $outdir/ortholog/$pair && python $Bin/orthlogPrediction/orthologPredict.py preprocess --hal $hal --refGenome $species[0] --targetGenomes $species[1] --extend 0 && python $Bin/orthlogPrediction/orthologPredict.py projection --hal $hal --refGenome $species[0] --targetGenomes $species[1] --extend 0\n";
        $cmdStep2 .= "sh $outdir/ortholog/$pair/ortholog/$pair.gff.intersect.ortholog.overlap.mafft.sh\n";
        $cmdStep3 .= "cd $outdir/ortholog/$pair && python $Bin/orthlogPrediction/orthologPredict.py combine --hal $hal --refGenome $species[0] --targetGenomes $species[1] --extend 0 && python $Bin/orthlogPrediction/orthologPredict.py final --hal $hal --refGenome $species[0] --targetGenomes $species[1] --extend 0\n";
    }
    open OUT1,">$outdir/ortholog/orthologStep1.sh" || die "fail creat orthologStep1.sh"; 
    print OUT1 $cmdStep1;
    close OUT1;
    open OUT2,">$outdir/ortholog/orthologStep2.sh" || die "fail creat orthologStep2.sh";
    print OUT2 $cmdStep2;
    close OUT2;
    open OUT3,">$outdir/ortholog/orthologStep3.sh" || die "fail creat orthologStep3.sh";
    print OUT3 $cmdStep3;
    close OUT3;
}
############################################
###merge ortholog table based on reference##
############################################
if($step eq "merge"){
    testmkdir("$outdir/merge");
    my %genes;
    foreach my $spe(@order){
        open GFF,"<$outdir/data/gff/${spe}.gff" or die"!";
        while(<GFF>){
            chomp;
            my @t=split;
            if($t[2] eq "mRNA"){
                if($t[-1]=~/ID=([^;]+);/){
                    $genes{$spe}{$1}=1;
                }
            }
        }
        close GFF;
    }
    foreach my $pairwise(@combine){
#        storeGene("$outdir/ortholog/${pairwise}/ortholog/${pairwise}.gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab",\%genes);
    }
    foreach my $species(@order){
#        my @subfiles = glob("$outdir/ortholog/*/ortholog/${species}-*.gff.intersect.ortholog.overlap.tab.best");
        my @subfiles = glob("$outdir/ortholog/*/ortholog/${species}-*.gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.best");
        my (%ortholog,@taxa);
        &readOrtholog(\@subfiles,\%ortholog,\@taxa,\%genes,$species);
        my $out;
#        foreach my $ref_gene(sort keys %ortholog){
        foreach my $ref_gene(sort keys %{$genes{$species}}){
            #print "$ref_gene\n";
            $ref_gene="$species-$ref_gene";
            $out .= "$species\t$ref_gene\t";
            foreach (@taxa){
                my $abb_name=$_;
                if (exists $ortholog{$ref_gene}{$abb_name}){
                    $out .= "$ortholog{$ref_gene}{$abb_name}\t";
                }
                else{
                    $out .= "-\t";
                }
            }
            $out .= "\n";
        }
    open OUT, ">$outdir/merge/${species}_ortholog_group" or die"!";
    print OUT "reference\t${species}\t".join ("\t", @taxa)."\n";
    print OUT $out;
    close OUT;
    }
    my @orthologTab = glob("$outdir/merge/*_ortholog_group");
    my @container;
    my @intial_table;
    my $flag=0;
    my $out_dir="$outdir/merge";
    foreach my $ref(@order){
        if($flag==0){
            push(@intial_table,"$outdir/merge/${ref}_ortholog_group");
            push(@container,$ref);
            $flag=1;
        }
        else{
            mergeOrthologTable(\@intial_table,$ref,"$outdir/merge/${ref}_ortholog_group",$out_dir,\@container);
        }
    }
    my $fileName = join("-",@order);
    fillIn("$outdir/merge/${fileName}.ortholog.table",\%genes,"$outdir/merge/merged.ortholog.table");
#    `cp $outdir/merge/${fileName}.ortholog.table $outdir/merge/merged.ortholog.table`;
}


############################################
###homologous goup construction           ##
############################################
if($step eq "homologGroup"){
    testmkdir("$outdir/homologGroup");
    `cat $outdir/ortholog/*/ortholog/*.gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab|cut -f 1,2,4|sort|uniq >$outdir/homologGroup/homolog.forHC`;
    `$Bin/geneFamily/hcluster_sg -w 0 -c $outdir/homologGroup/homolog.forHC >$outdir/homologGroup/homologGroup.cluster`;
}


############################################
###gene family construction               ##
############################################
if($step eq "genefamily"){
    testmkdir("$outdir/genefamily/input");
    testmkdir("$outdir/genefamily/output");
    my @pepfiles = glob("$outdir/data/pep/*.pep");
    my $files = join(" ",@pepfiles);
    `cat $files >$outdir/genefamily/input/all.pep`;
    my @allgene;
    getAll("$outdir/genefamily/input/all.pep",\@allgene);
    `ln -s $outdir/genefamily/input/all.pep $outdir/genefamily/output/all.pep`;
    `formatdb -i $outdir/genefamily/input/all.pep -p T -o T`;
    `perl $Bin/geneFamily/fastaDeal.pl -cutf 50 $outdir/genefamily/output/all.pep -outdir $outdir/genefamily/output`;
    my @subfiles = glob("$outdir/genefamily/output/all.pep.cut/all.pep*");
    my $blast_shell = "$outdir/genefamily/geneFamilyStep1.sh";
    open OUT1,">$blast_shell" || die "fail open $blast_shell\n";
    foreach my $subfile (@subfiles){
        print OUT1 "blastall -p blastp -m8 -e 1e-7 -F F -d $outdir/genefamily/input/all.pep -i $subfile -o $subfile.blast.m8 \n";
    }
    close OUT1;
    my $combine = "$outdir/genefamily/geneFamilyStep2.sh";
    open OUT2,">$combine" || die "fail open $combine\n";
    foreach my $subfile (@subfiles){
        print OUT2"cat $subfile.blast.m8 >> $outdir/genefamily/output/all_vs_all.blast.m8\n";   
    }
    close OUT2;
    create_cate_file("$outdir/data/category", \@allgene, "$outdir/data/category.gene");
    my $cluster = "$outdir/genefamily/geneFamilyStep3.sh";
    open OUT,">$cluster" || die "fail open $cluster\n";
    print OUT "perl $Bin/geneFamily/solar.pl -a prot2prot -f m8 -z $outdir/genefamily/output/all_vs_all.blast.m8 > $outdir/genefamily/output/all_vs_all.blast.m8.solar.raw;\n";
    print OUT "perl $Bin/geneFamily/filter_solar_paired_redundance.pl $outdir/genefamily/output/all_vs_all.blast.m8.solar.raw >$outdir/genefamily/output/all_vs_all.blast.m8.solar.noPaired;\n";
    print OUT "perl $Bin/geneFamily/filter_solar_align_rate.pl $outdir/genefamily/input/all.pep $outdir/genefamily/output/all_vs_all.blast.m8.solar.noPaired 0.33 > $outdir/genefamily/output/all_vs_all.blast.m8.solar;\n";
    print OUT "perl $Bin/geneFamily/bitScore_to_hclusterScore.pl $outdir/genefamily/output/all_vs_all.blast.m8.solar >$outdir/genefamily/output/all_vs_all.blast.m8.solar.forHC 2>$outdir/genefamily/output/all_vs_all.blast.m8.solar.forHC.warn;\n";
    print OUT "$Bin/geneFamily/hcluster_sg -w 10 -s 0.34 -m 500 -b 0.1 -C $outdir/data/category.gene $outdir/genefamily/output/all_vs_all.blast.m8.solar.forHC >$outdir/genefamily/output/all_vs_all.blast.m8.solar.forHC.hcluster\n";
    close OUT;
}

############################################
###merge gene family and ortholog         ##
############################################
if($step eq "final"){
    testmkdir("$outdir/final");
    my %geneFamily;
    hclster2Table("$outdir/genefamily/output/all_vs_all.blast.m8.solar.forHC.hcluster",\%geneFamily);
#    hclster2Table("$outdir/homologousGroup/output/all_vs_all.blast.m8.solar.forHC.hcluster",\%geneFamily);
    my %homologGroup;
    hclster2Table("$outdir/homologGroup/homologGroup.cluster",\%homologGroup);
    orthologFamily("$outdir/merge/merged.ortholog.table",\%homologGroup,"$outdir/final/final_homologGroup_ortholog.tsv");
    orthologFamily("$outdir/merge/merged.ortholog.table",\%geneFamily,"$outdir/final/final_genefamily_ortholog.tsv");

}

#########################
####### FUNCTIONS #######
#########################
###########################################################################
sub orthologFamily{
    my $ortholog = shift;
    my $family = shift;
    my $out = shift;
    open IN,"<$ortholog" or die"!";
    open OUT,">$out" or die"!";
    my $header = <IN>;
    print OUT "familyID\t$header";
    while(<IN>){
        chomp;
        my $line=$_;
        my @t = split;
        shift @t;
        my %judge;
        foreach my $gene(@t){
            if(exists $family->{$gene}){
                $judge{$family->{$gene}}=1;
            }
        }
        my @family = keys%judge;
        if(@family == 1){
            print OUT "$family[0]\t$line\n";
        }
        else{
            print OUT "single\t$line\n";
        }
    }
    close OUT;
    close IN;
}
sub hclster2Table{
    my $input=shift;
    my $hash=shift;
    open IN,"<$input" or die"!";
    while(<IN>){
        chomp;
        my @t = split;
        my @genes = split(/,/,$t[-1]);
        foreach my $gene(@genes){
            $hash->{$gene} = $t[0];
        }
    }
    close IN;
}
sub create_cate_file {
        my ($cate_file, $a_geneids, $out_file) = @_;

        my %Category;
        open IN,"$cate_file" || die $!;
        while (<IN>) {
                chomp;
                my ($spec,$cate) = (split /\s+/)[0,1];
                $Category{$spec} = $cate;
        }
        close IN;

        open OUT,">$out_file" || die $!;
        foreach my $gene (@$a_geneids) {
                my $spec = $1 if($gene =~ /^([^-]+)-/);
                print OUT "$gene\t$Category{$spec}\n";
        }
        close OUT;
}

sub getAll{
   my ($uniprot_file, $a_geneids) = @_; 
   open IN, "$uniprot_file" or die "Can't open $uniprot_file\n";
   $/=">"; <IN>; $/="\n";
   while (<IN>){
    chomp;
    my $head = $1 if(/^(\S+)/);
    $/=">";
    <IN>;
    $/="\n";
    push @$a_geneids, $head;
   }
   close IN;
}
sub parse_input{
    my $file = shift;
    eval{ checkIfExists($file) };
    die if($@);
    my $info = shift;
    my $rank = shift;
    open IN,"<$file" or die"!";
#    print "$file\n";
    while(<IN>){
        chomp;
        my @arr=split;
        eval{ checkIfExists($arr[0]) };
        eval{ checkIfExists($arr[1]) };
        $info->{$arr[0]}{"genome"} = $arr[1];
        $info->{$arr[0]}{"annotation"} = $arr[2];
        push(@$rank,$arr[0]);
    }
    close IN;

}

sub checkIfExists{
        my $file = shift @_;
        if(! -f $file){
                die "File $file does not exist";
        }
}


###############################################################################################################
#combine;modified from answer of Borodin
#https://stackoverflow.com/questions/10299961/in-perl-how-can-i-generate-all-possible-combinations-of-a-list
##############################################################################################################

sub combine {

  my ($list, $n) = @_;
  die "Insufficient list members" if $n > @$list;

  return map [$_], @$list if $n <= 1;

  my @comb;

  for (my $i = 0; $i+$n <= @$list; ++$i) {
    my $val  = $list->[$i];
    my @rest = @$list[$i+1..$#$list];
    push @comb, [$val, @$_] for combine(\@rest, $n-1);
  }

  return @comb;
}
###########################################################################################
sub readOrtholog{
        my $list=shift;       
        my $hshp = shift;
        my $hsha = shift;
        my $addGene = shift;
        my $ref = shift;
        my $specie="";
        foreach my $ln(@$list){
                my $taxa = basename($ln);
                my $tem=(split(/\./,$taxa))[0];
                $specie=(split(/-/,$tem))[1];
                push @$hsha,$specie;
                open ORTH,$ln or die "Fail to open $ln\n";
                while (<ORTH>)
                {
                        next if(/^#/);
                        my ($id1,$id2) = (split)[0,1];
                        $hshp->{$id1}{$specie} = $id2;
                }
                close ORTH;
                #####
#                foreach my $gene(keys%{$addGene->{$ref}}){
#                    if($addGene->{$ref}{$gene}==1){
#                        $hshp->{"$ref-$gene"}{$specie} = "-";
#                    }
#                    else{
#                        next;
#                    }
#                }
                #####
        }
}
sub mergeOrthologTable{
        my ($left,$taxon2,$ortholog_table2,$dest,$priority) = @_;
        testmkdir($dest);
        my $ortholog_table1 = $left->[0];
        eval{ checkIfExists($ortholog_table1) };
        die if($@);
        eval{ checkIfExists($ortholog_table2) };
        die if($@);
        my (@left,%ortholog1,%all_genes_left,%ortholog2,@order);
        fetchLeftTable($ortholog_table1,\%ortholog1,\%all_genes_left,\@order);

        fetchRightTable($ortholog_table2,\%ortholog2);
        foreach my $reference(keys %ortholog2){
            foreach my $gene(keys %{$ortholog2{$reference}}){
                my $switch=1;
                foreach my $rank(@$priority){
                    my @get1 = keys %{$ortholog1{$rank}};
                    foreach my $left_gene(@get1){
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
            if($switch==1){
               my @add=values %{$ortholog2{$reference}{$gene}};
                    my $flag=1;
                    foreach my $add_gene(@add){
                        if(exists $all_genes_left{$add_gene}){
                            $flag=0;
                            next;
                        }
                        else{
                            next;
                        }
                    }
                    if($flag==0){
                    next;
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
        }
}
    push(@$priority,$taxon2);
    my $file_name = join("-",@$priority);
    my $first_line = join("\t",@order);
    open OUT,"> $dest/$file_name.ortholog.table" or die "Cannot open the file!";
    print OUT "reference\t$first_line\n";
        foreach my $ref_print(keys %ortholog1){
            foreach my $gene_print(keys%{$ortholog1{$ref_print}}){
            print OUT "$ref_print\t";
                foreach my $species(@order){
                  print OUT"$ortholog1{$ref_print}{$gene_print}{$species}\t";
                }
                print OUT"\n";
            }
        }

    close OUT;
    $left->[0]="$dest/$file_name.ortholog.table";
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
                        $ortholog_hash->{$ref}{$t[int($column{$ref})]}{$column_rev{$index}} = $gene;
                        $index +=1;
            }

        }
        close TABLE;



}

sub fetchLeftTable{
        my($table,$ortholog_hash,$all_gene,$order) = @_;
        open TABLE,"<$table" or die"Cannot open the ortholog table!";
        my $first_line = <TABLE>;
        chomp($first_line);
        @$order = split(/\s/,$first_line);
        splice @$order,0,1;
        my $count = 0;
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
            foreach my $gene(@t){
                if($gene ne "-"){
                    $all_gene->{$gene} =1;
                    }
                        $ortholog_hash->{$ref}{$t[int($column{$ref})]}{$column_rev{$index}} = $gene;
                        $index +=1;
                }
            }
        close TABLE;
        }
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
sub storeGene{
    my $file = shift;
    my $hash = shift;
#    print "$file\n";
    open IN,"<$file" or die"!";
    while(<IN>){
        chomp;
        my @arr=split;
        my @spe1=split(/-/,$arr[0]);
        my @spe2=split(/-/,$arr[1]);
        $hash ->{$spe1[0]}{$spe1[1]}=0;
        $hash ->{$spe2[0]}{$spe2[1]}=0;
    }
    close IN;
}

sub fillIn{
    my($table,$gff_h,$out)=@_;
    open TABLE,"<$table" or die"!";
    open NT,">$out" or die"!";
    open LIST,">$out.id" or die"!";
    my $first_line = <TABLE>;
    chomp($first_line);
    print NT "$first_line\n";
    my @order = split(/\s/,$first_line);
    splice @order,0,1;
    my $count = 0;
    my %column;
    my %column_rev;
    foreach my $species(@order){
        $column{$species} = $count;
        $column_rev{$count} = $species;
        $count +=1;
        }
    while(<TABLE>){
        chomp;
        print NT "$_\n";
        my @t =split;
        my $index=0;
        my $ref = shift @t;
        foreach my $gene(@t){
            my $id=(split(/-/,$gene))[1];
            #print "$id\t$column_rev{$index}\n";
            if($index==$column{$ref}){
                print LIST "$id\n";           
            }
            $gff_h->{$column_rev{$index}}{$id}=0;
            $index +=1;
        }
    }
    close TABLE;
    foreach my $spe(keys%$gff_h){
        foreach my $genes(keys %{$gff_h->{$spe}}){
            next if($gff_h->{$spe}{$genes}==0);
            print NT "$spe";
            for(my $i=0;$i<$#order+1;$i++){
                if($i==$column{$spe}){
                    print NT "\t$spe-$genes";
                    print LIST "$genes\n";
                }else{
                    print NT "\t-";
                }
            }
            print NT "\n";

        }
    }
    close NT;
    close LIST;

}
