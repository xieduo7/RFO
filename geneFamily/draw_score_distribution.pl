#!/usr/bin/perl
use strict;

##draw the distribution figure for hcluster score

my $file = shift;

my $distribute_svg = "perl /ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_evolution/Tools/distribute_svg/distribute_svg.pl";
my $distribute_pre = "perl /ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_evolution/Animal/distribute_pre.pl";
my $svg2xxx = "perl /ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_evolution/Tools/svg_kit/svg2xxx.pl";


`more $file | awk '\$1!=\$2' > $file.score`;
`awk '{print \$3}' $file.score  | $distribute_pre --frequency --minborder 0 --maxborder 100 --binsize 1  --header rect --color blue --mark frequency  -Ystart 0 -Yend 10 -Ystep 2  -X \"Hcluster score\" -Y \"Percent in all\" > $file.score.lst`;
`$distribute_svg $file.score.lst $file.score.svg`;
`$svg2xxx $file.score.svg`;
`rm $file.score`;
