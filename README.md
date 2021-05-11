# RFO (Reference-Free Ortholog)

This repository provides a pipeline that takes as input a HAL-format multiple whole genome alignment and GFF annotation files, and generates a reference-free ortholog table.

# Requirements
After installing these software, you need to add paths of these software into your `PATH` global variable.
1. [Kent toolkit](http://hgdownload.soe.ucsc.edu/admin/exe/).
2. [bedtools](http://bedtools.readthedocs.io/en/latest/).
3. [perl 5](https://www.perl.org/)
4. [python 2](https://www.python.org)
5. [HAL toolkit](https://github.com/glennhickey/hal). To install the HAL toolkit, you must also have the [sonLib](https://github.com/benedictpaten/sonLib) repository in the same parent directory. Compile sonLib first, then compile hal. Once hal is compiled, you need to have the binaries on your path.
6. [MAFFT](https://mafft.cbrc.jp/alignment/software/)

# Main options
`--step`: Assign which step to run.
`--hal`: Input HAL alignment file. (REQUIRED).
`--outdir`: Output directory. Defaults to current dir.
# Input file table
A template for the input table is below.
```
Homo_sapiens Homo_sapiens.fa Homo_sapiens.gff
Pan_troglodytes Pan_troglodytes.fa Pan_troglodytes.gff
Gorilla_gorilla Gorilla_gorilla.fa Gorilla_gorilla.gff
Pongo_abelii Pongo_abelii.fa Pongo_abelii.gff
Nomascus_leucogenys Nomascus_leucogenys.fa Nomascus_leucogenys.gff
Macaca_mulatta Macaca_mulatta.fa Macaca_mulatta.gff
```
Your GFF file should be formatted as below.
```
chr12 protein_coding mRNA 82358566 82479021 . + . ID=ENSG00000127720;
chr12 protein_coding CDS 82358566 82358824 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82386803 82386967 . + 2 Parent=ENSG00000127720;
chr12 protein_coding CDS 82389816 82389922 . + 2 Parent=ENSG00000127720;
chr12 protein_coding CDS 82398795 82399394 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82402983 82403130 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82430893 82430987 . + 2 Parent=ENSG00000127720;
chr12 protein_coding CDS 82434695 82434724 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82438718 82438791 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82456727 82456820 . + 1 Parent=ENSG00000127720;
chr12 protein_coding CDS 82476644 82476718 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82477281 82477352 . + 0 Parent=ENSG00000127720;
chr12 protein_coding CDS 82478932 82479021 . + 0 Parent=ENSG00000127720;
chr10 protein_coding mRNA 45303439 45304422 . - . ID=ENSG00000256574;
chr10 protein_coding CDS 45303439 45304422 . - 0 Parent=ENSG00000256574;
```
# Running the pipeline

1.Data preparation
```
perl RFO_pipeline.pl --step datapre --hal example.hal --outdir output input.tsv
```

2.Extract orthologs

```
perl RFO_pipeline.pl --step ortholog --hal example.hal --outdir output input.tsv
cd ortholog
sh orthologStep1.sh
sh orthologStep2.sh
sh orthologStep3.sh
```


3.Get the Reference-free Ortholog table
```
perl RFO_pipeline.pl --step merge --hal example.hal --outdir output input.tsv
```

This pipeline will output the reference-free ortholog table to `merge/`.


# Citation

Reference-free whole-genome alignment promote identification of ortholog and lineage specific gene
