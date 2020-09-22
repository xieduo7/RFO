#!/user/bin/env python
#! coding:utf-8
import argparse
import tools.dataPrepare 
import subprocess
import os
import shutil
import tools.basic
import tools.filterData
###
#Parse Parameters
###

parser = argparse.ArgumentParser(
    description='A pipeline to predict ortholog from cactus alignment',
    epilog='DuoXie(xieduo@genomics.cn)',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=True)
#common_parser = argparse.ArgumentParser(add_help=False)
Input = parser.add_argument_group('Input file options')
#Input.add_argument("--refGff",help="The annotation gff of reference genome.Required",required=True)
parser.add_argument(
    "option",
    help="choose a process to run,"
    + "preprocess/filter"
    )
Input.add_argument(
    "--hal",
    help="HAL alignment file. Required",
    required=True)

Input.add_argument(
    "--refGenome",
    help="Reference genome sequence name."
    + "Must be present in HAL. Required",
    required=True
    )
Input.add_argument(
    "--targetGenomes",
    help="list of genomes to use. Required",
    required=True
    )
Input.add_argument(
    "--extend",
    help="the length of extending in extract alignment. Required",
    required=True
    )
#Input.add_argument(
#    "option",
#    help="choose a process to run,"
#    + "preprocess/filter"
#    )
args = parser.parse_args()

###
cwd = os.getcwd()
genomeFiles = cwd+"/genome_files/"
annotation = cwd+"/annotation/"
if args.option == "preprocess":
###
#Data Preparation 
###
#    nameConvert = {'Chimp': 'Clint_Chimp', 'Gorilla': 'Susie_Gorilla', 'Orangutan': 'Susie_Orangutan'}
    tools.dataPrepare.faSize(genomeFiles)
#    tools.dataPrepare.faSize(annotation)
    tools.dataPrepare.faTwoBit(genomeFiles)
 #   tools.dataPrepare.gffToPsl(annotation+args.refGenome+".gff",
  #                              annotation+args.refGenome+".size",
  #                              genomeFiles+args.refGenome+".size")
  #  tools.dataPrepare.gffRemoveRNA(annotation+args.targetGenomes+".gff",
  #                          annotation+args.targetGenomes+".no-mRNA.gff")
    tools.dataPrepare.sizeToBed(genomeFiles+args.refGenome+".size",args.refGenome)
    refName = args.refGenome
    targetName = args.targetGenomes
    #if args.targetGenomes in nameConvert:
  #      targetName = nameConvert[args.targetGenomes]
 #   else:
  #      targetName = args.targetGenomes
  #  if args.refGenome in nameConvert:
  #      refName = nameConvert[args.refGenome]
  #  else:
  #      refName = args.refGenome
###
#chaining
###

    tools.basic.testMkdir("chaining")
    halliftoverCmd = ['halLiftover','--outPSL',
                       args.hal,refName,genomeFiles+args.refGenome+".bed",
                       targetName,"chaining/"+args.refGenome
                       +"-"+args.targetGenomes+".psl"]
    halliftover = subprocess.Popen(halliftoverCmd)
    halliftover.wait()
    pslPosTarget = subprocess.Popen(['pslPosTarget',"chaining/"
                                    +args.refGenome+"-"
                                    +args.targetGenomes
                                    +".psl","chaining/"
                                    +args.refGenome+"-"+args.targetGenomes+".pos.psl"])
    pslPosTarget.wait()    
    axtChain = subprocess.Popen(
                            ['axtChain', 
                            '-psl', 
                            '-verbose=0', 
                            '-linearGap=medium',
                            cwd +"/chaining/"
                            + args.refGenome
                            + "-"+args.targetGenomes
                            + ".pos.psl",
                            genomeFiles+args.targetGenomes
                            +".fa.2bit",
                            genomeFiles+args.refGenome
                            +".fa.2bit",
                            cwd+"/chaining/"+args.refGenome
                            +"-"+args.targetGenomes+".pos.psl.chain"]
                            )
    axtChain.wait()
if args.option == "projection":
    tools.dataPrepare.targetGffExtend(annotation+args.refGenome+".gff",
                                    genomeFiles+args.refGenome+".size",args.extend,
                                    annotation+args.refGenome+".extend.gff")
    tools.dataPrepare.gffToPslBed(annotation+args.refGenome+".extend.gff",
                        genomeFiles+args.refGenome+".size")
    tools.dataPrepare.targetGffExtend(annotation+args.targetGenomes+".gff",
                            genomeFiles+args.targetGenomes+".size",
                            args.extend,
                            annotation+args.targetGenomes+".extend.gff")
    tools.dataPrepare.gffToSize(annotation+args.refGenome+".extend.gff")
    tools.dataPrepare.gffToSize(annotation+args.targetGenomes+".extend.gff")
    pslMap = subprocess.Popen(
                            ['pslMap',
                            '-chainMapFile',
                            annotation+args.refGenome+".extend.gff.psl",
                            cwd+"/chaining/"+args.refGenome
                            +"-"+args.targetGenomes+".pos.psl.chain",
                            annotation+args.targetGenomes+".gff.psl"])
    pslMap.wait()
    pslMapPostChain = subprocess.Popen(
                            ['pslMapPostChain',
                            annotation+args.targetGenomes
                            +".gff.psl",annotation
                            +args.targetGenomes
                            +".gff.post.psl"])
    pslMapPostChain.wait()
    sortPsl = open(annotation+args.targetGenomes+".gff.post.sort.psl","w")
    sort = subprocess.Popen(
                            ['sort','-k14,14',
                            '-k16,16n',
                            annotation+args.targetGenomes
                            +".gff.post.psl"],
                            stdout = sortPsl)
    sort.wait()
    sortPsl.close()
    pslGff = open(annotation+args.targetGenomes+".gff.post.psl.sort.gff","w")
    pslToGff = subprocess.Popen(
                            ['perl',
                            os.path.split(os.path.realpath(__file__))[0]
                            + "/extend/blat2gff.pl",
                            "-version",
                            "2",
                            annotation+args.targetGenomes
                            + ".gff.post.sort.psl"],
                            stdout = pslGff)
    pslToGff.wait()
    pslGff.close()
    tools.basic.testMkdir("ortholog")
    interect = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect","w")
    gffIntersect = subprocess.Popen(
                            ['bedtools',
                            'intersect',
                            '-wo',
                            '-s',
                            '-a',
                            annotation
                            + args.targetGenomes+".extend.gff",
                            '-b',
                            annotation+args.targetGenomes
                            +".gff.post.psl.sort.gff"],
                            stdout = interect)
    gffIntersect.wait()
    interect.close()
    tools.basic.getOrthologPre(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect")
    tools.basic.getOrtholog(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog")


###
#make shell
###
    tools.basic.testShell(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.mafft.sh")
    makeShell = subprocess.Popen(
                            ["perl",
                            os.path.split(os.path.realpath(__file__))[0]
                            + "/extend/run_mafft_ortholog.pl",
                            cwd+"/ortholog/"+args.refGenome
                            + "-"+args.targetGenomes
                            + ".gff.intersect.ortholog.overlap",
                            annotation+args.refGenome
                            +".fa",
                            annotation+args.targetGenomes
                            +".fa",
                            cwd+"/ortholog/"])
    makeShell.wait()

if args.option == "combine":
    tools.filterData.combine(
                        cwd
                        + "/ortholog/"
                        + args.refGenome
                        + "-"
                        + args.targetGenomes
                        + ".gff.intersect.ortholog.overlap.cut"

    )
    tools.basic.orthologToBed(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog")
    tools.basic.calculateOverlap(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed",cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog")
    tools.filterData.rbhPrepare(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap",
                                cwd+ "/ortholog/"+ args.refGenome+ "-"+ args.targetGenomes+ ".gff.intersect.ortholog.overlap.cut.id.ident",
                                annotation+args.refGenome+".extend.gff.size",
                                annotation+args.targetGenomes+".extend.gff.size",
                                args.refGenome,
                               args.targetGenomes )
#    getRbh = subprocess.Popen(
#            ["perl",
#            os.path.split(os.path.realpath(__file__))[0]+"/extend/getRbh_v2.pl",
#            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab",
#            args.refGenome,
#            args.targetGenomes,
#            cwd+"/ortholog/"])
#    getRbh.wait()
if args.option == "final":
############class ortholog into 1:1 1:n n:1 and n:n
    classFile = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class","w")
    orthologClass = subprocess.Popen(
        ["perl",
           os.path.split(os.path.realpath(__file__))[0]+"/extend/orthologClass.pl",
           cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab"],
        stdout = classFile)
    orthologClass.wait()
    classFile.close()
#######add gene synteny information
    geneOrder = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.order","w")
    addOrder = subprocess.Popen(
        ["perl",
            os.path.split(os.path.realpath(__file__))[0]+"/extend/geneOrder_add_orthologV2.pl",
            annotation+args.refGenome+".gff",
            annotation+args.targetGenomes+".gff",
            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab"],
            stdout = geneOrder)
    addOrder.wait()
    geneOrder.close()
########gene order ortholog
    geneOrt = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.order.genesynteny","w")
    getGenesyn = subprocess.Popen(
        ["perl",
            os.path.split(os.path.realpath(__file__))[0]+"/extend/syntenic_ort.pl",
            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.order",
            annotation+args.refGenome+".gff",
            annotation+args.targetGenomes+".gff"],
            stdout = geneOrt)
    getGenesyn.wait()
    geneOrt.close()
#########combine
    combine = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny","w")
    combineSyn = subprocess.Popen(
        ["perl",
            os.path.split(os.path.realpath(__file__))[0]+"/extend/gene_synteny_filter.pl",
            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class",
            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.order.genesynteny"],
            stdout = combine)
    combineSyn.wait()
    combine.close()
########1:1
    prim = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny.prime","w")
    getPrim = subprocess.Popen(
        ["perl",
            os.path.split(os.path.realpath(__file__))[0]+"/extend/gene_synteny_get.pl",
            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny"],
    stdout = prim)
    getPrim.wait()
    prim.close()
#########rbh
    rbh = open(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny.prime.rbh","w")
    getRbh = subprocess.Popen(
        ["perl",
            os.path.split(os.path.realpath(__file__))[0]+"/extend/getRbh_v10.pl",
            cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny.prime",
            args.refGenome,
            args.targetGenomes,
            cwd+"/ortholog/"],
            stdout = rbh)
    getRbh.wait()
    rbh.close()
#####combine and generate the final ortholog table
    tools.filterData.makeOrth(cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny.prime",cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.class.class.genesynteny.prime.rbh",cwd+"/ortholog/"+args.refGenome+"-"+args.targetGenomes+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.best",cwd+"/ortholog/"+args.targetGenomes+"-"+args.refGenome+".gff.intersect.ortholog.sorted.bed.sorted.merge.bed.overlap.tab.best")
    
