import os
import re
import subprocess
def testMkdir(path):
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)

def testShell(shell):
    if os.path.isfile(shell):
        os.remove(shell)
    else:
        pass
def getOrthologPre(intersect):
    orth = open(intersect+".ortholog","w")
    for l in open(intersect,"r"):
        l = l.rstrip("\n")
        arr=l.split("\t")
        refer = arr[17].rstrip(";").split("=")[1]
        target = arr[8].rstrip(";").split("=")[1].split(";")[0]
        orth.write(refer+"\t"+target+"\t"+arr[3]+"\t"+arr[4]+"\t"+arr[12]+"\t"+arr[13]+"\t"+arr[18]+"\n") 
    orth.close()

def orthologToBed(ortholog):
    f = open(ortholog+".bed","w")
    for l in open(ortholog):
        l = l.rstrip("\n")
        arr = l.split("\t")
        start = int(arr[4])-1
        f.write(arr[0]+"-"+arr[1]+"\t"+str(start)+"\t"+arr[5]+"\n")
    f.close()
    sortedBed = open(ortholog+".sorted.bed","w")
    sortBed = subprocess.Popen(['sortBed','-i',ortholog+".bed"],stdout = sortedBed)
    sortBed.wait()
    sortedBed.close()

def calculateOverlap(bed,ortholog):
    mergedBed = open(bed+".sorted.merge.bed","w")
    mergeBed = subprocess.Popen(['bedtools','merge','-i',bed],stdout = mergedBed)
    mergeBed.wait()
    mergedBed.close()
    overlapLen = {}
    for l in open(bed+".sorted.merge.bed"):
        arr = l.split("\t")
        geneA = arr[0].split("-")[0]
        geneB = arr[0].split("-")[1]
        start = int(arr[1])
        end = int(arr[2])
        length = end - start
        #print end
    #    print start
        #print length
#        length = int(arr[2]) - int(arr[1])
        if geneA in overlapLen.keys():
            if geneB in overlapLen[geneA].keys():
                overlapLen[geneA][geneB] += length
            else:
                overlapLen[geneA][geneB] = length
        else:
            overlapLen[geneA]={}
            overlapLen[geneA][geneB] = length
    orth = open(bed+".sorted.merge.bed.overlap","w")
    for i in overlapLen.keys():
        for j in overlapLen[i].keys():
            orth.write(i+"\t"+j+"\t"+str(overlapLen[i][j])+"\n")
    orth.close()

                
def getOrtholog(intersect):
    f = open(intersect,"r")
    lines = f.readlines()
    lines2 = {}.fromkeys(lines).keys()
    ortholog = {}
    for l in lines2:
        l=l.rstrip("\n")
        arr=l.split("\t")
        refer = arr[0]
#        target = arr[8].rstrip(";").split("=")[1]
        target = arr[1]
        length = int(arr[-1])
        if refer in ortholog.keys():
            if target in ortholog[refer].keys():
                ortholog[refer][target] += length
            else:
                ortholog[refer][target] = length
        else:
            ortholog[refer] = {}
            ortholog[refer][target] = length
    f.close()
    orth = open(intersect+".overlap","w")
    for i in ortholog.keys():
        for j in ortholog[i].keys():
            orth.write(i+"\t"+j+"\t"+str(ortholog[i][j])+"\n")
    orth.close()
