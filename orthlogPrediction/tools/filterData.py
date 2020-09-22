import os
import basic
import subprocess
import dataPrepare


def makeOrth(result,rbh,destfileA,destfileB):
    fa = open(destfileA,"w+")
    fb = open(destfileB,"w+")
    header = "#Referrence\tSubject\tR_ratio\tS_ratio\tIdentity\tSocre\tState\tLevel"
    fa.write(header+"\n")
    fb.write(header+"\n")
    for x in open(result):
        x = x.rstrip("\n")
        arr = x.split()    
        if(arr[-1] == "anc_copy" or arr[-1] == "ortholog_one2one"):
        ############g1  g2  align1  align2  ID  
            output1 = arr[0]+"\t"+ arr[1]+"\t"+arr[5]+"\t"+arr[6]+"\t"+arr[4]+"\t"+arr[3]+"\t"+arr[-1]+"\tL0"
            fa.write(output1+"\n")
            output2 = arr[1]+"\t"+ arr[0]+"\t"+arr[6]+"\t"+arr[5]+"\t"+arr[4]+"\t"+arr[3]+"\t"+arr[-1]+"\tL0"
            fb.write(output2+"\n")
    for y in open(rbh):
        y = y.rstrip("\n")
        art = y.split()
        output3 = art[0]+"\t"+ art[1]+"\t"+art[2]+"\t"+art[3]+"\t"+art[4]+"\t"+art[5]+"\t"+art[6]+"\t"+art[-1]
        fa.write(output3+"\n")
        output4 = art[1]+"\t"+ art[0]+"\t"+art[3]+"\t"+art[2]+"\t"+art[4]+"\t"+art[5]+"\t"+art[6]+"\t"+art[-1]
        fb.write(output4+"\n")
    fa.close()
    fb.close()

            
        


def combine(path):
    basic.testShell(path+".id.ident")
    f = open(path+".id.ident","w+")
    for folderName, subfolders, filenames in os.walk(path):
            for filename in filenames:
                if filename.find("mafft_ginsi.fa")!=-1:
                    ident = subprocess.Popen(
                                      ['perl',
                                      os.path.split(os.path.realpath(__file__))[0]
                                      + '/../extend/identity.pl',
                                        folderName + "/" + filename],
                                        stdout = f)
                    ident.wait()
    f.close()
#            if os.path.splitext(filename)[-1] == ".fa"
def rbhPrepare(overlap,identity,refSize,targetSize,reference,target):
    refSizes = dataPrepare.readSize(refSize)
    targetSizes = dataPrepare.readSize(targetSize)
    ident = {}
    overlapper = open(overlap+".tab","w")
    for x in open(identity):
        x = x.rstrip("\n")
        arr = x.split()
        if arr[0] in ident.keys():
            ident[arr[0]][arr[1]] = [arr[7],arr[8],arr[9]]
        else:
            ident[arr[0]] = {}
            ident[arr[0]][arr[1]] = [arr[7],arr[8],arr[9]]
    for y in open(overlap):
        y = y.rstrip("\n")
        over = y.split()
        if int(refSizes[over[0]]) > int(targetSizes[over[1]]):
            size = int(refSizes[over[0]])
        else:
            size = int(targetSizes[over[1]])
        percent = int(int(over[2])*100/size)
    #    percent = int(int(over[2])*100/((int(refSizes[over[0]])+int(targetSizes[over[1]]))/2))
        if over[0] in ident.keys() and over[1] in ident[over[0]].keys():
#            output = reference+"_"+over[0]+"\t"+target+"_"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+ident[over[0]][over[1]][0]+"\t"+ident[over[0]][over[1]][1]+"\t"+ident[over[0]][over[1]][2]
            output = reference+"-"+over[0]+"\t"+target+"-"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+ident[over[0]][over[1]][0]+"\t"+ident[over[0]][over[1]][1]+"\t"+ident[over[0]][over[1]][2]
#            output = over[0]+"\t"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+ident[over[0]][over[1]][0]+"\t"+ident[over[0]][over[1]][1]+"\t"+ident[over[0]][over[1]][2]
#            output += target+"_"+over[1]+"\t"+reference+"_"+over[0]+"\t"+over[2]+"\t"+str(percent)+"\t"+ident[over[0]][over[1]][0]+"\t"+ident[over[0]][over[1]][2]+"\t"+ident[over[0]][over[1]][1]
#            output = over[0]+"\t"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+ident[over[0]][over[1]][0]+"\t"+ident[over[0]][over[1][1]]+"\t"+ident[over[0]][over[1]][2]
            overlapper.write(output+"\n")
        else:
#            output = reference+"_"+over[0]+"\t"+target+"_"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+"0\t0\t0"
            output = reference+"-"+over[0]+"\t"+target+"-"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+"0\t0\t0"
#            output = over[0]+"\t"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+"0\t0\t0"
 #           output += target+"_"+over[0]+reference+"_"+"\t"+over[1]+"\t"+over[2]+"\t"+str(percent)+"\t"+"0\t0\t0"
            overlapper.write(output+"\n")
    overlapper.close()
    
    
def getRbh(path):
    f = open(path+".v2best","w+")
#    fIdent = open(path,"r")
#    ls = fIdent.readlines()
    genePair = {}
    filterOne = {}
    filterTwo = {}
    for l in open(path,"r"):
       #arr = l.split("\t")
       arr = l.split()
       filterOne[arr[0]] = []
       if arr[1] in genePair.keys():
           genePair[arr[1]].append([arr[0],arr[7],arr[8]])
       else:
           genePair[arr[1]]=[]
           genePair[arr[1]].append([arr[0],arr[7],arr[8]])
 #      if arr[1] in genePair.keys():
  #         genePair[arr[1]][arr[0]] = [arr[7],arr[8]] 
  #     else:
  #         genePair[arr[1]] = {}
  #         genePair[arr[1]][arr[0]] = [arr[7],arr[8]]
#    f.close()
    for gene1 in genePair.keys():
#        sortedGene = sorted(genePair.items(), key=lambda d: (d[1][1],d[1][2]),reverse = True)
        sortedGene = sorted(genePair[gene1], key=lambda x: (x[1],x[2]),reverse = True)
        print sortedGene
#        filterOne[sortedGene[0][0]]= []
        filterOne[sortedGene[0][0]].append([gene1,sortedGene[0][1],sortedGene[0][2]])
    for gene2 in filterOne.keys():

#        print filterOne[gene2][2]
        if filterOne[gene2]:
           sortedGeneTwo = sorted(filterOne[gene2], key=lambda x: (x[1],x[2]),reverse = True)
        #filterTwo[sortedGeneTwo[gene2][0][0]] = [gene2,sortedGeneTwo[gene2][0][1],sortedGeneTwo[gene2][0][2]]
           filterTwo[sortedGeneTwo[0][0]] = [gene2,sortedGeneTwo[0][1],sortedGeneTwo[0][2]]
        else:
            pass
    for gene in filterTwo.keys():
#        f.write(gene+"\t"+filterTwo[gene][0][0]+"\t"+str(filterTwo[gene][0][1])+"\t"+str(filterTwo[gene][0][2])+"\n")
        f.write(gene+"\t"+filterTwo[gene][0]+"\t"+str(filterTwo[gene][1])+"\t"+str(filterTwo[gene][2])+"\n")
    f.close()    
    
