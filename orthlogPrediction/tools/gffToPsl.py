import re

def readSize(ChromSizes):
        f=open(ChromSizes,"r")
        lines=f.readlines()
        size={}
        for l in lines:
                l=l.rstrip("\n")
                arr=l.split("\t")
                size[arr[0]]=arr[1]
        return size
        f.close()
def gffToPsl(gff,qChromSize,rChromSize):
        qChromSizes=readSize(qChromSize)
#        print  qChromSizes
        rChromSizes=readSize(rChromSize)
        f=open(gff,"r")
        psl=open(gff+".psl",'w')
        lines=f.readlines()
        count=0
        for l in lines:
                l=l.rstrip("\n")
                arr=l.split("\t")
                count+=1
                if l.find("mRNA")!=-1:
                        blockSizes=list()
                        blockCount=0
                        tBaseInsert=0
                        tStarts=list()
                        qStarts=list()
                        blockEnd=list()
                        tName=arr[0]
                        qName=re.search(r'ID=([^;]+);',arr[-1]).group(1)
                        strand=arr[6]
                        tStart=int(arr[3])-1
                        tEnd=arr[4]
                else:
                        tStarts.append(str(int(arr[3])-1)+",")
                        blockSizes.append(int(arr[4])-int(arr[3]))
                        blockEnd.append(arr[4])
                        blockCount+=1
                        if len(blockEnd)>1:
                                tBaseInsert+=int(arr[3])-int(blockEnd[-2])+1
                        if len(blockSizes)>1:
                                qStarts.append(str(sum(blockSizes[:-1]))+",")
                        else:
                                qStarts.append("0,")
                        if count==len(lines) or lines[count].find("mRNA")!=-1:
                                blockSizes = [str(x) for x in blockSizes]
#                                print str(qChromSizes[qName])
                                output=str(qChromSizes[qName])+"\t"+str(0)+"\t"  \
                                        +str(0)+"\t"+str(0)+"\t"+str(0)+"\t"+str(0)+"\t" \
                                        +str(blockCount-1)+"\t" \
                                       +str(tBaseInsert)+"\t"+strand+"\t"+qName+"\t"+str(qChromSizes[qName])+"\t"   \
                                       +str(0)+"\t"   \
                                       +str(qChromSizes[qName])+"\t"+tName+"\t"+str(rChromSizes[tName])+"\t" \
                                       +str(tStart)+"\t"+str(tEnd)+"\t"+str(blockCount)+"\t"+",".join(blockSizes)+","+ \
                                       "\t"+"".join(qStarts)+"\t"+"".join(tStarts)+"\n"
                                psl.write(output)
                                       
        f.close()
        psl.close()
