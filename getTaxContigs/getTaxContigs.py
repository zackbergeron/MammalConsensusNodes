from Bio import SeqIO
import glob
import sys

taxaFile = sys.argv[1]
alignmentFolder = sys.argv[2]
outputFolder = sys.argv[3]

taxaList = {}
with open(taxaFile) as taxaHandle:
    for line in taxaHandle:
        tax = line.strip()
        taxaList[tax] = open(outputFolder+"/"+tax+".fas", "w")

files = glob.glob(alignmentFolder+"/*.fasta")

for f in files:
    fnew = ".".join(f.split("/")[-1].split(".")[:-1])
    seqs = SeqIO.parse(f, "fasta")
    for seq in seqs:
        if seq.id in taxaList:
            print (">"+fnew+"\n"+str(seq.seq).upper().replace("-","N"), file=taxaList[seq.id])

for f in taxaList:
    taxaList[f].close()