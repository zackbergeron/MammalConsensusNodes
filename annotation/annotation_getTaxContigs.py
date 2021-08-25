from Bio import SeqIO
import glob
import sys

taxon = sys.argv[1]
taxonFile = open(taxon+".fasta", "w")
alignmentFolder = sys.argv[2]



files = glob.glob(alignmentFolder+"/*.fasta")

for f in files:
    fnew = ".".join(f.split("/")[-1].split(".")[:-1])
    seqs = SeqIO.parse(f, "fasta")
    for seq in seqs:
        if seq.id == taxon:
            print (">"+fnew+"\n"+str(seq.seq).upper().replace("-","N"), file=taxonFile)

taxonFile.close()