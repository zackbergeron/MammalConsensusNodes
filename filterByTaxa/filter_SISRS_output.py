import sys
import glob
from Bio import SeqIO

filtertable = sys.argv[1]
fastafolder = sys.argv[2]
#make this dir!
outputfolder = sys.argv[3]
minseqdata = float(sys.argv[4])
mintaxa = int(sys.argv[5])
groupfilter = int(sys.argv[6])

taxa_to_keep = {}
with open(filtertable) as thandle:
	for line in thandle:
		splitline = line.strip().split(",")
		taxa_to_keep[splitline[1]] = splitline[0]

with open("logfile2.txt", "w") as loghandle:
	print("log started", file = loghandle)

counter = 0
counter1 = 0
counter2 = 0
for f in glob.glob(fastafolder+"/*.fasta"):
	grouplist = []
	seqs_for_output = {}
	seqs = SeqIO.parse(f, "fasta")
	for seq in seqs:
		if float(len(str(seq.seq).upper().replace("N", "")))/len(seq.seq) > minseqdata:
			seqs_for_output[seq.id] = seq.seq
			grouplist.append(taxa_to_keep[seq.id])
	if len(set(grouplist)) >= groupfilter and len(seqs_for_output) >= mintaxa:
		with open(outputfolder+"/"+f.split("/")[-1], "w") as outhandle:
			for seq in sorted(seqs_for_output):
				print (">"+seq, file=outhandle)
				print (seqs_for_output[seq], file=outhandle)
		counter += 1
	counter1 += 1
	if counter1 - counter2 >= 10000:
		with open("logfile2.txt", "a") as loghandle:
			print ("processed", counter1, "loci, wrote", counter, "loci", file=loghandle)
		counter2 = counter1
print ("done")