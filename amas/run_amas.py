import sys
import glob
import subprocess

fasta_folder = sys.argv[1]
total_cores = sys.argv[2]
files = glob.glob(fasta_folder+"/*.fasta")
amas = sys.argv[3]

fileList = []
count = 0
first = True
for f in files:
	fileList.append(f)
	count += 1
	if count == 1000:
		cmd = "python "+amas+" summary -c "+total_cores+" -o amas_output_temp.txt -f fasta -d dna -i "+" ".join(fileList)
		subprocess.call(cmd, shell=True)
		if first:
			subprocess.call("sed -e '$a\\' amas_output_temp.txt > amas_total_results.txt", shell=True)
			first = False
		else:
			subprocess.call("sed -e 1d $f amas_output_temp.txt | sed -e '$a\\' >> amas_total_results.txt", shell=True)
		count = 0
		fileList = []
if len(fileList) > 0:
	cmd = "python "+amas+" summary -c "+total_cores+" -o amas_output_temp.txt -f fasta -d dna -i "+" ".join(fileList)
	subprocess.call(cmd, shell=True)
	subprocess.call("sed -e 1d $f amas_output_temp.txt | sed -e '$a\\' >> amas_total_results.txt", shell=True)