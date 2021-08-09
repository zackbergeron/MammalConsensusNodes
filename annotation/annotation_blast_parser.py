# this script takes BLAST outfmt 6 table and produces BED file,
# while also filtering out overlapping and / or closely matching sequences
# this BED file can then be annotated using bedtools intersect with the GFF file of interest

import sys

infile = sys.argv[1]

# parameters - adjust as needed
maxovlp = 5 #(greater than 5 bp overlap is considered overlapping)
maxgap = -10000 #(gap over 10Kbp is considered to be too large to join pieces, only one will be selected)
maxscoreratio = 0.05 #(bitscores that are only 5% different are considered to be too close to accurately discriminate)

# function to parse the BLASTN format 6 and get the needed fields out
def get_hit(row):
	qstart = int(row[6])
	qend = int(row[7])
	tstart = int(row[8])
	tend = int(row[9])
	bitscore = float(row[11])
	return [qstart, qend, tstart, tend, bitscore]

# compute overlap (positive) or gap (negative) between two BLASTN HSPs
def getOverlap(a, b):
	a0=min(a)-1 #start one pos less to detect 1 bp overlap
	a1=max(a)
	b0=min(b)-1 #start one pos less to detect 1 bp overlap
	b1=max(b)
	return min(a1, b1) - max(a0, b0)

# get all retained HSPs, select the best one, and produce the BED formatted output
def output_best_hit(hit_data, query_name):
	# if any HSPs retained
	if len(hit_data) > 0:
		# only 1 Target retained
		if len(hit_data) == 1:
			# additionally check that it's not empty
			if len(list(hit_data.values())[0]) > 0:
				# get final coordinates
				tstart, tend, strand = get_final_coordinates(list(hit_data.values())[0])
				# output in BED format
				return list(hit_data)[0]+"\t"+tstart+"\t"+tend+"\t"+query_name+"\t0\t"+strand
		else:
			# multiple Targets, pick the best one
			best_hit = ""
			best_hit_score = 0
			num_best_hits = 0
			for hitkey, hitvalue in hit_data.items():
				if len(hitvalue) > 0:
					for hitregion in hitvalue:
						if best_hit_score == 0:
							best_hit = hitkey
							best_hit_score = hitregion[4]
							num_best_hits = 1
						else:
							hitratio = hitregion[4]/best_hit_score
							if 1 - maxscoreratio <= hitratio <= 1 + maxscoreratio:
								best_hit_score = max(hitregion[4], best_hit_score)
								num_best_hits = 2
							else:
								if hitratio > 1 + maxscoreratio:
									best_hit = hitkey
									best_hit_score = hitregion[4]
									num_best_hits = 1

			if num_best_hits == 1:
				# one Target considerably better, output
				tstart, tend, strand = get_final_coordinates(hit_data[best_hit])
				return best_hit+"\t"+tstart+"\t"+tend+"\t"+query_name+"\t0\t"+strand
			else:
				# several closely scored Targets, output None (conservative approach)
				return None
	else:
		return None

# given a collection of HSPs (1 or more), get the overall range and strand
# thus if the HSPs are too disjunct, the overall range will be too large...
def get_final_coordinates(hit_data_value):
	# initialize min and max coordinates on Query and on Target
	minq = 0
	maxq = 0
	mint = 0
	maxt = 0
	strand = None
	for item in hit_data_value:
		if item[0] < item[1] and item[2] < item[3]:
			#forward
			if minq == 0 or item[0] < minq:
				minq = item[0]
			if item[1] > maxq:
				maxq = item[1]
			if mint == 0 or item[2] < mint:
				mint = item[2]
			if item[3] > maxt:
				maxt = item[3]
			strand = "+"

		elif item[0] > item[1] and item[2] < item[3]:
			#reverse on Query
			if minq == 0 or item[1] < minq:
				minq = item[1]
			if item[0] > maxq:
				maxq = item[0]
			if mint == 0 or item[2] < mint:
				mint = item[2]
			if item[3] > maxt:
				maxt = item[3]
			strand = "-"

		elif item[0] < item[1] and item[2] > item[3]:
			#reverse on Target
			if minq == 0 or item[0] < minq:
				minq = item[0]
			if item[1] > maxq:
				maxq = item[1]
			if mint == 0 or item[3] < mint:
				mint = item[3]
			if item[2] > maxt:
				maxt = item[2]
			strand = "-"
			
	#format and output - only Target coordinates and strand are needed
	return str(mint), str(maxt), strand

# main procedure
with open(infile) as inhandle:
	# stores last read Query name (SISRS locus name)
	previous_q = None
	# stores last read Target / Subject name (Chromosome name)
	previous_t = None
	# stores HSPs for a given Query
	previous_hit = {}
	for line in inhandle:
		split_line = line.strip().split("\t")
		# if this is beginning of the file parsing
		if previous_q == None:
			previous_q = split_line[0]
			previous_t = split_line[1]
			previous_hit[previous_t] = [get_hit(split_line)]
		# any subsequent line
		else:
			#same Query as before
			if split_line[0] == previous_q:
				# same Target as before
				if split_line[1] == previous_t:
					# get current line needed fields
					hit1 = get_hit(split_line)
					# iterate over previously recorded HSPs and check for overlaps
					for hit2 in previous_hit[previous_t]:
						#compute overlaps on Query and on Target
						qovlp = getOverlap(hit1[:2],hit2[:2])
						tovlp = getOverlap(hit1[2:4],hit2[2:4])
						# if two HSPs overlap on Target (same chromosome region?),
						# or have too large gap on Target (too disjunct, causes problems counting annotations if the region is too large)
						# or two HSPs overlap on Query (same SISRS contig region)
						if tovlp < maxgap or tovlp > maxovlp or qovlp > maxovlp:
							# compute ratio of the new HSP score over one of the previous HSP
							hitratio = hit1[4]/hit2[4]
							# if new HSP's score is considerably worse, do nothing (skipping)
							# if new HSP's score is considerably better, replace the old with the new
							if hitratio > 1 + maxscoreratio:
								previous_hit[previous_t][previous_hit[previous_t].index(hit2)] = hit1
								# the list might get duplicates of the same record, but that's ok
							# otherwise ignore the new and remove the old
							elif 1 - maxscoreratio <= hitratio <= 1 + maxscoreratio:
								del previous_hit[previous_t][previous_hit[previous_t].index(hit2)]
								break
						else:
							#otherwise add the HSP to the current Target
							previous_hit[previous_t].append(hit1)
							break
				# different target
				else:
					# add as is
					previous_t = split_line[1]
					previous_hit[previous_t] = [get_hit(split_line)]

			# new query - first check previous output
			else:
				line_output = output_best_hit(previous_hit, previous_q)
				if line_output != None:
					print (line_output)
				# then, default all the previous vars
				previous_q = split_line[0]
				previous_t = split_line[1]
				previous_hit = {}
				previous_hit[previous_t] = [get_hit(split_line)]
	#final output - output the last results
	final_output = output_best_hit(previous_hit, previous_q)
	if final_output != None:
		print (final_output)