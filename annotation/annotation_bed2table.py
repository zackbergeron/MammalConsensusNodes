# this script takes annotated BED file and 
# produces a count table of different marker types of interest for each (SISRS) locus
# currently the whole annotation table is stored in RAM, watch out when the files are large

import sys

infilename = sys.argv[1]
option = sys.argv[2]

# function to determine output feature type given the input features
# for example, combination of gene and exon but no CDS will yield UTR,
# while gene with no exon will yield intron
def determineType(inplist):
	#prep input / output
	typelist = [x.split("@")[0] for x in inplist]
	returnSet = set([])
	
	#features to ignore for deciding on what's not annotated vs 'other' features
	#the entire genome is annotated with region, so not helpful
	#match is unclear, also ignoring
	#mRNA and transcript are unnecessary given gene exon and CDS features, so ignoring as well
	ignored_features = ["region", "match","mRNA","transcript"]
	#feature types, that are not in ignored and not in special, will be set to "other" in the output
	special_features = ["pseudogene","gene","exon","CDS","lnc_RNA"]

	#parse special features
	if "pseudogene" in typelist:
		returnSet.add("pseudogene")
	if "gene" in typelist:
		if "exon" in typelist:
			if "CDS" in typelist:
				returnSet.add ("CDS")
			else:
				returnSet.add ("UTR")
		else:
			returnSet.add ("intron")
	if "lnc_RNA" in typelist:
		returnSet.add ("lnc_RNA")
	#figure out if there are any features that qualify for 'other'
	for ft in typelist:
		if ft not in ignored_features and ft not in special_features:
			returnSet.add ("other")
	#if not then set as unannotated
	if len(returnSet) == 0:
		returnSet.add ("unannotated")
	return returnSet

# main process
d = {}
with open(infilename) as infhandle:
	for line in infhandle:
		splitline = line.strip().split('\t')
		# get SISRS locus name
		locus = splitline[3]
		
		# get feature string: featureType@Start@End@Strand
		feature = splitline[8]+"@"+"@".join(splitline[9:11])+"@"+splitline[12]
		
		# get SISRS locus annotation string: Chromo@Start@End@Strand
		# normally each SISRS locus only mapped once so there should be
		# only 1 annotationR entry per SISRS locus
		annotationR = splitline[0]+"@"+"@".join(splitline[1:3])+"@"+splitline[5]

		# add data to dictionary - this will be kept in RAM
		if locus in d:
			if annotationR in d[locus]:
				d[locus][annotationR].add(feature)
			else:
				d[locus][annotationR] = set([feature])
		else:
			d[locus] = {annotationR: set([feature])}

# prepare for output
print ("locus,pseudogene,CDS,UTR,intron,lnc_RNA,other,unannotated")

# for each sisrs contig
for locus in d:
	# for each different mapping of sisrs contig (ideally, only 1)
	for annotationR in d[locus]:
		#initialize the dict to store counts for this mapping of the SISRS locus (normally the ONLY mapping)
		output_feature_dict = {"pseudogene":0, "CDS":0, "UTR":0, "intron":0, "lnc_RNA":0, "other":0, "unannotated":0}

		#get SISRS locus info
		locusstart = int(annotationR.split("@")[1])
		locusend = int(annotationR.split("@")[2])
		locusstrand = annotationR.split("@")[3]

		#parse all the features for this SISRS contig and mapping
		features = d[locus][annotationR]
		featuretypes = [x.split("@")[0] for x in features]
		featurestart = [int(x.split("@")[1])-1 for x in features]
		featureend = [int(x.split("@")[2]) for x in features]
		featurestrand = [x.split("@")[3] for x in features]

		#collect all features and their ranges:
		# for every SISRS locus position mark which features span it
		sitecounter = 0
		posdict = {}
		for site in range(locusstart, locusend):
			posdict[sitecounter] = []
			for featureN in range(len(features)):
				if featurestart[featureN] <= site and featureend[featureN] >= site:
					comboName = featuretypes[featureN]+"@"+featurestrand[featureN]
					posdict[sitecounter].append(comboName)
			sitecounter += 1
		
		# iterate over the obtained posdict
		# to figure out feature type for a given region,
		# considering the original features that overlapped there
		final_annotation = {}
		poscounter = 0
		currentType = None
		currentCounter = poscounter
		for seqpos in posdict:
			if currentType == None:
				currentType = determineType(posdict[seqpos])
			else:
				if currentType != determineType(posdict[seqpos]):
					final_annotation[str(currentCounter)+"@"+str(poscounter)] = currentType
					currentType = determineType(posdict[seqpos])
					currentCounter = poscounter
			poscounter+=1
		final_annotation[str(currentCounter)+"@"+str(poscounter)] = currentType

		#generate output - two options here:
		# get counts of each feature type for a locus or get sequence length proportions of each feature type for a locus
		
		if option == "l":
			total_annotation_bases = 0
		for fk, fv in final_annotation.items():
			for ft in fv:
				if option == "c":
					output_feature_dict[ft] += 1
				elif option == "l":
					annotation_length = int(fk.split("@")[1])-int(fk.split("@")[0])
					output_feature_dict[ft] += annotation_length
					total_annotation_bases += annotation_length
		if option == "c":
			print (locus+","+",".join([str(x) for x in output_feature_dict.values()]))
		elif option == "l":
			print (locus+","+",".join([str(round(x/total_annotation_bases,3)) for x in output_feature_dict.values()]))
