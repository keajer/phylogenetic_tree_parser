sequenceFile = "sequences.txt"
nucleotides = []

# these lists will contain the nucleotides that match the pattern
patternAAA = []
patternACC = []
patternAAC = []
patternACA = []
patternACT = []

# read in the 3 sequences from the text file
f = open(sequenceFile)
sequences = f.readlines()

# go through each sequence and separate
for index, sequence in enumerate(sequences):
	sequence = sequence.strip()
	sequences[index] = list(sequence)	# split string into characters
	
sequenceCharCount = len(sequences[0])	# get the number of characters in each sequence

# now go through and get the nucleotides vertically
for i in range(0, sequenceCharCount):
	nucleotides.append(list(sequences[0][i] + sequences[1][i] + sequences[2][i]))

# go thorugh nucleotides and sort them into respective patterns
for nucleotide in nucleotides:
	if nucleotide[0] == nucleotide[1] and nucleotide[0] == nucleotide[2]:		# if AAA pattern
		patternAAA.append(nucleotide[0] + nucleotide[1] + nucleotide[2])
		
	elif nucleotide[0] != nucleotide[1] and nucleotide[1] == nucleotide[2]:	# if ACC pattern
		patternACC.append(nucleotide[0] + nucleotide[1] + nucleotide[2])
		
	elif nucleotide[0] == nucleotide[1] and nucleotide[1] != nucleotide[2]:	# if AAC pattern
		patternAAC.append(nucleotide[0] + nucleotide[1] + nucleotide[2])
		
	elif nucleotide[0] == nucleotide[2] and nucleotide[1] != nucleotide[2]:	# if ACA pattern
		patternACA.append(nucleotide[0] + nucleotide[1] + nucleotide[2])
		
	else:	# if nothing else, then it is ACT pattern
		patternACT.append(nucleotide[0] + nucleotide[1] + nucleotide[2])
	
# now we can easily get counts of each pattern
countAAA = len(patternAAA)
countACC = len(patternACC)
countAAC = len(patternAAC)
countACA = len(patternACA)
countACT = len(patternACT)

# display counts
print("Pattern Counts\n\nAAA: " + str(countAAA) + \
      "\nACC: " + str(countACC) + \
	  "\nAAC: " + str(countAAC) + \
	  "\nACA: " + str(countACA) + \
	  "\nACT: " + str(countACT) + \
	  "\nTotal: " + str(sequenceCharCount))