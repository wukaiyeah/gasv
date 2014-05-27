#! usr/bin/python
# FilterFiles.py allows you to filter GASV output from a number of different
# sources, by placing maximum and minimum number of reads required per cluster
# from each source. 

# USAGE: python FilterInputFile (outputDir)

# The input file should have a line for each gasv input file, with 3 columns, 
# containing the path to the file, the, the max number of reads (N = No limit)
# in a cluster from that file, and the min number. 

# Example input file
"""
# GASVInputFile	MaxReads	MinReads
/path/to/dir/normal.gasv.in	5	0
/path/to/dir/tumor.gasv.in	N	5
"""
# Where N = No limit

###
# GASV OPTIONS
###
gasvOptions = ""

import sys
from subprocess import call

gasvLoc = "../bin/GASV.jar"
# Unique IDs.
IDs = ['a','b','c','d','e','f','g','h','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

def withinLimits(i, reads, limits):
	count = 0
	for read in reads:
		if read.startswith(IDs[i]): count += 1
	if limits[1] != 'N':
		up = int(limits[1])
		if count > up: return False
	down = int(limits[2])
	if count < down: return False
	return True

# Read In Command Line input
if len(sys.argv) not in [2,3]:
	print "USAGE:", sys.argv[0], "InputFile (outputDir)"

filename = sys.argv[1]
outputDir = "./"
if len(sys.argv) == 3:
	outputDir = sys.argv[2]

# Reads the input file to get the list of batch files and the
# limits for each one
files = []
with open(filename) as f:
	for line in f:
		if line.startswith("#"): continue
		line = line.strip().split("\t")
		if len(line) != 3: print "Invalid input format"
		files.append((line[0], line[1], line[2]))

# Parses the batch files and to get the list of individual files
# Creates a single gasv input file containing all files
indiv_files = []
outputFile = filename + ".gasv.in"
with open(outputFile, 'w') as out:
	for i,f in enumerate(files):
		filename = f[0]
		with open(filename) as f:
			for line in f:
				if line.startswith("#"): continue
				line = line.strip().split("\t")
				indiv_files.append((line[0], i))
				outLine = "\t".join([line[0]+".temp",] + line[1:])
				out.write(outLine+"\n")

# For each individual file, we create a temporary file. For each read we append
# an ID based on which source
for f in indiv_files:
	filename = f[0]
	tempFile = filename + ".temp"
	index = f[1]
	with open(filename) as f, open(tempFile, 'w') as out :
		for line in f:
			out.write(IDs[index] + line)

# Runs GASV 
call (["java", "-jar", gasvLoc, "--output", "reads", "--outputdir", outputDir, "--batch", outputFile])

# Filter gasv input file
outputFile = outputDir + outputFile + ".clusters"
fullOutputFile = outputFile + ".filtered"
with open(outputFile) as f, open(fullOutputFile,'w') as out:
	for line in f:
		# Check whether the cluster meets the criteria for all 
		reads = line.strip().split("\t")[-1].split(",")
		MeetsCriteria = True
		for i, limits in enumerate(files):
			if not withinLimits(i, reads, limits):
				MeetsCriteria = False
		if MeetsCriteria:
			out.write("\t".join(line.strip().split("\t")[:-1]) + "\n")

# Delete all temp files
# call(["rm", outputFile])  # Comment out to keep the not filtered file
for f in indiv_files:
	filename = f[0]
	tempFile = filename + ".temp"
	call (["rm", tempFile])
