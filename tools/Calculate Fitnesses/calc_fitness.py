# This script requires BioPython, which in turn has a good number of dependencies (some optional but very helpful).
# How to install BioPython and a list of its dependencies can be found here: http://biopython.org/DIST/docs/install/Installation.html
# K. McCoy









##### ARGUMENTS #####

def print_usage():
	print "\n" + "You are missing one or more required flags. A complete list of flags accepted by calc_fitness is as follows:" + "\n\n"
	print "\033[1m" + "Required" + "\033[0m" + "\n"
	print "-ref" + "\t\t" + "The name of the reference genome file, in GenBank format." + "\n"
	print "-t1" + "\t\t" + "The name of the bowtie mapfile from time 1." + "\n"
	print "-t2" + "\t\t" + "The name of the bowtie mapfile from time 2." + "\n"
	print "-out" + "\t\t" + "Name of a file to enter the .csv output." + "\n"
	print "\n"
	print "\033[1m" + "Optional" + "\033[0m" + "\n"
	print "-expansion" + "\t\t" + "Expansion factor (default: 250)" + "\n"
	print "-reads1" + "\t\t" + "The number of reads to be used to calculate the correction factor for time 0." + "\n\t\t" + "(default counted from bowtie output)" + "\n"
	print "-reads2" + "\t\t" + "The number of reads to be used to calculate the correction factor for time 6." + "\n\t\t" + "(default counted from bowtie output)" + "\n"
	print "-cutoff" + "\t\t" + "Discard any positions where the average of counted transcripts at time 0 and time 1 is below this number (default 0)" + "\n"
	print "-cutoff2" + "\t\t" + "Discard any positions within the normalization genes where the average of counted transcripts at time 0 and time 1 is below this number (default 10)" + "\n"
	print "-strand" + "\t\t" + "Use only the specified strand (+ or -) when counting transcripts (default: both)" + "\n"
	print "-normalize" + "\t" + "A file that contains a list of genes that should have a fitness of 1 - used for normalization and bottleneck calculations." + "\n"
	print "-b" + "\t" + "Calculate bottleneck value (the percentage of insertions randomly lost) from all genes (rather than only normalization genes)" + "\n"
	print "-maxweight" + "\t" + "The maximum weight a transposon gene can have in normalization calculations" + "\n"
	print "-multiply" + "\t" + "Multiply all fitness scores by a certain value (e.g., the fitness of a knockout). You should normalize the data." + "\n"
	print "-ef" + "\t\t" + "Exclude insertions that occur in the first N amount (%) of gene--becuase may not affect gene function." + "\n"
	print "-el" + "\t\t" + "Exclude insertions in the last N amount (%) of the gene--considering truncation may not affect gene function." + "\n"
	print "-wig" + "\t\t" + "Create a wiggle file for viewing in a genome browser. Provide a filename." + "\n"
	print "-uncol" + "\t\t" + "Use if reads were uncollapsed when mapped." + "\n"
	print "\n"

import argparse 
parser = argparse.ArgumentParser()
parser.add_argument("-ref", action="store", dest="ref_genome")
parser.add_argument("-t1", action="store", dest="mapfile1")
parser.add_argument("-t2", action="store", dest="mapfile2")
parser.add_argument("-out", action="store", dest="outfile")
parser.add_argument("-out2", action="store", dest="outfile2")
parser.add_argument("-expansion", action="store", dest="expansion_factor")
parser.add_argument("-reads1", action="store", dest="reads1")
parser.add_argument("-reads2", action="store", dest="reads2")
parser.add_argument("-cutoff", action="store", dest="cutoff")
parser.add_argument("-cutoff2", action="store", dest="cutoff2")
parser.add_argument("-strand", action="store", dest="usestrand")
parser.add_argument("-normalize", action="store", dest="normalize")
parser.add_argument("-b", action="store", dest="bottleall")
parser.add_argument("-maxweight", action="store", dest="max_weight")
parser.add_argument("-multiply", action="store", dest="multiply")
parser.add_argument("-ef", action="store", dest="exclude_first")
parser.add_argument("-el", action="store", dest="exclude_last")
parser.add_argument("-wig", action="store", dest="wig")
parser.add_argument("-uncol", action="store", dest="uncol")
arguments = parser.parse_args()

if (not arguments.ref_genome or not arguments.mapfile1 or not arguments.mapfile2 or not arguments.outfile):
	print_usage() 
	quit()

# Sets the default value of the expansion factor to 250, which is a trivial placeholder number.
	
if (not arguments.expansion_factor):
	arguments.expansion_factor = 250

# 75 is similarly trivial

if (not arguments.max_weight):
	arguments.max_weight = 75
	
# Sets the default value of cutoff to 0; cutoff exists to discard positions with a low number of counted transcripts, because fitnesses calculated from them may not be very accurate, by the same reasoning that studies with low sample sizes are innacurate. 
	
if (not arguments.cutoff):
	arguments.cutoff = 0
	
# Sets the default value of cutoff2 to 10; cutoff2 exists to discard positions within normalization genes with a low number of counted transcripts, because fitnesses calculated from them similarly may not be very accurate.
# This only has an effect if it's larger than cutoff, since the normalization step references a list of insertions already affected by cutoff.
	
if (not arguments.cutoff2):
	arguments.cutoff2 = 10

if (not arguments.usestrand):
	arguments.usestrand = "both"
	
	
	
	
	
	
##### PARSING THE REFERENCE GENOME #####

def get_time():
	import datetime
	return datetime.datetime.now().time()
print "\n" + "Starting: " + str(get_time()) + "\n"

from Bio import SeqIO
import os.path
handle = open(arguments.ref_genome, "rU")
for record in SeqIO.parse(handle, "genbank"):
    refname = record.id
    features = record.features
handle.close()

# Makes a dictionary out of each feature that's a gene - with its gene name, start location, end location, and strand as keys to their values. Then makes a list out of all those dictionaries for ease of accessing later on.

feature_list = []
for feature in features:
	if feature.type == "gene":
		gene = feature.qualifiers["locus_tag"]
		strand = feature.location.strand
		start = float(feature.location.start)
		end = float(feature.location.end)
		
# Exclude_first and exclude_last are used here to exclude whatever percentage of the genes you like from calculations; e.g. a value of 0.1 for exclude_last would exclude the last 10% of all genes!
# This can be useful because insertions at the very start or end of genes often don't actually break its function.
		
		if (arguments.exclude_first):
			start += (end - start) * float(arguments.exclude_first)
		if (arguments.exclude_last):
			end -= (end - start) * float(arguments.exclude_last)
		feature_dictionary = {"gene": gene, "start": start, "end": end, "strand": strand}
		feature_list.append(feature_dictionary)

print "Done generating feature lookup: " + str(get_time()) + "\n"










##### PARSING THE MAPFILES #####

with open(arguments.mapfile1) as file:
	r1 = file.readlines()
with open(arguments.mapfile2) as file:
	r2 = file.readlines()

# When called, goes through each line of the mapfile to find the strand (+/Watson or -/Crick), count, and position of the read. It may be helpful to look at how the mapfiles are formatted to understand how this code finds them. 

def read_mapfile(reads):
	plus_total = 0
	minus_total = 0
	plus_counts = {"total": 0, "sites": 0}
	minus_counts = {"total": 0, "sites": 0}	
	for read in reads:
		if (arguments.uncol):
			strand = read.split()[2]
			count = 1
			position = float(read.split()[4])
			if arguments.usestrand != "both" and strand != arguments.usestrand:
				continue
			if (strand == "+"):		
				sequence_length = len(read.split()[5])
				position += (sequence_length - 2)
				plus_counts["total"] += count
				plus_counts["sites"] += 1
				if position in plus_counts:
					plus_counts[position] += count
				else:
					plus_counts[position] = count
			else:
				minus_counts["total"] += count
				minus_counts["sites"] += 1
				if position in minus_counts:
					minus_counts[position] += count
				else:
					minus_counts[position] = count
		else:
			if "-" in read.split()[0]:
				strand = read.split()[1]
				count = float(read.split()[0].split("-")[1])		
				position = float(read.split()[3])
			else:
				continue

# If for some reason you want to skip all reads from one of the strands - for example, if you wanted to compare the two strands - that's done here.

			if arguments.usestrand != "both" and strand != arguments.usestrand:
				continue
			
# Makes dictionaries for the + & - strands, with each insert position as a key and the number of insertions there as its corresponding value.

			if (strand == "+"):		
				sequence_length = len(read.split()[4])
				
# The -2 in "(sequence_length -2)" comes from a fake "TA" in the read; see how the libraries are constructed for further on this 

				position += (sequence_length - 2)
				plus_counts["total"] += count
				plus_counts["sites"] += 1
				if position in plus_counts:
					plus_counts[position] += count
				else:
					plus_counts[position] = count
			else:
				minus_counts["total"] += count
				minus_counts["sites"] += 1
				if position in minus_counts:
					minus_counts[position] += count
				else:
					minus_counts[position] = count
	return (plus_counts, minus_counts)

# Calls read_mapfile(reads) to parse arguments.reads1 and arguments.reads2 (your reads from t1 and t2).





(plus_ref_1, minus_ref_1) = read_mapfile(r1)
print "Read first file: " + str(get_time()) + "\n"
(plus_ref_2, minus_ref_2) = read_mapfile(r2)
print "Read second file: " + str(get_time()) + "\n"

# The lines below are just printed for reference. The number of sites is the length of a given dictionary of sites - 1 because its last key, "total", isn't actually a site.

print "Reads:" + "\n"
print "1: + " + str(plus_ref_1["total"]) + " - " + str(minus_ref_1["total"]) + "\n"
print "2: + " + str(plus_ref_2["total"]) + " - " + str(minus_ref_2["total"]) + "\n"
print "Sites:" + "\n"
print "1: + " + str(plus_ref_1["sites"]) + " - " + str(minus_ref_1["sites"]) + "\n"
print "2: + " + str(plus_ref_2["sites"]) + " - " + str(minus_ref_2["sites"]) + "\n"










##### FITNESS CALCULATIONS #####

# If reads1 and reads2 weren't specified in the command line, sets them as the total number of reads (found in read_mapfile())

if not arguments.reads1:
	arguments.reads1 = plus_ref_1["total"] + minus_ref_1["total"]
if not arguments.reads2:
	arguments.reads2 = plus_ref_2["total"] + minus_ref_2["total"]

# Calculates the correction factors for reads from t1 and t2; cfactor1 and cfactor2 are the number of reads from t1 and t2 respectively divided by total, which is the average number of reads between the two.
# This is used later on to correct for pipetting errors, or any other error that would cause unequal amounts of DNA from t1 and t2 to be sequenced so that an unequal amount of reads is produced

total = (float(arguments.reads1) + float(arguments.reads2))/2
cfactor1 = float(arguments.reads1)/total
cfactor2 = float(arguments.reads2)/total
print "Cfactor 1: " + str(cfactor1) + "\n"
print "Cfactor 2: " + str(cfactor2) + "\n"
import math
import csv
results = [["position", "strand", "count_1", "count_2", "ratio", "mt_freq_t1", "mt_freq_t2", "pop_freq_t1", "pop_freq_t2", "gene", "D", "W", "nW"]]
genic = 0
total_inserts = 0
with open(arguments.ref_genome, "r") as file:
	firstline = file.readline()
genomelength = firstline.split()[2]
i = 0
while i < float(genomelength):

# At each possible location for an insertion in the genome, counts the number of actual insertions at t1 and which strand(s) the corresponding reads came from.

	c1 = 0
	if i in plus_ref_1:
		c1 = float(plus_ref_1[i])
		strand = "+/"
		if i in minus_ref_1:
			c1 += float(minus_ref_1[i])
			strand = "b/"
	elif i in minus_ref_1:
		c1 = float(minus_ref_1[i])
		strand = "-/"

# If there were no insertions at a certain location at t1 just continues to the next location; there can't be any comparison to make between t1 and t2 if there are no t1 insertions!

	else:
		i += 1
		continue
		
# At each location where there was an insertion at t1, counts the number of insertions at t2 and which strand(s) the corresponding reads came from.

	c2 = 0
	if i in plus_ref_2:
		c2 = float(plus_ref_2[i])
		if i in minus_ref_2:
			c2 += float(minus_ref_2[i])
			strand += "b"
		else:
			strand += "+"
	elif i in minus_ref_2:
		c2 = float(minus_ref_2[i])
		strand += "-"

# Corrects with cfactor1 and cfactor2

 	c1 /= cfactor1
 	if c2 != 0:
 		c2 /= cfactor2
 		ratio = c2/c1
 	else:
 		c2 = 0
 		ratio = 0

# Passes by all insertions with a number of reads smaller than the cutoff, as they may lead to inaccurate fitness calculations.

	if (c1 + c2)/2 < float(arguments.cutoff):
		i+= 1
		continue
		
# Calculates each insertion's frequency within the populations at t1 and t2.

 	mt_freq_t1 = c1/total
 	mt_freq_t2 = c2/total
	pop_freq_t1 = 1 - mt_freq_t1
	pop_freq_t2 = 1 - mt_freq_t2
	
# Calculates each insertion's fitness! This is from the fitness equation log((frequency of mutation @ time 2 / frequency of mutation @ time 1)*expansion factor)/log((frequency of population without the mutation @ time 2 / frequency of population without the mutation @ time 1)*expansion factor)

	w = 0
	if mt_freq_t2 != 0:
		top_w = math.log(mt_freq_t2*(float(arguments.expansion_factor)/mt_freq_t1))
		bot_w = math.log(pop_freq_t2*(float(arguments.expansion_factor)/pop_freq_t1))
		w = top_w/bot_w
		
# Checks which gene locus the insertion falls within, and records that.
	
	gene = ''
	for feature_dictionary in feature_list:
		if feature_dictionary["start"] <= i and feature_dictionary["end"] >= i:
			gene = "".join(feature_dictionary["gene"])
			genic += 1
			break
	total_inserts += 1
		
# Writes all relevant information on each insertion and its fitness to a cvs file: the location of the insertion, its strand, c1, c2, etc. (the variable names are self-explanatiory)
# w is written twice, because the second w will be normalized if normalization is called for, thus becoming nW.

	row = [i, strand, c1, c2, ratio, mt_freq_t1, mt_freq_t2, pop_freq_t1, pop_freq_t2, gene, arguments.expansion_factor, w, w]
	results.append(row)
	i += 1
with open(arguments.outfile, "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(results)

print "Done comparing mapfiles " + str(get_time()) + "\n"
print "Genic: " + str(genic) + "\n"
print "Total: " + str(total_inserts) + "\n"










##### BOTTLENECK VALUE CALCULATION #####

#the bottleneck value is calculated here if done from all genes - otherwise it's done in the normalization section if only taken from normalization genes

if (arguments.bottleall):
	overall_blank_count = 0
	for list in results:
		if (list[2] != 0 and list[3] == 0):
			overall_blank_count += 1
	overall_original_count = len(results)
		
	pc_blank_normals = float(overall_blank_count) / float(overall_original_count)
	
	with open(arguments.outfile2, "w") as f:
		f.write("bottleneck_value: " + str(pc_blank_normals) + "\n" + "total: " + str(total) + "\n" + "refname: " + refname)










##### NORMALIZATION #####

# If making a WIG file is requested in the arguments, starts a string to be added to and then written to the WIG file with a typical WIG file header.
# The header is just in a typical WIG file format; if you'd like to look into this more UCSC has notes on formatting WIG files on their site.

if (arguments.wig):
	wigstring = "track type=wiggle_0 name=" + arguments.wig + "\n" + "variableStep chrom=" + refname + "\n"

# Takes normalization genes (which should all be predicted or known to have fitness values of exactly 1.0, like transposons for example) and uses them to normalize the fitnesses of all insertion locations

if (arguments.normalize):
	with open(arguments.normalize) as file:
		transposon_genes = file.read().splitlines()
	print "Normalize genes loaded" + "\n" 
	blank_ws = 0
	sum = 0
	count = 0
	weights = []
	scores = []
	for list in results:
		if list[9] != '' and list[9] in transposon_genes:
			c1 = list[2]
			c2 = list[3]
			score = list[11]
			avg = (c1 + c2)/2
			
# Skips over those insertion locations with too few insertions - their fitness values are less accurate because they're based on such small insertion numbers.
			
			if float(c1) >= float(arguments.cutoff2):
			
# Sets a max weight, to prevent insertion location scores with huge weights from unbalancing the normalization.
			
				if (avg >= float(arguments.max_weight)):
					avg = float(arguments.max_weight)
                
# Tallies how many w values are 0 within the blank_ws value; you might get many transposon genes with a w value of 0 if a bottleneck occurs, for example, which is especially common with in vivo experiments. This is used later by aggregate.py
# For example, when studying a nasal infection in a mouse model, what bacteria "sticks" and is able to survive and what bacteria is swallowed and killed or otherwise flushed out tends to be a matter of chance not fitness; all mutants with an insertion in a specific transposon gene could be flushed out by chance!

				if score == 0:
					blank_ws += 1

				sum += score
				count += 1
				weights.append(avg)
				scores.append(score)
				
				print str(list[9]) + " " + str(score) + " " + str(c1)

# Counts and removes all "blank" fitness values of normalization genes - those that = 0 - because they most likely don't really have a fitness value of 0, and you just happened to not get any reads from that location at t2. 
    
	blank_count = 0
	original_count = len(scores)
	curr_count = original_count
	i = 0
	while i < curr_count:
		w_value = scores[i]
		if w_value == 0:
			blank_count += 1
			weights.pop(i)
			scores.pop(i)
			i -= 1
			curr_count = len(scores)
		i += 1

# If no normalization genes can pass the cutoff, normalization cannot occur, so this ends the script advises the user to try again and lower cutoff and/or cutoff2.
	
	if len(scores) == 0:
		print 'ERROR: The normalization genes do not have enough reads to pass cutoff and/or cutoff2; please lower one or both of those arguments.' + "\n"
		quit()
	
	pc_blank_normals = float(blank_count) / float(original_count)
	print "# blank out of " + str(original_count) + ": " + str(pc_blank_normals) + "\n"
	if (not arguments.bottleall):
		with open(arguments.outfile2, "w") as f:
			f.write("bottleneck_value: " + str(pc_blank_normals) + "\n" + "total: " + str(total) + "\n" + "refname: " + refname)
    
	average = sum / count
	i = 0
	weighted_sum = 0
	weight_sum = 0
	while i < len(weights):
		weighted_sum += weights[i]*scores[i]
		weight_sum += weights[i]
		i += 1
	weighted_average = weighted_sum/weight_sum
    
	print "Normalization step:" + "\n"
	print "Regular average: " + str(average) + "\n"
	print "Weighted Average: " + str(weighted_average) + "\n"
	print "Total Insertions: " + str(count) + "\n"

	old_ws = 0
	new_ws = 0
	wcount = 0
	
	for list in results:
		if list[11] == 'W':
			continue
		new_w = float(list[11])/weighted_average
		
# Sometimes you want to multiply all the fitness values by a constant; this does that.
# For example you might multiply all the values by a constant for a genetic interaction screen - where Tn-Seq is performed as usual except there's one background knockout all the mutants share.
		
		if arguments.multiply:
			new_w *= float(arguments.multiply)
	
		if float(list[11]) > 0:
			old_ws += float(list[11])
			new_ws += new_w
			wcount += 1
			
		list[12] = new_w
		
		if (arguments.wig):
			wigstring += str(list[0]) + " " + str(new_w) + "\n"
      	   
	old_w_mean = old_ws / wcount
	new_w_mean = new_ws / wcount
	print "Old W Average: " + str(old_w_mean) + "\n"
	print "New W Average: " + str(new_w_mean) + "\n"

with open(arguments.outfile, "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(results)
    	
if (arguments.wig):
	if (arguments.normalize):
		with open(arguments.wig, "wb") as wigfile:
			wigfile.write(wigstring)
	else:
		for list in results:
			wigstring += str(list[0]) + " " + str(list[11]) + "\n"
		with open(arguments.wig, "wb") as wigfile:
				wigfile.write(wigstring)
				
				
				
				
				
				































































































#                                                                                                     `````````````                                                         
#                                                                                                     `````````````                                                         
#                                                                                                     ``@@@@@@@@@``                                                         
#                                                                                                     ``@@@@@@@@@```                                                        
#                                                                                                     ``@@@@@@@@@``                                                         
#                                                                                                     ``@@@@@@@@@``                                                         
#                                                                                                     ``@@@@@@@@@``                                                         
#                                                                                                     ``@@@@@@@@@``                                                         
#                                                                                                    ```@@@@@@@@#``                                                         
#                                                                                                    ```@@@@@@@@#``                                                         
#                                                                                                    ```@@@@@@@@+``                                                         
#                                                                                                    ```@@@@@@@@'``                                                         
#                                                                                                    ```@@@@@@@@;``                                                         
#                                                                                                    ```@@@@@@@@;``                                                         
#                                                                                                    ```@@@@@@@@:``                                                         
#                                                                                                    ```@@@@@@@@,``                                                         
#                                                                                                    ``.@@@@@@@@.``                                                         
#                                                                                                    ``.@@@@@@@@```                                                         
#                                                                                                    ``.@@@@@@@@```                                                         
#                                                                                                    ``.@@@@@@@@```                                                         
#                                                                                                    ``.@@@@@@@@``                                                          
#                                                                                                    ``,@@@@@@@@``                                                          
#                                                                                                    ``,@@@@@@@@``                                                          
#                                                                                                    ``.@@@@@@@@``                                                          
#                                                                                                    ```@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``'@@@@@@@@``                                                          
#                                                                                                    ``;@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``:@@@@@@@@``                                                          
#                                                                                                    ``;@@@@@@@#``                                                          
#                                                                                                  ````+@@@@@@@#`````                                                       
#                                                                                               ```````#@@@@@@@#``````                                                      
#                                                                                               `````.,@@@@@@@@@...````                                                     
#                                                                                               ``@@@@@@@@@@@@@@@@@@;``                                                     
#                                                                                               ``@@@@@@@@@@@@@@@@@@;``                                                     
#                                                                                               ```````````````````````                                                     
#                                                                                                `````````````````````                                                      
#                                                                                                   ``````.```````                                                          
#                                                                                                    ````@.''```                                                            
#                                                                                                    ```#  `;```                                                            
#                                                                                                    ``.+    @```                                                           
#                                                                                                   ```@ ````,+```                                                          
#                                                                                                  ```;;````` @```                                                          
#                                                                                                  ```@ ``````,@```                                                         
#                                                                                                 ```,+```..```@```                                                         
#                                                                                                 ```@ ``....```@```                                                        
#                                                                                                ```+' ``....```#'``                                                        
#                                                                                                ```@```......`` @```                                                       
#                                                                                               ```'+```......```'@```                                                      
#                                                                                               ```@ ``........```@```                                                      
#                                                                                              ```'#```........````@```                                                     
#                                                                                              ```@ ``..........```#,``                                                     
#                                                                                             ```'#```...........`` @```                                                    
#                                                                                             ```@``.............```.+```                                                   
#                                                                                            ```:#```.............`` #```                                                   
#                                                 ```````                                    ```@ ```.......#......``.@```                                                  
#                                               ``````````                                  ```:@```#`......@......```@```                                                  
#                                             ``````#@@@``                                  ```@ `.`:.......@.......`` @```                                                 
#                                             ```.#@###@``                                 ```:@``..`+`....`@.......```@,``                                                 
#                                            ```'@####@```                                 ```@````..@@@@@@@@#,`..#```` @```                                                
#                                           ```#####@@```                                  ``;@ ,`.,@@.      `@@..#..```''``                                                
#                                           ``:####@#````                                 ```@``@`@@            @@:...`` @```                                               
#                                          ```@#####````                                  ``,@``.@,              ,@`...``:@```                                              
#                                          ``.####@```                                   ```@.` @`                 @....``@```                                              
#                                          ``####@```                                    ``,@  @.`                  @`.````@```                                             
#                                          ``@##@````                                   ```@, @:        ;#          `@..```@.``                                             
#                                         ```@##````                                    ``.@`,@         @@,          #...`` @```                                            
#                                         ```@#@```                                    ```@, #         `@@@           @`.```;'``                                            
#                                         ```##:``                                     ``.@ +,         .@@@           ,'..`` @```                                           
#                                         ``.##```                                    ```@, @          `@@@            @`.```,+```                                          
#                                      `````@##```                                    ```@`'.           @@@            :...```@```                ``````````                
#   ````````````````````````````````````````##@```                                   ```@:`@            @@@             #...`` #```          `````````````````              
#  ```````````````````````````````````````.###@```                                   ```@ `,        .@@@++'++#@@'`      #`..```#```     ````````````'@@@@@.````             
# `````+@####################@@@@@@@@@@@@#####@```                                  ```#;`,.     `@#...,.,,,,,,..;@,    @....`` @``````````````+@@@########@.``             
# `+@##########################################,````                    ```````````````@```@   +@,.,,,,,,,,,,,,,,,,,@   @....```#`````````'@@@##############@```            
# `@###########################################@``````````````````````````````````````+'``.,'.#.,,,,,,,,,,,,,,,,,,,,.++@......`` @````+@@@#######@+``````'###'``            
# ``:@@########@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@##@``````````````````````````......,`,,.,@ ```.##.,.,,,,,,,,,,,,,,,,,,,,.##......`` :.+@@#######@@:``````````###@```           
# ````````````````````````````````,#########@###@@@@#################################@'```...@.,,,,,,,,,,,,,,,,,,,,,,,#.........`'@######@@+```````````````@##```           
#  ```````````````````````````````@#########@#########################################```.....@:,,,,,,,,,,,,,,,,,,..;@..........`@####@@:````````       ```@##@``           
#                             `````@@####@@@##########@@@@@@@@@@@@@@@@@@@@@@@@@@@@#+@+```......@#.,,,,,,,,,,,,,,,,.##..........`` #@#`````````           ```##@```          
#                               ``.#@######@####:```````````````````````````````````@ ``.......@:#,,,,,,,,,,,,,,,;@@`............`` @``````              ```@##:``          
#                               ``:########@###@```````````````````````````````````#;```......+..`##,.,,,,,,,,.#@#..'............`` @````                 ``;##@```         
#                               ```@@####@@##@'````                            ````@ ``.......'.....@@#+;:;'#@@;`...#`............`` @```                 ```@##```         
#                                ````````````````                              ```@,```.............'..:''':.@`......:............```@.``                  ``@###```        
#                                 ``````````````                               ``.@```..............#........'`....................```@```                 ``.##@``````     
#                                                                             ```@.``..............`#........,.....................```@.``                 ```@#+,``````    
#                                                                             ``.@``.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,.................```@```               ````+#####@````   
#                                                                            ```@````......,........................,....,,.,.,,,,,..` @,``               ```;@@######````  
#                                                                            ``.@```.......,......................`......,...........```@```             ```+#########@```` 
#                                                                           ```@```...........`@@+...............+@@`....,...........```@:``             ``;@#########@'``` 
#                                                                           ``.@ ``............@@@@@`.........,@@@@@.....,............```@```            ``@#########@#@@```
#                                                                          ```@.```............@@@@@@@`.....'@@@@@@@.....,............```#'``            ``@###@###@#@#@@@``
#                                                                          ``.@```.............@@@@@@@@@..+#@@@@@@@@.....,.............`` #```           ``@#@@@@##@#@#@@@``
#                                                                         ```@````.............@@@@@@@@@@@@@@@@@@@@@.....,.............```'#```          ``'#@@@@##@#@@@@@,`
#                                                                         ``.@ ``.........,....@@@@@@@@@',##@@@@@@@@`....,..............```@```          ```@@@@@##@@'#@@@@`
#                                                                        ```@.```.........,....@@@@@@@#`...`#@@@@@@@`....,..............```.@```          ``#@@@@##;```@@@@`
#                                                                        ``.@ ``.....,,,,,,,.,.@@@@@#.,,,,,,,.#@@@@@,,,,,,,,,,,,:,,,,,,,,.``@```          ``#@@@@#.````@@@@`
#                                                                       ```@. ``...............@@@;......,......#@@@`...........,.........```@```         ``#@@@;``````@@@@`
#                                                                       ```@```................@,........,........+#`...........,.........```@.``         ``#@@@;``````@@@@`
#                                                                      ```@.``...........................,.........`............,..........```@```        ``#@@@'``  ``#@@;`
#                                                                      ``.@ ``................,..........,......................,..........```#:``        ``#@@@'``  ```#@``
#                                                                     ```@,``............................,......................,...........`` @```       ``+@@@'``   ````` 
#                                                                     ``.@```............................,......................,...........```#+``       ``;@@@+``     ``  
#                                                                    ```@,``.............................,......................,............```@```      ``'@@@+``         
#                                                                    ``.@```.........,...................,......................,............```'#```     ``;@@@+``         
#                                                                   ```@:```.........,...................,......................,.............`` @```     ``:@@@+``         
#                                                                   ```@`..,,,,,,,,,,,,,,,,,,..,.........,......................,.............```'@```    ``;@@@#``         
#                                                                   ``+'```...,.................,....,,,,,,,,,,,,,,,,,,,,,,,,,,,,,..........,...``@```    ``;@@@#``         
#                                                                  ```@ ``....,.................,.......................................,.......``;@```   ``:@@@#``         
#                                                                  ``'#```....,.................,.......................................,.......```@```   ``:@@@@``         
#                                                                 ```@```.....,.................,.......................................,........``;#``   ``:@@@@``         
#                                                                 ```@ ``.....,.................,.......................................,........`` @```  ``:@@@@``         
#                                                                 ``@````...............................................................,........`` @.``  ``;@@@@``         
#                                                                 ``@  ```..............................................................,........`` .#``  ``'@@@@``         
#                                                                 ``#  ``````.`.```````````````..````````````````````..`````````````````.``````````  @``  ``'@@@@``         
#                                                                 ``.   ``````````````````````````````````````````````````````````````````````````  .;``  ``'@@@@``         
#                                                                 ``@;`      ``              `             ` `  `  ````   `  `````    ` `        `,+@```  ``+@@@@``         
#                                                                 `````:;'++##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@+.````   ``+@@@@``         
#                                                                  ```````````````````````````+##@``````````````````@#@```````````````````````````````    ``+@@@@``         
#                                                                    `````````````````````````@##@``````````````````@##;````````````````````````````      ``+@@@@``         
#                                                                                         ````###,````          ````+##@```                               ``+@@@@``         
#                                                                                          ``,###```              ``.@##```                               ``'@@@@``         
#                                                                                          ``###@``               ```@##```                               ``'@@@@``         
#                                                                                         ```@##@``                ``@##+``                               ``'@@@@``         
#                                                                                         ```###.``                ``:##@``                               ``'@@@@``         
#                                                                                         ``:###```                ```##@```                              ``'@@@@``         
#                                                                                         ``@##@``                 ```@##```                              ``'@@@@``         
#                                                                                        ```@##'``                  ``@###``                              ``'@@@@``         
#                                                                                        ```@##```                  ```##@``                              ``'@@@@``         
#                                                                                        ``,###```                  ```@#@```                             ``'@@@@``         
#                                                                                        ``####``                    ``@##.``                             ``'@@@@``         
#                                                                                        ``@##@``                    ``;##@``                             ``'@@@@``         
#                                                                                 `````````@##@``                    ```##@````                           ``;@@@@``         
#                                                                            ``````````````@##;``                    ```###``````````````                 ``;@@@@``         
#                                                                          `````````.,;.```###```                     ``@##:``````````````                ``;@@@@``         
#                                                                         `````#@#########@@##```                     ``###@@@@@@###@#@'```               ``;@@@@``         
#                                                                        ```@@###############@``                      ``,################``               ``;@@@@``         
#                                                                        ``'@################+``                      ```###############+``               ``;@@@@``         
#                                                                         ``````````````````````                       ``###########@#,````               ``.@@@@``         
#                                                                         `````````````````````                        ```````````````````                ```@@@.`          
#                                                                                                                       ````````````````                   ```````          
#                                                                                                                                                                           
#
