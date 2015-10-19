# Blueberry Scripts


<b> Blueberry.pl</b>

  <u>Description</u>: scans the genome with a fixed size sliding window, calculating a) organismal fitness and b) essentiality based on insertion representation calculations for the genomic region included in the window<br />
  <u>Usage</u>: Input .csv file(s) (outputted from calc_fitness.pl) with reference genome. Outputs following files into a single directory: 
  
  - <b>unmatched.txt</b> Genome indexes of insertion sites from Tn-Seq data that could not be matched with a TA site in the given genome at that exact index
  - <b>allTAsites.txt</b> All TA sites that occur in the genome and whether an insertion occurred there
  - <b>nullDist.txt</b> Statistics (size, mean, variance, etc) on null distribution libraries used in the essentials analysis
  - <b>essentialWindows.csv</b> All information about each window (including fitness and p-value for essentiality)
  - <b>essentialWindows.wig</b> Only p-values for each window in wig file for Integrated Genome Viewer (IGV)
  - <b>fitWindows.csv</b> Only fitness-related information for each window
  - <b>fitWindows.wig</b> Only fitness for each window in wig file for Integrated Genome Viewer (IGV)
  - <b>groupedWindows.txt</b> Grouped consecutive windows with same fitness value 
  - <b>log.txt</b> If the --log option was selected, all standard output goes here
  
<b>merge.pl</b>

  Description: Merges and sorts multiple comma-separated (.csv) files<br />
  Usage: Input .csv files, outputs merge sorted file
  
<b>countTA.pl</b>

  Description: Counts number of TA sites in the genome<br />
  Usage: input genome, outputs number of TA sites

<b>customWindow.pl</b>

  Description: Given a tab delimited text file of genome region name, start position, and end position, customWindow.pl 
	outputs Tn-Seq stats for those regions: # of mutants, # of unique insertions, # of TA sites, ratio of insertions
	to TA sites, aggregate fitness value, p-value for essentiality <br />
  Usage: input GBK file, fasta file for genome sequence, results files, custom region file	

<b>getCoordsGBK.py<\b>

  Description: Get gene names and coordinates from a GBK file.<br /> 
	Output file can be inputted into Gviz as an annotation track for visualization.<br />
  Usage: python getCoordsGBK.py <input GBK file> <name of output file>

<b>singleFit.pl<\b>
  
  Description: Gets a single fitness value for each insertion based on multiple mutants (lanes)<br />
  Usage: perl ../singleFit.pl <all results.csv files>



  
