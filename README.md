# Tn-Seq Analysis Scripts


<b> slidingWindow.pl</b>

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
  
