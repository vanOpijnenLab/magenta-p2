getTracks<-function(fromCoord,toCoord){
  #Load libraries
  library(seqinr)
  library(Biostrings)
  library(Gviz)
  
  genomeFile="NC_003028.fa"
  geneFile="tigr4_genes.txt"
  singleFitFile="singleFit.txt"
  

  #so Gviz doesn't look for NC_003028 in the UCSC site
  options(ucscChromosomeNames = FALSE)
  
  #GenomeAxisTrack: Number coordinates for genome (no file input)
  gTrack <- GenomeAxisTrack()
  
  #SequenceTrack: Just fasta file for NC_003028. Modify it so first line at > reads ">NC_003028"
  fcol <- c(A = "darkorange", C = "yellow", T = "darkred",G = "darkgreen")
  sTrack<-SequenceTrack(genomeFile,name="Sequence",fontcolor=fcol)
  
  #AnnotationTrack1: from NC_003028.gbk. Got genes and coordinates from genbank file using script getCoordsGBK.py
  
  # is.numeric(genes[,2])
  # Make sure coordinate columns 2 and 3 are numeric. If not then reassign as numerically converted column
  #geneFile="tigr4_new.txt"
  geneDF<-as.data.frame(read.table(file = geneFile,header=TRUE))
  anTrack<-AnnotationTrack(start=geneDF$start,end=geneDF$end,chromosome="NC_003028",strand=geneDF$strand,id=geneDF$id,showFeatureId=TRUE, name="Genes",fill="gray",fontcolor.item="black")
  
  #DataTrack: single aggregate fitness per insertion site
  #singleFitFile="singleFit.txt"
  insertFit<-as.data.frame(read.table(file=singleFitFile))
  colnames(insertFit)<-c("seqnames","start","end","fit")
  insertFit$seqnames="NC_003028"
  insertData<-makeGRangesFromDataFrame(insertFit,keep.extra.columns=TRUE,ignore.strand=TRUE,seqinfo=NULL,seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack<-DataTrack(insertData,chromosome="NC_003028",ylim=c(-1,1),type="histogram",
                      col.baseline="black",baseline=0,fill="darkred",
                      name="Insertions")
  #background.panel="#FFFFE0",
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  
  fitTrack<-AnnotationTrack(insertFit,chromosome="NC_003028",
                            showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                            fontcolor.item="black",col.frame="transparent",
                            background.title="transparent",col="transparent",name="Fitness Value")
  
  Wavg<-function(fit,fromCoord,toCoord){
    sum=0; count=0;i=1;
    while (fit[i,2]<=fromCoord) i<-i+1
    while (fit[i,2]<toCoord){
      sum<-sum+fit[i,4]; count<-count+1;
      print(fit[i,4]);
      i<-i+1;
    }
    avg<-(sum/count);
    return(avg);
  }
  
  
  
  ht1 <- HighlightTrack(trackList = list(gTrack,sTrack,anTrack,aggTrack,fitTrack), start =fromCoord,end= toCoord,chromosome = "NC_003028")
  #plotTracks(list(gTrack,sTrack,anTrack,aggTrack,fitTrack),
             #background.title="darkred",fontsize=17,chromosome="NC_003028",
             #from=fromCoord,to=toCoord)
  plotTracks(ht1,chromosome="NC_003028",background.title="darkred",fontsize=17,from=fromCoord-15,to=toCoord+15)
  
  confirm=paste0("Tracks plotted for genomic coordinates ",fromCoord-10," to ",toCoord+10,".")
  return(confirm)
  
  
}