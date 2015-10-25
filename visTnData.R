getTracks<-function(fromCoord,toCoord,x,y){
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
  insertFit$fit=insertFit$fit-1;
  insertFit$seqnames="NC_003028"
  insertData<-makeGRangesFromDataFrame(insertFit,keep.extra.columns=TRUE,ignore.strand=TRUE,seqinfo=NULL,seqnames.field="seqnames",start.field="start",starts.in.df.are.0based=FALSE)
  aggTrack<-DataTrack(insertData,chromosome="NC_003028",ylim=c(-1,1),type="histogram",
                      col.baseline="black",baseline=0,fill="darkred",
                      name="Insertions")
  #background.panel="#FFFFE0",
  
  #AnnotationTrack2: Fitness values that corresponding to DataTrack1
  output<-Wavg(fit,fromCoord,toCoord)
  avg<-round(as.numeric(output[1]),digits=4)
  count<-output[2]
  avgText<-paste0("Average= ",avg,"  |  [Fitness= ", avg+1, "]")
  avgTrack<-AnnotationTrack(start=fromCoord,end=toCoord,chromosome="NC_003028",id=avgText,showFeatureId=TRUE,fill="darkred")
  
  #Track for individual fitness values when the window is small enough
  fitTrack<-AnnotationTrack(insertFit,chromosome="NC_003028",
                            showFeatureId=TRUE,showId=FALSE,id=insertFit$fit,fill="transparent",
                            fontcolor.item="black",col.frame="transparent",
                            background.title="transparent",col="transparent",name="Fitness Value")
  
#Keep track pollution low:
  if (count>40 || (toCoord-fromCoord>10000)){
    ht1 <- HighlightTrack(trackList = list(gTrack,sTrack,aggTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "NC_003028")
  }
  else{
    ht1 <- HighlightTrack(trackList = list(gTrack,sTrack,anTrack,aggTrack,fitTrack,avgTrack), start =fromCoord,end= toCoord,chromosome = "NC_003028")
    
  }
  #plotTracks(list(gTrack,sTrack,anTrack,aggTrack,fitTrack),
             #background.title="darkred",fontsize=17,chromosome="NC_003028",
             #from=fromCoord,to=toCoord)
  
  #TO CHANGE SURROUNDING WINDOWS ALTER:
  if (missing(x))x=50
  if (missing(y))y=50

  plotTracks(ht1,chromosome="NC_003028",background.title="darkred",fontsize=17,from=fromCoord-x,to=toCoord+y)
  
  confirm=paste0("Tracks plotted for genomic coordinates ",fromCoord," to ",toCoord," with a window (",fromCoord-x,",",toCoord+y,")")
  return(confirm)
  
  
}

Wavg<-function(fit,fromCoord,toCoord){
  sum=0; count=0;i=1;
  while (fit[i,2]<=fromCoord) i<-i+1
  while (fit[i,2]<toCoord){
    sum<-sum+fit[i,4]; count<-count+1;
    i<-i+1;
  }
  avg<-(sum/count)-1;
  output<-list(avg,count)
  return(output);
}