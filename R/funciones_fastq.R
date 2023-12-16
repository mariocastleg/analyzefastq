#' Read a FASTQ file
#' 
#'This function reads a FASTQQ file and returns the sequences in it and their ASCII scores in a dataframe
#'
#' @param file Input FASTQ file
#' @return A dataframe with 2 columns: sequences and scores
#'@export
read_fastq<-function(file){   
  archivo<-readLines(file)
  sequences<-archivo[seq(2,length(archivo),4)]
  score<-scores<-archivo[seq(4,length(archivo),4)]
  return(data.frame(sequences,score))
}

#' Length of reads
#' 
#'This function calculates the length of the reads in the FASTQ file
#'
#' @param sequences_vector Vector with all the sequences in the FASTQ file
#' @return Prints a summary of the values and returns the mean
#'@export
reads_length<-function(sequences_vector){
  longitudes<-sapply(sequences_vector,FUN=function(x)nchar(x))
  print(summary(longitudes))
  return(mean(longitudes))
}

#' Number of G and C
#' 
#'This function calculates the percentage of GC in a given sequence
#'
#' @param secuencia Input DNA sequence
#' @return GC percentage value
#' @export
count_GC<-function(secuencia){
  finds <- sum(strsplit(secuencia,'')[[1]]%in% c('C','G'))
  percent<-finds*100/nchar(secuencia)
  return(percent)
}
#' Calculation of GC content for a set of sequences
#' 
#'This function calculates GC content for a set of sequences in a vector
#'
#' @param sequences_vector Vector with all the sequences in the FASTQ file
#' @return A vector with all the GC percentage values ready for plotting
#' @export
GC_content<-function(sequences_vector){ 
  gc_content<-sapply(sequences_vector,FUN=count_GC) 
  return(gc_content)  
}

#' Quality Scores from ASCII values
#' 
#'This function calculates quality scores given the sequences with ASCII characters
#'
#' @param sequences_vector Vector with all the ASCII scores for the sequences in the FASTQ file
#' @param read_length Length of reads (100 by default)
#' @return A dataframe with all the quality scores for each position
#' @export
scoring<-function(ascii_vector,read_length=100){
  calidades<-sapply(ascii_vector,FUN=function(x)utf8ToInt(x))
  quality_df <- as.data.frame(t(calidades))
  colnames(quality_df)<-as.factor(1:read_length)  
  return(quality_df)
}

#' First quartile for Quality scores dataframe
#' 
#'This function calculates the first quartile of the quality scores in each position of the readsGC content for a set of sequences in a vector
#'
#' @param quality_df Dataframe with all the quality scores; positions as the column names, values in names
#' @return A vector with all the first quartile values for the scores in each position
#' @export
cuartil<-function(quality_df){  
  cuartiles<-apply(quality_df,MARGIN=2,FUN=function(x)quantile(x,probs=0.25))
  return(cuartiles)
}
