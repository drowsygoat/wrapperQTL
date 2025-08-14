
# Will's funtion, 2b documented

peak_range<-function(
                     qtl_data, # overall_results_table                     
                     METHOD="DistanceFromPeak", #Method used to obtain ranges: DistanceFromPeak, AdjacentPoints.
                     SIGTHRESHOLD=0.01, #P-value significance threshold. Stop searching for peaks if the p-value of the next peak is above SIGTHRESHOLD.
                     SEARCHDISTANCE=3e+6, #Maximum distance from local peak (if DistanceFromPeak method), or maximum distance between adjacent significant points (if AdjacentPoints method).
                     LOCALPEAKTHRESHOLD=0.01, #Relative threshold of local peak (0.01 represents 2 orders of magnitude).
                     MINGAPSIZE=0){ #Minimum size of gaps between ranges. Adjacent ranges with a gap smaller than or equal to MINGAPSIZE will be merged together.
  
  #1. Find peak position (lowest p-value).
  #2. Set local peak p-value threshold (e.g., 0.01 = two orders of magnitude below local peak) to limit local peak height.
  #3. DistanceFromPeak: Find the farthest significant points (below local peak p-value threshold) within the defined 
  #distance from local peak. This method limits the extension of the region, without limiting the gap between adjacent significant points..
  #3. AdjacentPoints: Find the farthest significant points (below local peak p-value threshold) upstream/downstream 
  #the last picked point upstream/downstream the local peak. This method limits the gap between adjacent significant points, not 
  #the extension of the total regions.
  #4. Save peak and range.
  #5. Remove range from data.
  #6. Go back to step 1.
  
  library(tidyverse)
  library(data.table)

  CHRVECTOR <- qtl_data$chr_snp #Vector of chromosomes.
  POSVECTOR <- qtl_data$location #Vector of chromosomal positions.
  PVALUEVECTOR <- qtl_data$FDR #Vector of p-values (FDR).


  if(missing(CHRVECTOR)){
    CHRVECTOR<-rep("ChromosomeC",length(POSVECTOR))
  }
  
  INPUTDATA<-data.frame(Chromosome=CHRVECTOR,
                        Position=POSVECTOR,
                        Pvalue=PVALUEVECTOR) 
  
  #Create output data frame.
  OUTPUT<-data.frame(Chromosome=NA,
                     PeakPosition=NA,
                     PeakPvalue=NA,
                     RegionStart=NA,
                     RegionEnd=NA,
                     NPeaks=NA,
                     NSignificant=NA)
  
  for(C in unique(INPUTDATA$Chromosome)){
    
    #Get data for chromosome C.
    DATA_CHRC<-INPUTDATA %>% 
      filter(Chromosome==C) %>% 
      arrange(Position)
    
    #Create output data frame.
    OUTPUT_CHRC<-data.frame(Chromosome=NA,
                            PeakPosition=NA,
                            PeakPvalue=NA,
                            RegionStart=NA,
                            RegionEnd=NA)
    
    #Continue searching while there are significant positions that haven't been assigned to a range.
    while(nrow(DATA_CHRC)>0){ 
      #Find peak value.
      PEAKVALUE<-min(DATA_CHRC$Pvalue) 
      #No meaningful peaks left.
      if(PEAKVALUE>SIGTHRESHOLD){break} 
      #Local peak dynamic threshold.
      X<-PEAKVALUE/LOCALPEAKTHRESHOLD 
      #if(X>SIGTHRESHOLD){X<-SIGTHRESHOLD}
      #Find all peak indices.
      PEAKINDICES<-which(DATA_CHRC$Pvalue==PEAKVALUE) 
      
      if(METHOD=="DistanceFromPeak"){
        #Find ranges around each peak.
        RANGES_CHRC<-lapply(PEAKINDICES,function(idx){
          PEAKPOSITION<-DATA_CHRC$Position[idx]
          PEAKPVALUE<-DATA_CHRC$Pvalue[idx]
          DATASUBSET<-DATA_CHRC %>%
            filter(Position>=(PEAKPOSITION-SEARCHDISTANCE),
                   Position<=(PEAKPOSITION+SEARCHDISTANCE))
          #Search upstream peak.
          LEFTLIMIT<-min(DATASUBSET$Position[DATASUBSET$Pvalue<=X])
          #Search downstream peak.
          RIGHTLIMIT<-max(DATASUBSET$Position[DATASUBSET$Pvalue<=X])
          return(c(PEAKPOSITION,PEAKPVALUE,LEFTLIMIT,RIGHTLIMIT))
        })
        RANGES_CHRC<-do.call(rbind,RANGES_CHRC)
        RANGES_CHRC<-as.data.frame(RANGES_CHRC) %>%
          mutate(Chromosome=C)
        colnames(RANGES_CHRC)<-c("PeakPosition","PeakPvalue","RegionStart","RegionEnd","Chromosome")
      }
      
      if(METHOD=="AdjacentPoints"){
        #Find ranges around each peak.
        RANGES_CHRC<-lapply(PEAKINDICES,function(idx){
          PEAKPOSITION<-DATA_CHRC$Position[idx]
          PEAKPVALUE<-DATA_CHRC$Pvalue[idx]
          DATASUBSET<-DATA_CHRC %>%
            filter(Position>=(PEAKPOSITION-SEARCHDISTANCE),
                   Position<=(PEAKPOSITION+SEARCHDISTANCE))
          #Search upstream peak.
          CONTINUE<-TRUE
          LEFTLIMIT<-PEAKPOSITION
          while(CONTINUE==TRUE){
            DATASUBSET<-DATA_CHRC %>%
              filter(Position>=(LEFTLIMIT-SEARCHDISTANCE),
                     Position<LEFTLIMIT)
            if(nrow(DATASUBSET)>0){
              if(sum(DATASUBSET$Pvalue<=X)>0){
                LEFTLIMIT<-min(DATASUBSET$Position[DATASUBSET$Pvalue<=X])
              }else{
                CONTINUE<-FALSE
              }
            }else{
              CONTINUE<-FALSE
            }
          }
          #Search downstream peak.
          CONTINUE<-TRUE
          RIGHTLIMIT<-PEAKPOSITION
          while(CONTINUE==TRUE){
            DATASUBSET<-DATA_CHRC %>%
              filter(Position>RIGHTLIMIT,
                     Position<=(RIGHTLIMIT+SEARCHDISTANCE))
            if(nrow(DATASUBSET)>0){
              if(sum(DATASUBSET$Pvalue<=X)>0){
                RIGHTLIMIT<-max(DATASUBSET$Position[DATASUBSET$Pvalue<=X])
              }else{
                CONTINUE<-FALSE
              }
            }else{
              CONTINUE<-FALSE
            }
          }
          return(c(PEAKPOSITION,PEAKPVALUE,LEFTLIMIT,RIGHTLIMIT))
        })
        RANGES_CHRC<-do.call(rbind,RANGES_CHRC)
        RANGES_CHRC<-as.data.frame(RANGES_CHRC) %>%
          mutate(Chromosome=C)
        colnames(RANGES_CHRC)<-c("PeakPosition","PeakPvalue","RegionStart","RegionEnd","Chromosome")
      }
      
      #Save ranges.
      OUTPUT_CHRC<-rbind(OUTPUT_CHRC,RANGES_CHRC)
      #Remove ranges from the data.
      for(i in 1:nrow(RANGES_CHRC)){
        DATA_CHRC<-DATA_CHRC %>%
          filter(!(Position>=RANGES_CHRC$RegionStart[i] & Position<=RANGES_CHRC$RegionEnd[i]))
      }
    }
    
    #Merge overlapping regions.
    #If several peaks with the same p-value exist within a range, repeat range for each peak.
    if(nrow(OUTPUT_CHRC[!is.na(OUTPUT_CHRC$PeakPosition),])>0){
      OUTPUT_CHRC<-OUTPUT_CHRC %>% 
        filter(!is.na(Chromosome)) %>%
        arrange(RegionStart) %>% 
        group_by(Index=cumsum(cummax(lag(RegionEnd+MINGAPSIZE/2,default=data.table::first(RegionEnd+MINGAPSIZE/2)))<RegionStart-MINGAPSIZE/2)) %>% 
        dplyr::summarise(RegionStart=min(RegionStart),
                         RegionEnd=max(RegionEnd),
                         PeakPosition=PeakPosition[PeakPvalue==min(PeakPvalue)],
                         PeakPvalue=min(PeakPvalue),
                         NPeaks=n()) %>%
        mutate(Chromosome=C) %>%
        as.data.frame
      OUTPUT_CHRC<-OUTPUT_CHRC[,c("Chromosome","PeakPosition","PeakPvalue","RegionStart","RegionEnd","NPeaks")]
      
      #Get number of significant points per region.
      DATA_CHRC<-INPUTDATA %>% 
        filter(Chromosome==C,
               Pvalue<=SIGTHRESHOLD)
      OUTPUT_CHRC$NSignificant<-NA
      for(i in 1:nrow(OUTPUT_CHRC)){
        OUTPUT_CHRC$NSignificant[i]<-DATA_CHRC %>%
          filter(Position>=OUTPUT_CHRC$RegionStart[i],
                 Position<=OUTPUT_CHRC$RegionEnd[i]) %>%
          nrow()
      }
      
      #Save ranges.
      OUTPUT<-rbind(OUTPUT,OUTPUT_CHRC)
    }
    
  }
  
  #Remove NA entries.
  OUTPUT<-OUTPUT %>%
    filter(!is.na(PeakPosition)) %>%
    arrange(Chromosome,PeakPvalue)
  OUTPUT$RegionLength<-OUTPUT$RegionEnd-OUTPUT$RegionStart+1
  
  return(OUTPUT)
  
}

# #Test data.
# TESTDATA<-data.frame(Chromosome=rep(1,100),
#                      Position=sample(1:1000,100,replace=FALSE),
#                      Pvalue=runif(100,min=1e-10,max=1e-4))
# 
# p<-ggplot(TESTDATA)+
#   geom_point(aes(x=Position,y=-log10(Pvalue)),color="blue",size=0.3)+
#   geom_hline(yintercept=-log10(0.01),linetype="dashed")+
#   labs(x="Position",y=expression("-Log"[10]~"Pvalue"))+
#   theme_bw()
# p
# 
# RANGES1<-peak_range(CHRVECTOR=TESTDATA$Chromosome,
#                    POSVECTOR=TESTDATA$Position,
#                    PVALUEVECTOR=TESTDATA$Pvalue,
#                    METHOD="DistanceFromPeak",
#                    SIGTHRESHOLD=0.01,
#                    SEARCHDISTANCE=20,
#                    LOCALPEAKTHRESHOLD=0.1,
#                    MINGAPSIZE=0)
# 
# RANGES2<-peak_range(CHRVECTOR=TESTDATA$Chromosome,
#                     POSVECTOR=TESTDATA$Position,
#                     PVALUEVECTOR=TESTDATA$Pvalue,
#                     METHOD="DistanceFromPeak",
#                     SIGTHRESHOLD=0.01,
#                     SEARCHDISTANCE=20,
#                     LOCALPEAKTHRESHOLD=0.1,
#                     MINGAPSIZE=10)
# 
# RANGES3<-peak_range(CHRVECTOR=TESTDATA$Chromosome,
#                    POSVECTOR=TESTDATA$Position,
#                    PVALUEVECTOR=TESTDATA$Pvalue,
#                    METHOD="AdjacentPoints",
#                    SIGTHRESHOLD=0.01,
#                    SEARCHDISTANCE=20,
#                    LOCALPEAKTHRESHOLD=0.1,
#                    MINGAPSIZE=0)
# 
# RANGES4<-peak_range(CHRVECTOR=TESTDATA$Chromosome,
#                     POSVECTOR=TESTDATA$Position,
#                     PVALUEVECTOR=TESTDATA$Pvalue,
#                     METHOD="AdjacentPoints",
#                     SIGTHRESHOLD=0.01,
#                     SEARCHDISTANCE=20,
#                     LOCALPEAKTHRESHOLD=0.1,
#                     MINGAPSIZE=10)
# 
# RANGES1
# RANGES2
# RANGES3
# RANGES4
# 
# p1<-ggplot(TESTDATA)+
#   geom_point(aes(x=Position,y=-log10(Pvalue)),color="blue",size=0.3)+
#   geom_rect(data=RANGES1[RANGES1$RegionLength>1,],aes(xmin=RegionStart,xmax=RegionEnd,ymin=-Inf,ymax=Inf),
#             fill="gray",
#             alpha=0.3,
#             inherit.aes=FALSE)+
#   geom_point(data=RANGES1[RANGES1$RegionLength>1,],aes(x=PeakPosition,y=-log10(PeakPvalue)),color="red",size=0.5)+
#   geom_hline(yintercept=-log10(0.01),linetype="dashed")+
#   labs(x="Position",y=expression("-Log"[10]~"PvalueLRT"))+
#   theme_bw()
# p1
# 
# p2<-ggplot(TESTDATA)+
#   geom_point(aes(x=Position,y=-log10(Pvalue)),color="blue",size=0.3)+
#   geom_rect(data=RANGES2[RANGES2$RegionLength>1,],aes(xmin=RegionStart,xmax=RegionEnd,ymin=-Inf,ymax=Inf),
#             fill="gray",
#             alpha=0.3,
#             inherit.aes=FALSE)+
#   geom_point(data=RANGES2[RANGES2$RegionLength>1,],aes(x=PeakPosition,y=-log10(PeakPvalue)),color="red",size=0.5)+
#   geom_hline(yintercept=-log10(0.01),linetype="dashed")+
#   labs(x="Position",y=expression("-Log"[10]~"PvalueLRT"))+
#   theme_bw()
# p2
# 
# p3<-ggplot(TESTDATA)+
#   geom_point(aes(x=Position,y=-log10(Pvalue)),color="blue",size=0.3)+
#   geom_rect(data=RANGES3[RANGES3$RegionLength>1,],aes(xmin=RegionStart,xmax=RegionEnd,ymin=-Inf,ymax=Inf),
#             fill="gray",
#             alpha=0.3,
#             inherit.aes=FALSE)+
#   geom_point(data=RANGES3[RANGES3$RegionLength>1,],aes(x=PeakPosition,y=-log10(PeakPvalue)),color="red",size=0.5)+
#   geom_hline(yintercept=-log10(0.01),linetype="dashed")+
#   labs(x="Position",y=expression("-Log"[10]~"PvalueLRT"))+
#   theme_bw()
# p3
# 
# p4<-ggplot(TESTDATA)+
#   geom_point(aes(x=Position,y=-log10(Pvalue)),color="blue",size=0.3)+
#   geom_rect(data=RANGES4[RANGES4$RegionLength>1,],aes(xmin=RegionStart,xmax=RegionEnd,ymin=-Inf,ymax=Inf),
#             fill="gray",
#             alpha=0.3,
#             inherit.aes=FALSE)+
#   geom_point(data=RANGES4[RANGES4$RegionLength>1,],aes(x=PeakPosition,y=-log10(PeakPvalue)),color="red",size=0.5)+
#   geom_hline(yintercept=-log10(0.01),linetype="dashed")+
#   labs(x="Position",y=expression("-Log"[10]~"PvalueLRT"))+
#   theme_bw()
# p4