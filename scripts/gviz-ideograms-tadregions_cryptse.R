# script: locus maps for TAD regions of se.oi
library(Gviz)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
genesens <- genes(edb)
seqlevels(genesens) <- paste0("chr",seqlevels(genesens))

#==========================
# coordinates of interest
#==========================
for(i in 1:6){
  se.oi <- se.oi.mrgn.dmpen[5]
  tad.oi <- unlist(strsplit(ssef.auto[names(ssef.auto)==se.oi]$tadregion,";"))
  # note: be sure list is length 1, if >1 then form a loop
  coor.obj <- c(gsub(".*:|-.*","",tad.oi),gsub(".*-","",tad.oi),gsub(":.*","",tad.oi)); 
  names(coor.obj) <- c("start","end","chr")
  tadoi.gr <- makeGRangesFromDataFrame(data.frame(start=coor.obj[1],end=coor.obj[2],chr=coor.obj[3]))
  #===============================
  # prepare the granges objects
  #===============================
  # gene coordinates and ideograms
  {
    # genome axis
    gtrack <- GenomeAxisTrack()
    
    # chr ideogram
    itrack <- IdeogramTrack(genome = "hg19", chromosome = as.character(coor.obj[3]))
    
  }
  
  # NEW: ChIP-Seq Coverage, H3K27ac marks OL
  h3k27ac.track <- AlignmentsTrack(range=subsetByOverlaps(chgr.h3k27ac.allgr$redranges.4crypt,tadoi.gr),name="H3K27ac coverage",
                                   type = "coverage")
  
  # cohen2017 objects
  {
    # GVEL
    c <- list.cohen2017.vels
    gvel.all <- makeGRangesFromDataFrame(data.frame(start=c$si2_gainVEL.csv$start,
                                                    end=c$si2_gainVEL.csv$end,
                                                    chr=c$si2_gainVEL.csv$chr,
                                                    stringsAsFactors = F))
    #gvel.all <- subsetByOverlaps(gvel.all,seoi.gr)
    at.gvel <- AnnotationTrack(reduce(subsetByOverlaps(gvel.all,tadoi.gr)),name="GVEL",col="chartreuse2",fill="chartreuse2")
    names(at.gvel) <- "GVEL"
    
    # LVEL
    lvel.all <- makeGRangesFromDataFrame(data.frame(start=c$si3_lostVEL.csv$start,
                                                    end=c$si3_lostVEL.csv$end,
                                                    chr=c$si3_lostVEL.csv$chr,
                                                    stringsAsFactors = F))
    #lvel.all <- subsetByOverlaps(lvel.all,seoi.gr)
    at.lvel <- AnnotationTrack(reduce(subsetByOverlaps(lvel.all,tadoi.gr)),name="LVEL",col="coral3",fill="coral3")
    names(at.lvel) <- "LVEL"
    
  }
  
  # gene constructs
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                      symbol=unique(subsetByOverlaps(genesens,tadoi.gr)$gene_name),
                                      name = "ENSEMBL",stacking="squish")
  
  # cpg islands
  cpgIslands <- UcscTrack(genome = "hg19", chromosome = as.character(coor.obj[3]),
                          track = "cpgIslandExt", from = as.numeric(coor.obj[1]), to = as.numeric(coor.obj[2]),
                          trackType = "AnnotationTrack", start = "chromStart",
                          end = "chromEnd", id = "name", shape = "box",
                          fill = "darkblue", name = "CpG Islands")
  
  # modular regions:
  # for TAD regions show the ranges of overlapping super-enhancers
  at.se<- AnnotationTrack(subsetByOverlaps(ssef.auto,tadoi.gr),name="SE's",col="orange",fill="orange")
  names(at.se) <- "SE's"
  
  
  #=========================
  # 1. ideograms with Gviz
  #=========================
  
  # Super-enhancer plot
  #gviz.seplotlist <- list()
 
  # write plot to pdf file
  pdftitle <- paste0("tad-ideo_",gsub(":","_",tad.oi),".pdf")
  {
    pdf(file=pdftitle,width=7,height=9)
    plotTracks(list(itrack,gtrack,at.se,biomTrack,h3k27ac.track,at.gvel,at.lvel,cpgIslands),geneSymbols=T,
               type="coverage",
               from = as.numeric(coor.obj[1]), to = as.numeric(coor.obj[2]),
               background.panel="#FFFEDB",background.title="darkslategray4",
               main=paste0("SE ID ",se.oi,"\nTAD ID ",tad.oi))
    
    #plotTracks(list(itrack,gtrack,dTrack,biomTrack,at.gvel,at.lvel,cpgIslands,gcContent),geneSymbols=T,
    #           from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]),
    #           type = c("a", "p", "confint"),
    #           legend=T,background.panel="#FFFEDB",background.title="darkslategray4",
    #           main=paste0("SE ID ",se.oi))
    
  }
  dev.off()
}

#=======
# misc
#=======
