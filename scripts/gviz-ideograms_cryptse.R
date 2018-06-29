# script: locus maps for TAD regions of se.oi

library(Gviz)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
genesens <- genes(edb)
seqlevels(genesens) <- paste0("chr",seqlevels(genesens))
genesens <- genesens[mcols(genesens)$gene_biotype=="protein_coding"]

# dmr stuff
dmrc.gr <- dmr.crx.infolist$dmr.crx.ranges
dmrm.gr <- dmr.match.infolist$dmr.match.ranges
dmrc.gr <- makeGRangesFromDataFrame(data.frame(start=gsub(".*:|-.*","",dmrc.gr),
                                    end=gsub(".*-","",dmrc.gr),
                                    chr=gsub(":.*","",dmrc.gr)))
dmrm.gr <- makeGRangesFromDataFrame(data.frame(start=gsub(".*:|-.*","",dmrm.gr),
                                    end=gsub(".*-","",dmrm.gr),
                                    chr=gsub(":.*","",dmrm.gr)))
# methy cpg
cggr <- granges(g3)

#==========================
# coordinates of interest
#==========================
for(i in 1:6){
  se.oi <- se.oi.mrgn.dmpen[i]
  coor.obj <- c(gsub(".*:|-.*","",se.oi),gsub(".*-","",se.oi),gsub(":.*","",se.oi)); names(coor.obj) <- c("start","end","chr")
  seoi.gr <- makeGRangesFromDataFrame(data.frame(start=coor.obj[1],end=coor.obj[2],chr=coor.obj[3]))
  
  #========================
  # prepare the data track
  #========================
  cgoi <- names(subsetByOverlaps(cggr,seoi.gr))
  cggr.oi <- cggr[names(cggr) %in% cgoi]
  #identical(cgoi,names(cggr.oi))
  bi <- getBeta(g3[cgoi,])
  colnames(bi) <- paste0(g3$patient_id,"_",g3$hist)
  dfcg <- data.frame(chr=seqnames(cggr.oi),start=start(cggr.oi),end=end(cggr.oi))
  rownames(dfcg) <- names(cggr.oi)
  identical(rownames(dfcg),rownames(bi))
  dfcg <- cbind(dfcg,bi)
  dfcg.gr <- makeGRangesFromDataFrame(dfcg,keep.extra.columns = T)
  dTrack <- DataTrack(dfcg.gr,groups=gsub(".*_","",colnames(bi)),
                      type=c("a","p","confint"))
   
  #===============================
  # prepare the granges objects
  #===============================
  # genome axis
  gtrack <- GenomeAxisTrack()
  # chr ideogram
  itrack <- IdeogramTrack(genome = "hg19", chromosome = as.character(coor.obj[3]))
  
  # NEW: ChIP-Seq Coverage, H3K27ac marks OL
  h3k27ac.track <- AlignmentsTrack(range=subsetByOverlaps(chgr.h3k27ac.allgr$redranges.4crypt,seoi.gr),stacking="hide")
  
  # cohen objects
  # GVEL
  c <- list.cohen2017.vels
  gvel.all <- makeGRangesFromDataFrame(data.frame(start=c$si2_gainVEL.csv$start,
                                                  end=c$si2_gainVEL.csv$end,
                                                  chr=c$si2_gainVEL.csv$chr,
                                                  stringsAsFactors = F))
  #gvel.all <- subsetByOverlaps(gvel.all,seoi.gr)
  at.gvel <- AnnotationTrack(subsetByOverlaps(gvel.all,seoi.gr),name="GVEL",col="chartreuse2",fill="chartreuse2")
  names(at.gvel) <- "GVEL"
  
  # LVEL
  lvel.all <- makeGRangesFromDataFrame(data.frame(start=c$si3_lostVEL.csv$start,
                                                  end=c$si3_lostVEL.csv$end,
                                                  chr=c$si3_lostVEL.csv$chr,
                                                  stringsAsFactors = F))
  #lvel.all <- subsetByOverlaps(lvel.all,seoi.gr)
  at.lvel <- AnnotationTrack(subsetByOverlaps(lvel.all,seoi.gr),name="LVEL",col="coral3",fill="coral3")
  names(at.lvel) <- "LVEL"
  
  # dmr data: match
  at.dmr.match <- AnnotationTrack(subsetByOverlaps(dmrm.gr,seoi.gr),name="DMR: Paired",col="purple",fill="purple")
  names(at.dmr.match) <- "DMR: Paired"
  # dmr data: crx
  at.dmr.crx <- AnnotationTrack(subsetByOverlaps(dmrc.gr,seoi.gr),name="DMR: Cross",col="purple",fill="purple")
  names(at.dmr.crx) <- "DMR: Cross"

  # gene constructs
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                      symbol=unique(subsetByOverlaps(genesens,seoi.gr)$gene_name),
                                      name = "ENSEMBL",stacking="full")
  
  # gc content object
  gcContent <- UcscTrack(genome="hg19", chromosome=as.character(coor.obj[3]),
                         track="GC Percent", table="gc5Base", from=as.numeric(coor.obj[1]), to=as.numeric(coor.obj[2]),
                         trackType="DataTrack", start="start", end="end", data="score",
                         type="hist", window=-1, windowSize=1500, fill.histogram="blue",
                         col.histogram="blue", ylim=c(0, 100), name="GC Percent")
  # cpg islands
  cpgIslands <- UcscTrack(genome = "hg19", chromosome = as.character(coor.obj[3]),
                          track = "cpgIslandExt", from = as.numeric(coor.obj[1]), to = as.numeric(coor.obj[2]),
                          trackType = "AnnotationTrack", start = "chromStart",
                          end = "chromEnd", id = "name", shape = "box",
                          fill = "darkblue", name = "CpG Islands")
  
  
  #=========================
  # 1. ideograms with Gviz
  #=========================
  
  # Super-enhancer plot
  #gviz.seplotlist <- list()
  se.plottracks.oi <- list(itrack,gtrack,dTrack,h3k27ac.track,at.dmr.match,at.dmr.crx,biomTrack,at.gvel,at.lvel,cpgIslands,gcContent)
  #names(se.plottracks.oi) <- c("gviz.tracklist","coor")
  #gviz.seplotlist[[length(gviz.seplotlist)+1]] <- se.plottracks.oi
  #names(gviz.seplotlist)[length(gviz.seplotlist)] <- paste0(se.oi)
  
  # write plot to pdf file
  #pdftitle <- paste0("se-ideo_",gsub(":","_",se.oi),".pdf")
  
  {
    #pdf(file=pdftitle,width=7,height=12)
    pdf(file="seoi1",width=7,height=12)
    plotTracks(se.plottracks.oi,geneSymbols=T,
               from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]),
               type = c("a", "p", "confint"),
               legend=T,background.panel="#FFFEDB",background.title="darkslategray4",
               main=paste0("SE ID ",se.oi))
    dev.off()
    
    #plotTracks(list(itrack,gtrack,dTrack,biomTrack,at.gvel,at.lvel,cpgIslands,gcContent),geneSymbols=T,
    #           from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]),
    #           type = c("a", "p", "confint"),
    #           legend=T,background.panel="#FFFEDB",background.title="darkslategray4",
    #           main=paste0("SE ID ",se.oi))
    
  }
  dev.off()
}

save(gviz.seplotlist,file="gviz_seplotlist_robj.rda")

#=======
# misc
#=======
# examples
{
  plotTracks(list(biomTrack),geneSymbols=T,
             from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]))
  
  plotTracks(list(at.gvel,at.lvel),geneSymbols=T,
             from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]),
             names=c("ensembl","gvel","lvel"))
  
  
  plotTracks(list(biomTrack,at.gvel,at.lvel),geneSymbols=T,
             from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]))
  
  plotTracks(list(biomTrack,at.gvel,at.lvel,gcContent,cpgIslands),geneSymbols=T,
             from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]))
}


#====================
# new ideogram stuff
#====================

# filt = 5 se.oi with overlapping dmrs
#se.oi.mrgn.dmpen <- se.oi.mrgn.dmpen[se.oi.mrgn.dmpen!="chr9:36254473-36338542"]

# regions to plot
# i = 1
for(i in 1:5){
  se.oi <- se.oi.mrgn.dmpen[i]
  gviz.img.titlei.pdf <- paste0("ideo-seoi_wh3k27ac_",gsub(":","_",se.oi),".pdf"); 
  gviz.img.titlei.jpeg <- paste0("ideo-seoi_wh3k27ac_",gsub(":","_",se.oi),".jpg")
  coor.obj <- c(gsub(".*:|-.*","",se.oi),gsub(".*-","",se.oi),gsub(":.*","",se.oi)); names(coor.obj) <- c("start","end","chr")
  seoi.gr <- makeGRangesFromDataFrame(data.frame(start=coor.obj[1],end=coor.obj[2],chr=coor.obj[3]))
  
  # 1. DATA TRACK
  {
    cgoi <- names(subsetByOverlaps(cggr,seoi.gr))
    cggr.oi <- cggr[names(cggr) %in% cgoi]
    #identical(cgoi,names(cggr.oi))
    bi <- getBeta(g3[cgoi,])
    colnames(bi) <- paste0(g3$patient_id,"_",g3$hist)
    dfcg <- data.frame(chr=seqnames(cggr.oi),start=start(cggr.oi),end=end(cggr.oi))
    rownames(dfcg) <- names(cggr.oi)
    identical(rownames(dfcg),rownames(bi))
    dfcg <- cbind(dfcg,bi)
    dfcg.gr <- makeGRangesFromDataFrame(dfcg,keep.extra.columns = T)
    dTrack <- DataTrack(dfcg.gr,groups=gsub(".*_","",colnames(bi)),
                        type=c("a","p","confint"))
  }
  
  # genome axis
  gtrack <- GenomeAxisTrack()
  # chr ideogram
  itrack <- IdeogramTrack(genome = "hg19", chromosome = as.character(coor.obj[3]))
  
  # NEW: ChIP-Seq Coverage, H3K27ac marks OL
  h3k27ac.track <- AlignmentsTrack(range=subsetByOverlaps(chgr.h3k27ac.allgr$redranges.4crypt,seoi.gr),name="H3K27ac\ncoverage",
                                   type = "coverage")
  #h3k27ac.track <- DataTrack(range=subsetByOverlaps(chgr.h3k27ac.allgr$redranges.4crypt,seoi.gr),
  #                           type="l",name="H3K27ac coverage")
  #h3k27ac.track <- dTrack4 <- DataTrack(range = lfi[1], genome = "hg19",
  #                                      type = "l", name = "Coverage", window = -1, 
  #                                      chromosome = "chr1")
  
  # cohen objects
  {
    # GVEL
    c <- list.cohen2017.vels
    gvel.all <- makeGRangesFromDataFrame(data.frame(start=c$si2_gainVEL.csv$start,
                                                    end=c$si2_gainVEL.csv$end,
                                                    chr=c$si2_gainVEL.csv$chr,
                                                    stringsAsFactors = F))
    #gvel.all <- subsetByOverlaps(gvel.all,seoi.gr)
    at.gvel <- AnnotationTrack(subsetByOverlaps(gvel.all,seoi.gr),name="GVEL",col="chartreuse2",fill="chartreuse2")
    names(at.gvel) <- "GVEL"
    
    # LVEL
    lvel.all <- makeGRangesFromDataFrame(data.frame(start=c$si3_lostVEL.csv$start,
                                                    end=c$si3_lostVEL.csv$end,
                                                    chr=c$si3_lostVEL.csv$chr,
                                                    stringsAsFactors = F))
    #lvel.all <- subsetByOverlaps(lvel.all,seoi.gr)
    at.lvel <- AnnotationTrack(subsetByOverlaps(lvel.all,seoi.gr),name="LVEL",col="coral3",fill="coral3")
    names(at.lvel) <- "LVEL"
  }
  
  
  # dmr data: match
  at.dmr.match <- AnnotationTrack(subsetByOverlaps(dmrm.gr,seoi.gr),name="DMR: Paired",col="purple",fill="purple")
  names(at.dmr.match) <- "DMR: Paired"
  # dmr data: crx
  at.dmr.crx <- AnnotationTrack(subsetByOverlaps(dmrc.gr,seoi.gr),name="DMR: Cross",col="purple",fill="purple")
  names(at.dmr.crx) <- "DMR: Cross"
  
  # gene constructs
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                      symbol=unique(subsetByOverlaps(genesens,seoi.gr)$gene_name),
                                      name = "ENSEMBL",stacking="full")
  
  # gc content object
  gcContent <- UcscTrack(genome="hg19", chromosome=as.character(coor.obj[3]),
                         track="GC Percent", table="gc5Base", from=as.numeric(coor.obj[1]), to=as.numeric(coor.obj[2]),
                         trackType="DataTrack", start="start", end="end", data="score",
                         type="hist", window=-1, windowSize=1500, fill.histogram="blue",
                         col.histogram="blue", ylim=c(0, 100), name="GC Percent")
  # cpg islands
  cpgIslands <- UcscTrack(genome = "hg19", chromosome = as.character(coor.obj[3]),
                          track = "cpgIslandExt", from = as.numeric(coor.obj[1]), to = as.numeric(coor.obj[2]),
                          trackType = "AnnotationTrack", start = "chromStart",
                          end = "chromEnd", id = "name", shape = "box",
                          fill = "darkblue", name = "CpG Islands")
  
  # Super-enhancer plot
  se.plottracks.oi <- list(itrack,gtrack,dTrack,at.dmr.match,at.dmr.crx,biomTrack,h3k27ac.track,at.gvel,at.lvel,cpgIslands,gcContent)
  
  jpeg(file=gviz.img.titlei.jpeg,width=6,height=12,units="in",res=500)
  {
    plotTracks(se.plottracks.oi,geneSymbols=T,
               from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]),
               type = c("a", "p", "confint","coverage"),
               legend=T,background.panel="#FFFEDB",background.title="darkslategray4",
               main=paste0("SE ID ",se.oi))
    dev.off()
  }
  
  pdf(file=gviz.img.titlei.pdf,width=7,height=12)
  {
    plotTracks(se.plottracks.oi,geneSymbols=T,
               from=as.numeric(coor.obj[1]),to=as.numeric(coor.obj[2]),
               type = c("a", "p", "confint","coverage"),
               legend=T,background.panel="#FFFEDB",background.title="darkslategray4",
               main=paste0("SE ID ",se.oi))
    dev.off()
  }
  
}



