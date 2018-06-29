# script: differential signal and enrichment of markers of interest

# 0. Enrichment of CpG Coverage and Coverage Bias at Crypt Enhancers
{
  cggr <- granges(g3); save(cggr,file="grobj-cpg_hm450filt-platform.rda")
  
  # Super-enhancer coverage
  ssef.auto <- ssef[!seqnames(ssef) %in% c("chrX","chry")]; length(ssef.auto)
  ssef.auto <- ssef.auto[!ssef.auto$ncpg450==0]; length(ssef.auto) # 1206
  save(ssef.auto,file="ssef-info-filtered_autosomal-w-cpgcoverage.rda")
  
  # Null coverage - random ranges of similar size to super-enhancers
  mean(end(ssef)-start(ssef)) # mean = 53387.64
  {
    library(BSgenome.Hsapiens.UCSC.hg19)
    seqlengths(NFKB) = seqlengths(Hsapiens)
    
    chrlengths = seqlengths(Hsapiens)
    
    uchr <- as.character(unique(seqnames(ssef.auto)))
    freqchr <- as.data.frame(table(seqnames(ssef.auto))); freqchr <- freqchr[!freqchr$Var1 %in% c("chrX","chrY"),]
    
    nullgr.df <- as.data.frame(matrix(nrow=0,ncol=3)); colnames(nullgr.df) <- c("chr","start","end")
    
    for(i in 1:nrow(freqchr)){
      chri <- as.character(freqchr[i,1])
      lenchri <- chrlengths[names(chrlengths)==chri]
      
      for(j in 1:freqchr[i,2]){
        
        ssefj <- ssef[seqnames(ssef)==chri] # check that chr is same
        ssefj <- ssefj[sample(length(ssefj),1)]
        lenj <- end(ssefj)-start(ssefj) # get random sample length
        
        indj <- sample((as.numeric(lenchri)-lenj),1) # subtract length so only valid ranges retained
        
        nullgr.df <- rbind(nullgr.df,
                           data.frame(chr=freqchr[i,1],
                                      start=indj,
                                      end=indj+lenj))
        message(i,";",j)
      }
      message(i)
    }
    
    nullgr <- makeGRangesFromDataFrame(nullgr.df)
    length(nullgr); length(reduce(nullgr))
  }
  nullgr <- reduce(nullgr); summary(end(nullgr)-start(nullgr)) # mean = 55626
  binom.test(length(subsetByOverlaps(cggr,ssef.auto)),length(cggr),
             p=(length(subsetByOverlaps(cggr,nullgr))/length(cggr))) # p-value < 2.2e-16 => ssef is under-enriched
}

#====================================================
# 1. Enrichment of differentially methylated CpGs
#====================================================
# 1A. Matched sample test
{
  msamp.all <- colnames(betaval.ssef$`chr1:967724-1016615`$betaval.match)
  msamp.t <- msamp.all[substr(msamp.all,14,15)=="01"]
  msamp.nm <- msamp.all[substr(msamp.all,14,15)=="11"]
  
  g3mt <- g3[,g3$patient_id %in% substr(msamp.all,9,12)]
  table(g3m$hist); summary(as.data.frame(table(g3m$patient_id))) # 38 samples, each group
  
  g3mt <- g3[,g3$patient_id %in% substr(msamp.all,9,12) & g3$hist=="01"]
  g3mn <- g3[,g3$patient_id %in% substr(msamp.all,9,12) & g3$hist=="11"]
  g3mt <- g3mt[,order(match(g3mt$patient_id,g3mn$patient_id))]; identical(g3mt$patient_id,g3mn$patient_id)
  
  b.g3mt <- getBeta(g3mt); colnames(b.g3mt) <- g3mt$patient_id
  b.g3mn <- getBeta(g3mn); colnames(b.g3mn) <- g3mn$patient_id
  identical(colnames(b.g3mt),colnames(b.g3mn))
  
  tdfm <- as.data.frame(matrix(nrow=0,ncol=4)); colnames(tdf) <- c("cpg.id","t.test.p","t.test.est")
  for(i in 1:nrow(b.g3mn)){
    ti <- t.test(b.g3mt[i,],b.g3mn[i,],paired=T)
    tdfm <- rbind(tdfm,data.frame(cpg.id=rownames(b.g3mn)[i],
                                  t.test.p=as.numeric(ti$p.value),
                                  t.test.est=as.numeric(ti$estimate),
                                  stringsAsFactors = F))
    message(i)
  }
  #tdfm.1to165102 <- tdfm
  #save(tdfm.1to165102,file="tdfm_dmp_1to165102.rda")
  
  # 1A cont. Enrichment of DMPs (matched t-test) among sse.auto
  tdfm.all$padj <- p.adjust(tdfm.all$t.test.p,method="BH")
  dim(tdfm.all[tdfm.all$padj<1e-10,]) # 18481
  tdfm.dmp <- dmpm <- tdfm.all[tdfm.all$padj<1e-10,]$cpg.id
  # enrichment test - among sse's
  bndf.dmpm <- as.data.frame(matrix(nrow=0,ncol=3))
  colnames(bndf.dmpm) <- c("se.id","ratiodif","bt.pval")
  r1 <- length(subsetByOverlaps(cggr[names(cggr) %in% dmpm],ssef.auto))/length(subsetByOverlaps(cggr,ssef.auto))
  for(i in 1:length(ssef.auto)){
    ssef.autoi <- ssef.auto[i]
    r2i <- length(subsetByOverlaps(cggr[names(cggr) %in% dmpm],ssef.autoi))/length(subsetByOverlaps(cggr,ssef.autoi))
    bti <- binom.test(length(subsetByOverlaps(cggr[names(cggr) %in% dmpm],ssef.autoi)),
                      length(subsetByOverlaps(cggr,ssef.autoi)),
                      p=r1)
    
    bndf.dmpm <- rbind(bndf.dmpm,data.frame(se.id=names(ssef.autoi),
                                            ratiodif=r1-r2i,
                                            bt.pval=as.numeric(bti$p.value)))
    
    message(i)
  }
  save(bndf.dmpm,file="binomdf-dmp-matched_coad-nm-vs-t.rda")
  summary(bndf.dmpm$bt.pval)
  dmpm.en.seid <- as.character(bndf.dmpm[bndf.dmpm$bt.pval<0.01 & bndf.dmpm$ratiodif<0,]$se.id) # filter for enrichment of dmps
  save(dmpm.en.seid,file="seid-enriched01-dmp-matched_coad-nm-vs-t.rda")
}

# 1B. Cross-sectional comparison
{
  g3c <- g3[,!(g3$hist=="01" & g3$patient_id %in% g3[,g3$hist=="11"]$patient_id)]
  table(g3c$hist); summary(as.data.frame(table(g3c$patient_id)))
  #01  11 
  #255  38
  
  fdfc <- dmpFinder(getBeta(g3c),pheno=g3c$hist,type="categorical")
  save(fdfc,file="ftestdf-crxtest-dmp_coad-nm-vs-t.rda")
  
  # enrichment test - vs. background
  dmpc <- rownames(fdfc[fdfc$qval<1e-5,]); length(dmpc) # [1] 93047
  binom.test(length(subsetByOverlaps(cggr[names(cggr) %in% dmpc],ssef.auto)),
             length(subsetByOverlaps(cggr,ssef.auto)),
             p=(length(subsetByOverlaps(cggr[names(cggr) %in% dmpc],nullgr))/length(subsetByOverlaps(cggr,nullgr)))) # [1] 184129, underenriched
  binom.test(length(subsetByOverlaps(cggr[names(cggr) %in% dmpc],ssef.auto)),length(subsetByOverlaps(cggr,ssef.auto)),
             p=(length(dmpc)/length(cggr))) # p-value < 2.2e-16
  # enrichment test - among sse's
  bndf.dmpc <- as.data.frame(matrix(nrow=0,ncol=3))
  colnames(bndf.dmpc) <- c("se.id","ratiodif","bt.pval")
  r1 <- length(subsetByOverlaps(cggr[names(cggr) %in% dmpc],ssef.auto))/length(subsetByOverlaps(cggr,ssef.auto))
  
  for(i in 1:length(ssef.auto)){
    ssef.autoi <- ssef.auto[i]
    r2i <- length(subsetByOverlaps(cggr[names(cggr) %in% dmpc],ssef.autoi))/length(subsetByOverlaps(cggr,ssef.autoi))
    bti <- binom.test(length(subsetByOverlaps(cggr[names(cggr) %in% dmpc],ssef.autoi)),
                      length(subsetByOverlaps(cggr,ssef.autoi)),
                      p=r1)
    
    bndf.dmpc <- rbind(bndf.dmpc,data.frame(se.id=names(ssef.autoi),
                                            ratiodif=r1-r2i,
                                            bt.pval=as.numeric(bti$p.value)))
    
    message(i)
  }
  save(bndf.dmpc,file="binomdf-dmp-crx_coad-nm-vs-t.rda")
  summary(bndf.dmpc$bt.pval)
  
  dmpc.en.seid <- as.character(bndf.dmpc[bndf.dmpc$bt.pval<0.01 & bndf.dmpc$ratiodif<0,]$se.id) # filter for enrichment of dmps
  save(dmpc.en.seid,file="seid-enriched01-dmp-crx_coad-nm-vs-t.rda")
}

# NOTES: generate normalized smooth plots, overlaying region methylation at the se.id's of interest here...
{
  length(unique(dmpm.en.seid,dmpc.en.seid)) # 31
  length(intersect(dmpm.en.seid,dmpc.en.seid)) # 15
  
}

#========================================
# 2. Enrichment of DEGs (TAD-assoc. mRNA)
#========================================
library(methyIntegratoR)
data(gene_promoters_grlist_hg19); pgr <- genepromoters.grlist$txdb19.promotersgr

# 2.0 get TAD-assoc. gene id's
{
  
  tad.regions <- ssef.auto$tadregion; tad.regions <- tad.regions[!tad.regions==":-"]
  tadgr <- makeGRangesFromDataFrame(data.frame(chr=gsub(":.*","",tad.regions),
                                               start=gsub(".*:|-.*","",tad.regions),
                                               end=gsub(".*-","",tad.regions),
                                               stringsAsFactors = F))
  names(tadgr) <- paste0(seqnames(tadgr),":",start(tadgr),"-",end(tadgr))
  save(tadgr,file="grobj-tadregions-sseautofilt.rda")
  
  testgenes <- unique(unlist(strsplit(ssef.auto$tadgeneid,";"))) # len = 10641
}

# 2A. matched test
{
  exm.mrna <- ex[gsub("\\|.*","",rownames(ex)) %in% testgenes & !duplicated(gsub("\\|.*","",rownames(ex))),
                 substr(colnames(ex),14,15) %in% c("01","11")]
  dim(exm.mrna)
  
  nmid <- intersect(substr(colnames(exm.mrna[,substr(colnames(exm.mrna),14,15)=="11"]),9,12),
                    substr(colnames(exm.mrna[,substr(colnames(exm.mrna),14,15)=="01"]),9,12))
  exm.mrna <- exm.mrna[,substr(colnames(exm.mrna),9,12) %in% nmid]
  summary(as.data.frame(table(substr(colnames(exm.mrna),9,12))))
  table(substr(colnames(exm.mrna),14,15))
  #01 11 
  #32 32
  save(exm.mrna,file="expr-mrna-match_coad-nm-vs-t.rda")
  
  tdfm.mrna <- as.data.frame(matrix(nrow=0,ncol=5))
  colnames(tdfm.mrna) <- c("gene.id","tad.id","sse.id","xdifsamp.tn","xdifpop.tn","ttpunadj")
  
  for(i in 1:nrow(exm.mrna)){
    
    genei <- gsub("\\|.*","",rownames(exm.mrna)[i])
    exti <- as.numeric(exm.mrna[i,substr(colnames(exm.mrna),14,15)=="01"])
    exni <- as.numeric(exm.mrna[i,substr(colnames(exm.mrna),14,15)=="11"]) 
    #identical(substr(names(exti),9,12),substr(names(exni),9,12))
    
    popdif.tn <- mean(exti,na.rm=T)-mean(exni,na.rm=T)
    sampdif.tn <- mean(exti-exni,na.rm=T)
    
    tpi <- as.numeric(t.test(exti,exni,paired=T)$p.value)
    tdfm.mrna <- rbind(tdfm.mrna,data.frame(gene.id=genei,
                                              tad.id="NA",
                                              sse.id="NA",
                                              xdifsamp.tn=sampdif.tn,
                                              xdifpop.tn=popdif.tn,
                                              ttpunadj=tpi,
                                              stringsAsFactors = F))
    message(i)
    }
    
    save(tdfm.mrna,file="df-ttest-exmrnatad_coad-t-vs-n-match.rda")
}

# 2B. crx test
{
  exc.mrna <- ex[gsub("\\|.*","",rownames(ex)) %in% testgenes & !duplicated(gsub("\\|.*","",rownames(ex))),]
  nmid <- substr(colnames(exc.mrna[,substr(colnames(exc.mrna),14,15)=="11"]),9,12)
  exc.mrna <- exc.mrna[,substr(colnames(exc.mrna),14,15) %in% c("01","11") & 
                         !(substr(colnames(exc.mrna),9,12) %in% nmid & substr(colnames(exc.mrna),14,15)=="01")]
  dim(exc.mrna)
  table(substr(colnames(exc.mrna),14,15))
  #01  11 
  #347  51
  summary(as.data.frame(table(substr(colnames(exc.mrna),9,12))))
  save(exc.mrna,file="expr-mrna-crx_coad-nm-vs-t.rda")
  
  tdfc.mrna <- as.data.frame(matrix(nrow=0,ncol=5))
  colnames(tdfc.mrna) <- c("gene.id","tad.id","sse.id","xdifsamp.tn","xdifpop.tn","ttpunadj")
  
  for(i in 1:nrow(exc.mrna)){
    
    genei <- gsub("\\|.*","",rownames(exc.mrna)[i])
    exti <- as.numeric(exc.mrna[i,substr(colnames(exc.mrna),14,15)=="01"])
    exni <- as.numeric(exc.mrna[i,substr(colnames(exc.mrna),14,15)=="11"]) 
    #identical(substr(names(exti),9,12),substr(names(exni),9,12))
    
    popdif.tn <- mean(exti,na.rm=T)-mean(exni,na.rm=T)
    sampdif.tn <- mean(exti-exni,na.rm=T)
    
    tpi <- as.numeric(t.test(exti,exni,paired=F)$p.value) # not paired t-test!
    tdfc.mrna <- rbind(tdfc.mrna,data.frame(gene.id=genei,
                                            tad.id="NA",
                                            sse.id="NA",
                                            xdifsamp.tn=sampdif.tn,
                                            xdifpop.tn=popdif.tn,
                                            ttpunadj=tpi,
                                            stringsAsFactors = F))
    message(i)
  }
  
  save(tdfc.mrna,file="df-ttest-exmrnatad_coad-t-vs-n-crx.rda")
}

# 2C enrichment tests and overlapping results
tdfm.mrna$padj <- p.adjust(tdfm.mrna$ttpunadj,method="BH")
tdfc.mrna$padj <- p.adjust(tdfc.mrna$ttpunadj,method="BH")

bm <- nrow(tdfm.mrna[tdfm.mrna$padj<0.01,])/nrow(tdfm.mrna) # [1] 0.481166
bc <- nrow(tdfc.mrna[tdfc.mrna$padj<0.01,])/nrow(tdfc.mrna) # [1] 0.7098129

exdf.ssetad <- as.data.frame(matrix(nrow=0,ncol=8))
colnames(exdf.ssetad) <- c("sse.id","ntadgene","ndeg.crx","ndeg.match","frac.crx","frac.match",
                           "binom.match.punadj","binom.crx.punadj")
padj.sig <- 0.01

for(i in 1:length(ssef.auto)){
  ssei <- ssef.auto[i]
  ntadgenei <- ssef.auto[i]$ntadgenes
  
  if(ntadgenei>0){
    genesi <- unlist(strsplit(ssef.auto[i]$tadgeneid,";"))
    
    mi <- tdfm.mrna[tdfm.mrna$gene.id %in% genesi,]
    ci <- tdfc.mrna[tdfc.mrna$gene.id %in% genesi,]
    
    ndeg.crxi <- nrow(ci[ci$padj<padj.sig,])
    ndeg.matchi <- nrow(mi[mi$padj<padj.sig,])
    
    if(ndeg.crxi>0){
      bncpi <- as.numeric(binom.test(nrow(ci[ci$padj<padj.sig,]),nrow(ci),p=bc)$p.value)
    } else{
      bncpi <- "NA"
    }
    
    if(ndeg.matchi>0){
      bnmpi <- as.numeric(binom.test(nrow(mi[mi$padj<padj.sig,]),nrow(mi),p=bm)$p.value)
    } else{
      bnmpi <- "NA"
    }
    
    exdf.ssetad <- rbind(exdf.ssetad,
                         data.frame(sse.id=names(ssei),
                                    ntadgene=ntadgenei,
                                    ndeg.crx=ndeg.crxi,
                                    ndeg.match=ndeg.matchi,
                                    frac.crx=(nrow(ci[ci$padj<padj.sig,])/nrow(ci)),
                                    frac.match=(nrow(mi[mi$padj<padj.sig,])/nrow(mi)),
                                    binom.match.punadj=bnmpi,
                                    binom.crx.punadj=bncpi,
                                    stringsAsFactors = F))
  }
  message(i)
}


bn.exdf.ssetad <- exdf.ssetad
head(bn.exdf.ssetad)

bn.exdf.ssetad$bnm.padj <- p.adjust(as.numeric(bn.exdf.ssetad$binom.match.punadj),method="BH")
bn.exdf.ssetad$bnc.padj <- p.adjust(as.numeric(bn.exdf.ssetad$binom.crx.punadj),method="BH")

save(bn.exdf.ssetad,file="bn-df-degmrna_coad-n-vs-t-both.rda")

b <- bn.exdf.ssetad; dim(b)
b <- b[!(is.na(b$bnm.padj)|is.na(b$bnc.padj)),]; dim(b)
b$bnm.padj <- p.adjust(as.numeric(b$binom.match.punadj),method="BH")
b$bnc.padj <- p.adjust(as.numeric(b$binom.crx.punadj),method="BH")

pbh.sig <- 0.1
nrow(b[b$bnm.padj<pbh.sig,]) # [1] 5
nrow(b[b$bnc.padj<pbh.sig,]) # [1] 0
nrow(b[b$bnc.padj<pbh.sig|b$bnm.padj<pbh.sig,]) # [1] 5
nrow(b[b$bnc.padj<pbh.sig & b$bnm.padj<pbh.sig,]) # [1] 0

length(intersect(b[b$bnc.padj<pbh.sig & b$bnm.padj<pbh.sig,]$sse.id,
                 se.oi.mrgn.dmpen)) # [1] 0
length(intersect(b[b$bnc.padj<pbh.sig|b$bnm.padj<pbh.sig,]$sse.id,
                 se.oi.mrgn.dmpen))


# 2D. Rerun tests on log2-transformed data
{
  tdfm.mrna2 <- as.data.frame(matrix(nrow=0,ncol=5))
  colnames(tdfm.mrna) <- c("gene.id","tad.id","sse.id","xdifsamp.tn","xdifpop.tn","ttpunadj")
  exm.mrna2 <- as.matrix(exm.mrna); class(exm.mrna2) <- "numeric"
  exm.mrna2 <- log2(exm.mrna2+0.1)
  
  for(i in 1:nrow(exm.mrna)){
    
    genei <- gsub("\\|.*","",rownames(exm.mrna2)[i])
    exti <- as.numeric(exm.mrna2[i,substr(colnames(exm.mrna2),14,15)=="01"])
    exni <- as.numeric(exm.mrna2[i,substr(colnames(exm.mrna2),14,15)=="11"]) 
    #identical(substr(names(exti),9,12),substr(names(exni),9,12))
    
    popdif.tn <- mean(exti,na.rm=T)-mean(exni,na.rm=T)
    sampdif.tn <- mean(exti-exni,na.rm=T)
    
    tpi <- as.numeric(t.test(exti,exni,paired=T)$p.value)
    tdfm.mrna2 <- rbind(tdfm.mrna2,data.frame(gene.id=genei,
                                              tad.id="NA",
                                              sse.id="NA",
                                              xdifsamp.tn=sampdif.tn,
                                              xdifpop.tn=popdif.tn,
                                              ttpunadj=tpi,
                                              stringsAsFactors = F))
    message(i)
  }
  
  save(tdfm.mrna2,file="df-ttest-exmrnatad-log2dat_coad-t-vs-n-match.rda")
}

padj.sig <- 0.01
t1 <- as.data.frame(tdfm.mrna)
t2 <- tdfm.mrna2
t1$padj <- p.adjust(t1$ttpunadj,method="BH")
t2$padj <- p.adjust(t2$ttpunadj,method="BH")
nrow(t2[t2$padj<padj.sig,]) # [1] 4484
length(intersect(t1[t1$padj<padj.sig,1],
                 t2[t2$padj<padj.sig,1])) # [1] 3712/4484 = 83%
