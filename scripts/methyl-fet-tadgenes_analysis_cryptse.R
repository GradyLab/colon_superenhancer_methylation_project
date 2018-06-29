# script: FET analysis of TAD region cg classes

# prepare the methylation patient data
{
  nmid <- g3[,g3$hist=="11"]$patient_id
  gc <- g3[,!(g3$patient_id %in% nmid & g3$hist=="01")]
  summary(as.data.frame(table(gc$patient_id)))
  
  gm <- g3[,g3$patient_id %in% nmid]
  summary(as.data.frame(table(gm$patient_id)))
  table(gm$hist)
}

#================================================================
# Prepare methylation data at gene promoters in TAD regions
#================================================================
{
  library(methyIntegratoR)
  data(gene_promoters_grlist_hg19); pgr <- genepromoters.grlist$txdb19.promotersgr
  length(pgr)
  pgr <- pgr[!is.na(pgr$symbol)]; length(pgr)
}

# 0. get the tad gene ids
tgid <- unique(unlist(strsplit(ssef.auto$tadgeneid,";")))

# 1. get the tad gene cpgs, and store in a list
{
  cggr <- granges(g3)
  tgcg.list <- list()
  for(i in 1:length(tgid)){
    tgcg.list[[i]] <- names(subsetByOverlaps(cggr,pgr[pgr$symbol %in% tgid[i]]))
    names(tgcg.list)[i] <- tgid[i]
    message(i)
  }
}

# 2. get the methylation scores for each promoter region
{
  tg.score.list <- list()
  for(i in 1:length(tgcg.list)){
    
    genei <- names(tgcg.list)[i]
    tgi <- tgcg.list[[i]]
    
    if(length(tgi)>=3){
      btc <- getBeta(gc[tgi,gc$hist=="01"])
      bnc <- getBeta(gc[tgi,gc$hist=="11"])
      
      btm <- getBeta(gm[tgi,gm$hist=="01"])
      bnm <- getBeta(gm[tgi,gm$hist=="11"])
      
      crx.score <- list(nrow(btc), length(rownames(btc)[rowMeans(btc)>rowMeans(bnc)]))
      names(crx.score) <- c("ncg","nhtc")
      
      match.score <- list(nrow(btm),length(rownames(btm)[rowMeans(btm)>rowMeans(bnm)]))
      names(match.score) <- c("ncg","nhtm")
      
      tg.score.list[[i]] <- list(crx.score,match.score); names(tg.score.list[[i]]) <- c("crx","match")
      names(tg.score.list)[i] <- genei
    }
    
    message(i)
  }
}

# 3. get background comparison of all SE scores and their assoc. TAD promoter scores
{
  dftad.seauto <- as.data.frame(matrix(nrow=0,ncol=12))
  colnames(dftad.seauto) <- c("se.id",
                              "gene.id",
                              "se.m.score",
                              "gene.m.prscore",
                              "se.m.perc",
                              "gene.m.prperc",
                              "fet.m.pval",
                              "se.c.score",
                              "gene.c.prscore",
                              "se.c.perc",
                              "gene.c.prperc",
                              "fet.c.pval")
  
  for(i in 1:length(ssef.auto)){
    segri <- ssef.auto[i]
    
    genelisti <- unique(unlist(strsplit(segri$tadgeneid,";")))
    
    if(length(genelisti)>0){
      for(j in 1:length(genelisti)){
        
        if(genelisti[j] %in% names(tg.score.list)){
          # se methyl
          secgi <- names(subsetByOverlaps(cggr,segri))
          
          bse.tc <- getBeta(gc[secgi,gc$hist=="01"]); rtc <- rownames(bse.tc)
          bse.nc <- getBeta(gc[secgi,gc$hist=="11"])
          bse.tm <- getBeta(gm[secgi,gm$hist=="01"]); rtm <- rownames(bse.tm)
          bse.nm <- getBeta(gm[secgi,gm$hist=="11"])
          
          sec.ncg <- length(rownames(bse.tc))
          sec.tcg <- length(rtc[rowMeans(bse.tc)>rowMeans(bse.nc)])
          sem.ncg <- length(rownames(bse.tm))
          sem.tcg <- length(rtm[rowMeans(bse.tm)>rowMeans(bse.nm)])
          
          
          # pr methyl
          prc.ncg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$crx$ncg
          prc.tcg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$crx$nhtc
          prm.ncg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$match$ncg
          prm.tcg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$match$nhtm
          
          fetc <- fisher.test(matrix(c((sec.ncg-sec.tcg),sec.tcg,
                                       (prc.ncg-prc.tcg),prc.tcg),nrow=2))
          
          fetm <- fisher.test(matrix(c((sem.ncg-sem.tcg),sem.tcg,
                                       (prm.ncg-prm.tcg),prm.tcg),nrow=2))
          
          dftad.seauto <- rbind(dftad.seauto,data.frame(se.id=names(segri),
                                                        gene.id=genelisti[j],
                                                        se.m.score=sem.tcg,
                                                        gene.m.prscore=prm.tcg,
                                                        se.m.perc=sem.tcg/sem.ncg,
                                                        gene.m.prperc=prm.tcg/prm.ncg,
                                                        fet.m.pval=as.numeric(fetm$p.value),
                                                        se.c.score=sec.tcg,
                                                        gene.c.prscore=prc.tcg,
                                                        se.c.perc=sec.tcg/sec.ncg,
                                                        gene.c.prperc=prc.tcg/prc.ncg,
                                                        fet.c.pval=as.numeric(fetc$p.value),
                                                        stringsAsFactors=F))
        }
        
        message(i," : ",j)
      }
    }
    message(i)
  }
}


# 4. get comparisons of SE's of interest and their TAD promoter scores
{
  dftad.seauto.test <- as.data.frame(matrix(nrow=0,ncol=12))
  colnames(dftad.seauto.test) <- c("se.id",
                                   "gene.id",
                                   "se.m.score",
                                   "gene.m.prscore",
                                   "se.m.perc",
                                   "gene.m.prperc",
                                   "fet.m.pval",
                                   "se.c.score",
                                   "gene.c.prscore",
                                   "se.c.perc",
                                   "gene.c.prperc",
                                   "fet.c.pval")
  
  seoigr <- ssef.auto[names(ssef.auto) %in% se.oi.mrgn.dmpen]
  
  for(i in 1:length(seoigr)){
    segri <- seoigr[i]
    
    genelisti <- unique(unlist(strsplit(segri$tadgeneid,";")))
    
    if(length(genelisti)>0){
      for(j in 1:length(genelisti)){
        
        if(genelisti[j] %in% names(tg.score.list)){
          # se methyl
          secgi <- names(subsetByOverlaps(cggr,segri))
          
          bse.tc <- getBeta(gc[secgi,gc$hist=="01"]); rtc <- rownames(bse.tc)
          bse.nc <- getBeta(gc[secgi,gc$hist=="11"])
          bse.tm <- getBeta(gm[secgi,gm$hist=="01"]); rtm <- rownames(bse.tm)
          bse.nm <- getBeta(gm[secgi,gm$hist=="11"])
          
          sec.ncg <- length(rownames(bse.tc))
          sec.tcg <- length(rtc[rowMeans(bse.tc)>rowMeans(bse.nc)])
          sem.ncg <- length(rownames(bse.tm))
          sem.tcg <- length(rtm[rowMeans(bse.tm)>rowMeans(bse.nm)])
          
          
          # pr methyl
          prc.ncg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$crx$ncg
          prc.tcg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$crx$nhtc
          prm.ncg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$match$ncg
          prm.tcg <- tg.score.list[[which(names(tg.score.list)==genelisti[j])]]$match$nhtm
          
          fetc <- fisher.test(matrix(c((sec.ncg-sec.tcg),sec.tcg,
                                       (prc.ncg-prc.tcg),prc.tcg),nrow=2))
          
          fetm <- fisher.test(matrix(c((sem.ncg-sem.tcg),sem.tcg,
                                       (prm.ncg-prm.tcg),prm.tcg),nrow=2))
          
          dftad.seauto.test <- rbind(dftad.seauto.test,data.frame(se.id=names(segri),
                                                                  gene.id=genelisti[j],
                                                                  se.m.score=sem.tcg,
                                                                  gene.m.prscore=prm.tcg,
                                                                  se.m.perc=sem.tcg/sem.ncg,
                                                                  gene.m.prperc=prm.tcg/prm.ncg,
                                                                  fet.m.pval=as.numeric(fetm$p.value),
                                                                  se.c.score=sec.tcg,
                                                                  gene.c.prscore=prc.tcg,
                                                                  se.c.perc=sec.tcg/sec.ncg,
                                                                  gene.c.prperc=prc.tcg/prc.ncg,
                                                                  fet.c.pval=as.numeric(fetc$p.value),
                                                                  stringsAsFactors=F))
        }
        
        message(i," : ",j)
      }
    }
    message(i)
  }
  
  prscore.testlist <- list(dftad.seauto,dftad.seauto.test)
  names(prscore.testlist) <- c("bg.allse","test.seoi")
  save(prscore.testlist,file="prmethy-scoretest_tad-genes-seoi.rda")
  
}

# summaries and histograms
{
  head(prscore.testlist$bg.allse)
  
  summary(prscore.testlist$bg.allse$se.c.perc-prscore.testlist$bg.allse$gene.c.prperc)
  hist(prscore.testlist$bg.allse$se.c.perc-prscore.testlist$bg.allse$gene.c.prperc)
  summary(prscore.testlist$bg.allse$se.m.perc-prscore.testlist$bg.allse$gene.m.prperc)
  hist(prscore.testlist$bg.allse$se.m.perc-prscore.testlist$bg.allse$gene.m.prperc)
  
  summary(prscore.testlist$test.seoi$se.c.perc-prscore.testlist$test.seoi$gene.c.prperc)
  hist(prscore.testlist$test.seoi$se.c.perc-prscore.testlist$test.seoi$gene.c.prperc)
  
}

# plot all differences
{
  plot(density(prscore.testlist$bg.allse$se.c.perc-prscore.testlist$bg.allse$gene.c.prperc),col="red",
       ylim=c(0,2),xlab="SE - TGP score (cross-sectional, percent T-methyl. CpGs)",
       main="Distributions of region differences in\nMethyl. Score (percent T-methyl. CpGs, SE - TGP)")
  lines(density(prscore.testlist$test.seoi$se.c.perc-prscore.testlist$test.seoi$gene.c.prperc),col="blue")
  legend("top",horiz="T",legend=c("All SE","SE of interest"),lty=c(1,1),col=c("red","blue"))
  
}

# plot only significant differences in each
{
  pbh.sig <- 0.05
  whichm.seall <- p.adjust(prscore.testlist$bg.allse$fet.m.pval,method="BH")
  whichm.seall <- which(whichm.seall<pbh.sig)
  whichc.seall <- p.adjust(prscore.testlist$bg.allse$fet.c.pval,method="BH")
  whichc.seall <- which(whichc.seall<pbh.sig)
  
  whichm.seoi <- p.adjust(prscore.testlist$test.seoi$fet.m.pval,method="BH")
  whichm.seoi <- which(whichm.seoi<pbh.sig)
  whichc.seoi <- p.adjust(prscore.testlist$test.seoi$fet.c.pval,method="BH")
  whichc.seoi <- which(whichc.seoi<pbh.sig)
}

# density plots, composite
{
  dev.off()
  par(mfrow=c(2,2),oma=c(2,2,4,1))
  # plot 1. all densities, crx
  main1 <- paste0("All Pairs (Cross-section)\n All SE (N = ",length(prscore.testlist$bg.allse$se.c.perc)," pairs) vs. SE.OI (N = ",length(prscore.testlist$test.seoi$se.c.perc)," pairs)")
  plot(density(prscore.testlist$bg.allse$se.c.perc-prscore.testlist$bg.allse$gene.c.prperc),col="red",
       ylim=c(0,2),xlab="",main=main1)
  lines(density(prscore.testlist$test.seoi$se.c.perc-prscore.testlist$test.seoi$gene.c.prperc),col="blue")
  # plot 2. all densities, match
  main2 <- paste0("All Pairs (Matched)\n All SE (N = ",length(prscore.testlist$bg.allse$se.m.perc)," pairs) vs. SE.OI (N = ",length(prscore.testlist$test.seoi$se.m.perc)," pairs)")
  plot(density(prscore.testlist$bg.allse$se.m.perc-prscore.testlist$bg.allse$gene.m.prperc),col="red",
       ylim=c(0,2),xlab="",main=main2)
  lines(density(prscore.testlist$test.seoi$se.m.perc-prscore.testlist$test.seoi$gene.m.prperc),col="blue")
  # plot 3. sig densities, crx
  main3 <- paste0("Sig. Pairs (P-bh < ",pbh.sig,", Cross-section)\n All SE (N = ",length(prscore.testlist$bg.allse$se.c.perc[whichc.seall])," pairs) vs. SE.OI (N = ",length(prscore.testlist$test.seoi$gene.c.prperc[whichc.seoi])," pairs)")
  plot(density(prscore.testlist$bg.allse$se.c.perc[whichc.seall]-prscore.testlist$bg.allse$gene.c.prperc[whichc.seall]),col="red",
       ylim=c(0,2),xlab="",main=main3)
  lines(density(prscore.testlist$test.seoi$se.c.perc[whichc.seoi]-prscore.testlist$test.seoi$gene.c.prperc[whichc.seoi]),col="blue")
  # plot 4. sig densities, matched
  main4 <- paste0("Sig. Pairs (P-bh < ",pbh.sig,", Matched)\n All SE (N = ",length(prscore.testlist$bg.allse$se.m.perc[whichm.seall])," pairs) vs. SE.OI (N = ",length(prscore.testlist$test.seoi$gene.m.prperc[whichm.seoi])," pairs)")
  plot(density(prscore.testlist$bg.allse$se.c.perc[whichm.seall]-prscore.testlist$bg.allse$gene.c.prperc[whichm.seall]),col="red",
       ylim=c(0,2),xlab="",main=main4)
  lines(density(prscore.testlist$test.seoi$se.c.perc[whichm.seoi]-prscore.testlist$test.seoi$gene.c.prperc[whichm.seoi]),col="blue")
  
  mtext("Distributions of Methyl. Percent\nDifferences (% CpGs where T > N, SE - TAD gene Promoter/TGP score)\nP-values from Fisher Exact Tests",side=3,outer=T)
  mtext("SE - TGP (% T > N methyl CpGS)",outer=T,side=1)
  legend("top",horiz="T",legend=c("All SE","SE of interest"),lty=c(1,1),col=c("red","blue"),bty="n",cex=0.8)
  
}

# info and get regions of highest difference in seoi
{
  dfoi <- prscore.testlist$test.seoi
  dim(dfoi[(dfoi$se.c.perc-dfoi$gene.c.prperc)>0.5,]) # [1] 15 12, pairs
  dim(dfoi[(dfoi$se.m.perc-dfoi$gene.m.prperc)>0.5,]) # [1] 24 12, pairs
  
  # cross-section test
  unique(dfoi[(dfoi$se.c.perc-dfoi$gene.c.prperc)>0.5,]$se.id)
  #[1] "chr1:47895497-47915572"   "chr7:27134948-27144668"   "chr9:36254473-36338542"  
  #[4] "chr9:124396261-124462692" "chr11:65584824-65625749" 
  dfoi[(dfoi$se.c.perc-dfoi$gene.c.prperc)>0.5,]$gene.id
  #[1] "SLC5A9"       "FOXE3"        "HOXA-AS3"     "HOTTIP"       "HOXA10-HOXA9" "HOXA11-AS"   
  #[7] "HOXA10"       "MIR196B"      "LINC00950"    "PSMD5-AS1"    "PSMD5"        "MIR4690"     
  #[13] "CATSPER1"     "CTSW"         "TSGA10IP" 
  
  # matched test
  unique(dfoi[(dfoi$se.m.perc-dfoi$gene.m.prperc)>0.5,]$se.id)
  #[1] "chr1:47895497-47915572"   "chr7:27134948-27144668"   "chr9:36254473-36338542"  
  #[4] "chr9:124396261-124462692" "chr11:65584824-65625749"
  dfoi[(dfoi$se.m.perc-dfoi$gene.m.prperc)>0.5,]$gene.id
  #[1] "SLC5A9"       "CMPK1"        "STIL"         "HOXA-AS3"     "HOTTIP"       "HOXA10-HOXA9"
  #[7] "HOXA11-AS"    "HOXA4"        "HOXA7"        "HOXA10"       "HOXA13"       "CLTA"        
  #[13] "FAM221B"      "TMEM8B"       "HINT2"        "LINC00950"    "MELK"         "LOC100288842"
  #[19] "PSMD5-AS1"    "PSMD5"        "KAT5"         "CATSPER1"     "CTSW"         "TSGA10IP" 
  
  # combined
  intersect(unique(dfoi[(dfoi$se.c.perc-dfoi$gene.c.prperc)>0.5,]$se.id),unique(dfoi[(dfoi$se.m.perc-dfoi$gene.m.prperc)>0.5,]$se.id))
  #[1] "chr1:47895497-47915572"   "chr7:27134948-27144668"   "chr9:36254473-36338542"  
  #[4] "chr9:124396261-124462692" "chr11:65584824-65625749" 
  intersect(dfoi[(dfoi$se.c.perc-dfoi$gene.c.prperc)>0.5,]$gene.id,dfoi[(dfoi$se.m.perc-dfoi$gene.m.prperc)>0.5,]$gene.id)
  #[1] "SLC5A9"       "HOXA-AS3"     "HOTTIP"       "HOXA10-HOXA9" "HOXA11-AS"    "HOXA10"      
  #[7] "LINC00950"    "PSMD5-AS1"    "PSMD5"        "CATSPER1"     "CTSW"         "TSGA10IP" 
  
  
}

#======================================================================
# prepare methylation data at intergenic, OpenSea CpGs in TAD regions
#======================================================================



