# script: analyze and prepare expression data for super-enhancers of interest

ecm <- as.matrix(exc.mrna);class(ecm) <- "numeric"; summary(rowMeans(ecm))
emm <- as.matrix(exm.mrna);class(emm) <- "numeric"; summary(rowMeans(emm))

#==================================================
# tad-assoc gene expr, multiple scales, for se oi
#==================================================

exc.dflist <- list()
exm.dflist <- list()

#se.oi <- se.oi.mrgn.dmpen[6]
se.oi.list <- names(ssef.auto)
for(k in 1:length(se.oi.list)){
  se.oi <- se.oi.list[k]
  seoi.gr <- ssef.auto[names(ssef.auto)==se.oi]
  se.genes <- unlist(strsplit(mcols(seoi.gr)$tadgeneid,";"))
  segenes.available <- intersect(gsub("\\|.*","",rownames(ecm)),se.genes)
  
  # cross-sectional data
  {
    if(length(segenes.available)>1){
      # rsem and log2 scale data
      tadgexc.rsem <- ecm[gsub("\\|.*","",rownames(ecm)) %in% segenes.available,]
      tadgexc.l2rsem <- log2(ecm[gsub("\\|.*","",rownames(ecm)) %in% segenes.available,]+0.1)
      
      tgec.rmt.rsem <- rowMeans(tadgexc.rsem[,substr(colnames(tadgexc.rsem),14,15)=="01"])
      tgec.rmn.rsem <- rowMeans(tadgexc.rsem[,substr(colnames(tadgexc.rsem),14,15)=="11"])
      summary(tgec.rmt.rsem)
      summary(tgec.rmn.rsem)
      tgec.rmt.rsem-tgec.rmn.rsem
      
      tgec.rmt.l2rsem <- rowMeans(log2(tadgexc.rsem[,substr(colnames(tadgexc.rsem),14,15)=="01"]+0.1))
      tgec.rmn.l2rsem <- rowMeans(log2(tadgexc.rsem[,substr(colnames(tadgexc.rsem),14,15)=="11"]+0.1))
      summary(tgec.rmt.l2rsem)
      summary(tgec.rmn.l2rsem)
      tgec.rmt.l2rsem-tgec.rmn.l2rsem
      
      # log2fc data, tumor only
      ec.xn <- tgec.rmn.rsem
      ec.l2fc.t <- log2((tadgexc.rsem[,substr(colnames(tadgexc.rsem),14,15)=="01"]+0.1)/(ec.xn+0.1))
      summary(rowMeans(ec.l2fc.t))
      
      testdfc <- as.data.frame(matrix(nrow=0,ncol=11))
      colnames(testdfc) <- c("se.id","gene.id","tpunadj","wcpunadj",
                             "xt.rsem","xn.rsem","xt.l2rsem","xn.l2rsem",
                             "xtn.rsemdif","xtn.l2rsemdif",
                             "xt.l2fc")
      for(i in 1:length(segenes.available)){
        genei <- segenes.available[i]
        
        tdati <- tadgexc.l2rsem[gsub("\\|.*","",rownames(tadgexc.l2rsem))==genei,
                                substr(colnames(tadgexc.l2rsem),14,15)=="01"]
        ndati <- tadgexc.l2rsem[gsub("\\|.*","",rownames(tadgexc.l2rsem))==genei,
                                substr(colnames(tadgexc.l2rsem),14,15)=="11"]
        
        if(sd(tdati) > 0 & sd(ndati) > 0){
          tip <- as.numeric(t.test(tdati,ndati)$p.value)
          wcp <- as.numeric(wilcox.test(tdati,ndati)$p.value)
          
          testdfc <- rbind(testdfc,data.frame(se.id=se.oi,
                                              gene.id=genei,
                                              tpunadj=as.numeric(tip),
                                              wcpunadj=as.numeric(wcp),
                                              xt.rsem=mean(as.numeric(2^tdati)),
                                              xn.rsem=mean(as.numeric(2^tdati)),
                                              xt.l2rsem=mean(as.numeric(tdati)),
                                              xn.l2rsem=mean(as.numeric(ndati)),
                                              xtn.rsemdif=mean(as.numeric(2^tdati))-mean(as.numeric(2^ndati)),
                                              xtn.l2rsemdif=mean(as.numeric(tdati))-mean(as.numeric(ndati)),
                                              xt.l2fc=mean(as.numeric(ec.l2fc.t[gsub("\\|.*","",rownames(ec.l2fc.t))==genei,])),
                                              stringsAsFactors = F))
        }
        
        message(k,":",i)
      }
      
      exc.dflist[[length(exc.dflist)+1]] <- testdfc
      names(exc.dflist)[length(exc.dflist)] <- se.oi
      
      message(k)
    }
    }
    
}

save(exc.dflist,file="expr-crx-analysis_seALL-dflist.rda")


# matched/paired data
for(k in 1:length(se.oi.list)){
  {
    #se.oi <- se.oi.mrgn.dmpen[6]
    se.oi <- se.oi.list[k]
    seoi.gr <- ssef.auto[names(ssef.auto)==se.oi]
    se.genes <- unlist(strsplit(mcols(seoi.gr)$tadgeneid,";"))
    segenes.available <- intersect(gsub("\\|.*","",rownames(ecm)),se.genes)
    
    if(length(segenes.available)>1){
      # rsem and log2 scale data
      tadgexm.rsem <- emm[gsub("\\|.*","",rownames(emm)) %in% segenes.available,]
      tadgexm.l2rsem <- log2(emm[gsub("\\|.*","",rownames(emm)) %in% segenes.available,]+0.1)
      
      
      # log2fc data, tumor only
      em.xn <- tgec.rmn.rsem
      em.l2fc.t <- log2((tadgexm.rsem[,substr(colnames(tadgexm.rsem),14,15)=="01"]+0.1)/(em.xn+0.1))
      summary(rowMeans(em.l2fc.t))
      
      testdfm <- as.data.frame(matrix(nrow=0,ncol=11))
      colnames(testdfm) <- c("se.id","gene.id","tpunadj","wcpunadj",
                             "xt.rsem","xn.rsem","xt.l2rsem","xn.l2rsem",
                             "xtn.rsemdif","xtn.l2rsemdif",
                             "xt.l2fc")
      for(i in 1:length(segenes.available)){
        genei <- segenes.available[i]
        
        tdati <- tadgexm.l2rsem[gsub("\\|.*","",rownames(tadgexm.l2rsem))==genei,
                                substr(colnames(tadgexm.l2rsem),14,15)=="01"]
        ndati <- tadgexm.l2rsem[gsub("\\|.*","",rownames(tadgexm.l2rsem))==genei,
                                substr(colnames(tadgexm.l2rsem),14,15)=="11"]
        
        #identical(substr(colnames(tdati),9,12),substr(colnames(ndati),9,12))
        
        tip <- as.numeric(t.test(tdati,ndati,paired=T)$p.value)
        wcp <- as.numeric(wilcox.test(tdati,ndati,paired=T)$p.value)
        
        testdfm <- rbind(testdfm,data.frame(se.id=se.oi,
                                            gene.id=genei,
                                            tpunadj=as.numeric(tip),
                                            wcpunadj=as.numeric(wcp),
                                            xt.rsem=mean(as.numeric(2^tdati)),
                                            xn.rsem=mean(as.numeric(2^tdati)),
                                            xt.l2rsem=mean(as.numeric(tdati)),
                                            xn.l2rsem=mean(as.numeric(ndati)),
                                            xtn.rsemdif=mean(as.numeric(2^tdati))-mean(as.numeric(2^ndati)),
                                            xtn.l2rsemdif=mean(as.numeric(tdati))-mean(as.numeric(ndati)),
                                            xt.l2fc=mean(as.numeric(em.l2fc.t[gsub("\\|.*","",rownames(em.l2fc.t))==genei,])),
                                            stringsAsFactors = F))
        message(k,":",i)
      }
      
      exm.dflist[[length(exm.dflist)+1]] <- testdfm
      names(exm.dflist)[length(exm.dflist)] <- se.oi
    }
    }
    message(k)
}

save(exm.dflist,file="expr-matched-analysis_seALL-dflist.rda")

#======================================
# get list of genes of interest by se
#======================================
# Notes: in crx, significant genes (0.01) with <= -1 mean log2FC tumor expr
# in matched, get significnat genes with <= -2 mean log2FC tumor expr

segoi.list <- list()
#for(i in 1:6){
for(i in 1:length(exc.dflist)){
  sei <- names(exc.dflist)[i]
  ci <- exc.dflist[[i]]
  ci <- ci[ci$tpunadj<0.01 & ci$xt.l2fc<= -1,]$gene.id
  
  mi <- exm.dflist[[i]]
  mi <- mi[mi$tpunadj<0.01 & mi$xt.l2fc <= -2,]$gene.id
  
  segoi.list[[i]] <- unique(ci,mi)
  names(segoi.list)[i] <- sei
  message(i)
}

save(segoi.list,file="seALL-goi-expfilt-alltest_genelist.rda")

#==================
# prepare goi expr
#==================
goilist <- c(); for(i in 1:6){goilist <- c(goilist,segoi.list[[i]])}

goiexc <- exc.mrna[gsub("\\|.*","",rownames(exc.mrna)) %in% goilist,]
goiexm <- exm.mrna[gsub("\\|.*","",rownames(exm.mrna)) %in% goilist,]

goiex.list <- list(goiexc,goiexm)
names(goiex.list) <- c("goiex.rsem.crx","goiex.rsem.match")
save(goiex.list,file="rsem-ex-list_se-goi-filt.rda")

#===============================================
# prepare gene promoter and se methylation data
#===============================================
prcg <- repman450[repman450$group %in% c("TSS1500","TSS200","1stExon"),]
prcg.oi <- prcg[prcg$name %in% as.character(unlist(segoi.list)),]$cpg
bvalpr.segoi <- getBeta(g3[rownames(g3) %in% prcg.oi,])
save(bvalpr.segoi,file="betval-methyl-coad-t-n-all_segoi-promotercg.rda")

ssef.auto.oi <- ssef.auto[names(ssef.auto) %in% se.oi.mrgn.dmpen]
save(ssef.auto.oi,file="ssef-grlist-info_seoi.rda")








