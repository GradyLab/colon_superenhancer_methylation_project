# script: region methylation testing, flank and border cpgs vs. se regions

seoi <- se.oi.mrgn.dmpen

#================================================
# get flanking region methylation, and compare
#================================================

rgnmethydf.seoi <- as.data.frame(matrix(nrow=0,ncol=9))
colnames(rgnmethydf.seoi) <- c("se.id","ncpg.se",
                                        "ncpg.fr",
                                        "ncpgse.crx.tmethyl",
                                        "ncpgfr.crx.tmethyl",
                                        "ncpgse.match.tmethyl",
                                        "ncpgfr.match.tmethyl",
                                        "fet.crx.pval",
                                        "fet.match.pval",
                                        stringsAsFactors = F)

for(i in 1:length(seoi)){
  fi.name <- ssef.auto[names(ssef.auto)==seoi[i]]$fr.id
  fr.crx <- betaval.fr[[which(names(betaval.fr)==fi.name)]]$beta.crx
  fr.match <- betaval.fr[[which(names(betaval.fr)==fi.name)]]$beta.match
  
  frm.ti <- rowMeans(fr.match[,substr(colnames(fr.match),14,15)=="01"])
  frm.ni <- rowMeans(fr.match[,substr(colnames(fr.match),14,15)=="11"])
  frm.nmethyt <- length(frm.ti[frm.ti>frm.ni])
  
  frc.ti <- rowMeans(fr.crx[,substr(colnames(fr.crx),14,15)=="01"])
  frc.ni <- rowMeans(fr.crx[,substr(colnames(fr.crx),14,15)=="11"])
  frc.nmethyt <- length(frc.ti[frc.ti>frc.ni])
  
  
  sei.name <- seoi[i]
  se.crx <- betaval.ssef[[which(names(betaval.ssef)==sei.name)]]$betaval.crx
  se.match <- betaval.ssef[[which(names(betaval.ssef)==sei.name)]]$betaval.match
  
  sem.ti <- rowMeans(se.match[,substr(colnames(se.match),14,15)=="01"])
  sem.ni <- rowMeans(se.match[,substr(colnames(se.match),14,15)=="11"])
  sem.nmethyt <- length(sem.ti[sem.ti>sem.ni])

  sec.ti <- rowMeans(se.crx[,substr(colnames(se.crx),14,15)=="01"])
  sec.ni <- rowMeans(se.crx[,substr(colnames(se.crx),14,15)=="11"])
  sec.nmethyt <- length(sec.ti[sec.ti>sec.ni])
  
  
  ncg.sei <- as.numeric(ssef.auto[names(ssef.auto)==sei.name]$ncpg450)
  ncg.fri <- as.numeric(ssef.auto[names(ssef.auto)==sei.name]$frol.ncpg450filt)
  
  # fet 
  matchfet <- as.numeric(fisher.test(matrix(c(ncg.sei-sem.nmethyt,
                                              ncg.fri-frm.nmethyt,
                                              sem.nmethyt,
                                              frm.nmethyt),nrow=2))$p.value)
  crxfet <- as.numeric(fisher.test(matrix(c(ncg.sei-sec.nmethyt,
                                              ncg.fri-frc.nmethyt,
                                              sec.nmethyt,
                                              frc.nmethyt),nrow=2))$p.value)
  
  rgnmethydf.seoi <- rbind(rgnmethydf.seoi,data.frame(se.id=sei.name,
                                                      ncpg.se=ncg.sei,
                                                      ncpg.fr=ncg.fri,
                                                      ncpgse.crx.tmethyl=sec.nmethyt,
                                                      ncpgfr.crx.tmethyl=frc.nmethyt,
                                                      ncpgse.match.tmethyl=sem.nmethyt,
                                                      ncpgfr.match.tmethyl=frm.nmethyt,
                                                      fet.crx.pval=crxfet,
                                                      fet.match.pval=matchfet,
                                                      stringsAsFactors = F))
  
  message(i)
  
}

save(rgnmethydf.seoi,file="df-rgnmethy-fet_coad-t-vs-n.rda")
write.csv(rgnmethydf.seoi,file="df-rgnmethy-fet_coad-t-vs-n.csv",row.names=F)

#==========================
# region smooth overlays
#========================

i=6

jtitle=paste0("smoothplots_",gsub(":","_",seoi[i]),".jpg")
jpeg(jtitle,6,8,units="in",res=400)
{
  # flanking region
  {
    fi.name <- ssef.auto[names(ssef.auto)==seoi[i]]$fr.id
    fr.crx <- betaval.fr[[which(names(betaval.fr)==fi.name)]]$beta.crx
    fr.match <- betaval.fr[[which(names(betaval.fr)==fi.name)]]$beta.match
    frcg <- rownames(fr.match)
    
    frm.ti <- rowMeans(fr.match[,substr(colnames(fr.match),14,15)=="01"])
    frm.ni <- rowMeans(fr.match[,substr(colnames(fr.match),14,15)=="11"])
    frm.nmethyt <- length(frm.ti[frm.ti>frm.ni])
    
    frc.ti <- rowMeans(fr.crx[,substr(colnames(fr.crx),14,15)=="01"])
    frc.ni <- rowMeans(fr.crx[,substr(colnames(fr.crx),14,15)=="11"])
    frc.nmethyt <- length(frc.ti[frc.ti>frc.ni])
  }
  
  
  # super-enhancer region
  {
    sei.name <- seoi[i]
    se.crx <- betaval.ssef[[which(names(betaval.ssef)==sei.name)]]$betaval.crx
    se.match <- betaval.ssef[[which(names(betaval.ssef)==sei.name)]]$betaval.match
    secg <- rownames(se.crx)
    
    sem.ti <- rowMeans(se.match[,substr(colnames(se.match),14,15)=="01"])
    sem.ni <- rowMeans(se.match[,substr(colnames(se.match),14,15)=="11"])
    sem.nmethyt <- length(sem.ti[sem.ti>sem.ni])
    
    sec.ti <- rowMeans(se.crx[,substr(colnames(se.crx),14,15)=="01"])
    sec.ni <- rowMeans(se.crx[,substr(colnames(se.crx),14,15)=="11"])
    sec.nmethyt <- length(sec.ti[sec.ti>sec.ni])
  }
  
  cggri <- cggr[names(cggr) %in% c(frcg,secg)]
  sic.methyl <- rbind(fr.crx,se.crx); sic.methyl <- sic.methyl[order(match(rownames(sic.methyl),names(cggri))),]
  sim.methyl <- rbind(fr.match,se.match); sim.methyl <- sim.methyl[order(match(rownames(sim.methyl),names(cggri))),]
  
  identical(rownames(sic.methyl),names(cggri)); identical(rownames(sim.methyl),names(cggri))
  
  par(mfrow=c(2,1))
  # crx data
  {
    plot(predict(loess(rowMeans(sic.methyl[,substr(colnames(sic.methyl),14,15)=="01"])~start(cggri))),
         type="l",col="purple",ylab="Methyl. Prediction",main=paste0("Cross-sectional Data\nSE.ID:",sei.name),xlab="Pos")
    lines(predict(loess(rowMeans(sic.methyl[,substr(colnames(sic.methyl),14,15)=="11"])~start(cggri))),type="l",col="green")
    abline(v=(which(rownames(sic.methyl)==rownames(se.crx)[1])+1),col="blue")
    abline(v=(which(rownames(sic.methyl)==rownames(se.crx)[nrow(se.crx)])+1),col="blue")
    abline(h=0.5)
    
    legend("topright",legend=c("N","T","FR"),lty=c(1,1,1),lwd=c(1,1,1),col=c("green","purple","blue"),h="T",bty="n")
    
  }
  
  # matched data
  {
    plot(predict(loess(rowMeans(sim.methyl[,substr(colnames(sim.methyl),14,15)=="01"])~start(cggri))),
         type="l",col="purple",ylab="Methyl. Prediction",main=paste0("Matched Data\nSE.ID:",sei.name),xlab="Pos")
    lines(predict(loess(rowMeans(sim.methyl[,substr(colnames(sim.methyl),14,15)=="11"])~start(cggri))),type="l",col="green")
    abline(v=(which(rownames(sim.methyl)==rownames(se.crx)[1])+1),col="blue")
    abline(v=(which(rownames(sim.methyl)==rownames(se.crx)[nrow(se.crx)])+1),col="blue")
    abline(h=0.5)
    
    legend("topright",legend=c("N","T","FR"),lty=c(1,1,1),lwd=c(1,1,1),col=c("green","purple","blue"),h="T",bty="n")
  }
}

dev.off()
