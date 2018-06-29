# script_analysis: get methylation region summaries for se.oi

#======================
# methyl region tallies
#======================

dfsem <- as.data.frame(matrix(nrow=0,ncol=8)); 
colnames(dfsem) <- c("se.id","n.cpg","n.methyl.tm","n.methyl.tc","xdif.tnm","xdif.tnc","mdif.tnm","mdif.tnc")
se.oi.list <- unique(c(dmpc.en.seid,dmpm.en.seid))
for(i in 1:length(se.oi.list)){
  sei <- se.oi.list[i]
  mi <- betaval.ssef[[which(names(betaval.ssef)==sei)]]$betaval.match
  mi.dif <- rowMeans(mi[,substr(colnames(mi),14,15)=="01"])-rowMeans(mi[,substr(colnames(mi),14,15)=="11"])
  mi.tn <- length(mi.dif[mi.dif>0])
  mi.xdif <- mean(mi.dif);mi.mdif <- median(mi.dif)
  
  ci <- betaval.ssef[[which(names(betaval.ssef)==sei)]]$betaval.crx
  ci.dif <- rowMeans(ci[,substr(colnames(ci),14,15)=="01"])-rowMeans(ci[,substr(colnames(ci),14,15)=="11"])
  ci.tn <- length(ci.dif[ci.dif>0])
  ci.xdif <- mean(ci.dif);ci.mdif <- median(ci.dif)
  
  dfsem <- rbind(dfsem,data.frame(se.id=sei,
                                  n.cpg=nrow(mi),
                                  n.methyl.tm=mi.tn,
                                  n.methyl.tc=ci.tn,
                                  xdif.tnm=mi.xdif,
                                  xdif.tnc=ci.xdif,
                                  mdif.tnm=mi.mdif,
                                  mdif.tnc=mi.mdif,
                                  stringsAsFactors = F))
  message(i)
  
}
 
save(dfsem,file="df-methyl-region_se-oi_coad-t-vs-n-mc.rda") 

nrow(dfsem[dfsem$xdif.tnm>0 & dfsem$xdif.tnc>0,]) # 34
nrow(dfsem[dfsem$xdif.tnm>0.1 & dfsem$xdif.tnc>0,]) # 7
nrow(dfsem[dfsem$xdif.tnm>0.1 & dfsem$xdif.tnc>0.1,]) # 6

nrow(dfsem[dfsem$mdif.tnm>0 & dfsem$mdif.tnc>0,]) # 30
nrow(dfsem[dfsem$mdif.tnm>0.1 & dfsem$mdif.tnc>0,]) # 5
nrow(dfsem[dfsem$mdif.tnm>0.1 & dfsem$mdif.tnc>0.1,]) # 5

write.csv(dfsem[dfsem$xdif.tnm>0.1 & dfsem$xdif.tnc>0.1,],file="df-methyrgn_se-oi-thyp.csv")

se.oi.mrgn.dmpen <- dfsem[dfsem$xdif.tnm>0.1 & dfsem$xdif.tnc>0.1,]$se.id
save(se.oi.mrgn.dmpen,file="seoi-list-mrgn-dmpen_coad-t-vs-n.rda")

#============================================
# region methylation and methyl differences
#============================================

for(i in 1:length(se.oi.mrgn.dmpen)){
  seoii <- se.oi.mrgn.dmpen[i]
  bpic <- betaval.ssef[[which(names(betaval.ssef)==seoii)]]$betaval.crx
  bpim <- betaval.ssef[[which(names(betaval.ssef)==seoii)]]$betaval.match
  
  titlei <- paste0(seoii,".jpg")
  
  jpeg(titlei,4,8,units="in",res=400)
  #par(mfrow=c(4,1),oma=c(4,2,2,1))
  par(mfrow=c(4,1))
  
  boxplot(t(bpic[,substr(colnames(bpic),14,15)=="01"]),las=2,ylab="Tumor Methyl.",xaxt="n",main=seoii)
  boxplot(t(bpic[,substr(colnames(bpic),14,15)=="11"]),las=2,ylab="Normal Methyl.",xaxt="n")
  barplot(rowMeans(bpic[,substr(colnames(bpic),14,15)=="01"])-rowMeans(bpic[,substr(colnames(bpic),14,15)=="11"]),
          ylab="T-N (means) Cross-sectional",xlab="",xaxt="n")
  barplot(rowMeans(bpim[,substr(colnames(bpim),14,15)=="01"])-rowMeans(bpim[,substr(colnames(bpim),14,15)=="11"]),
          ylab="T-N (means) Matched",las=2,cex.axis = 0.5)
  
  dev.off()
}
