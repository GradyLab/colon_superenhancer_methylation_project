# script: DEG enrichment testing
# note: prepare exm- and exc-dflist objects as in script "script_analysis_seExpr.R"

#===================
# MATCHED EXPR DEGs
#===================
# aggregate matched expr data

seall.me <- exm.dflist[[1]]

for(i in 2:length(exm.dflist)){
  seall.me <- rbind(seall.me,exm.dflist[[i]])
  message(i)
}

dim(seall.me) # 19551    11

# deg set, under-expressed in T
nrow(seall.me[seall.me$se.id %in% names(ssef5),]) # 60 TAD-assoc. genes from SE OI N=5
length(unique(seall.me$se.id)) # 1136/1206 SEs with genes
length(unique(seall.me$gene.id)) # 8301
length(unique(seall.me[seall.me$se.id %in% names(ssef5),]$se.id)) # all 5 se oi represented
summary(as.data.frame(table(seall.me$se.id))) # mean = 17.21 genes/se (overall)
summary(as.data.frame(table(seall.me[seall.me$se.id %in% names(ssef5),]$se.id))) # mean = 12 genes/se.oi (N=5)

# distribution comparisons

# 1. compare all genes, at se oi to random set of 60 as background
summary(seall.me[seall.me$se.id %in% names(ssef5),]$xtn.l2rsemdif)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-3.4402 -0.7000 -0.2058 -0.2523  0.2895  3.3407
summary(seall.me[c(sample(nrow(seall.me),60)),]$xtn.l2rsemdif)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-3.22118 -0.33586 -0.09636 -0.13586  0.31882  1.94352

# 2. compare all se.oi degs to random set of n degs as background
seall.me$tpbh <- p.adjust(seall.me$tpunadj,method="BH")
seall.me.deg <- seall.me[!duplicated(seall.me$gene.id),]; nrow(seall.me.deg) # 8301
seall.me.deg <- seall.me.deg[seall.me.deg$tpbh<0.05,]
nrow(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]) # 27
nrand <- 27
summary(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]$xtn.l2rsemdif)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-3.4402 -1.0778 -0.3546 -0.3487  0.3171  3.3407
summary(seall.me.deg[sample(nrow(seall.me.deg),nrand),]$xtn.l2rsemdif)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-3.62172 -1.04408 -0.25955 -0.03924  0.70351  4.49938

par(mfrow=c(2,1),oma=c(3,3,3,1))

plot(density(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]$xtn.l2rsemdif),
     ylim=c(0,1),lwd=2,main="Distributions of DEG Expression:\nAll DEGs",ylab="",
     xlab="")
vali.vector <- list()
for(i in 1:1000){
  vali <- seall.me.deg[sample(nrow(seall.me.deg),nrand),]$xtn.l2rsemdif
  vali <- vali[!is.na(vali)]
  abline(v=median(vali),col=rgb(0.5,0.1,0.3,0.2))
  vali.vector[[i]] <- vali
}
for(i in 1:length(vali.vector)){
  lines(density(vali.vector[[i]]),col=rgb(0.8,0.1,0.1,0.3))
}
lines(density(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]$xtn.l2rsemdif),
      ylim=c(0,1),lwd=2,col="black")
abline(v=median(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]$xtn.l2rsemdif),
       col="blue",lwd=2)
legend("topright",bty="n",legend=c("SE.OI","SE.Rand1K","SE.OI.med","SE.Rand1K.med"),
       lwd=c(2,1,2,1),lty=c(1,1,1,1),col=c("black",rgb(0.8,0.1,0.1,0.8),"blue",rgb(0.5,0.1,0.3,0.8)))

med.val <- c()
for(i in 1:length(vali.vector)){
  med.val <- c(med.val,median(vali.vector[[i]]))
}
plot(density(med.val),col=rgb(0.5,0.1,0.3,0.8),lty=2,xlab="Expr. Diff (mean, T - N, log2[RSEM])",ylab="",
     xlim=c(-5,5),main="Only Median DEG Expr. Diff.")
abline(v=median(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]$xtn.l2rsemdif),
       col="blue",lwd=2)
legend("topright",bty="n",legend=c("SE.OI.med","SE.Rand1K.med"),
       lwd=c(1,1),lty=c(1,2),col=c("blue",rgb(0.5,0.1,0.3,0.8)),
       cex=0.8)

mtext("Relative DEG Density",side=2,outer=T)


med.emp <- median(seall.me.deg[seall.me.deg$se.id %in% names(ssef5),]$xtn.l2rsemdif) # -0.3282743
table(med.val>med.emp)
#FALSE  TRUE 
#344   656

legend("topleft",legend="64.1% Med.Rand. > SE.OI.Med",
       bty="n",cex =0.7)


#207/1000 => 79% of reps have greater medain


