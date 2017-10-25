## ----setup, include=FALSE, message=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 5.5, fig.height = 3.5, warning = FALSE)
library(knitr)

## ----env.set.up, warning = F, message=FALSE------------------------------
library(devtools)
install_github("jtlovell/qtlTools")
library(qtlTools)

## ----loadmultitrait------------------------------------------------------
data(multitrait)

## ----fakemap-------------------------------------------------------------
map<-pullMap(multitrait)
map$bp<-0
for(i in unique(map$chr)){
  n<-sum(map$chr==i)
  p<-sin((1:n/n)*pi)
  map$bp[map$chr==i]<-cumsum(p*1000000)
}

## ----fakegff-------------------------------------------------------------
gff<-data.frame(chr = rep(paste0("scaffold_",1:5),each = 200),
   feature = rep("gene",1000),
   start = rep(seq(from = 0, to = max(map$bp), length = 200), 5),
   end = rep(seq(from = 0, to = max(map$bp), length = 200))+1000,
   strand = rep("+",1000),
   attribute = paste0("gene",1:1000,";","gene",1:1000,".1"), stringsAsFactors=F)

## ------------------------------------------------------------------------
geneCM<-findGenecM(cross = multitrait, marker.info = map, gff = gff,
   gffCols = c("chr","feature","start","end","strand","attribute"))

## ----plotrecom, fig.height = 4, fig.width = 6, fig.cap = "Physical and mapping positions of fake data for the 5 A. thaliana chromosomes."----
par(mfrow=c(2,3))
for(i in unique(map$chr)){
  with(geneCM[geneCM$chr==i,], plot(pos,bp, col="grey", 
                                ylab = "physical position (bp)",
                                xlab = "mapping position (cM)"))
  with(map[map$chr==i,], points(pos,bp, col=i, pch = 19, cex=.8))
}

## ---- fig.height = 3, fig.width = 6, fig.cap = "Scanone profile with confidence intervals."----
multitrait<-calc.genoprob(multitrait)
s1<-scanone(multitrait, method="hk", pheno.col=1)
perm<-scanone(multitrait, n.perm=100, method="hk",pheno.col=1, verbose=FALSE)
cis<-calcCis(cross = multitrait, s1.output=s1, perm.output=perm, drop=5)

par(mfrow = c(1,1))
plot(s1)
segmentsOnPeaks(multitrait, s1.output=s1, calcCisOutput = cis, int.y = 13.1)

## ------------------------------------------------------------------------
candGenes<-findGenesInterval(findGenecM.output = geneCM, calcCis.output = cis)
print(candGenes)

## ---- fig.height = 5, fig.width = 5, fig.cap = "correlation matrix of simulated gene expression data."----
cross<-subset(multitrait, ind = !is.na(pull.pheno(multitrait, 1)))
phe<-pull.pheno(cross, 1)

mult.fact<-exp(seq(from = 0, to = 50, length.out = 50))
facs<-sapply(1:50, function(x){
  scale(sapply(scale(phe), function(y) rnorm(n = 1, mean = y, sd = mult.fact[x])))
})
plot(sapply(1:50, function(x) cor(phe, facs[,x])),
     ylab = "cor. coef. (expression ~ chr5 QTL genotype)",
     xlab = "gene id")

## ------------------------------------------------------------------------
expression.covariates = facs
colnames(expression.covariates)<-paste0("gene",1:ncol(expression.covariates))
pheno.col = 1

qtl = makeqtl(cross, chr = max(s1)$chr, pos = max(s1)$pos, what = "prob")

test_additive<-covScanQTL(cross = cross,
                 pheno.col = 1,
                 qtl = qtl,
                 expression.covariates = expression.covariates,
                 qtl.method = "hk",
                 nperm = 20)
kable(head(test_additive), caption = "top candidate genes using a simple, additive covariate scan approach.")

## ------------------------------------------------------------------------
qtl2 = makeqtl(cross, chr = summary(s1)$chr[4:5], 
               pos = summary(s1)$pos[4:5], what = "prob")
test_epistasis<-covScanQTL(cross = cross,
                 pheno.col = 1,
                 qtl = qtl2,
                 which.epiqtl = 1,
                 focalqtl.index = 2,
                 expression.covariates = expression.covariates,
                 qtl.method = "hk",
                 nperm = 20)

kable(head(test_epistasis), caption = "top candidate genes using a covariate scan that incorporates epistasis")

## ------------------------------------------------------------------------
data(fake.f2)
cross<-fake.f2
phe<-pull.pheno(cross, 1)
mult.fact<-exp(seq(from = 0, to = 50, length.out = 50))
facs<-sapply(1:50, function(x){
  scale(sapply(scale(phe), function(y) rnorm(n = 1, mean = y, sd = mult.fact[x])))
})
facs<-data.frame(facs)
colnames(facs)<-paste0("gene",1:ncol(facs))
cross<-calc.genoprob(cross)
sex = data.frame(sex = pull.pheno(cross, pheno.col = "sex"))
s1<-scanone(cross, addcovar = sex)
qtl = makeqtl(cross, chr = summary(s1)$chr[1], pos = summary(s1)$pos[1], what = "prob")

test_QTLxE<-covScanQTL(cross = cross,
                 pheno.col = 1,
                 qtl = qtl,
                 addcovar = sex,
                 intcovar = sex,
                 expression.covariates = facs,
                 qtl.method = "hk",
                 nperm = 20)

kable(head(test_QTLxE), caption = "top candidate genes using a covariate scan that incorporates GxE")

