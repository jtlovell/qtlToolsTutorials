## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, include=FALSE------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE)

library(qtlToolsTutorials)
library(qtlTools)
data(markerData)

## ---- eval = F-----------------------------------------------------------
#  library(devtools)
#  install_github("jtlovell/qtlTools")
#  install_github("jtlovell/qtlToolsTutorials")
#  library(qtlToolsTutorials)
#  library(qtlTools)
#  
#  data(markerData)
#  data("makeGeneticMapTutorialTmpData")

## ---- echo = F, fig.cap = "distribution of reads mapping to the HAL2 P. hallii reference genome. Since it is an F2, We expect that most loci have a 1:2:1 AA:AB:BB ratio."----
kable(cbind(markerData[1:10,1:3],apply(markerData[1:10,4:10],2,round,2)), 
            caption = "example of 7 libraries and 10 markers from the markerData dataset. Each value represents the proportion of reads mapping to the reference genome vs. the alternative genome for a given marker 'HAL2' and library (column names)")
hist(as.matrix(markerData[,-c(1:3)])*100, breaks = 50,
     main = "Distribution of % reads mapping to the reference genome",
     xlab = "% reads mapping to reference")

## ------------------------------------------------------------------------
mar.chr = markerData$chr_HAL2
mar.pos = markerData$pos_HAL2
mar.id = markerData$HAL2
prop.ref = markerData[,-c(1:3)]

## ---- eval = F-----------------------------------------------------------
#  sw<-swGenotype(reference.prop = prop.ref,
#                 chr = mar.chr,
#                 pos = mar.pos,
#                 marker.id = mar.id)

## ---- echo = F, fig.cap = "Distribution of marker calls across the first two P. hallii chromosomes. We retain the hard calls (blue), which are the majority vote for any 20-marker window. "----
plot(sw$means[sw$means$chr %in% c("Chr01","Chr02"),4]*100, 
     xlab = "marker order (Chr 1-2)", 
     ylab = "% reads mapping to reference",
     pch = 16, cex = 1.5, col = "grey")
points(prop.ref[sw$means$chr %in% c("Chr01","Chr02"),1]*100, cex = .25)
points((sw$softCalls[sw$softCalls$chr %in% c("Chr01","Chr02"),4]*100)+2,
       pch = 16,cex = .5, col = "purple")
points((sw$hardCalls[sw$hardCalls$chr %in% c("Chr01","Chr02"),4]*100)-2,
       pch = 16,cex = .5, col = "dodgerblue")
legend("bottomright",
       c("mean by window","raw calls","softcalls","hardcalls (maj. vote)"),
       pch = c(16,16,16), pt.cex = c(1.5,.2,.5,.5), col = c("grey","black","purple","dodgerblue"))

## ---- eval = T-----------------------------------------------------------
hard.calls<-sw$hardCalls 
rownames(hard.calls)<-paste0(gsub("Chr0","",hard.calls$chr),"_",hard.calls$pos)
hard.calls<-hard.calls[,-c(1:3)]

## ---- eval = T-----------------------------------------------------------
hard.calls<-hard.calls[!duplicated(hard.calls),]
hard.calls<-t(hard.calls)

## ---- eval = T-----------------------------------------------------------
hard.calls[hard.calls==1]<-2
hard.calls[hard.calls==0.5]<-1
hard.calls[hard.calls==0]<-0

geno2cross(hard.calls)

## ---- eval = T-----------------------------------------------------------
geno2cross(hard.calls,crossfile = "~/Desktop/cross.csv",
           chr = sample(1:9, replace = T, size = ncol(hard.calls)),
           pos = runif(ncol(hard.calls), min = 0, max = 100))

## ---- eval = T-----------------------------------------------------------
cross<-read.cross("csv", file="~/Desktop/cross.csv",
                  genotypes=c(0,1,2), crosstype="f2")

## ------------------------------------------------------------------------
summary(cross)

## ---- echo = F, fig.cap = "genotype frequecies for each individual. Note that many individuals have 0 frequencies of either the AA or BB. These should be excluded."----
g <- pull.geno(cross)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))

## ------------------------------------------------------------------------
cross<-subset(cross,
              ind = gfreq[2,]<.8 & # <80% heterozygosity
                gfreq[1,]>0.02 & # >2% AA homozygote
                gfreq[3,]>0.02) # >2% BB homozygote

## ----eval = F------------------------------------------------------------
#  cross.sub<-dropSimilarMarkers(cross, blockSize = 500,
#                                rf.threshold = 0.02, runFullMatrix = T,
#                                byChr = F)

## ---- fig.cap = "recombination fraction among markers. Warmer colors have lower recombination fractions (more similar), darkblue colors indicate marker pairs that are uncorrelated. Since we randomized the position of markers, it is not surprising that there is no association between position and recombination fraction"----
summary(cross.sub)
par(mfrow = c(1,1))
plot.rf(cross.sub)

## ------------------------------------------------------------------------
lgmar<-formLinkageGroups(cross.sub, reorgMarkers=F, max.rf = .23,verbose = F)
head(lgmar)

## ------------------------------------------------------------------------
lgmar$true.chr<-splitText(rownames(lgmar))
lgmar$true.pos<-as.numeric(splitText(rownames(lgmar), num = 2))

## ---- eval = F-----------------------------------------------------------
#  marlist<-split(rownames(lgmar), as.factor(lgmar$true.chr))
#  cross2<-newLG(cross.sub, marlist, keep.rf = T)

## ---- fig.cap = "recombination fractions of markers, now the correlations are higher within chromosomes"----
plot.rf(cross2)

## ---- eval = F-----------------------------------------------------------
#  cross3<-tspOrder(cross = cross2,
#                   concorde_path = "/Users/John/Documents/concorde/TSP")

## ---- fig.cap = "recombination fractions of markers, now proximate markers have the lowest recombination fraction (highest correlation)"----
plot.rf(cross3)

## ---- eval = F-----------------------------------------------------------
#  cross4<-reDoCross(cross3, window = 3, min.distance = 2,
#                    initialEstMap = T, verbose = T)

## ---- fig.cap = "order of markers before (grey) and after (black) matching to physical orientiation of chromosomes."----
map<-pullMap(cross4)
map$phys.chr<-splitText(map$marker.name)
map$phys.pos<-as.numeric(splitText(map$marker.name, num = 2))/1e6
kable(head(map))
cross5<-matchMarkerOrder(cross4)

## ---- echo = F, fig.height = 8, fig.cap = "The frequency of the three genotypes, A/A (Red), B/B (Blue), A/B (Cyan), are plotted. Note some global bias towards the B/B genotype. This is not a problem, at least at this level of bias. On the right arm of Chr09 there is a near lack of B/B genotypes. Originally we dropped these distorted markers, which forced Chr09 to be split into two chromosomes. However, we now know that this distortion is real, caused by some genetic factor that when in the B/B state is selected against."----
gt<-geno.table(cross5, scanone.output = T)
par(mfrow = c(2,1), las = 1)
plot(gt, chr = 2, lod = 3:5, col = c("darkred","cyan", "darkblue"), type = "n",
     ylab = "genotype frequency", xlab = "Chr02 mapping position (cM)",
     main = "no segregation distortion on Chr02")
abline(h = colMeans(gt[,5:7]), col = c("darkred","cyan", "darkblue"), lty = 2)
plot(gt, chr = 2, lod = 3:5, col = c("darkred","cyan", "darkblue"), add = T)

plot(gt, chr = 9, lod = 3:5, col = c("darkred","cyan", "darkblue"), type = "n",
     ylab = "genotype frequency", xlab = "Chr09 mapping position (cM)",
     main = "Strong segregation distortion on Chr09")
abline(h = colMeans(gt[,5:7]), col = c("darkred","cyan", "darkblue"), lty = 2)
plot(gt, chr = 9, lod = 3:5, col = c("darkred","cyan", "darkblue"), add = T)

## ------------------------------------------------------------------------
el<-calc.errorlod(subset(cross5, chr = 8), error.prob=0.001, 
                  map.function="kosambi")
tel<- top.errorlod(el, cutoff = 2)

## ----echo = F, fig.cap = "position of crossover events (x) and the distribtution of each genotype on Chr08"----
kable(tel, caption = "The markers with an error LOD score >2 on Chr08")
badlibs<-which(as.character(getid(el)) %in% tel$id[tel$chr == 8])
par(mfrow = c(1,1))
plotGeno(el, ind=badlibs, cutoff=2, min.sep=2)

## ------------------------------------------------------------------------
cross.clean<-cross5
for(i in 1:nrow(tel)) {
  chr <- tel$chr[i]
  id <- tel$id[i]
  mar <- tel$marker[i]
  cross.clean$geno[[chr]]$data[cross5$pheno$id==id, mar] <- NA
}
cross6<-cross.clean
map1<-est.map(cross5, chr = 8, error.prob = 0.001, map.function = "kosambi")
map2<-est.map(cross6, chr = 8, error.prob = 0.001, map.function = "kosambi")
chrlen(map1) # original chr08 length
chrlen(map2) # new chr08 length

