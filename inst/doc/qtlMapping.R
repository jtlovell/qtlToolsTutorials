## ---- include=FALSE, message=FALSE---------------------------------------
suppressPackageStartupMessages(library(knitr))
knitr::opts_chunk$set(echo = TRUE, 
                      message=FALSE, 
                      error=FALSE, 
                      warning = FALSE, 
                      fig.width = 5.5, fig.height = 3.5)

library(devtools)
install_github("jtlovell/qtlTools")
install_github("jtlovell/qtlToolsTutorials")

suppressPackageStartupMessages(library(qtlTools))
suppressPackageStartupMessages(library(qtlToolsTutorials))
suppressPackageStartupMessages(library(qtlDesign))
suppressPackageStartupMessages(library(ggplot2))

## ------------------------------------------------------------------------
data("completeQTL_tutorial_data")

## ------------------------------------------------------------------------
summary(cr)

## ------------------------------------------------------------------------
phe.mat<-pull.pheno(cr)
geno.mat<-pull.geno(cr)
map<-pullMap(cr)

## ------------------------------------------------------------------------
colnames(geno.mat)<-paste0(map$chr,"_", round(map$pos,4))
rownames(geno.mat)<-getid(cr)

## ---- echo = F-----------------------------------------------------------
kable(geno.mat[1:5,1:5], row.names = T, caption = "The first 5 markers and F2 lines in the genotype matrix. Genotypes are coded as 1 = A/A, 2 = A/B, 3 = B/B.")

## ------------------------------------------------------------------------
geno2cross(geno.mat, crossfile = "~/Downloads/cross4tut.csv")

cross<-read.cross("csv", file="~/Downloads/cross4tut.csv",
                  genotypes=c(1,2,3), crosstype="f2")

## ---- fig.cap = "The scale of missing data in the cross. Note that the level of entropy in our cross does not exceed 20%. Not too shaby."----
plot.info(cross)

## ------------------------------------------------------------------------
cross<-sim.geno(cross, step = 1, stepwidth = "max", n.draws = 64,
                map.function = "kosambi", # use for plants
                error.prob = 0.0001) # set this for your genotyping platform. 

## ------------------------------------------------------------------------
cross<-calc.genoprob(cross, step = 1, stepwidth = "max", 
                map.function = "kosambi", # use for plants
                error.prob = 0.0001) # set this for your genotyping platform. 

## ------------------------------------------------------------------------
identical(phe.mat$id, as.numeric(getid(cross)))

## ------------------------------------------------------------------------
cross$pheno<-phe.mat

## ---- fig.cap = "comparison of one-way scans using imputations (black) and genotype probabilites (blue)"----
s1.imp<-scanone(cross, method = "imp", pheno.col = "phenotype")
s1.hk<-scanone(cross, method = "hk", pheno.col = "phenotype")
plot(s1.imp, s1.hk, lty = c(1,3), col = c("black", "lightblue"))

## ---- fig.cap = "comparison of one-way scans for just Chr09 using imputations (black) and genotype probabilites (blue)"----
plot(s1.imp, s1.hk, lty = c(1,1), col = c("black", rgb(0,0,1,.5)), chr = 9)

## ---- fig.cap = "permutation distribution and scanone with permutation threshold", fig.height = 5----
par(mfrow = c(2,1))
s1.perm<-scanone(cross, method = "hk", pheno.col = "phenotype", n.perm = 100)
plot(s1.perm)
plot(s1.hk)
add.threshold(s1.hk, perms = s1.perm, col = "red", lty = 2)
summary(s1.hk, perms = s1.perm, alpha = 0.05, pvalues = T)

## ------------------------------------------------------------------------
covar = data.frame(covar = as.numeric(as.factor(pull.pheno(cross, "Treatment"))))

## ---- fig.cap = "scanone plots for GxE (black), additive covariate effect (red) or no covariate (blue)"----
s1.no.covar<-scanone(cross, method = "hk", pheno.col = "phenotype")
s1.add.covar<-scanone(cross, method = "hk", pheno.col = "phenotype", addcovar = covar)
s1.int.covar<-scanone(cross, method = "hk", pheno.col = "phenotype", addcovar = covar, intcovar = covar)
plot(s1.int.covar, s1.add.covar,s1.no.covar)

## ---- fig.cap = "interaction between treatment and genotype. Note that the genotype effect only slightly depends on the treatment. "----
set.seed(42)
perm.no.covar<-scanone(cross, method = "hk", pheno.col = "phenotype", 
                       n.perm = 100, perm.strata = covar[,1], addcovar = NULL, verbose = F)
set.seed(42)
perm.add.covar<-scanone(cross, method = "hk", pheno.col = "phenotype", 
                       n.perm = 100, perm.strata = covar[,1], addcovar = covar, verbose = F)
set.seed(42)
perm.int.covar<-scanone(cross, method = "hk", pheno.col = "phenotype", 
                       n.perm = 100, perm.strata = covar[,1],
                       addcovar = covar, intcovar = covar, verbose = F)
kable(summary(s1.add.covar-s1.no.covar, 
        perms = perm.add.covar-perm.no.covar,
        pvalues = T), caption = "scanone results for the difference between no covariate andadditive covariate. Note that the Chr03 QTL is significant.")

kable(summary(s1.int.covar-s1.add.covar, 
        perms = perm.int.covar-perm.add.covar,
        pvalues = T), caption = "scanone results for the difference between interactive covariate and additive covariate. Note that the Chr03 QTL marginally significant.")

effectplot(cross, pheno.col = "phenotype", mname1 = "3_1.0175",
           mark2 = covar[,1], geno2 = c("dry","wet",""))

## ---- fig.cap = "scan looking for QTL in addition to the one on the top of Chr03"----
qtl = makeqtl(cross, chr = 3, pos = 11, what = "prob")
scan4secondqtl_add<-addqtl(cross, pheno.col = "phenotype", covar = covar, 
                           formula = "y~Q1+Q1*covar+Q2+covar", qtl = qtl,
                           method = "hk")
scan4secondqtl_int<-addqtl(cross, pheno.col = "phenotype", covar = covar, 
                           formula = "y~Q1+Q1*covar+Q2*covar+covar", qtl = qtl,
                           method = "hk")
plot(scan4secondqtl_int, scan4secondqtl_add)

