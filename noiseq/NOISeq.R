### R code from vignette source 'NOISeq.Rnw'
### Encoding: UTF-8
rm(list=ls())

###################################################
### code chunk number 1: options
###################################################
options(digits=3, width=95)


###################################################
### code chunk number 2: data
###################################################
library(NOISeq)
data(Marioni)


###################################################
### code chunk number 3: NOISeq.Rnw:88-89
###################################################
head(mycounts)


###################################################
### code chunk number 4: factors
###################################################
myfactors = data.frame(Tissue=c("Kidney","Liver","Kidney","Liver","Liver","Kidney","Liver",
                                "Kidney","Liver","Kidney"),
                       TissueRun = c("Kidney_1","Liver_1","Kidney_1","Liver_1","Liver_1",
                                     "Kidney_1","Liver_1","Kidney_2","Liver_2","Kidney_2"),
                       Run = c(rep("R1", 7), rep("R2", 3)))
myfactors


###################################################
### code chunk number 5: NOISeq.Rnw:119-123
###################################################
head(mylength)
head(mygc)
head(mybiotypes)
head(mychroms)


###################################################
### code chunk number 6: readData
###################################################
mydata <- readData(data=mycounts,length=mylength, gc=mygc, biotype=mybiotypes,
                   chromosome=mychroms, factors=myfactors)
mydata


###################################################
### code chunk number 7: NOISeq.Rnw:152-156
###################################################
str(mydata)
head(assayData(mydata)$exprs)
head(pData(mydata))
head(featureData(mydata)@data)


###################################################
### code chunk number 8: readData2
###################################################
mydata <- readData(data=mycounts,chromosome=mychroms, factors=myfactors)


###################################################
### code chunk number 9: readData3
###################################################
mydata <- addData(mydata, length=mylength, biotype=mybiotypes, gc = mygc)


###################################################
### code chunk number 10: dat
###################################################
myexplodata <- dat(mydata, type = "biodetection")
explo.plot(myexplodata, plottype = "persample")


###################################################
### code chunk number 11: nicedata
###################################################
mynicedata <- dat2save(myexplodata)


###################################################
### code chunk number 12: fig_biodetection
###################################################
mybiodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
explo.plot(mybiodetection, samples=c(1,2), plottype = "persample")


###################################################
### code chunk number 13: fig_biodetection2
###################################################
par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
explo.plot(mybiodetection, samples=c(1,2), toplot = "protein_coding", plottype = "comparison")


###################################################
### code chunk number 14: fig_boxplot1
###################################################
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")


###################################################
### code chunk number 15: fig_sat1
###################################################
mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:2, yleftlim = NULL, yrightlim = NULL)


###################################################
### code chunk number 16: fig_sat2
###################################################
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:4)


###################################################
### code chunk number 17: fig_boxplot2
###################################################
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")


###################################################
### code chunk number 18: fig_boxplot3
###################################################
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")


###################################################
### code chunk number 19: fig_length
###################################################
mylengthbias = dat(mydata, factor = "Tissue", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")


###################################################
### code chunk number 20: showmodels
###################################################
show(mylengthbias)


###################################################
### code chunk number 21: fig_GC
###################################################
myGCbias = dat(mydata, factor = "Tissue", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")


###################################################
### code chunk number 22: fig_countdistr
###################################################
mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd)


###################################################
### code chunk number 23: randomBatchEffect
###################################################
set.seed(123)
mycounts2 = mycounts
mycounts2[,1:4] = mycounts2[,1:4] + runif(nrow(mycounts2)*4, 3, 5)
myfactors = data.frame(myfactors, "batch" = c(rep(1,4), rep(2,6)))
mydata2 = readData(mycounts2, factors = myfactors)


###################################################
### code chunk number 24: fig_PCA
###################################################
myPCA = dat(mydata2, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Tissue")
explo.plot(myPCA, factor = "batch")


###################################################
### code chunk number 25: QCreportExample
###################################################
QCreport(mydata, samples = NULL, factor = "Tissue", norm = FALSE)


###################################################
### code chunk number 26: normalization
###################################################
myRPKM = rpkm(assayData(mydata)$exprs, long = mylength, k = 0, lc = 1)
myUQUA = uqua(assayData(mydata)$exprs, long = mylength, lc = 0.5, k = 0)
myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)
head(myRPKM[,1:4])


###################################################
### code chunk number 27: filtering
###################################################
myfilt = filtered.data(mycounts, factor = myfactors$Tissue, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")


###################################################
### code chunk number 28: fig_knownBatch
###################################################
mydata2corr1 = ARSyNseq(mydata2, factor = "batch", batch = TRUE, norm = "rpkm",  logtransf = FALSE)
myPCA = dat(mydata2corr1, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Tissue")
explo.plot(myPCA, factor = "batch")


###################################################
### code chunk number 29: fig_unknownBatch
###################################################
mydata2corr2 = ARSyNseq(mydata2, factor = "Tissue", batch = FALSE, norm = "rpkm",  logtransf = FALSE)
myPCA = dat(mydata2corr2, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Tissue")
explo.plot(myPCA, factor = "batch")


###################################################
### code chunk number 30: results
###################################################
mynoiseq = noiseq(mydata, k = 0.5, norm = "rpkm", factor="Tissue", pnr = 0.2, 
                  nss = 5, v = 0.02, lc = 1, replicates = "technical")
head(mynoiseq@results[[1]])


###################################################
### code chunk number 31: NOISeq.Rnw:797-799
###################################################
mynoiseq.tmm = noiseq(mydata, k = 0.5, norm = "tmm", factor="TissueRun", 
                      conditions = c("Kidney_1","Liver_1"), lc = 0, replicates = "technical")


###################################################
### code chunk number 32: NOISeq.Rnw:821-823
###################################################
myresults <- noiseq(mydata, factor = "Tissue", k = NULL, norm="n", pnr = 0.2, 
                    nss = 5, v = 0.02, lc = 1, replicates = "no")


###################################################
### code chunk number 33: NOISeq.Rnw:875-877
###################################################
mynoiseqbio = noiseqbio(mydata, k = 0.5, norm = "rpkm", factor="Tissue", lc = 1, r = 20, adj = 1.5, plot = FALSE,
                        a0per = 0.9, random.seed = 12345, filter = 2)


###################################################
### code chunk number 34: NOISeq.Rnw:922-923
###################################################
head(mynoiseq@results[[1]])


###################################################
### code chunk number 35: NOISeq.Rnw:943-946
###################################################
mynoiseq.deg = degenes(mynoiseq, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(mynoiseq, q = 0.8, M = "down")


###################################################
### code chunk number 36: fig_summ_expr
###################################################
DE.plot(mynoiseq, q = 0.9, graphic = "expr", log.scale = TRUE)


###################################################
### code chunk number 37: fig_summ_MD
###################################################
DE.plot(mynoiseq, q = 0.8, graphic = "MD")


###################################################
### code chunk number 38: fig_manhattan
###################################################
DE.plot(mynoiseq, chromosomes = c(1,2), log.scale = TRUE,
        join = FALSE, q = 0.8, graphic = "chrom")


###################################################
### code chunk number 39: fig_distrDEG
###################################################
DE.plot(mynoiseq, chromosomes = NULL, q = 0.8, graphic = "distr")


###################################################
### code chunk number 40: session
###################################################
sessionInfo()

