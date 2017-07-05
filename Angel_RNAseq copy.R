setwd("/Users/apple/Desktop/UC_Davis/Lab/Michelmore/RLKgenes/")
path = "/Users/apple/Desktop/UC_Davis/Lab/Michelmore/RLKgenes/RNASEQ_RLK_DEF.csv" #edit this line to reflect where the file is on your computer
rna = read.csv(paste(path))
View(rna)
dim(rna)
rna[,1] = as.character(rna[,1])
rna[1,1] = "Water"
rna[2,1] = "Nitrogen"
rna[3,1] = "Tissue"
rna[4,1] = "Date"
rna[5,1] = "Rep"
des  = rna[1:5,]
#install.packages("tidyverse") 
des = as.tibble(t(des))
View(des)  
colnames(des) = des[1,]
des = des[-1,]
treatment = des %>%
  unite(Treatment, Water, Nitrogen, Tissue, Date, Rep, sep = "_")
View(treatment)
rna = as.tibble(rna)
t.rna = as.tibble(t(rna))
View(t.rna)
t.rna = t.rna[-1,]
t.rna[,6] = as.character(t.rna[,6])
t.rna[,6] = treatment
View(t.rna)
rna2 = as.tibble(t(t.rna))
View(rna2)
rna2 = data.frame(rna[,1],rna2)
rna2 = as.tibble(rna2)
colnames(rna2) = rna2[6,]
rna2 = rna2[-c(1:6),]
colnames(rna2) = gsub(" ",".",colnames(rna2))
View(rna2)
cleaned_rna <- rna2 %>%
  gather(before.water.treatments.commenced_before.nitrogen.treatments.commenced_Leaf_7.17.12_1:LOW_HIGH_Leaf_8.16.12_4, key = "Treatments", value = "Expression") %>%
  separate(Treatments, into = c("Water","Nitrogen","Tissue","Date","Rep"), sep = "_")
View(cleaned_rna)
cleaned_rna["Expression"] = as.numeric(unlist(cleaned_rna["Expression"]))

##########differential expression analysis###########
edge.rna = rna[-c(1:5),]
View(edge.rna)
t.edge.rna = as.data.frame(edge.rna)
t.edge.rna[,1] = as.character(unlist(t.edge.rna[,1]))
edge.rna = as.data.frame(t(t.edge.rna))
colnames(edge.rna) = as.character(unlist(edge.rna[1,]))
rownames(edge.rna) = as.character(unlist(edge.rna[,1]))
edge.rna = edge.rna[-1,-1]
t.edge.rna = t(edge.rna)
t.edge.rna = as.data.frame(unlist(apply(t.edge.rna, 2, function(x) as.numeric(paste(as.factor(x))))))
rownames(t.edge.rna) = colnames(edge.rna)
View(t.edge.rna)

## make a histogram of counts for each of the sample  
pdf("fpkm distribution of each sample.pdf")
for(i in 1:68) {
  hist(t.edge.rna[,i], breaks = 500, main = colnames(t.edge.rna)[i])
}
#log transformation
for(i in 1:68) {
  hist(log(t.edge.rna[,i]), breaks = 500, main = colnames(t.edge.rna)[i])
}
dev.off()

#creating sample desctiption file
des = rna[1:6,]
View(des)
des = as.data.frame(t(des))
des = des[,c(6,1,2,3,4,5)]
des[,2] = substring(des[,2],1,6)
des[,3] = substring(des[,3],1,6)
colnames(des) = c("Sample", "Water","Nitrogen","Tissue","Date","Rep")
des = des %>%
  unite(group, Water, Nitrogen, Tissue, Date, Rep, sep = "_", remove = FALSE)
des = des[-1,]
colnames(t.edge.rna) = des$group

library(edgeR)

dge.data <- DGEList(counts=t.edge.rna, group=des$group)
dim(dge.data) 
dge.data <- calcNormFactors(dge.data, method = "TMM")
dge.data$samples

#library(tidyverse)
WN  = des %>%
  unite(WN, Water, Nitrogen, sep = "_") %>%
  select(WN, Tissue, Date)
WN

pdf("BCVplot_by_Tissue_Date_and_Treatment.pdf")
plotMDS(dge.data, method = "bcv", col = as.numeric(des$Tissue), main = "Whole dataset")
plotMDS(dge.data[,des$Tissue == "Leaf"], method = "bcv", col = as.numeric(des$Date[des$Tissue=="Leaf"]), main = "Leaf tissue, colored by date of sampling")
plotMDS(dge.data[,des$Tissue == "Leaf"], method = "bcv", col = as.numeric(as.factor(WN$WN[WN$Tissue == "Leaf"])), main = "Leaf, colored by treatment")
plotMDS(dge.data[,des$Tissue == "Leaf" & des$Date == "7.29.12"], method = "bcv", col = as.numeric(as.factor(WN$WN[WN$Tissue == "Leaf" & WN$Date == "7.29.12"])), main = "Leaf, 7.29.12, colored by treatment")
plotMDS(dge.data[,des$Tissue == "Leaf" & des$Date == "8.1.12"], method = "bcv", col = as.numeric(as.factor(WN$WN[WN$Tissue == "Leaf" & WN$Date == "8.1.12"])), main = "Leaf, 8.1.12, colored by treatment")
plotMDS(dge.data[,des$Tissue == "Leaf" & des$Date == "8.16.12"], method = "bcv", col = as.numeric(as.factor(WN$WN[WN$Tissue == "Leaf" & WN$Date == "8.16.12"])), main = "Leaf, 8.16.12, colored by treatment")
plotMDS(dge.data[,des$Tissue == "Head"], method = "bcv", col = as.numeric(as.factor(WN$WN[WN$Tissue == "Head"])), main = "Head, 8.16.12, colored by treatment")
dev.off()

####From this point on we need to subset the data because the "before" treatment for water and nitrogen,
####the first and last dates, and tissue type are co-variant. This will result in a less-than-full-rank 
####design matrix, causing the analysis to terminate. 

####First of all, we are going to take out just the leaf tissue, and only at the time points after the
####water and nitrogen treatments have been implemented. 

#t.edge.rna = read.csv("expression_alone.csv", row.names = 1)
#des = read.csv("sample_description.csv", row.names = 1)
leaf.data = t.edge.rna[,des$Tissue == "Leaf"]
leaf.data = leaf.data[,substr(colnames(leaf.data),1,3) != "bef"]
dim(leaf.data)
View(leaf.data)
des.leaf = des[des$Tissue == "Leaf",]
des.leaf = des.leaf[des.leaf$Date != "7.17.12",]
dim(des.leaf)
View(des.leaf)
dge.leaf <- DGEList(counts=leaf.data, group=des.leaf$group)
dim(dge.leaf) 
dge.leaf <- calcNormFactors(dge.leaf, method = "TMM")
dge.leaf$samples

##Building design matrix
design.leaf  = model.matrix(~Water+Nitrogen+Date, data = des.leaf)
rownames(design.leaf)  = des.leaf$sample
design.leaf #somehow this design matrix is (as usual) not very correct. Need to correct manually.
dim(design.leaf)
design.leaf = design.leaf[,-c(2,4,6)] #this nees to be adjusted every time, becasue the design matrix generated is wrong differently every time...

##Fiding and plotting dispersion
#First the overall dispersion
dge.leaf <- estimateGLMCommonDisp(dge.leaf,design.leaf,verbose = TRUE)
#Then a trended dispersion based on count level
dge.leaf <- estimateGLMTrendedDisp(dge.leaf,design.leaf)
#And lastly we calculate the gene-wise dispersion, using the prior estimates to "squeeze" the dispersion towards the common dispersion.
dge.leaf <- estimateGLMTagwiseDisp(dge.leaf,design.leaf) 
#We can examine this with a plot
pdf("dispersion_leaf.pdf")
plotBCV(dge.leaf)
dev.off()

fit.leaf <- glmFit(dge.leaf, design.leaf)

#To find genes that are differentially expressed on 8.16.12
gt.lrt <- glmLRT(fit.leaf,coef = "Date8.16.12")
topTags(gt.lrt, n=242)
DE.leaf = decideTestsDGE(gt.lrt,p.value=0.01)
summary(DE.leaf)
DEgene.hvst <- topTags(gt.lrt,n = Inf)$table[topTags(gt.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.hvst)
write.csv(DEgene.hvst, "DEgenes.NH2O.harvest.csv")
#to find genes thar are differentially expressed due to low N treatment
ln.lrt <- glmLRT(fit.leaf,coef = "NitrogenLOW")
topTags(ln.lrt, n=242)
DE.leaf.ln = decideTestsDGE(ln.lrt,p.value=0.01)
summary(DE.leaf.ln)
DEgene.ln <- topTags(ln.lrt,n = Inf)$table[topTags(ln.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.ln)
write.csv(DEgenes.ln, "DEgenes.NH2O.LowN.csv")
#to find genes differentially expressed due to low water treatment
lw.lrt <- glmLRT(fit.leaf,coef = "WaterLOW")
topTags(lw.lrt, n=242)
DE.leaf.lw = decideTestsDGE(lw.lrt,p.value=0.01)
summary(DE.leaf.lw) #no upregulated genes 
DEgene.lw <- topTags(lw.lrt,n = Inf)$table[topTags(lw.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.lw)
write.csv(downreg.leaf.lw, "DEgenes.NH2O.LowH2O.csv")

#here is a plotting function for gene-treatment interaction
plotDE <- function(genes, dge, sample.description) {
  require(ggplot2)
  require(reshape2)
  tmp.data <- t(log2(cpm(dge[genes,])+1))
  tmp.data <- merge(tmp.data,sample.description,by.x="row.names",by.y="sample")
  tmp.data <- melt(tmp.data,value.name="log2_cpm",variable.name="gene")
  pl <- ggplot(tmp.data,aes(x=gt,y=log2_cpm,fill=trt))
  pl <- pl + facet_wrap( ~ gene)
  pl <- pl + ylab("log2(cpm)") + xlab("genotype")
  pl <- pl + geom_boxplot()
  pl + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
}
#The function didn't work for me lol.
#plotDE(rownames(DEgene.hvst)[1:9],dge.leaf,des.leaf)

#Now onto treatment interactions
design.leaf.inter <- model.matrix(~Water + Nitrogen + Date + Water*Nitrogen + Water*Date + Nitrogen*Date,data = des.leaf)
design.leaf.inter
design.leaf.inter = design.leaf.inter[,-c(2,4,6,9,10,11,13,14,15,17,19,20,21,23)] #adjust again
design.leaf.inter
rownames(design.leaf.inter) <- des.leaf[,1]

fit.inter = glmFit(dge.leaf, design.leaf.inter)

#to find genes that are differentially expressed due to water*nitrogen interaction:
wn.lrt = glmLRT(fit.inter, coef = "WaterLOW:NitrogenLOW")
topTags(wn.lrt)
summary(decideTestsDGE(wn.lrt,p.value=0.01))

#interaction of water and 8.16.12:

w816.lrt = glmLRT(fit.inter, coef = "WaterLOW:Date8.16.12")
topTags(w816.lrt)
summary(decideTestsDGE(w816.lrt,p.value=0.01))

#interaction of nitrogen and 8.16.12:

n816.lrt = glmLRT(fit.inter, coef = "NitrogenLOW:Date8.16.12")
topTags(n816.lrt)
summary(decideTestsDGE(n816.lrt,p.value=0.01))

#interaction of water and 8.1.12:

w81.lrt = glmLRT(fit.inter, coef = "WaterLOW:Date8.1.12")
topTags(w81.lrt)
summary(decideTestsDGE(w81.lrt,p.value=0.01))

#interaction of nitrogen and 8.1.12:

n81.lrt = glmLRT(fit.inter, coef = "NitrogenLOW:Date8.16.12")
topTags(n81.lrt)
summary(decideTestsDGE(n81.lrt,p.value=0.01))

#in conclsion, none of the RLK genes' expression level responded to treatment interactions or treatment*date interactions. 

####now we are going to look at the pair-wise expression changes between heads and leaves at harvest.

harvest.data = t.edge.rna[,des$Date == "8.16.12"]
dim(harvest.data)
View(harvest.data)
des.harvest = des[des$Date =="8.16.12",]
dim(des.harvest)
View(des.harvest)
dge.harvest <- DGEList(counts=harvest.data, group=des.harvest$group)
dim(dge.harvest) 
dge.harvest <- calcNormFactors(dge.harvest, method = "TMM")
dge.harvest$samples

##Building design matrix
design.harvest  = model.matrix(~Water+Nitrogen+Tissue, data = des.harvest)
rownames(design.harvest)  = des.harvest$sample
design.harvest #somehow this design matrix is (as usual) not very correct. Need to correct manually.
dim(design.harvest)
design.harvest = design.harvest[,-c(2,4)] #this nees to be adjusted every time, becasue the design matrix generated is wrong differently every time...

##Fiding and plotting dispersion
#First the overall dispersion
dge.harvest <- estimateGLMCommonDisp(dge.harvest,design.harvest,verbose = TRUE)
#Then a trended dispersion based on count level
dge.harvest <- estimateGLMTrendedDisp(dge.harvest,design.harvest)
#And lastly we calculate the gene-wise dispersion, using the prior estimates to "squeeze" the dispersion towards the common dispersion.
dge.harvest <- estimateGLMTagwiseDisp(dge.harvest,design.harvest) 
#We can examine this with a plot
pdf("dispersion_harvest.pdf")
plotBCV(dge.harvest)
dev.off()

##Find differentially expressed genes by tissue type at harvest
fit.harvest = glmFit(dge.harvest,design.harvest)

tissue.lrt = glmLRT(fit.harvest, coef = "TissueLeaf")
topTags(tissue.lrt)
summary(decideTestsDGE(tissue.lrt, p.value = 0.01))
DEgene.tissue <- topTags(tissue.lrt,n = Inf)$table[topTags(tissue.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.tissue)
write.csv(DEgene.tissue, "DEgenes.tissue.harvest.csv")

h2o.lrt = glmLRT(fit.harvest, coef = "WaterLOW")
topTags(h2o.lrt)
summary(decideTestsDGE(h2o.lrt, p.value = 0.01))
DEgene.h2o <- topTags(h2o.lrt,n = Inf)$table[topTags(h2o.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.h2o)
write.csv(DEgene.h2o, "DEgenes.water.harvest.csv")

n.lrt = glmLRT(fit.harvest, coef = "NitrogenLOW")
topTags(n.lrt)
summary(decideTestsDGE(n.lrt, p.value = 0.01))
DEgene.n <- topTags(n.lrt,n = Inf)$table[topTags(n.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.n)
write.csv(DEgene.n, "DEgenes.nitrogen.harvest.csv")

##interaction at harvest
design.harvest.inter <- model.matrix(~Water + Nitrogen + Tissue + Water*Nitrogen + Water*Tissue + Nitrogen*Tissue, data = des.harvest)
design.harvest.inter
design.harvest.inter = design.harvest.inter[,-c(2,4,7,8,9,11,13)]
fit.harvest.inter = glmFit(dge.harvest,design.harvest.inter)

wn.harv.lrt = glmLRT(fit.harvest.inter, coef = "WaterLOW:NitrogenLOW")
topTags(wn.harv.lrt)
summary(decideTestsDGE(wn.harv.lrt, p.value = 0.01))
DEgene.wn.harv <- topTags(wn.harv.lrt,n = Inf)$table[topTags(wn.harv.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.wn.harv)
write.csv(DEgene.wn.harv, "DEgenes.H2ON.harvest.csv")

wl.harv.lrt = glmLRT(fit.harvest.inter, coef = "WaterLOW:TissueLeaf")
topTags(wl.harv.lrt)
summary(decideTestsDGE(wl.harv.lrt, p.value = 0.01))
DEgene.wl.harv <- topTags(wl.harv.lrt,n = Inf)$table[topTags(wl.harv.lrt,n = Inf)$table$FDR<0.01,]
dim(DEgene.wl.harv)
write.csv(DEgene.wl.harv, "DEgenes.H2Oleaf.harvest.csv")

nl.harv.lrt = glmLRT(fit.harvest.inter, coef = "NitrogenLOW:TissueLeaf")
topTags(nl.harv.lrt)
summary(decideTestsDGE(nl.harv.lrt, p.value = 0.01))
#no interaction genes here. 


