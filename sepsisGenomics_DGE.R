#Author Akram Mohammed
#November 2018
#Kamaleswaran Lab

#Set working directory for download

wd = "/Users/mukhyala/gmu_phd/classes/702_BiologicalDataAnalysis/septic_shock_degs"
setwd(paste(wd,"GSE66099","data",sep=.Platform$file.sep))
library(affy)
library(limma)

#data1 <- read.csv("groupbyGene_SepticShock_2.csv", comment.char="#")
pd <-read.AnnotatedDataFrame("SepticShock_ClassLabel_SS.csv",sep=",", header=T)
#table <-read.csv("groupbyGene_SepticShock_2_Transpose.csv",header=T)
table=read.csv('rma2.txt',sep="\t")
#tableSub <- subset(table, select = c(X0:X180))
grptab = aggregate(. ~ SYMBOL, data = table[ -c(1,2,3)], FUN = mean)
tableSub = grptab[,rownames(pd)]

#genes <- table$index
genes= grptab$SYMBOL

#design <- model.matrix(~ 0+factor(c(data1$Class)))
design <- model.matrix(~ 0+factor(c(pd$Mortality)))
colnames(design) <- c("S", "NS")
contrast.matrix <- makeContrasts(diff=S-NS, levels=design)

fit <- lmFit(tableSub, design = design) # What should the design be?
fit2 <- contrasts.fit(fit, contrast.matrix)
options(digits=3)
fit3 <- eBayes(fit2)

writefile = topTable(fit3, number=Inf, genelist=genes, adjust.method = "BH", sort.by="logFC")

write.csv(writefile, file="DEG_Septic_Shock_2_FDR.csv")
