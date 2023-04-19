library(RColorBrewer)


pd <-read.AnnotatedDataFrame("SepticShock_ClassLabel_SS.csv",sep=",", header=T)

makePie <- function(values) {
  lev = unique(values)
  genericTable = table(values)
  percentLabels = round(100*genericTable/sum(genericTable), 1)
  pieLabels = paste(percentLabels, "%", lev, sep="")
  pie(genericTable, cex=0.5,col=rainbow(length(lev)), labels=pieLabels)
}

makePie(pd$Gender)

coul <- brewer.pal(8, "Set2") 
par(mar=c(2,18,1,1))
barplot(table(pd$Organism), horiz=TRUE, las=2, col=coul)
barplot(table(pd$Source.of.infection), horiz=TRUE, las=2, col=coul)
barplot(table(pd$Age.Category), horiz=TRUE, las=2, col=coul)
stripchart(unlist(tableSub[grep("^CDC20$",genes),]) ~ ssFactor, method="jitter", col=c("darkgreen","red"), vertical=TRUE, ylab="CDC20 expresion", xlab="Mortality")

top50 = topTable(fit3, number=50, genelist=genes, adjust.method = "BH", sort.by="p", fc=1.5)
biomarkers50 = as.matrix(tableSub[genes %in% top50$ID,])
rownames(biomarkers50) = top50$ID
