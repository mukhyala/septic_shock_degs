library(RColorBrewer)
library(ComplexHeatmap)
library(EnhancedVolcano)

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

dge = topTable(fit3, number=50, genelist=genes, adjust.method = "fdr", sort.by="P")
dTableSub = tableSub[dge$ID,]
mat = t(scale((dTableSub)))
rownames(mat)=NULL
row_ha = rowAnnotation(Mortality=pd$Mortality, Pathogen = pd$Organism.Yes.No, GramPositive=pd$GRAM.POS_yes1_no0,ImmunoSuppression=pd$Immunosuppression, Gender=pd$Gender, InfectionSource=pd$Source.of.infection)
Heatmap(mat,right_annotation = row_ha)

EnhancedVolcano::EnhancedVolcano(stats_df, 
                                 lab = stats_df$ID,
                                 x = "logFC",
                                 y = "adj.P.Val",
                                 pCutoff = 0.05, 
                                 FCcutoff=0.5, 
                                 ylim=c(0,3.5), 
                                 xlim=c(-1.5,1.5)
                                 )
