library(enrichR)
dge = topTable(fit3, number=50, genelist=genes, adjust.method = "BH", sort.by="P")
dgeList=dge$P.Value
names(dgeList)=dge$ID
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021","KEGG_2021_Human","Reactome_2022")
enriched <- enrichr(names(dgeList), dbs)
plotEnrich(enriched$GO_Molecular_Function_2015, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched$GO_Molecular_Function_2021, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched$GO_Cellular_Component_2021, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched$GO_Biological_Process_2021, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched$KEGG_2021_Human, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched$Reactome_2022, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

