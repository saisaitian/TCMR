library(dplyr)
data("data_logFC")

load("D:/github/TCMR/data/TCGA_BRCA_signature.Rdata")

idx <- intersect(res_RNA2$row,rownames(data_logFC))

res3 <- res_RNA2[res_RNA2$row%in%idx,]

DERNA=res3%>%subset(abs(log2FoldChange)>2&padj<0.05)

data_logFC2 <- data_logFC[idx,]

drug1 <- data_logFC2$Puerarin
drug2 <- data_logFC2$Anhydroicaritin

names(drug1) <- rownames(data_logFC2)
names(drug2) <- rownames(data_logFC2)

down_gene1 <- drug1[order(-drug1)]%>%
  tail(500)%>%
  names()

up_gene1 <- drug1[order(-drug1)]%>%
  head(500)%>%
  names()


down_gene2 <- drug2[order(-drug2)]%>%
  tail(500)%>%
  names()

up_gene2 <- drug2[order(-drug2)]%>%
  head(500)%>%
  names()

# disease_siguature ----------------------------------------------------------


down_disease <- DERNA%>%
  filter(log2FoldChange<0)%>%
  select(row)

up_disease <- DERNA%>%
  filter(log2FoldChange>0)%>%
  select(row)



id1 <- length(Reduce(intersect,list(up_disease$row,down_gene1,down_gene2)))

id2 <- length(Reduce(intersect,list(down_disease$row,up_gene1,up_gene2)))

id3 <- length(Reduce(intersect,list(up_disease$row,up_gene1,up_gene2)))

id4 <- length(Reduce(intersect,list(down_disease$row,down_gene1,down_gene2)))


synergyscore <- id1+id2-id3-id4

