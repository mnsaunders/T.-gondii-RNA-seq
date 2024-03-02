RHvsLPSinfect1<-filter(RHvsLPSinfect, padj<0.05)
MockvsLPSinfect1<-filter(MockvsLPSinfect, padj<0.05)

BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
library("AnnotationDbi")
BiocManager::install("clusterProfiler")
install.packages("Bioconductor")
library("BiocManager")
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
library("clusterProfiler")
library("tidyverse")
library("RColorBrewer")

remotes::install_github('YuLab-SMU/ggtree')
remotes::install_github('YuLab-SMU/gosemsim')
remotes::install_github("YuLab-SMU/enrichplot")
remotes::install_github("YuLab-SMU/clusterProfiler")

RH<-RHvsLPSinfect1 %>% remove_rownames %>% column_to_rownames(var="...1")
Mock<-MockvsLPSinfect1 %>% remove_rownames %>% column_to_rownames(var="...1")

RH1<-rownames(RH)
Mock1<-rownames(Mock)

GO_RH<-enrichGO(gene=RH1, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
GO_Mock<-enrichGO(gene=Mock1, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")

RH<-as.data.frame(GO_RH)
RHplot<-plot(barplot(GO_RH, showCategory = 20))
RHplot+labs(title = "RH vs LPS showing neutrophils")
Mockplot<-plot(barplot(GO_Mock, showCategory = 20))
Mockplot+labs(title = "Mock vs LPS showing neutrophils")+scale_color_continuous(aes='p.adjust', low="yellow", high="darkgreen")





LPS<-LPSplusinfections1 %>% remove_rownames %>% column_to_rownames(var="...1")

LPS1<-rownames(LPS)

GO_LPS<-enrichGO(gene=LPS1, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")

neutroplot<-plot(barplot(GO_LPS, showCategory = 20))
neutroplot+labs(title = "LPS vs RHLPS showing neutrophils")