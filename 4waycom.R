library('DESeq2')
library('tidyverse')
library("ggrepel")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("apeglm")
library("apeglm")
library("pheatmap")
library("dplyr")



data<-X4wayLPSRHMOCKRHLPS %>% remove_rownames %>% column_to_rownames(var="Geneid")
Counts<-filter(data, LPS_1 >0, LPS_2 >0,LPS_3 >0, LPS_4 >0)

condition <- factor(c(rep("RHLPS",4),rep("RH",4),rep("Mock",4),rep("LPS",4)))
coldata<-data.frame(row.names=colnames(Counts), condition)

dds<-DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

deseq<-DESeq(dds)
RM<-resultsNames(deseq)
results_names<-RM[-1]
results_names
res<-results(deseq)
resM <- results(deseq, name = "condition_Mock_vs_LPS")
resR <- results(deseq, name = "condition_RH_vs_LPS")
resRH <- results(deseq, name = "condition_RHLPS_vs_LPS")

resRD<-as.data.frame(resR)%>%rownames_to_column(var="gene_id")
resRD$condition<-"condition_RH_vs_LPS"
resMD<-as.data.frame(resM)%>%rownames_to_column(var="gene_id")
resMD$condition<-"condition_Mock_vs_LPS"
resRHD<-as.data.frame(resRH)%>%rownames_to_column(var="gene_id")
resRHD$condition<-"condition_RHLPS_vs_LPS"



comdata<-resRD %>%
  left_join(resRHD, by="gene_id") %>% filter( padj.x<1e-15, padj.y<1e-15) 
comdataup<-subset(res, padj<1e-40)
comdataup<-subset(comdataup, log2FoldChange>1.5)
comdataup<-as.data.frame(comdataup)
comdataup<-rownames_to_column(comdataup)

comdatadown<-subset(res, padj<1e-10)
comdatadown<-subset(comdatadown, log2FoldChange<1.5)
comdatadown<-as.data.frame(comdatadown)
comdatadown<-rownames_to_column(comdatadown)

RH<-RHvsLPSinfect %>% column_to_rownames(var="...1")
Mock<-MockvsLPSinfect %>% column_to_rownames(var="...1")

RH1<-rownames(RH)
Mock1<-rownames(Mock)

GO_RH<-enrichGO(gene=RH1, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
GO_Mock<-enrichGO(gene=Mock1, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")

RHdataframe<-as.data.frame(GO_RH)
RHplot<-plot(barplot(GO_RH, showCategory = 20))
RHplot+labs(title = "RH vs LPS showing neutrophils")
Mockplot<-plot(barplot(GO_Mock, showCategory = 20))
Mockplot+labs(title = "Mock vs LPS showing neutrophils")


comdata<-comdata%>%
  dplyr::select(gene_id, log2FoldChange.x,log2FoldChange.y)
colnames(comdata)[colnames(comdata) == 'log2FoldChange.x'] <- 'T. gondii verses LPS Affects in Infected Neutrophils' 
colnames(comdata)[colnames(comdata) == 'log2FoldChange.y'] <- 'T. gondii Affects in Neutrophils When Combined with LPS' 
rownames(comdata) <- NULL
comdata<- comdata%>%column_to_rownames('gene_id')
pheatmap(comdata, fontsize = 8, angle_col = 45)


normalized_counts<-counts(deseq, normalized=T)
Transfromed_counts<-log10(normalized_counts+1)
df<-as.data.frame(Transfromed_counts)
df<-rownames_to_column(df)
up<-df%>%left_join(comdataup, by='rowname')
up<-na.omit(up)
up<-up[,-18:-23]
up<-remove_rownames(up)
up<-column_to_rownames(up)

pheatmap(up, fontsize = 5)

down<-df%>%left_join(comdatadown, by='rowname')
down<-na.omit(down)
down<-down[,-18:-23]
down<-remove_rownames(down)
down<-column_to_rownames(down)

pheatmap(down, fontsize = 5)



vsdata<-vst(deseq, blind=FALSE)
pcaO<-plotPCA(vsdata, intgroup="condition")
pcaO + geom_text_repel(aes(label = name), size=3.7)+
  geom_point(size=3)+
  scale_colour_manual(values = c('#B7D2DF','lightblue3','lightblue4','#232323'))+ theme_bw()+
  theme(legend.position = "none", text = element_text(size=14, family="roboto"))
  

plotDispEsts(deseq)
res<- results(deseq, contrast=c("condition", "RH", "LPS", "Mock"))
res
nona<-na.omit(res)
nona
nona$diffexp<-'insignificant'
nona$diffexp[nona$log2FoldChange>0 & nona$padj<0.05] <-'up'
nona$diffexp[nona$log2FoldChange<0 & nona$padj<0.05] <-'down'
write.csv(nona, "RHvsLPSinfect.csv")
RHvsLPSinfect1<-filter(RHvsLPSinfect, padj<0.05)
names(RHvsLPSinfect1)[names(RHvsLPSinfect1) == '...1'] <- 'gene_id'

genedata_RHLPS<-RHvsLPSinfect1 %>%
  left_join(Human_annotation_table, by = "gene_id")

ggplot(genedat, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexp))+
  geom_point()+ geom_text_repel(aes(label = gene), check_overlap = TRUE, size=3, fontface = "bold")+
  geom_vline(xintercept = c(0,0), col="grey")+ theme_bw()+
  scale_color_manual(values=c("red3","grey","green3"))+
  labs(x="Gene Regulation", title= "RH vs LPS neutrophils")

up<-filter(genedat, log2FoldChange> 6)
down<-filter(genedat, log2FoldChange< -6)
top<-rbind(up,down)
top<-filter(top, gene!="NA")


ggplot(top, aes(y=log2FoldChange, x=gene, fill=diffexp))+
  geom_bar(stat="identity")+ geom_hline(yintercept = c(0,0), col="grey")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),axis.text.x=element_text(angle=80,hjust=1,vjust=1, size = 5.5, face="bold.italic"))+
  labs(y="Gene Regulation", title= "RH vs LPS neutrophils")


tx_lengths <- read.csv("mart_export1.txt")
tx_lengths %>% group_by(Gene.name) %>%
  summarise(mean_tx_len = mean(tx_lengths$Transcript.length..including.UTRs.and.CDS.)) %>%
              ungroup() -> tx_lengths


fc4tpm <- left_join(X4wayLPSRHMOCKRHLPS, tx_lengths, by = c("Geneid" = "Gene.name")) %>% na.omit()
mean_tx_len<-fc4tpm$mean_tx_len
rate <- Counts / mean_tx_len
tpm<-rate / sum(rate) * 1e6
logtpm<-log10(tpm)
logtpm[sapply(logtpm, is.infinite)] <- NA
na.omit(logtpm)

logsub<-subset(logtpm, !(RH_2 > -2.3 & RH_2 < 2))
pheatmap(logsub, fontsize = 8, angle_col = 45)

