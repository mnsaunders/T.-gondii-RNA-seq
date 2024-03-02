library('DESeq2')
library('tidyverse')
library("ggrepel")

data<-RHvsLPScom %>% remove_rownames %>% column_to_rownames(var="Geneid")
Counts<-filter(data, LPS_1 >0, LPS_2 >0,LPS_3 >0, LPS_4 >0)

condition <- factor(c("RH","RH","RH","RH","LPS","LPS","LPS","LPS"))
coldata<-data.frame(row.names=colnames(Counts), condition)

dds<-DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

deseq<-DESeq(dds)

vsdata<-vst(deseq, blind=FALSE)
pcaO<-plotPCA(vsdata, intgroup="condition")
pcaO + geom_text_repel(aes(label = name), size=2.5)+
  geom_point(size=3)+
  theme_bw()

plotDispEsts(deseq)
res<- results(deseq, contrast=c("condition", "RH", "LPS"))
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

ggplot(genedata_RHLPS, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexp))+
  geom_point()+ geom_text_repel(aes(label = gene), check_overlap = TRUE, size=3, fontface = "bold", show.legend = FALSE)+
  geom_vline(xintercept = c(0,0), col="grey")+ theme_bw()+
  scale_color_manual(values=c("red3","grey","navy"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size=13, family="roboto"), legend.title = element_blank(), legend.position="bottom")+
  labs(x="log2FoldChange", y="-log(pvalue)", title= "T. gondii Affects Gene Regulation in Neutrophils")

up<-filter(genedat, log2FoldChange> 6)
down<-filter(genedat, log2FoldChange< -6)
top<-rbind(up,down)
top<-filter(top, gene!="NA")


ggplot(top, aes(y=log2FoldChange, x=gene, fill=diffexp))+
  geom_bar(stat="identity")+ geom_hline(yintercept = c(0,0), col="grey")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),axis.text.x=element_text(angle=80,hjust=1,vjust=1, size = 5.5, face="bold.italic"))+
  labs(y="Gene Regulation", title= "RH vs LPS neutrophils")

