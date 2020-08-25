#### 1. 加载包 ####
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(magrittr)
library(dplyr)
library(forcats)
library(ggplot2)
library(Cairo)

#### 2. 读取数据->生成ipr矩阵 ####
getwd()
setwd("C:/Users/99039/Desktop/interproscan")
#将interproscan软件输出结果test_proteins.fasta.tsv保存为csv文件test_proteins.fasta.csv
#test_proteins.fasta.tsv共15列，每一列题头依次为：Protein_Accession	Sequence_MD5_digest	Sequence_length	Analysis	Signature_Accession	Signature_Description	Start_location	Stop_location	Score	Status	Date	InterPro_annotations_accession	InterPro_annotations_description	GO_annotations	Pathways_annotations
#将空格补充为NA
interpro <- read.csv ("test_proteins.fasta.csv",header = T)
#去掉Protein_Accession列含有空格的行
interpro<-subset(interpro,Protein_Accession!="NA")
#去掉InterPro_annotations_accession列含有空格的行
interpro<-subset(interpro,InterPro_annotations_accession!="NA")
#去掉InterPro_annotations_description列含有空格的行
interpro<-subset(interpro,GO_annotations!="NA")

ipr = interpro %>% select(Protein_Accession,InterPro_annotations_accession,GO_annotations) %>% group_by(Protein_Accession,InterPro_annotations_accession) %>%
  summarise(GO_annotations=GO_annotations[1]) %>% group_by(GO_annotations) %>% summarise(Count=n())%>%
  arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))
ipr


#### 3. 作图->存为pdf ####
png_path="GO_annotations_Histogram.pdf"
CairoPDF(png_path, width = 12, height = 7)
ggplot(ipr) +
  geom_bar(aes(x=GO_annotations,y=Count, fill=Count), stat='identity') + 
  coord_flip() +
  scale_fill_gradient(expression(Count),low="red", high = "blue") +
  xlab("") +
  ylab("Gene count") +
  scale_y_continuous(expand=c(0, 0))+
  theme(
    axis.text.x=element_text(color="black",size=rel(1.5)),
    axis.text.y=element_text(color="black", size=rel(1.6)),
    axis.title.x = element_text(color="black", size=rel(1.6)),
    legend.text=element_text(color="black",size=rel(1.0)),
    legend.title = element_text(color="black",size=rel(1.1))
    # legend.position=c(0,1),legend.justification=c(-1,0)
    # legend.position="top",
  )
dev.off()
