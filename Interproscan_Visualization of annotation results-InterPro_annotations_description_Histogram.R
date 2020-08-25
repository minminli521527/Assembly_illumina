#### 1. ���ذ� ####
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(magrittr)
library(dplyr)
library(forcats)
library(ggplot2)
library(Cairo)

#### 2. ��ȡ����->����ipr���� ####
getwd()
setwd("C:/Users/99039/Desktop/interproscan")
#��interproscan����������test_proteins.fasta.tsv����Ϊcsv�ļ�test_proteins.fasta.csv
#test_proteins.fasta.tsv��15�У�ÿһ����ͷ����Ϊ��Protein_Accession	Sequence_MD5_digest	Sequence_length	Analysis	Signature_Accession	Signature_Description	Start_location	Stop_location	Score	Status	Date	InterPro_annotations_accession	InterPro_annotations_description	GO_annotations	Pathways_annotations
#���ո񲹳�ΪNA
interpro <- read.csv ("interproscan.csv",header = T)
#ȥ��Protein_Accession�к��пո����
interpro<-subset(interpro,Protein_Accession!="NA")
#ȥ��InterPro_annotations_accession�к��пո����
interpro<-subset(interpro,InterPro_annotations_accession!="NA")
#ȥ��InterPro_annotations_description�к��пո����
interpro<-subset(interpro,InterPro_annotations_description!="NA")

ipr = interpro %>% select(Protein_Accession,InterPro_annotations_accession,InterPro_annotations_description) %>% group_by(Protein_Accession,InterPro_annotations_accession) %>%
  summarise(InterPro_annotations_description=InterPro_annotations_description[1]) %>% group_by(InterPro_annotations_description) %>% summarise(Count=n())%>%
  arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))
ipr

#ipr = interpro %>% select(Protein_Accession,InterPro_annotations_accession,InterPro_annotations_description) %>% group_by(Protein_Accession,InterPro_annotations_accession) %>%
#summarise(InterPro_annotations_description=InterPro_annotations_description[1]) %>% group_by(InterPro_annotations_accession, InterPro_annotations_description) %>% summarise(Count=n())%>%
#arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))
#ipr

#### 3. ��ͼ->��Ϊpdf ####
#pdf("InterPro_annotations_description.pdf")
png_path="InterPro_annotations_description_Histogram.pdf"
#CairoPNG(png_path, width = 12, height = 7, units='in', dpi=600)
CairoPDF(png_path, width = 12, height = 7)
ggplot(ipr) +
  geom_bar(aes(x=InterPro_annotations_description,y=Count, fill=Count), stat='identity') + 
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