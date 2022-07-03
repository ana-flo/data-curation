#This script curates the annotation from datasets taken from the Tirosh collection that is available online here:www.weizmann.ac.il/sites/3CA 
#They have already run infer CNV to call malignant cell and they call cell types. The curation renames the cells according to the ontology and 
#and add patient and sample information
#also have to add a sample barcode list for the 10x format

library(data.table)
library(tidyverse)

#curate annotation from Olalelakn et al 2021
#link to publication : annotation in Table 1 https://www.sciencedirect.com/science/article/pii/S2211124721005076#tbl1
#neoadjuvant therapy is not specified
#sequencing is dropseq
#GSE147082

df.annotation <- read.table("H:/data/10x datasets/Tirosh - collection/Olalekan et al 2021 ovarian/Cells.txt", header=TRUE, sep=" ")
colnames(df.annotation) <- c("CellID", "PatientID","SampleID" ,"cell.type_no.tumor", "malignant")
df.annotation$cell.type_no.tumor <- recode(df.annotation$cell.type_no.tumor, "B_cell"="B cell", "DendriticCell_Plasmacytoid"="plasmacytoid dendritic cell", "Endothelial"="endothelial cell", 
                                          "Epithelial"="epithelial cell", "ESC"="embryonic stem cell", "Fibroblast"="fibroblast", "Monocyte"="monocyte","MSC"="mesenchymal stem cell","Plasma"=
                                            "plasma cell", "T_cell"="T cell", "Undecide1_astro"="unclassified", "Undecided2_CMP_BoneMarrow"="unlcassified")
df.annotation$cell.type <- df.annotation$cell.type_no.tumor
df.annotation$Dataset_accession <- "GSE147082"
df.annotation$Dataset_processing <- "Tirosh collection"

for (i in 1: nrow(df.annotation)){
  if (df.annotation$malignant[i] == "yes" ) df.annotation$cell.type[i] <- "tumor cell"
}

df.annotation$Sample_tissue_of_origin <- "omentum"
df.annotation$Sample_type <- "metastasis"
df.annotation$Cell_sroting <- "none"
df.annotation$Dissociation_temperatue <- "37C"
df.annotation$Sequencing_platform <- "DROPseq"
df.annotation$Indication <- df.annotation$PatientID
df.annotation$Treatment <- df.annotation$PatientID
df.annotation$Age <- df.annotation$PatientID
df.annotation$Race <- df.annotation$PatientID
df.annotation$Stage <- df.annotation$PatientID
df.annotation$Stage_PMN <- df.annotation$PatientID
df.annotation$Histological_grade <- df.annotation$PatientID
df.annotation$Tumor_origin <- df.annotation$PatientID

df.annotation$Indication <- recode(df.annotation$Indication,"PT-1"="serous ovarian cancer" ,"PT-2"="high-grade serous ovarian cancer","PT-3"="high-grade serous ovarian cancer","PT-4"="high-grade serous ovarian cancer"
                                   ,"PT-5"="high-grade serous ovarian cancer","PT-6"="uterine malignant mixed mullerian tumor")

df.annotation$Treatment <- recode(df.annotation$Treatment,"PT-1"="neoadjuvant therapy" ,"PT-2"="none","PT-3"="neoadjuvant therapy","PT-4"="none"
                                   ,"PT-5"="neoadjuvant therapy","PT-6"="neoadjuvant therapy")
df.annotation$Sample_timepoint <- df.annotation$Treatment
df.annotation$Sample_timepoint <- recode(df.annotation$Sample_timepoint, "none"="baseline", "neoadjuvant therapy"="post")

df.annotation$Age <- recode(df.annotation$Age,"PT-1"="62" ,"PT-2"="56","PT-3"="66","PT-4"="46" ,"PT-5"="71","PT-6"="66")
df.annotation$Race <- recode(df.annotation$Race,"PT-1"="white" ,"PT-2"="white","PT-3"="black","PT-4"="asian" ,"PT-5"="black","PT-6"="asian")
df.annotation$Stage <- recode(df.annotation$Stage,"PT-1"="IVb" ,"PT-2"="IIIc","PT-3"="IIIc","PT-4"="IIIc" ,"PT-5"="IIIc","PT-6"="IIIc")
df.annotation$Stage_PMN <- recode(df.annotation$Stage_PMN,"PT-1"="ypT3a Nx M1" ,"PT-2"="pT3c Nx Mx","PT-3"="ypT3c N1a","PT-4"="pT3c Nx" ,"PT-5"="pT3c N1 M1","PT-6"="ypT3c Nx")
df.annotation$Histological_grade <- recode(df.annotation$Histological_grade,"PT-1"="not applicable" ,"PT-2"="high grade","PT-3"="high grade","PT-4"="high grade" ,"PT-5"="high grade","PT-6"="high grade")
df.annotation$Tumor_origin <- recode(df.annotation$Tumor_origin,"PT-1"="unknown" ,"PT-2"="left ovary","PT-3"="left fallopian tube (STIC)","PT-4"="left fallopian tube (STIC)" ,"PT-5"="left fallopian tube (STIC)","PT-6"="fallopian tube")

write.table(df.annotation, "H:/data/10x datasets/Tirosh - collection/Olalekan et al 2021 ovarian/MetaData-GSE147082.txt", sep="\t", row.names = FALSE, quote = FALSE)