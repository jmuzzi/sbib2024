#####################################################
############### CURSO R - SBIB ######################
######## Por Jessica Magno e Joao Muzzi #############

## PRATICA: Manipulacao de dados

#-- Definindo o diretorio
setwd("Documentos/SBIB_2024/")

#-- Download dos dados genomicos
url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-OV.htseq_counts.tsv.gz"
download.file(url = url, destfile = "tcga_ov/gexp.txt")

gexp_ov <- read.delim("tcga_ov/gexp.txt", dec = ".")
head(gexp_ov)[1:5]

#-- Download dos dados clinicos (metadados)
url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-OV.GDC_phenotype.tsv.gz"
download.file(url = url, destfile = "tcga_ov/clinical.txt")

clinical_ov <- read.delim("tcga_ov/clinical.txt")
head(clinical_ov)[1:5]

#-- Download dos dados de sobrevida (survival)
url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-OV.survival.tsv"
download.file(url = url, destfile = "tcga_ov/survival.txt")

survival_ov <- read.delim("tcga_ov/survival.txt")
head(survival_ov)

rm(url)

## Ajustando os nomes das amostras

# Manipulando os nomes de genes na matriz de contagem
rownames(gexp_ov)
rownames(gexp_ov) <- gexp_ov$Ensembl_ID
head(gexp_ov)[1:5]
identical(rownames(gexp_ov), gexp_ov$Ensembl_ID) #TRUE
gexp_ov <- gexp_ov[,-1]
head(gexp_ov)[1:5]

ens_gexp <- rownames(gexp_ov)
ens_gexp <- sub("\\..*", "", ens_gexp)
rownames(gexp_ov) <- ens_gexp
head(gexp_ov)[1:5]

# Manipulando os nomes das amostras na matriz de contagem
colnames(gexp_ov)
amostras_gexp <- colnames(gexp_ov)
amostras_gexp <- gsub("\\.", "-", amostras_gexp)
colnames(gexp_ov) <- amostras_gexp
head(gexp_ov)[1:5]

rm(ens_gexp, amostras_gexp)

# Filtrando dados clinicos e de sobrevida para manter as mesmas amostras
# que estao na matriz de contagem

survival_ov <- survival_ov[survival_ov$sample %in% colnames(gexp_ov),]

setdiff(colnames(gexp_ov),survival_ov$sample) # "TCGA-04-1357-01A"

which(colnames(gexp_ov) == "TCGA-04-1357-01A") #102
gexp_ov <- gexp_ov[,-102]

setdiff(colnames(gexp_ov),survival_ov$sample)

clinical_ov <- clinical_ov[clinical_ov$submitter_id.samples %in% colnames(gexp_ov),]
# Temos as mesmas 378 amostras nos dados de contagem, clinico e de sobrevida

# Deixando as amostras na mesma ordem entre os diferentes dados
rownames(clinical_ov)
rownames(clinical_ov) <- clinical_ov$submitter_id.samples
clinical_ov <- clinical_ov[colnames(gexp_ov),]
identical(rownames(clinical_ov), colnames(gexp_ov)) #TRUE

rownames(survival_ov)
rownames(survival_ov) <- survival_ov$sample
survival_ov <- survival_ov[colnames(gexp_ov),]
identical(rownames(survival_ov), colnames(gexp_ov)) #TRUE

save(gexp_ov, file = "data_ov/gexp_ov.RData")
save(clinical_ov, file = "data_ov/clinical_ov.RData")
save(survival_ov, file = "data_ov/survival_ov.RData")

############################################################

#-- Carregando os dados salvos
load("data_ov/clinical_ov.RData")
load("data_ov/gexp_ov.RData")
load("data_ov/survival_ov.RData")

clinical_ov$clinical_stage <- ifelse(clinical_ov$clinical_stage == "",
                                     NA,
                                     clinical_ov$clinical_stage)
summary(clinical_ov$clinical_stage)
clinical_ov$clinical_stage <- as.factor(clinical_ov$clinical_stage)
summary(clinical_ov$clinical_stage)

# Dados de expressao para o gene BRCA 1 (ENSG00000012048)
which(rownames(gexp_ov) == "ENSG00000012048") #310
BRCA1_ov <- gexp_ov[310,]
BRCA1_ov <- t(BRCA1_ov)
identical(rownames(BRCA1_ov), rownames(clinical_ov)) #TRUE

clinical_ov <- cbind(clinical_ov, BRCA1_ov)
colnames(clinical_ov)[99]
colnames(clinical_ov)[99] <- "BRCA1"

# Dados de expressao para o gene BRCA 2 (ENSG00000139618)
which(rownames(gexp_ov) == "ENSG00000139618") #7852
BRCA2_ov <- gexp_ov[7852,]
BRCA2_ov <- t(BRCA2_ov)
identical(rownames(BRCA2_ov), rownames(clinical_ov)) #TRUE

clinical_ov <- cbind(clinical_ov, BRCA2_ov)
colnames(clinical_ov)[100]
colnames(clinical_ov)[100] <- "BRCA2"

library(ggplot2)
library(cowplot)


# Criando um boxplot
stage_brca1 <- ggplot(clinical_ov[complete.cases(clinical_ov$clinical_stage),],
                   aes(x = clinical_stage, y = BRCA1, fill = clinical_stage)) +
  geom_boxplot() +
  labs(title = "BRCA1 expression in Ovarian Cancer",
       fill = "Clinical Stage") +
  xlab("Clinical Stage") +
  ylab("BRCA1 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
stage_brca1

stage_brca2 <- ggplot(clinical_ov[complete.cases(clinical_ov$clinical_stage),],
                      aes(x = clinical_stage, y = BRCA2, fill = clinical_stage)) +
  geom_boxplot() +
  labs(title = "BRCA2 expression in Ovarian Cancer",
       fill = "Clinical Stage") +
  xlab("Clinical Stage") +
  ylab("BRCA2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
stage_brca2

# Unindo os plots em uma unica figura
brca1_brca2 <- plot_grid(stage_brca1, stage_brca2, labels = "AUTO",
                         ncol = 2)

# Exemplo de exportacao da figura em pdf
ggsave2(filename = "BRCA1_BRCA2_ov.pdf", brca1_brca2,
        width = 11, height = 6, units = "in")
