#####################################################
############### CURSO R - SBIB ######################
######## Por Jessica Magno e Joao Muzzi #############

## PRATICA: Parte 2

library(tidyverse)

## Baixando os dados preparados na parte 1
load("data_ov/clinical_ov.RData")
load("data_ov/gexp_ov.RData")

class(gexp_ov)

### Extraindo dados de expressão para KRT7 (ENSG00000135480)
### https://pubmed.ncbi.nlm.nih.gov/33416175/
KRT7 <- as.data.frame(t(gexp_ov["ENSG00000135480", ]))
colnames(KRT7) <- "KRT7_exp"

summary(KRT7)

# Puxar a mediana de expressão para definir pacientes com alta e baixa expressão 
med_KRT7 <- median(KRT7$KRT7_exp)

KRT7$KRT7_level <- factor(
  ifelse(
    KRT7$KRT7_exp > med_KRT7,
    "KRT7_high",
    "KRT7_low"),
  levels = c("KRT7_high", "KRT7_low")
    )

summary(KRT7)
# KRT7_exp          KRT7_level 
# Min.   : 7.864   KRT7_high:189  
# 1st Qu.:14.082   KRT7_low :189  
# Median :14.846                  
# Mean   :14.761                  
# 3rd Qu.:15.623                  
# Max.   :17.821              

identical(rownames(clinical_ov), rownames(KRT7)) # T

clinical_ov$KRT7_level <- KRT7$KRT7_level

summary(as.factor(clinical_ov$clinical_stage))
#        Stage IC  Stage IIA  Stage IIB  Stage IIC Stage IIIA Stage IIIB 
# 3          1          3          3         17          7         13 
# Stage IIIC   Stage IV 
# 274         57 

### Simplificar as categorias de estágios tumorais
clinical_ov$Stage <- ifelse(
  grepl("IV", clinical_ov$clinical_stage),
  "IV",
  ifelse(
    grepl("III", clinical_ov$clinical_stage),
    "III",
    ifelse(
      grepl("II", clinical_ov$clinical_stage),
      "II",
      ifelse(
        grepl("I", clinical_ov$clinical_stage),
        "I",
        clinical_ov$clinical_stage
      )
    )
  )
)
summary(as.factor(clinical_ov$Stage))
#     I  II III  IV 
# 3   1  23 294  57 

### com table criamos uma tabela de contingência
### com prop.table calculamos a proporção dessa tabela (em percentual)
tab <- prop.table(
  table(
    clinical_ov$Stage, clinical_ov$KRT7_level), 
  margin = 1)*100

tab

# Transformar a tabela em um formato longo
tab_long <- reshape2::melt(tab)

# Exibir o formato longo da tabela
print(tab_long)

# Criar o gráfico de barras
g <- tab_long %>%
  dplyr::filter(
    Var1 != ""
  ) %>%
  ggplot(aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Tumor Stage", y = "Percentage", fill = "KRT7 Level") +
  ggtitle("Distribution of KRT7 Levels by Stage") +
  theme_minimal()
plotly::ggplotly(g)

### https://xenabrowser.net/datapages/?dataset=TCGA-OV.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
### Os dados brutos estão transformados (log2(counts+1))
### voltar para contagens brutas
gexp_ov <- round(2^gexp_ov - 1)

hist(as.numeric(gexp_ov["ENSG00000135480", ]))
hist(KRT7$KRT7_exp)

### Com o pacote AnnotationHub vamos baixar os dados atualizados
### de nomenclatura de genes para converter ENSEMBL ID em Symbol
### Também colhemos informação sobre o tipo do gene, nesse caso
### vamos selecionar somente os codificantes
library(AnnotationHub)
vignette("AnnotationHub")
### gerando rowAnnotation
ah <- AnnotationHub()
query(ah, "EnsDb")
edb <- query(ah, pattern = c("Homo sapiens", "EnsDb", 111))[[1]] #https://www.ensembl.org/index.html?redirect=no
gns <- genes(edb, return.type= "DataFrame")
gns <- gns[gns$gene_biotype=="protein_coding",]
#gns$length <- width(gns)
rowAnnotation <- as.data.frame(gns)

write_csv(rowAnnotation, "./data_ov/RowAnnotation.csv")

rowAnnotation <- as.data.frame(read_csv("data_ov/RowAnnotation.csv"))

gexp_ov <- gexp_ov %>%
  dplyr::filter(
    rownames(gexp_ov) %in% rowAnnotation$gene_id,
    rowMeans(gexp_ov) > 1
  )
rownames(rowAnnotation) <- rowAnnotation$gene_id
rowAnnotation <- rowAnnotation[rownames(gexp_ov), c("gene_id", "symbol")]

identical(rownames(rowAnnotation), rownames(gexp_ov)) # T

### Agora faremos a expressão diferencial
### Vamos utilizar o pacote DESeq2 o qual foi feito para trabalhar
### com dados de RNA-seq e contagens brutas. Recomendamos a leitura
### da vinheta para aprofundar. 
### Outro pacote bastante famoso é o limma
library(DESeq2)
vignette("DESeq2")

### criando objeto DESeq com as contagens brutas e dados clinicos
dds <- DESeqDataSetFromMatrix(countData = gexp_ov,
                              colData = clinical_ov,
                              design = ~ KRT7_level)

dds <- DESeq(dds)

# Extraindo resultados
res <- results(dds, contrast = c("KRT7_level", "KRT7_high","KRT7_low" ))

## Criando data frame com Symbols para fácil identificação dos genes
identical(rownames(res),rownames(rowAnnotation)) #TRUE
res.df <- cbind(rowAnnotation, as.data.frame(res)) #18222genes

write_csv(res.df, "./data_ov/res_df.csv")
res.df <- read_csv("./data_ov/res_df.csv")

### explicação das colunas dos resultados do DESeq2
?`DESeq2-package`
# For results: a DESeqResults object, which is a simple subclass of DataFrame. 
# This object contains the results columns: baseMean, log2FoldChange,
# lfcSE, stat, pvalue and padj, and also includes metadata columns
# of variable information. The lfcSE gives the standard error of the 
# log2FoldChange. For the Wald test, stat is the Wald statistic: 
#   the log2FoldChange divided by lfcSE, which is compared to 
# a standard Normal distribution to generate a two-tailed pvalue.

res.df <- res.df[order(res.df$padj),]

res.df$'D.E.genes' <- as.factor(
  ifelse(res.df$padj < (0.05) & res.df$log2FoldChange > 0,
         "up",
         ifelse(res.df$padj < (0.05) & res.df$log2FoldChange < 0, 
                "down", 
                "NS")))
summary(res.df$D.E.genes)

# D.E.genes   
# down: 1571  
# NS  :13168  
# up  : 3483  

# Salvando os resultados
writexl::write_xlsx(
  res.df, 
  path = "./Dif_exp_genes.xls") 

### Tendo uma lista de genes diferencialmente expressos iremos fazer 
### o enriquecimento de vias
### Para aprofundar veja site abaixo
### https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/#interpreting-gsea-results
library(fgsea)
library(msigdbr)
### para aprofundar nos gene sets do MSigDb
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
### Artigo original sobre Hallmarks of Cancer
# https://pubmed.ncbi.nlm.nih.gov/26771021/

### Aqui criamos uma lista com cada um dos 50 Hallmarks e seus
### respectivos genes
h_gene_sets = msigdbr(species = "Homo sapiens", 
                      category = "H") ## H for Hallmarks
msigdbr_list = split(x = h_gene_sets$ensembl_gene, 
                     f = h_gene_sets$gs_name)

### Vamos gerar uma lista de DEGs para o enriquecimento
degs <- res.df %>%
  dplyr::filter(
    padj < 0.05
  )
gen <- degs$log2FoldChange
names(gen) <- degs$gene_id
gen <- gen[order(gen)]
gen

### Finalmente fazemos o GSEA
df.gsea <- fgsea(pathways = msigdbr_list, 
                 stats= gen, 
                 nPermSimple=1000)
View(df.gsea)
### Definindo o sentido do enriquecimento
df.gsea$effect <- ifelse(df.gsea$NES>1,">1","<1")
### Preparando o data frame para plot
df.gsea <- df.gsea[order(df.gsea$padj),]

# filtrando Hallmarks significativos
df.gsea <- df.gsea %>%
  dplyr::filter(
  padj < 0.05
)
df.gsea$pathway <- gsub("_", " ", df.gsea$pathway)
df.gsea$pathway <- gsub("HALLMARK", "", df.gsea$pathway)
df.gsea$pathway <- factor(df.gsea$pathway,
                          levels = rev(df.gsea$pathway))

### Plotting GSEA tables
g <- ggplot(df.gsea) +
  geom_point(aes(y=pathway, x = -log10(padj), 
                 fill = effect, size = abs(NES)), 
             show.legend = c(size=T, colour=T), alpha=0.6, shape = 21) +
  scale_fill_manual(values = c("blue","red"), aesthetics="fill", drop=F)+
  labs(size="Absolute\nNES",fill="Enrichment\n(NES)", color="Process Category",
       title= "GSEA Hallmark" ) +
  scale_y_discrete(expand = expansion(mult = c(0.05,0.05))) +
  scale_size_continuous(range = c(3, 10)) + 
  ylab("Hallmarks gene sets") + 
  xlab("-log (P-adj)")+
  theme(
    plot.title=element_text(hjust=0),
    legend.position="right",
    title = element_text(size=12), 
    axis.text.x = element_text(vjust=0.5, hjust=1),
    axis.ticks = element_line(linetype=1, color="grey"),
    legend.title=element_text(size=10, face = "bold"),
    aspect.ratio=2,
    panel.border = element_rect(colour = "grey", fill=NA, size=1))
g
ggsave(g, filename = "~/ACC/CD34/Results/Fig2_GSEA_ACC.pdf.gsea", 
       height = 11, width = 7, units="in", device="pdf.gsea")



