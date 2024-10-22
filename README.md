---
title: "Génomique fonctionnelle: Analyse Transcriptomique dite Bulk"
author: Martin Braud, Thomas Stervinou
email: thomas.stervinou@univ-nantes.fr
geometry: margin=1cm
output: 
 html_document:
  toc: true
  toc_float: true
  toc_collapsed: true
  toc_depth: 6
  theme: united
editor_options: 
  markdown: 
    wrap: 100
---

# Objectifs du TP

**Les informations pour l'installation de RStudio, R et les différents paquets sont décrits dans le
[README.txt]{.underline}**

> Le but du TP est d'être capable d'analyser des données de transcriptomique à partir d'une matrice
> de comptage

-   Prise en main de R avec Rstudio

-   Etre capable de lire un fichier de matrice

-   Etre capable d'identifier avec divers outils la qualité des projets

-   Etre capable d'utiliser les outils d'analyse d'expression différentielle de gènes (DGE)

# 1. Introduction à l'analyse Transcriptomique

> L'analyse à pour but de quantifier et comparer l'expression des gènes de différents échantillons.
> Ce processus peut être décrit en plusieurs étapes:

-   Extraction et traitement des ARN des tissus/cellules

-   Préparation des librairies de séquençages / identification de chaque échantillon

-   Séquençage des échantillons

-   Extraction des séquences (FASTQ)

-   Alignements des séquences (BAM), Annotation et comptage du nombre de molécules uniques associées
    à un gène (matrice de comptage)

*- Extraction des informations de la matrice de comptage et analyse différentielle*

## 1.1 Workflow d'analyse transcriptomique

![](Workflow_Transcripto.svg)

## 1.2 Annotation des reads

> Les fragments d'ADNc séquencés (reads) sont annoté lors de l'alignement. Dans notre cas présent
> l'alignement s'effectue avec BWA sur le fichier du transcriptome humain en format FASTA.

<div>

**\>NR_152111** 1 GAGGAGAGAGCAGAGTATACCGCAGACATCATTTCTACTACAGTGGCGGAGCCGTACAGG
ACCTGTTTCACTGCAGGGGGATCCAAAACAAGCCCCGTGGAGCAGCAGCCAGA.......

**\>NR_152112** 1 GAGGAGAGAGCAGAGTATACCGCAGACATCATTTCTACTACAGTGGCGGAGCCGTACAGG
ACCTGTTTCACTGCAGGGGGATCCAAAACAAGCCCCGTGGAGCAGCAGCCAGA.......

**\>NR_137427** 1 AGTTCTCACGCTAGGTCTCCTCCGGGCGCTTCCCTAGCCCGTTCGCCGCCTGAGAGGGAC
GCTGTTCCGCCGCGTGGAAGCTTCGAGTCTCGACTCCACTGTTGACCCCTAGA.......

</div>

> Les transcripts annotés sont en format RefSeq (ncbi). Un second fichier d'annotation permet
> d'obtenir les symboles des gènes.

<div>

...

**ATXN1** NM_000332,NM_001128164,NM_001357857,**NR_152111**,**NR_152112**,NR_152113,NR_152114

**KLRG1**
NM_001329099,NM_001329101,NM_001329102,NM_001329103,NM_005810,NR_137426,**NR_137427**,NR_137428 ...

...

</div>

# 2. Exercices Pratique : Analyse de deux jeux de données DGEseq {.tabset}

> Deux jeux de données sont disponible pour ce TP consistant en deux matrices ainsi qu'un fichier de
> méta données.

-   Le premier jeu, *matrice1.csv*, contient les valeurs de comptage pour chaque gènes pour 50
    échantillons de peau séquencé en 2021.

-   Le second jeu *matrice2.csv* consistant en la suite de la première analyse et contient les
    valeurs de comptage pour chaque gènes pour 60 échantillons de peau séquencé en 2023

-   Le fichier de meta-données contient les informations relatives aux deux jeux de données.

> Les deux jeux de données correspondent à deux groupes de patients sélectionnés à 1 an
> d’intervalle.

> Ces patients sont atteint de dermatite atopique. Deux prélèvements ont été effectués sur chaque
> patient, un de peau saine (PNL) et un de peau lésée (PL). De plus des patient contrôle ne
> souffrant pas de dermatite atopique on été inclus (CTRL)

## 2.1 Contrôle Qualité des données: Matrice-1 2021

### 2.1.1 Lecture et manipulation de fichiers de comptage

```{r matrice1}
setwd('/Chemin/vers/dossier/TP/')
#Charge la matrice dans une variable à partir d'un fichier CVS
mat1 <- read.csv('Data/matrice1.csv', row.names=1)
#Visualise la matrice sous forme d'un data frame
head(mat1)
```

```{r}
#Visulaliser les noms des echantillons en colonne : names()
nom_ech <- names(mat1)
nom_ech
#Compte le nombre des échantillons
length(nom_ech)
```

```{r}
#Visualiser les noms gènes identifiés en ligne : row.names
nb_gene <- row.names(mat1)
head(nb_gene)
#Compte le nombre de gène
length(nb_gene)
```

```{r}
#Chargement des meta données
metadata <- read.csv('Data/matrice_metadata.csv', row.names=1)
head(metadata)
```

```{r}
#Attribut du data frame metadata
attributes(metadata)
```

```{r}
#Selection par ligne et colonne d'une partie des metadonnées sous forme [ligne,colonne]
#Selection des 3 première lignes et quatre première colonnes
metadata[1:30,1:6]
```

```{r}
#Selection par nom de colonne et nom de ligne
metadata[c('A01_2021','B01_2021','A03_2021'),c('description','conditions')]
```

```{r}
#Identification des valeurs contenu dans une colonne : table
table(metadata$Gender)
```

```{r}
#Avec table() combien de variable différente trouve t'on pour les "conditions"
table(metadata$conditions)
```

> On identifie trois variable pour les conditions :

> -   12 échantillons contrôle (CTRL)
> -   49 échantillons peaux lésées (PL)
> -   49 échantillons peaux non lésées (PNL)

```{r}
#Selection d'un sous ensemble par valeur avec subset()
#Sous ensemble des échantillons contrôles (CTRL)
subset(metadata, subset = conditions=='CTRL')
```

```{r}
#Selection des metadata du premier séquençage NovaSeq_211109
metadata_2021 <- subset(metadata, batch=='NovaSeq_211109')
table(metadata_2021$batch)
```

### 2.1.2 Identification de la qualité des échantillons par le nombre de gènes exprimés et d'UMI (reads)

```{r}
#Comptage du nombre de gènes exprimé par écahntillons
genes_mat1 <- colSums(mat1 > 0)
genes_mat1
```

```{r, width=50}
#Visualisation graphique du nombre de gènes exprimés par échantillons
library(ggplot2)
library(ggrepel)
genes_df1 <- data.frame('gene_count'=genes_mat1, "samples"=names(genes_mat1))
genes_plot1 <- ggplot(data=genes_df1, aes(x=samples, y=gene_count, fill='2021')) + geom_bar(stat="identity") + ggtitle('Comptage des gènes exprimés par échantillons') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6), plot.title = element_text(hjust = 0.5, size=10))
show(genes_plot1)
```

> L'échantillon DA5_2PNL à beaucoup moins de gènes que les autres échantillons

```{r}
#Visualisation graphique du nombre de "UMI" exprimés par échantillons
umi_genes <- data.frame('umi_count'=colSums(mat1), samples=names(mat1))
plot_umis <- ggplot(data=umi_genes, aes(x=samples, y=umi_count, fill='2021')) + geom_bar(stat="identity") + ggtitle('Comptages des UMIs par échantillons') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6), plot.title = element_text(hjust = 0.5, size=5))
show(plot_umis)
```

```{r}
#Quels échantillons < 500k UMI ?
subset(umi_genes, umi_count < 500000)
```

> Les échantillons DA15_4PL et DA5_2PNL présentent un compte d'UMIs beaucoup plus faible que les
> autres échantillons

```{r}
#Identification de l'expression des gènes :  combien d'UMI par gènes?
counts_umi <- data.frame('umi_count'=rowSums(mat1))
umi_pergenes_plt <- ggplot(counts_umi, aes(x=umi_count)) + geom_histogram(bins=100, colour=4,fill='lightblue') + ggtitle("Distribution du nombre d'UMI par gène")
show(umi_pergenes_plt)
```

> Beaucoup de gène ne sont peu ou pas exprimé. Afin de déterminer un cut-off, visualisons les gènes
> qui ont moins de 500 UMI à travers tout les échantillons

```{r}
#Identification de l'expression des gènes :  combien d'UMI par gènes?
umi_pergenes_plt_500 <- ggplot(subset(counts_umi, umi_count <500), aes(x=umi_count)) + geom_histogram(bins=100, colour=4,fill='lightblue') + ggtitle("Distribution du nombre d'UMI par gène ayant moins de 500 UMI")
show(umi_pergenes_plt_500)
```

> D'après le graphique, on peut fixer un cut-off à 200 UMI

```{r}
mat1 <- subset(mat1, rowSums(mat1) > 200)
```

### 2.1.3 Identification de la distribution des échantillons en fonction de l'expression des gènes : ACP

> Par simplicité de représentation nous utilisons le package DESeq2 utilisé aussi plus tard pour
> l'expression différentielle.
>
> Les données peuvent être transformée pour visualisation en utilisant la transformation des
> variances (VST) ou rlog. Utile pour identifier des outliers
>
> VST : Variance Stabilizing Transformation $\neq$ Normalisation

```{r, warning=FALSE}
suppressMessages(library(DESeq2))
row.names(metadata_2021) <- metadata_2021$description
#Formatage de la matrice de comptage et des meta données associées
dds1 <- DESeqDataSetFromMatrix(countData = mat1,colData = metadata_2021 ,design = ~ conditions)
#Normalisation des valeurs de comptage
vsd1 <- vst(dds1,blind=FALSE)
#ACP
pcaData1 <- plotPCA(vsd1, intgroup=c("conditions", "Gender"), returnData=TRUE)
#Récupération du % de variance pour la PC1 et PC2
percentVar1 <- round(100 * attr(pcaData1, "percentVar"))
pcaplot1 <- ggplot(pcaData1, aes(PC1, PC2, color=conditions, shape=Gender)) + geom_point(size=2) + xlab(paste0("PC1: ",percentVar1[1],"% variance")) + ylab(paste0("PC2: ",percentVar1[2],"% variance")) + coord_fixed() + theme(plot.title=element_text(size=20, hjust=0.5), axis.title=element_text(size=15)) + ggtitle('PCA collapsed samples Conditions and Batch') + geom_text_repel(data=pcaData1,aes(label=name), size=3)
show(pcaplot1)
```

> Les patients sont identifiés ici par le numéro avant la condition. Si l'on veut retirer les
> échantillons DA15_4PL et DA5_2PNL ont doit retirer tout les échantillon lié à ce patient à savoir
> 4PL, 4PNL, 2PL et 2PNL

```{r, message=FALSE}
grep('*_4PNL|*_4PL', names(mat1), value=TRUE)
grep('*_2PL|*_2PNL', names(mat1), value=TRUE)
```

```{r}
#On retire alors les échantillons DA15_4PL, DA5_2PNL, DA13_4PNL, DA7_2PL
ech_garder <- !names(mat1) %in% c('DA15_4PL', 'DA5_2PNL', 'DA13_4PNL', 'DA7_2PL')
mat1 <- mat1[,ech_garder]
#De même pour les metadata
metadata_2021 <- metadata_2021[ech_garder,]
#On peut combiner nos deux éléments en créant de nouveau un objet DESeq
dds_res1 <- DESeqDataSetFromMatrix(countData = mat1,colData = metadata_2021 ,design = ~ conditions)
```

### 2.1.4 Sauvegarde de la matrice nettoyées

> On peut sauvegarder tout objets créé précédemment en utilisant la commande saveRDS()

```{r}
#Sauvegarde de l'object dds_res1

saveRDS(dds_res1,'Data/DDS1_Matrice_nettoyée.rds')

```

## 2.2 Contrôle Qualité des données: Matrice-2 2023

### 2.2.1 Lecture et manipulation de fichiers de comptage

```{r}
#Charge la matrice 2023 dans une variable à partir d'un fichier CVS
mat2 <- read.csv('Data/matrice2.csv', row.names=1)
head(mat2)
```

```{r}
#Visulaliser les noms des echantillons en colonne : names()
```

```{r}
#Visualiser les noms gènes identifiés en ligne : row.names
```

```{r}
#Selection des metadata du premier séquençage NovaSeq_211109
```

### 2.2.2 Identification de la qualité des échantillons par le nombre de gènes exprimés et d'UMI (reads)

```{r}
#Comptage du nombre de gènes exprimé par écahntillons
```

```{r, width=50}
#Visualisation graphique du nombre de gènes exprimés par échantillons
```

> *Que peux-ton dire sur le nombre de gènes exprimés par échantillon?*

```{r, fig.widht=10}
#Visualisation graphique du nombre de "UMI" exprimés par échantillons
```

> *Que peut-on dire sur le nombre d'UMI par échantillons?*

```{r}
#Identification de l'expression des gènes :  combien d'UMI par gènes?
```

```{r}
#Identification de l'expression des gènes :  combien d'UMI par gènes?
```

> *D'après le graphique, quelle valeur semble un bon cutoff pour le nombre d'UMI par gène?*

```{r}
```

### 2.2.3 Identification de la distribution des échantillons en fonction de l'expression des gènes : ACP

```{r, warning=FALSE}
```

### 2.2.4 Sauvegarde de la matrice nettoyées

```{r}
#Sauvegarde de l'object dds_res1

```

## 2.3 Intégration des deux jeux de données

### 2.3.1 Lecture des matrices sauvées

```{r}
dds_res1 <- readRDS('Data/DDS1_Matrice_nettoyée.rds')
dds_res2 <- readRDS('Data/DDS2_Matrice_nettoyée.rds')
```

### 2.3.2 Fusions des matrices

```{r}
#Intégration des deux jeux de données
# Extraction des matrices de comptage brutes
mat1 <- counts(dds_res1)
mat2 <- counts(dds_res2)
metadata_2021 <- colData(dds_res1)
metadata_2023 <- colData(dds_res2)

#Certains échantillons ont été dupliqué de la première analyse pour la deuxième afin de faciliter l'intégration 
mat3 <- merge(mat1,mat2,by='row.names', all=TRUE, suffixes=c('_2021','_2023'))
#Donne la colonne 'Row.names' comme nom de ligne de mat3 et supprime cette colonne
row.names(mat3) <- mat3$Row.names
mat3$Row.names <- NULL
#Remplace les NA par 0 dans la matrice de comptage
mat3[is.na(mat3)] <- 0
```

```{r}
#Integration des metadata
metadata_all <- rbind(metadata_2021,metadata_2023)
row.names(metadata_all) <- names(mat3)
```

```{r}
#Comptage du nombre de gènes exprimé par échantillons
genes_mat3 <- colSums(mat3 > 0)
genes_mat3
```

```{r, width=50}
#Visualisation graphique du nombre de gènes exprimés par échantillons
```

```{r, width=50}
#Visualisation graphique du nombre de "UMI" exprimés par échantillons

```

> En Intégrant deux jeux de données ayant été préparé à different moment, dans différentes
> conditions et par différentes personnes cela peut engendrer des différences d'expression entre les
> jeux de données qui impacteront les résultats

### 2.3.3 Identification de la distribution des échantillons en fonction de l'expression des gènes : ACP

```{r, warning=FALSE, message=FALSE}
#Formatage de la matrice de comptage et des meta données associées

#Normalisation des valeurs de comptage

#ACP

#Récupération du % de variance pour la PC1 et PC2

```

*Qu'observe t-on à propos de la répartition des deux jeux de données?*

### 2.3.4 Correction de l'effet batch

```{r}
#Correction de l'effet "Batch": Utilisation de ComBatSeq (https://bioconductor.org/packages/release/bioc/html/sva.html) 
library(sva)
dds_corr <- ComBat_seq(as.matrix(counts(dds3)), batch=dds3$batch, group=dds3$conditions)
dds4 <- DESeqDataSetFromMatrix(countData = dds_corr,colData = metadata_all,design = ~ batch + conditions)
vsd4 <- vst(dds4,blind=FALSE)
```

```{r, warning=FALSE}
library(cowplot)
pcaData4 <- plotPCA(vsd4, intgroup=c("replicates", "batch"), returnData=TRUE)
percentVar4 <- round(100 * attr(pcaData4, "percentVar"))
pcaplot4 <- ggplot(pcaData4, aes(PC1, PC2, color=batch, shape=replicates)) + geom_point(size=2) + xlab(paste0("PC1: ",percentVar1[1],"% variance")) + ylab(paste0("PC2: ",percentVar1[2],"% variance")) + coord_fixed() + theme(plot.title=element_text(size=20, hjust=0.5), axis.title=element_text(size=15)) + ggtitle('ACP échantillons intégrés effet Batch corrigé') + geom_text_repel(data=pcaData4,aes(label=name), size=3)
show(pcaplot4)
#Visualiser les deux graphiques avant après correction avec plot_grid()
```

-   **Visualiser les deux graphique avant et après correction avec plot_grid()**
-   **Visualiser les réplicats en couleur sur la PCA**

```{r}
#Réduction de chaque paire d'échantillon duplicat en un seul
dds4 <- collapseReplicates(dds4, dds4$description)
metadata_all <- colData(dds4)
vsd4 <- vst(dds4, blind=FALSE)
```

### 2.3.5 Correlation des covariables : Corrélation de Pearson

> Les métadonnées représentent différentes valeurs associées à chaque échantillons, certaine de ces
> valeurs peuvent influer sur l'expressions des gènes en plus des conditions que nous voulons
> analyser. Afin de d'identifier si certaines conditions peuvent avoir un impact nous allons
> comparer la corrélation entre les valeurs de métadonnées et celles obtenues par analyse en
> composante principale

```{r}
library(mixOmics)
#Suppression des échantillons contrôle nous soumis aux conditions
vsd_noctrl <- vsd4[,row.names(metadata_all)[!metadata_all$conditions == 'CTRL']]
metadata_noctrl <- colData(vsd_noctrl)
metadata_noctrl$conditions <- droplevels(metadata_noctrl$conditions)
#ACP sur ces données avec 10 composante maximum
pca_corr <- mixOmics::pca(t(assay(vsd_noctrl)),10)
```

```{r}
#Formatage des meta données pour la matrice de corrélation et ajout des CP
matrix_corr <- metadata_noctrl[,c('conditions','Gender','IL22.IL22RA2')]
conditions <- as.character(matrix_corr$conditions)
matrix_corr$conditions <- conditions
#Pour pouvoir être comparée, les données doivent être numérique et renseignée (pas de NA)
matrix_corr$conditions[matrix_corr$conditions=='PNL'] <- 1
matrix_corr$conditions[matrix_corr$conditions=='PL'] <- 2
matrix_corr$conditions <- as.numeric(matrix_corr$conditions)
matrix_corr$Gender[matrix_corr$Gender=='F'] <- 1
matrix_corr$Gender[matrix_corr$Gender=='M'] <- 2
matrix_corr$IL22.IL22RA2[is.na(matrix_corr$IL22.IL22RA2)] <- 0
matrix_corr$IL22.IL22RA2 <- as.numeric(matrix_corr$IL22.IL22RA2)
matrix_corr_pca <- cbind(matrix_corr,
		                    PC1 = pca_corr$variates$X[,1], PC2 = pca_corr$variates$X[,2],
		                    PC3 = pca_corr$variates$X[,3], PC4 = pca_corr$variates$X[,4],
		                    PC5 = pca_corr$variates$X[,5], PC6 = pca_corr$variates$X[,6],
		                    PC7 = pca_corr$variates$X[,7], PC8 = pca_corr$variates$X[,8],
		                    PC9 = pca_corr$variates$X[,9], PC10 = pca_corr$variates$X[,10]
	)
```

```{r ,warning=FALSE, message=FALSE}
library(Hmisc)
library(corrplot)
#Calcul de la matrice de corrélation et visualisation
corr_m = rcorr(as.matrix(matrix_corr_pca),type="pearson")
corrplot(corr_m$r, p.mat=corr_m$P, method='square',sig.level = c(0.001, 0.01, 0.05),insig = 'label_sig',pch.cex = 0.9, col=rev(COL2('RdBu')), diag=FALSE)

```

> **En regardant les 2 première composantes (expliquant 25% de la variance des échantillons) on voit
> une corrélation entre l'expression des gènes et ?**

> **En fonction des co-variables identifiée, nous allons redéfinir le modèle utilisé pour l'analyse
> différentielle. Ce modèle prends en compte les co-variables dont l'impact devrait être diminué.**

> De plus nous travaillions sur des échantillons appariés, un patient présente deux conditons PL et
> PNL. Il faut donc prendre en compte cette variable dans le design du modèle

```{r}
#Renseigner les covariables impliquées manquantes
dds4$IL22.IL22RA2[is.na(dds4$IL22.IL22RA2)] <- 0
dds4$Gender <- as.factor(dds4$Gender)
#Redéfinition du design du modèle en fonction des covariables:
design(dds4) <-  ~ Patient + IL22.IL22RA2 + Gender + conditions
```

#### 2.3.5.1 Création de l'objet pour le calcul de l'expression différentielle : DESeq()

> DESeq2 utilise la méthode de la médiane des ratios pour calculer le facteur de taille (size
> factor) et normaliser le comptage des reads (UMI) entre échantillons. Puisque l'analyse
> d'expression différentielle compare le comptage entre groupe d'échantillons pour un gène donné, il
> n'est pas nécessaire de prendre en compte la longueur de ce gène.

#### 2.3.5.2 Normalisation : Estimation du facteur de taille (Size Factor)

```{r}
subset_raw <- counts(dds4)[1:3,]
subset_raw
```

-   Pour chaque gènes un pseudo échantillon de référence est créé correspondant à la moyenne
    géométrique de tous les échantillons:

    Soit la moyenne géométrique du gène A pour n échantillons :

$$\sqrt{\prod_{i=1}^{n} compte^{A}}^{n_{ième}}  = \;\;\;\;\; \left({\prod_{i=1}^{n} compte^{A}}\right) ^ {1/n}$$

La moyenne géométrique du gène A = produit(counts de A) puissance (1/nombre d'échantillons)

```{r}
# Calcul moyenne géométrique
geom_mean <- rowProds(subset_raw, useNames=TRUE) ^ (1/length(colnames(subset_raw)))
geom_mean
```

-   Pour chaque échantillon on calcul le ratio par rapport à l'échantillon de référence

```{r}
subset_ratio <- subset_raw / geom_mean
subset_ratio
```

-   Le facteur de de normalisation ou facteur de taille correspond à la médiane des ratios des gènes
    de chaque échantillons

```{r}
# Ici nous n'avons qu'un seul gène
norm_factor <- apply(subset_ratio,2,median)
norm_factor
```

-   On normalise alors les comptes de chaque échantillon avec le facteur de normalisation calculé

```{r}
subset_norm <- t(apply(subset_raw,1, '/', norm_factor))
subset_norm
```

#### 2.3.5.3 Normaliser toute la matrice de compte

```{r}
count_raw <- counts(dds4)
count_prod <- count_raw
count_prod[count_prod==0] <- NA
#Produit des counts sans les valeures nulles
row_prod_count <- rowProds(count_prod, useNames=TRUE, na.rm=TRUE)
#Nombre ce valeurs restante par genes
num_val_rest <- length(colnames(count_raw)) - apply(count_prod,1, function(x){sum(is.na(x))})
geom_mean <- row_prod_count ^ (1/length(colnames(count_raw)))

count_ratio <- count_raw / geom_mean
norm_factor <- apply(count_ratio,2,median)
count_norm <- t(apply(count_raw,1, '/', norm_factor))
count_norm[1:3,1:10]
```

**De petite différence existe en fonction de la méthode utilisée pour le produit de valeurs nulle**

> Plusieurs étapes sont effectuées par la fonction DESeq() :

> Définition des différentes étapes de DESeq2 et les fonctions indépendante pour effectuer le
> calcul:

-   Estimation de la taille des facteurs : Normalisation =\> **estimateSizeFactors()**

-   Estimation de la dispersion : Représente la variance de l'expression d'un gène pour une moyenne
    donnée

    -   Estimation de la dispersion liées aux gènes =\> **estimateDispersions**

-   Relation liée à la dispersion de la moyenne

-   Estimation finale de la dispersion - Ajustement du modèle et test (distribution binomiale
    negative et Test de Wald) =\> **nbinomWaldTest()**

```{r}
dds4 <- DESeq(dds4)
```

> Certaines valeurs calculées à chaque étapes peuvent être extraites

```{r}
cat('Facteur: ')
head(sizeFactors(dds4))
#Visualiser les dispersions
cat('Dispersion: ')
head(dispersions(dds4))
```

## 2.4 Expression Différentielle {.tabset}

> Sûr ce jeux de données nous disposons de 3 conditions que nous voulons comparer à savoir: - CTRL :
> contrôle - PL : Peaux lésées - PNL : Peaux non lésée

> Nous allons effectuer les comparaison suivante afin d'identifier les gènes dont l'expression varie
> entre deux conditions:

-   PNL vs CTRL

-   PL vs CTRL

-   PL vs PNL

### 2.4.1 Comparaison entre CTRL et PNL

#### 2.4.1.1 Extraction des résultats de comparaison et visualisation

> La fonction results() effectue la comparaison entre la variable d'interêts PNL et la référence
> CTRL

```{r}
result_pnl_ctrl_05 <- results(dds4, contrast=c('conditions','PNL','CTRL'), alpha=0.05)
summary(result_pnl_ctrl_05)
```

```{r}
plotMA(result_pnl_ctrl_05, ylim=c(-2,2))
```

```{r}
head(result_pnl_ctrl_05)
```

> **log2FoldChange:** Definit par le FoldChnage comme le ratio d'expression PNL/CTRL

> **lfcSE:** L'erreur type (Standard Error) utilisé pour le calcul du test de Wald

> \*\*stat: Wald statistic = log2FoldChange / lfcSE\*\*\*

> **pvalue:** la valeur statistique de signifcance dut test de Wald\*\*

> **padj:** la p-value corrigée incluant le FDR (False Discovery Rate) avec la méthode
> Benjamini-Hochberg

```{r}
#Function d'annotation des gènes différentiellement exprimés
get_df <- function(results,n_genes, fc, pval, ptype){
	df = data.frame(results)
	df <- df[order(df[[ptype]]),]
	df$diffexp <- NA
	df$diffexp[df$log2FoldChange > fc & df[[ptype]] < pval] <- "UP"
	df$diffexp[df$log2FoldChange < -fc & df[[ptype]] < pval] <- "DOWN"
	df$topgenes <- NA
	deg_id <- which(!is.na(df$diffexp))
	if (length(deg_id) < n_genes) {
		df$topgenes[deg_id] <- row.names(df)[deg_id]
	} else {
		df$topgenes[deg_id[1:n_genes]] <- rownames(df)[deg_id[1:n_genes]]
	}
	return(df)
}
```

<br>

```{r}
df_pnl_ctrl <- get_df(result_pnl_ctrl_05,50, 0.6,0.05,'padj')
df_pnl_ctrl <- df_pnl_ctrl[!is.na(df_pnl_ctrl$padj),]
head(df_pnl_ctrl)
```

```{r}
#Nombre de gènes UP et DOWN
c('UP'=length(which(df_pnl_ctrl$diffexp=='UP')),
     'DOWN'=length(which(df_pnl_ctrl$diffexp=='DOWN')))
```

#### 2.4.1.2 Représentation en Volcano Plot

```{r, warning=FALSE}
genes_df_pnl_ctrl <- row.names(df_pnl_ctrl[!is.na(df_pnl_ctrl$diffexp),])
vlcoP_pnl_ctrl1 <- ggplot(data=df_pnl_ctrl, aes(x=log2FoldChange, y=-log10(padj), col=diffexp)) + ggtitle(paste('DE genes NLS vs CTRL (log2FC=1 ;padj=0.01)')) + geom_point() + theme_minimal() + theme(plot.title=element_text(size=8, hjust=0.5)) + geom_text_repel(aes(x=log2FoldChange, y=-log10(padj), label=topgenes), max.overlaps=Inf)
show(vlcoP_pnl_ctrl1)
```

#### 2.4.1.3 Représentation des gènes en Heatmap : Zscore; le nombre d'écart type par rapport à la moyenne:

$$z = (x – μ) / σ$$

-   x: une valeur donnée
-   μ: la moyenne de la population
-   σ: l'écart type de la population

> Chaque gènes représente une population, nous allons donc calculer le Z-score pour chaque
> échantillons en fonction de chaque gènes

```{r}
#Exemple sur ALAS2 sur les échantillons CTRL et PNL
alas2 <- assay(vsd4)['ALAS2',vsd4$conditions %in% c('PNL','CTRL')]
#Calcul de la moyenne
moyenne <- mean(alas2)
#Calcucul de l'écart type
etype <- sd(alas2)
#calcul du zscore
zscore <- (alas2 - moyenne) / etype
zscore
```

> La fonction scale() permet d'effectuer le calcul du z-core pour chaque colonne

```{r, fig.width=10, fig.height=13}
library(RColorBrewer)
library(pheatmap)
#Définition des couleurs pour les co-variables
ann_color=list("conditions"=c(CTRL=brewer.pal(8,'Greys')[5], PNL=brewer.pal(8,'GnBu')[6]),
		"Gender"=c(F=brewer.pal(8,'Purples')[5],M=brewer.pal(8,'Greens')[5])
		)
#Renverser la matrice pour que les colonnes soient des gènes
vsd_z <- t(scale(t(assay(vsd4)[genes_df_pnl_ctrl,vsd4$conditions %in% c('PNL','CTRL')])))
#Créé une matrice des co-variables à identifier
df_hp_sub <- as.data.frame(colData(dds4)[,c("conditions","Gender","IL22.IL22RA2")])
names(df_hp_sub) <- c("conditions","Gender","IL22.IL22RA2")
rownames(df_hp_sub) <- rownames(colData(dds4))
exp_gene_sub_hp <- pheatmap(vsd_z, cluster_rows=TRUE, show_rownames=TRUE,fontsize_row=8, cluster_cols=TRUE, annotation_col=df_hp_sub, annotation_colors=ann_color,heigh=100, width=100, main='PNL vs CTRL Zscore samples')
show(exp_gene_sub_hp)
```

#### 2.4.1.4 Représentation de la dispersion des échantillons en boxplot

```{r}
#Soit les 10 gènes les plus significatif:
top_genes <- row.names(df_pnl_ctrl)[1:20]
data_boxP <- data.frame("zscore"=vsd_z['FOS',], "conditions"=colData(dds4)[colnames(vsd_z),]$conditions, 'gene'='FOS')
```

```{r}
#Dessiner le boxplot pour ALAS2
boxplot1 <- ggplot(data_boxP, aes(x=conditions,y=zscore, color=conditions)) + geom_boxplot() + geom_point() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6))
show(boxplot1)
```

```{r}
#Formatage des données pour les 10 premiers gènes
data_boxP_all <- data.frame()
for(gene in top_genes){
  data_boxP_all <- rbind(data_boxP_all,data.frame("zscore"=vsd_z[gene,], "conditions"=colData(dds4)[colnames(vsd_z),]$conditions, 'genes'=gene))
}
```

```{r}
#Boxplot de tous les gènes par conditions
boxplot2 <- ggplot(data_boxP_all, aes(x=genes,y=zscore, color=conditions)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6))
show(boxplot2)
```

#### 2.4.1.5 Enrichissement des gènes différentillement exprimé (GO enrichment)

-   Enrichissmeent effectué uniquement sur les gènes différentiellement exprimés
-   Utilisation de tout les gènes exprimés comme "backgroud"

```{r, fig.width=10, fig.height=12}
library(gprofiler2)
#Get enriched pathways
background <- row.names(df_pnl_ctrl)
gostres <- gost(query = list('NLS_vs_CTRL_DEgenes'=genes_df_pnl_ctrl), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = background, 
                numeric_ns = "", sources = c('GO:MF','GO:BP','GO:CC','KEGG','REAC'), as_short_link = FALSE)
top40_id <- head(gostres$result[order(gostres$result$p_value),],40)$term_id
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
pp <- publish_gostplot(p, highlight_terms = top40_id[1:24])
```

```{r}
head(gostres$result[order(gostres$result$p_value),c('term_name','p_value')],20)
```

#### 2.4.1.6 Gene Set Enrichment Analysis : Analyse d'enrichissemment des gènes exprimé ordonné

-   Téléchargement de la base de donnée GMT :
    <https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/msigdb_v2023.2.Hs_files_to_download_locally.zip>
-   Création de la liste des gènes exprimés ordonnée

```{r}
path_gsea_db <- 'Data/c5.all.v2023.2.Hs.symbols.gmt'
#Cération de la liste de gènes ordonné en fonction de la pvalue et du signe de fold Change
list_ranked_genes <- sign(df_pnl_ctrl$log2FoldChange) * -log10(df_pnl_ctrl$pvalue)
names(list_ranked_genes) <- row.names(df_pnl_ctrl)
list_ranked_genes <- list_ranked_genes[order(list_ranked_genes,decreasing=TRUE)]
list_ranked_genes <- list_ranked_genes[!is.na(list_ranked_genes)]
list_ranked_genes[1:10]
```

-   Calcul de l'enrichissement des gènes

```{r}
library(fgsea)
library(stringr)
pathways_GO <- gmtPathways(path_gsea_db)
fgseaResGO <- fgsea(pathways_GO, list_ranked_genes)
#collP <- collapsePathways(fgseaResGO,pathways_GO,list_ranked_genes)
#fgseaResGO.main <- fgseaResGO[pathway %in% collP$mainPathways]
fgseaResGO.main <- fgseaResGO
```

```{r}
# Extraction des pathways significativement enrichis
fgseaResGO_subset_order <- fgseaResGO.main[padj<0.01][order(NES, decreasing=TRUE)]
fgseaResGO_subset_order[1:10,]
# Sélection des top 15 pathways positif et negatif
fgsea_toplot <- rbind(head(fgseaResGO_subset_order,15),tail(fgseaResGO_subset_order,15))
# Représentation des pathways
gsea_plot <- ggplot(fgsea_toplot, aes(reorder(pathway,NES), NES), labels=leadingEdge) + 
				geom_col(aes(fill=padj)) + coord_flip() + 
				labs(x='GO Pathways', y='Normalized Enrichment Score', title=str_wrap(paste('Top 15 Gene Ontology pathways enrichi par GSEA',sep=' '),60)) +
				theme(axis.text = element_text(size = 7, angle=30),plot.margin = margin(10, 10, 100, 10))
show(gsea_plot)
```

<br>

### 2.4.2 Comparaison entre CTRL et PL

#### 2.4.2.1 Extraction des résultats de comparaison et visualisation

```{r}
#Extraire les gènes différentiellement exprimé avec results()
```

```{r}
#PlotMA
```

<br>

```{r}
#Identification des gènes UP and DOWN

```

```{r}
#Nombre de gènes UP et DOWN
```

#### 2.4.2.2 Représentation en Volcano Plot

```{r, warning=FALSE}
#Volcano Plot
```

#### 2.4.2.3 Représentation Heatmap

```{r, fig.width=10, fig.height=13}
#Calcul du Z-score et Heatmap
```

#### 2.4.2.4 Représentation de la dispersion des échantillons en boxplot

```{r}
#Formatage des données pour les 20 premiers gènes

#Boxplot de tous les gènes par conditions

```

#### 2.4.2.5 Enrichissement des gènes différentillement exprimé (GO enrichment)

```{r, fig.width=10, fig.height=12}
library(gprofiler2)
#Enrichissement des pathways

```

#### 2.4.2.6 Gene Set Enrichment Analysis : Analyse d'enrichissemment des gènes exprimé ordonné

```{r}
#Cération de la liste de gènes ordonné en fonction de la pvalue et du signe de fold Change

```

```{r}
#FGSEA
```

```{r}
# Extraction des pathways significativement enrichis

# Sélection des top 15 pathways positif et negatif

# Représentation des pathways

```

### 2.4.3 Comparaison entre PL et PNL

#### 2.4.3.1 Extraction des résultats de comparaison et visualisation

```{r}
#Extraire les gènes différentiellement exprimé avec results()
```

```{r}
#PlotMA
```

<br>

```{r}
#Identification des gènes UP and DOWN

```

```{r}
#Nombre de gènes UP et DOWN
```

#### 2.4.3.2 Représentation en Volcano Plot

```{r, warning=FALSE}
#Volcano Plot
```

#### 2.4.3.3 Représentation Heatmap

```{r, fig.width=10, fig.height=13}
#Calcul du Z-score et Heatmap
```

#### 2.4.3.4 Représentation de la dispersion des échantillons en boxplot

```{r}
#Formatage des données pour les 20 premiers gènes

#Boxplot de tous les gènes par conditions

```

#### 2.4.3.5 Enrichissement des gènes différentillement exprimé (GO enrichment)

```{r, fig.width=10, fig.height=12}
library(gprofiler2)
#Enrichissement des pathways

```

#### 2.4.3.6 Gene Set Enrichment Analysis : Analyse d'enrichissemment des gènes exprimé ordonné

```{r}
#Cération de la liste de gènes ordonné en fonction de la pvalue et du signe de fold Change

```

```{r}
#FGSEA
```

```{r}
# Extraction des pathways significativement enrichis

# Sélection des top 15 pathways positif et negatif

# Représentation des pathways

```

####{ - }

### Références

<https://gitlab.univ-nantes.fr/bird_pipeline_registry/srp-pipeline/>
<https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html>
<https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#what-are-the-exact-steps-performed-by-deseq>
