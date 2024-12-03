# --------------------------------------------------- #
# Librerias ----------------------------------------- #
# --------------------------------------------------- #

rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("simpleaffy")

library(car)
library(ggplot2)
library(Rcpp)
library(Matrix)
library(RColorBrewer)
library(lattice)
library(dplyr)
library(org.Hs.eg.db)
library(affy)
library(GEOquery)
library(CLL)
library(affydata)
library(vsn)
library(gcrma)
library(latticeExtra)
library(oligo)
library(limma)
library(siggenes)
library(ggrepel)
library(apeglm)
library(hexbin)
library(gtools)
require(DESeq2)
library(genefilter)
library(corrplot)
library(viridis)
library(tidyverse)
library(ade4)
library(VennDiagram)
library(gplots)

setwd("~/Documentos/ProyectoGBM/expresion")
