# ##################################################################### # 
# --------------------------------------------------------------------- #
# Carga de datos ------------------------------------------------------ #
# --------------------------------------------------------------------- #
# ##################################################################### #
gset <- getGEO("GSE13276", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)[,-c(9,11)] # Quitando las replicas
View(ex)
# GBM=glioblastoma, TCi=tejido circundante, Co=control de materia blanca
colnames(ex) <- c("GBM1","GBM2","GBM3","GBM4","GBM5",
                  "TCi1","TCi2","TCi3","TCi4","TCi5",
                  "Co1","Co2","Co3")
ex <- cbind(ex,"Probe.Set.ID"=rownames(ex))

# ##################################################################### # 
# --------------------------------------------------------------------- #
# Simbolo de Gen ------------------------------------------------------ #
# --------------------------------------------------------------------- #
# ##################################################################### #
anot <- read.csv("02 HG-U133A_2.na36.annot.csv",skip=25,header=T)
anot <- data.frame(Probe.Set.ID=anot$Probe.Set.ID, 
                   gene=anot$Gene.Symbol, stringsAsFactors = F) # nombres Gene Symbol
dim(anot) #22277     2
length(unique(anot$Probe.Set.ID)) #22277
length(unique(anot$gene)) #13333
anot <- anot[!anot$gene=="---",] #quitar "---"
any(is.na(anot)==TRUE) #no hay NA's
dim(anot) #21055     2

ex_anot <- merge(ex, anot, by = "Probe.Set.ID")
dim(ex_anot) #21055    15
any(is.na(ex_anot)==TRUE) #no hay NA's
length(unique(ex_anot$Probe.Set.ID)) #21055
length(unique(ex_anot$gene)) #13332

ex <- aggregate(ex_anot[,-ncol(ex_anot)],
             by=list(ex_anot$gene), max) #maximo con la suma por fila
dim(ex) #13332    15
rownames(ex) <- ex[,1]
ex <- ex[,-c(1,2)]
dim(ex) #13332  13

for(i in 1:dim(ex)[2]){
  ex[,i] <- as.numeric(ex[,i])
}

# ##################################################################### # 
# --------------------------------------------------------------------- #
# Preprocesamiento ---------------------------------------------------- #
# --------------------------------------------------------------------- #
# ##################################################################### #
boxplot(ex, las=2, col=c(rep("red",5),rep("cyan",5),rep("green",3)))
corrplot.mixed(cor(ex))  #c1 posible molestia
plotDensities(ex, legend=F)

dd <- dist2(ex)
diag(dd) <- 0
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
legend <- list(top=list(fun=dendrogramGrob,args=list(x=dd.row,side="top")))
lp <- levelplot(dd[row.ord,row.ord],xlab="", ylab="",legend=legend, las=3,  labels=list(cex=0.2))
lp

meanSdPlot(as.matrix(ex), ranks=TRUE)

select <- (0==rowSums(ex<=0)) # filtro genes
length(select) #13332   no se hace filtro

# ##################################################################### # 
# --------------------------------------------------------------------- #
# Normalizacion ------------------------------------------------------- #
# --------------------------------------------------------------------- #
# ##################################################################### #
ex_norm = justvsn(as.matrix(ex))

boxplot(ex_norm, las=2, col=c(rep("red",5),rep("cyan",5),rep("green",3)))
corrplot.mixed(cor(ex_norm))
plotDensities(ex_norm, legend=F)

dd  <- dist2(ex_norm)
diag(dd) <- 0
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
legend <- list(top=list(fun=dendrogramGrob,args=list(x=dd.row,side="top")))
lp <- levelplot(dd[row.ord,row.ord],xlab="", ylab="",legend=legend, las=3,  labels=list(cex=0.2))
lp

meanSdPlot(as.matrix(ex_norm), ranks=TRUE)

# ##################################################################### # 
# --------------------------------------------------------------------- #
# Guardar y cargar datos finales -------------------------------------- #
# --------------------------------------------------------------------- #
# ##################################################################### #
# write.csv(ex_norm,"03 datF_GSE13276.csv")
ex_norm <- read.csv("03 datF_GSE13276.csv",row.names=1, 
                    colClasses=c("character",rep("numeric",13)))
apply(ex_norm, 2, class)

# ##################################################################### # 
# --------------------------------------------------------------------- #
# LIMMA --------------------------------------------------------------- #
# --------------------------------------------------------------------- #
# ##################################################################### #
design <-  cbind(GBM=c(1,1,1,1,1,0,0,0,0,0,0,0,0),
                 TCi=c(0,0,0,0,0,1,1,1,1,1,0,0,0),
                 Co =c(0,0,0,0,0,0,0,0,0,0,1,1,1)
)
design
dim(design) #13 3
colnames(design)

fit = lmFit(ex_norm, design)
contrast.matrix = makeContrasts(GBM-TCi, GBM-Co, levels=design)
fit = contrasts.fit(fit, contrast.matrix)

efit = eBayes(fit)

topTable(efit, coef="GBM - TCi", genelist=fit$genes[,7], sort.by="P")

results = decideTests(efit, adjust.method = "BH", p.value = 0.05)

x <- list(
  "GBM - TCi" = rownames(results[results[,1] != 0,1]), 
  "GBM - Co" = rownames(results[results[,2] != 0,2]), 
  "TOTAL" = rownames(results)
)

vennDiagram(results, 
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("orange", "cyan"),
            show.include = T)
results
dim(results[rowSums(results)!=0,])
results[rowSums(results)!=0,]

# ##################################################################### # 
# --------------------------------------------------------------------- #
# DEGs ---------------------------------------------------------------- #
# --------------------------------------------------------------------- #
# ##################################################################### #
write(rownames(results[rowSums(results)!=0,]), "04 DEGs_GSE13276.txt")

