library(Seurat)
library(SeuratData)
library(SeuratDisk)

#####upload data and edit####
setwd("~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test")
#dataset obtained from https://ftp.cngb.org/pub/CNSA/data3/CNP0001102/Single_Cell/CSE0000013/
covid <- readRDS(file = "~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test/covid.RDS")
View(covid@meta.data)
covid <- UpdateSeuratObject(covid)
# all cell types, only covid and flu 
covid.subset <- subset(covid, subset = Stage != "Ctrl")
covid.subset <- subset(covid.subset, subset = Stage != "IAV-Stage1")
covid.subset <- subset(covid.subset, subset = Stage != "IAV-Stage2")
covid.subset <- subset(covid.subset, subset = Stage != "IAV-Stage3")

#replace d# w/ days since sx onset, did this manually based on metadata authors provided
#replace cov 1 w/ correct days 
covid.subset@meta.data$batch <- as.character(covid.subset@meta.data$batch)
indices <- grepl("^COV-1-D1$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-1-D7"
indices <- grepl("^COV-1-D4$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-1-D10"
indices <- grepl("^COV-1-D16$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-1-D22"
#replace cov 2 w/ correct days 
indices <- grepl("^COV-2-D16$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-2-D22"
indices <- grepl("^COV-2-D10$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-2-D16"
indices <- grepl("^COV-2-D7$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-2-D13"
indices <- grepl("^COV-2-D4$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-2-D10"
indices <- grepl("^COV-2-D1$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-2-D7"
#replace cov 3 w/ correct days 
indices <- grepl("^COV-3-D1$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-3-D13"
indices <- grepl("^COV-3-D16$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-3-D28"
indices <- grepl("^COV-3-D4$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-3-D16"
#replace cov 4w/ correct days 
indices <- grepl("^COV-4-D4$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-4-D13"
indices <- grepl("^COV-4-D16$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-4-D25"
#replace cov 5 w/ correct days 
indices <- grepl("^COV-5-D1$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-5-D3"
indices <- grepl("^COV-5-D13$",covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-5-D15"
indices <- grepl("^COV-5-D7$", covid.subset@meta.data$batch)
covid.subset@meta.data$batch[indices] <- "COV-5-D9"

# Convert "batch" column back to factor if desired
covid.subset@meta.data$batch <- as.factor(covid.subset@meta.data$batch)
#View(covid.subset@meta.data)
covid.subset.time <- covid.subset

covid.subset.time[["time"]] <- covid.subset.time$batch
covid.subset.time@meta.data$time <- sub(".*-D", "", covid.subset.time@meta.data$time)
View(covid.subset.time@meta.data)

DimPlot(covid.subset.time, group.by = "batch")


####subset by time d3 and make anndata for seacels####
covid.subset.time.d3 <- subset(covid.subset.time, subset = time == "3")
covid.subset.time.d3[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d3@meta.data, is.factor)
covid.subset.time.d3@meta.data[i] <- lapply(covid.subset.time.d3@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d3,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d3),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid3.h5Seurat")
Convert("covid3.h5Seurat", dest = "h5ad")

####subset by time d7 and make anndata for seacels####
covid.subset.time.d7 <- subset(covid.subset.time, subset = time == "7")
covid.subset.time.d7[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d7@meta.data, is.factor)
covid.subset.time.d7@meta.data[i] <- lapply(covid.subset.time.d7@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d7,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d7),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid7.h5Seurat")
Convert("covid7.h5Seurat", dest = "h5ad")






####subset by time d9 and make anndata for seacels####
covid.subset.time.d9 <- subset(covid.subset.time, subset = time == "9")
covid.subset.time.d9[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d9@meta.data, is.factor)
covid.subset.time.d9@meta.data[i] <- lapply(covid.subset.time.d9@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d9,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d9),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid9.h5Seurat")
Convert("covid9.h5Seurat", dest = "h5ad")

####subset by time d10 and make anndata for seacels####
covid.subset.time.d10 <- subset(covid.subset.time, subset = time == "10")
covid.subset.time.d10[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d10@meta.data, is.factor)
covid.subset.time.d10@meta.data[i] <- lapply(covid.subset.time.d10@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d10,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d10),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid10.h5Seurat")
Convert("covid10.h5Seurat", dest = "h5ad")




####subset by time d13 and make anndata for seacels####
covid.subset.time.d13 <- subset(covid.subset.time, subset = time == "13")
covid.subset.time.d13[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d13@meta.data, is.factor)
covid.subset.time.d13@meta.data[i] <- lapply(covid.subset.time.d13@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d13,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d13),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid13.h5Seurat")
Convert("covid13.h5Seurat", dest = "h5ad")

####subset by time d16 and make anndata for seacels####
covid.subset.time.d16 <- subset(covid.subset.time, subset = time == "16")
covid.subset.time.d16[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d16@meta.data, is.factor)
covid.subset.time.d16@meta.data[i] <- lapply(covid.subset.time.d16@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d16,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d16),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid16.h5Seurat")
Convert("covid16.h5Seurat", dest = "h5ad")

####subset by time d15 and make anndata for seacels####
covid.subset.time.d15 <- subset(covid.subset.time, subset = time == "15")
covid.subset.time.d15[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d15@meta.data, is.factor)
covid.subset.time.d15@meta.data[i] <- lapply(covid.subset.time.d15@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d15,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d15),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid15.h5Seurat")
Convert("covid15.h5Seurat", dest = "h5ad")

####subset by time d22 and make anndata for seacels####
covid.subset.time.d22 <- subset(covid.subset.time, subset = time == "22")
covid.subset.time.d22[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d22@meta.data, is.factor)
covid.subset.time.d22@meta.data[i] <- lapply(covid.subset.time.d22@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d22,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d22),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid22.h5Seurat")
Convert
####subset by time d25 and make anndata for seacels####
covid.subset.time.d25 <- subset(covid.subset.time, subset = time == "25")
covid.subset.time.d25[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d25@meta.data, is.factor)
covid.subset.time.d25@meta.data[i] <- lapply(covid.subset.time.d25@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d25,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d25),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid25.h5Seurat")
Convert("covid25.h5Seurat", dest = "h5ad")
####subset by time d28 and make anndata for seacels####
covid.subset.time.d28 <- subset(covid.subset.time, subset = time == "28")
covid.subset.time.d28[["RNA"]] <- NULL
i <- sapply(covid.subset.time.d28@meta.data, is.factor)
covid.subset.time.d28@meta.data[i] <- lapply(covid.subset.time.d28@meta.data[i], as.character)
covid.update <- DietSeurat(covid.subset.time.d28,
                           counts = T,
                           data = T,
                           scale.data = F,
                           features = rownames(covid.subset.time.d28),
                           assays = "integrated",
                           dimreducs = c("pca", "umap"),
                           misc = T,
                           graphs = c("RNA_nn", "RNA_snn"))
covid.update <- UpdateSeuratObject(object = covid.update)
SaveH5Seurat(covid.update, filename = "covid28.h5Seurat")
Convert("covid28.h5Seurat", dest = "h5ad")

#### after seacells
#### merge seacells by time point (after creating seacells in python) ####
library(Seurat)
#for all time points
# Define the time points to iterate over
time_points <- c("3", "7", "9", "10", "13", "15", "16", "22", "25", "28")

# Iterate over the time points
for (time_point in time_points) {
  covid.subset.time.loop <- covid.subset.time
  # Read seacells file
  seacell_file <- paste0("seacells.d", time_point, ".cells.csv")
  seacell.assign <- as.data.frame(read.csv(seacell_file, row.names = 1))
  
  # Subset seurat object
  covid.subset.time.loop <- subset(covid.subset.time.loop, subset = time == time_point)
  
  # Add seacells as metadata
  covid.subset.time.loop[["seacells"]] <- seacell.assign
  
  # Aggregate expression by metacells
  Idents(covid.subset.time.loop) <- covid.subset.time.loop$cell_type
  
  DefaultAssay(covid.subset.time.loop) <- "RNA"
  covid.subset.time.loop[["integrated"]] <- NULL
  
  # Calculate average expression
  covid.pb <- AverageExpression(covid.subset.time.loop, group.by = "seacells",
                                return.seurat = TRUE, slot = "data")
  
  # Retrieve the data matrix
  covid.pb.data <- as.matrix(covid.pb@assays$RNA$data)
  colnames(covid.pb.data) <- paste("d", time_point, colnames(covid.pb.data), sep = "_")
  
  #create covid.pb.d#
  assign(paste("covid.pb.d", time_point, sep = ""), covid.pb.data)
}
colnames(covid.pb.d3)
colnames(covid.pb.d7)
colnames(covid.pb.d9)
colnames(covid.pb.d10)
colnames(covid.pb.d13)
colnames(covid.pb.d15)
colnames(covid.pb.d16)
colnames(covid.pb.d22)
colnames(covid.pb.d25)
colnames(covid.pb.d28)

# Merge the data matrices
merged.data <- cbind(covid.pb.d3, covid.pb.d7, covid.pb.d9, covid.pb.d10,
                     covid.pb.d13, covid.pb.d15, covid.pb.d16, covid.pb.d22,
                     covid.pb.d25, covid.pb.d28)
merged.data
View(head(merged.data))
write.csv(head(merged.data), "merged.data.all.cells.csv")
covid.masigpro <- merged.data

#find avg number of genes detected
covid.masigpro.seurat <- CreateSeuratObject(covid.masigpro)
mean(covid.masigpro.seurat$nFeature_RNA)

#find average variance
gene_variances.seacell <- apply(covid.masigpro, 1, sd)
average_variance <- mean(gene_variances.seacell)

####masigpro all metacells together #####

#for all cells...design matrix made in excel
exp.data <- read.csv("all.cells.masigpro.design.2.csv")
rownames(exp.data) <- exp.data[,1]
exp.data <- exp.data[,-1]
exp.data <- exp.data[,-3]
exp.data <- exp.data[,-(4:18)]
View(exp.data)

#define regression model
design <- make.design.matrix(exp.data, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

#masigpro for rna seq 
library(MASS)
NBp <- p.vector(covid.masigpro, design, counts=F, 
                min.obs = 6)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

see.genes <- see.genes(get$sig.genes, k= 14)

View(as.data.frame(see.genes$cut))


#### masigpro by cell type ####

# Create  empty lists to store data frames
exp.data.list <- list()
covid.masigpro.list <- list()
design.list <- list()

# Loop through cell types
for (type in cell.types) {
  # Read the CSV file
  exp.data <- read.csv("all.cells.masigpro.design.2.csv")
  
  # Set row names
  rownames(exp.data) <- exp.data[, 1]
  
  # Remove the first column
  exp.data <- exp.data[, -1]
  
  # Subset the data based on cell type
  exp.data <- subset(exp.data, cell.type == type)
  
  # Store the data frame in the list
  exp.data.list[[type]] <- exp.data
  
}

for (type in cell.types) {
  exp.data.type <- exp.data.list[[type]]
  
  # Select only numeric columns
  numeric_cols <- sapply(exp.data.type, is.numeric)
  exp.data.type <- exp.data.type[, numeric_cols]
  
  # Remove columns that sum to zero
  exp.data.type <- exp.data.type[, colSums(exp.data.type) != 0]
  
  # Remove specific columns (e.g., columns 3 and 4)
  exp.data.type <- exp.data.type[, -c(3, 4)]
  
  # Add column cell type with all values equal to 1
  exp.data.type[, type] <- 1
  
  # Subset covid.masigpro based on column names
  covid.masigpro.type <- covid.masigpro[, rownames(exp.data.type)]
  
  # Store the data frames in the lists
  covid.masigpro.list[[type]] <- covid.masigpro.type
  
  exp.data.list[[type]] <- exp.data.type
  
}

for (type in cell.types) {
  #define regression model
  exp.data.type <- exp.data.list[[type]]
  covid.masigpro.type <- covid.masigpro.list[[type]]
  design <- make.design.matrix(exp.data.type, degree = 2, time.col = 1,
                               repl.col = 2, group.cols = 3)
  design.list[[type]] <- design
}

#Naive T cells
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Naive T cells`, design.list$`Naive T cells`, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get.naive.t.cells <-get.siggenes(NBt, vars="all", rsq = 0.5)
get.naive.t.cells$summary
see.gene <-  see.genes(get.naive.t.cells$sig.genes , k = 19, legend = T,
                       cex.legend = .8, )
View(as.data.frame(see.gene$cut))
get.naive.t.cells$summary[get.naive.t.cells$summary %in% get.naive.t.cells$summary]

#Activated CD4 T Cells
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Activated CD4 T cells`, design.list$`Activated CD4 T cells`, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get.activated.cd4.t.cells <-get.siggenes(NBt, vars="all", rsq = 0.5)
get.activated.cd4.t.cells$summary
see.genes <-  see.genes(get.activated.cd4.t.cells$sig.genes , k = 5, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# NKs
library(MASS)
NBp <- p.vector(covid.masigpro.list$NKs, design.list$`NKs`, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get.NKs <-get.siggenes(NBt, vars="all", rsq = 0.5)
get.NKs$summary
see.genes <-  see.genes(get.NKs$sig.genes , k=6, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# `Naive B cells`
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Naive B cells`, design.list$`Naive B cells`, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.naive.b.cells<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.naive.b.cells$summary
see.genes <-  see.genes(get.naive.b.cells$sig.genes, 
                        k = 6, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# MAIT
library(MASS)
NBp <- p.vector(covid.masigpro.list$MAIT, design.list$MAIT, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get.MAIT<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.MAIT$summary
see.genes <-  see.genes(get.MAIT$sig.genes, 
                        k = 2, legend = F,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# cytotoxic cd8 t cells
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Cytotoxic CD8 T cells`,
                design.list$`Cytotoxic CD8 T cells`, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.Cytotoxic.CD8.T.cells<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.Cytotoxic.CD8.T.cells$summary
ssee.genes <-  see.genes(get.Cytotoxic.CD8.T.cells$sig.genes, 
                        k = 4, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# memory b cells
library(MASS)

design.list$`Memory B cells` <- make.design.matrix(exp.data.list$`Memory B cells`, degree = 2, time.col = 1,
                                                   repl.col = 2, group.cols = 3)
NBp <- p.vector(covid.masigpro.list$`Memory B cells`,
                design.list$`Memory B cells`, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.Memory.B.cells<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.Memory.B.cells$summary
see.genes <-  see.genes(get.Memory.B.cells$sig.genes, 
                        k = 9, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# cycling t cells
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Cycling T cells`,
                design.list$`Cycling T cells`, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.Cycling.T.cells<-get.siggenes(NBt, vars="all", rsq = 0.7)
get.Cycling.T.cells$summary
see.genes <-  see.genes(get.Cycling.T.cells$sig.genes, 
                        k = 2,legend = F,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# xcl+ nks
library(MASS)
NBp <- p.vector(covid.masigpro.list$`XCL+ NKs`,
                design.list$`XCL+ NKs`, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.XCL.NKs<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.XCL.NKs$summary
see.genes <-  see.genes(get.XCL.NKs$sig.genes, 
                        k = 4, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# plasma
library(MASS)
NBp <- p.vector(covid.masigpro.list$Plasma,
                design.list$Plasma, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.Plasma<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.Plasma$summary
see.genes <-  see.genes(get.Plasma$sig.genes, 
                        k = 7, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# monocytes
library(MASS)
NBp <- p.vector(covid.masigpro.list$Monocytes,
                design.list$Monocytes, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.Monocytes<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.Monocytes$summary
see.genes <-  see.genes(get.Monocytes$sig.genes, 
                        k = 7, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# stem cells
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Stem cells`,
                design.list$`Stem cells`, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.stem.cells<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.stem.cells$summary
see.genes <-  see.genes(get.stem.cells$sig.genes, 
                        k = 4, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# Megakaryocytes
library(MASS)
NBp <- p.vector(covid.masigpro.list$Megakaryocytes,
                design.list$Megakaryocytes, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.Megakaryocytes<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.Megakaryocytes$summary
see.genes <-  see.genes(get.Megakaryocytes$sig.genes, 
                        k = 4, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# cycling plasma
library(MASS)
NBp <- p.vector(covid.masigpro.list$`Cycling Plasma`,
                design.list$`Cycling Plasma`, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get.Cycling.Plasma<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.Cycling.Plasma$summary
see.genes <-  see.genes(get.Cycling.Plasma$sig.genes, 
                        k = 4, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

# DCs
library(MASS)
NBp <- p.vector(covid.masigpro.list$DCs,
                design.list$DCs, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.DCs<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.DCs$summary
see.genes <-  see.genes(get.DCs$sig.genes, 
                        k = 9, legend = T,
                        cex.legend = 1)
View(as.data.frame(see.genes$cut))

#summarize #unique DEGs across seleected cells types
unique_row_names <- unique(c(rownames(get.activated.cd4.t.cells$sig.genes$sig.pvalues),
       rownames(get.Cytotoxic.CD8.T.cells$sig.genes$sig.pvalues),
       rownames(get.naive.b.cells$sig.genes$sig.pvalues),
       rownames(get.naive.t.cells$sig.genes$sig.pvalues),
       rownames(get.NKs$sig.genes$sig.pvalues),
       rownames(get.Plasma$sig.genes$sig.pvalues),
       rownames(get.Memory.B.cells$sig.genes$sig.pvalues),
       rownames(get.XCL.NKs$sig.genes$sig.pvalues)
       ))

View(data.frame(unique_row_names))

# code below was used to make line plot for sig genes by cell type
# conditional statements checked whether trend was maxima or minima
# after that, checked get."insert cell type name"$sig.genes$coefficients, not in code but can add in later
# if coefficient p value significant ie not linear, determine trend visually
#for naive t cells####
library(ggpmisc)
library(gridExtra)
library(dplyr)
#get sig gene list for naive t cells
genes <- get.naive.t.cells$summary
#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.naive.t.cells$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data.list$`Naive T cells`$time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#determine whether the averages of the first and last time points are greater than all points between them == minima
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#determine whether the averages of the first and last time points are less than all points between them == maxima
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots[1:36])


#for naive b cells####
genes <- get.naive.b.cells$summary
by.type <-as.data.frame(t(get.naive.b.cells$sig.genes$sig.profiles[genes,]))
by.type$time <- exp.data.list$`Naive B cells`$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots[1:36])

get.naive.b.cells

#for activated cd4 t cells####
genes <- get.activated.cd4.t.cells$summary
by.type <-as.data.frame(t(get.activated.cd4.t.cells$sig.genes$sig.profiles[genes,]))
by.type$time <- exp.data.list$`Activated CD4 T cells`$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#determine if overfit from d3, eliminate these genes downstream
results1 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > 0
      &averages[2, col] ==0
      &averages[3, col] ==0
      &averages[4, col] ==0
      &averages[5, col] ==0
      &averages[6, col] ==0
      &averages[7, col] ==0
      &averages[8, col] ==0
      &averages[9, col] ==0
      &averages[10, col] ==0) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      ){
    print(paste("T for column", col))
  } 
}


plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
plots1 = plots[1:36]
plots2 = plots[37:72]
plots3 = plots[73:108]
plots4 = plots[109:144]
gridExtra::grid.arrange(grobs = plots1)
gridExtra::grid.arrange(grobs = plots2)
gridExtra::grid.arrange(grobs = plots3)
gridExtra::grid.arrange(grobs = plots4)

#for cycling t cells ####
genes <- get.Cycling.T.cells$summary
by.type <-as.data.frame(t(get.Cycling.T.cells$sig.genes$sig.profiles[genes,]))
by.type$time <- exp.data.list$`Cycling T cells`$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[9, col] > averages[2, col]
      &averages[9, col] > averages[3, col]
      &averages[9, col] > averages[4, col]
      &averages[9, col] > averages[5, col]
      &averages[9, col] > averages[6, col]
      &averages[9, col] > averages[7, col]
      &averages[9, col] > averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[9, col] < averages[2, col]
      &averages[9, col] < averages[3, col]
      &averages[9, col] < averages[4, col]
      &averages[9, col] < averages[5, col]
      &averages[9, col] < averages[6, col]
      &averages[9, col] < averages[7, col]
      &averages[9, col] < averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#determine if overfit from d3 or d28, eliminate these genes downstream
results1 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > 0
      &averages[2, col] ==0
      &averages[3, col] ==0
      &averages[4, col] ==0
      &averages[5, col] ==0
      &averages[6, col] ==0
      &averages[7, col] ==0
      &averages[8, col] ==0
      &averages[9, col] ==0
      &averages[10, col] ==0) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
  ){
    print(paste("T for column", col))
  } 
}

results1 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] ==0
      &averages[2, col] ==0
      &averages[3, col] ==0
      &averages[4, col] ==0
      &averages[5, col] ==0
      &averages[6, col] ==0
      &averages[7, col] ==0
      &averages[8, col] ==0
      &averages[9, col] ==0
      &averages[10, col] >0) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
  ){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("expression")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots)

plots1 = plots[1:36]
plots2 = plots[37:72]
gridExtra::grid.arrange(grobs = plots1)
gridExtra::grid.arrange(grobs = plots2)

#for cytotoxic cd8 ####
genes <- get.Cytotoxic.CD8.T.cells
by.type <-as.data.frame(t(get.Cytotoxic.CD8.T.cells$sig.genes$sig.profiles))
by.type$time <- exp.data.list$`Cytotoxic CD8 T cells`$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots)

plots1 = plots[1:36]
plots2 = plots[37:72]
gridExtra::grid.arrange(grobs = plots1)
gridExtra::grid.arrange(grobs = plots2)

#for DCs ####
genes <- get.DCs
by.type <-as.data.frame(t(get.DCs$sig.genes$sig.profiles))
by.type$time <- exp.data.list$DCs$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[5, col] > averages[2, col]
      &averages[5, col] > averages[3, col]
      &averages[5, col] > averages[4, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[5, col] < averages[2, col]
      &averages[5, col] < averages[3, col]
      &averages[5, col] < averages[4, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-40], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("expression")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots)

#for MAIT ####
genes <- get.MAIT
by.type <-as.data.frame(t(get.MAIT$sig.genes$sig.profiles))
by.type$time <- exp.data.list$MAIT$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("expression")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots)

plots1 = plots[1:36]
plots2 = plots[37:72]
gridExtra::grid.arrange(grobs = plots1)
gridExtra::grid.arrange(grobs = plots2)

#for mega ####
genes <- get.Megakaryocytes
by.type <-as.data.frame(t(get.Megakaryocytes$sig.genes$sig.profiles))
by.type$time <- exp.data.list$Megakaryocytes$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-54], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("expression")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots)

#for mem b ####
genes <- get.Memory.B.cells$summary
by.type <-as.data.frame(t(get.Memory.B.cells$sig.genes$sig.profiles))
by.type$time <- exp.data.list$`Memory B cells`$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[7, col] > averages[2, col]
      &averages[7, col] > averages[3, col]
      &averages[7, col] > averages[4, col]
      &averages[7, col] > averages[5, col]
      &averages[7, col] > averages[6, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[7, col] < averages[2, col]
      &averages[7, col] < averages[3, col]
      &averages[7, col] < averages[4, col]
      &averages[7, col] < averages[5, col]
      &averages[7, col] < averages[6, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#check early time point outliers

for (col in colnames(averages)) {
  if (averages[1, col] > 0
      &averages[2, col] ==0
      &averages[3, col] ==0
      &averages[4, col] ==0
      &averages[5, col] ==0
      &averages[6, col] ==0
      &averages[7, col] ==0
      &averages[8, col] ==0
      &averages[9, col] ==0
      &averages[10, col] ==0) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
  ){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("expression")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots[1:36])

#for NKs ####
genes <- get.NKs$summary
by.type <-as.data.frame(t(get.NKs$sig.genes$sig.profiles))
by.type$time <- exp.data.list$NKs$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[9, col] > averages[2, col]
      &averages[9, col] > averages[3, col]
      &averages[9, col] > averages[4, col]
      &averages[9, col] > averages[5, col]
      &averages[9, col] > averages[6, col]
      &averages[9, col] > averages[7, col]
      &averages[9, col] > averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[9, col] < averages[2, col]
      &averages[9, col] < averages[3, col]
      &averages[9, col] < averages[4, col]
      &averages[9, col] < averages[5, col]
      &averages[9, col] < averages[6, col]
      &averages[9, col] < averages[7, col]
      &averages[9, col] < averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots[1:36])

#for Plasma ####
genes <- get.Plasma
by.type <-as.data.frame(t(get.Plasma$sig.genes$sig.profiles))
by.type$time <- exp.data.list$Plasma$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
      &averages[1, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] > averages[2, col]
      &averages[10, col] > averages[3, col]
      &averages[10, col] > averages[4, col]
      &averages[10, col] > averages[5, col]
      &averages[10, col] > averages[6, col]
      &averages[10, col] > averages[7, col]
      &averages[10, col] > averages[8, col]
      &averages[10, col] > averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
      &averages[1, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[10, col] < averages[2, col]
      &averages[10, col] < averages[3, col]
      &averages[10, col] < averages[4, col]
      &averages[10, col] < averages[5, col]
      &averages[10, col] < averages[6, col]
      &averages[10, col] < averages[7, col]
      &averages[10, col] < averages[8, col]
      &averages[10, col] < averages[9, col]) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots[1:36])

#for stem 
genes <- get.stem.cells
by.type <-as.data.frame(t(get.stem.cells$sig.genes$sig.profiles))
by.type$time <- exp.data.list$`Stem cells`$time
View(by.type)

plots <- lapply(names(by.type)[-8], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("expression")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots)

#for xcl nks ####
genes <- get.XCL.NKs
by.type <-as.data.frame(t(get.XCL.NKs$sig.genes$sig.profiles))
by.type$time <- exp.data.list$`XCL+ NKs`$time
View(by.type)
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))

#for min
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] > averages[2, col]
      &averages[1, col] > averages[3, col]
      &averages[1, col] > averages[4, col]
      &averages[1, col] > averages[5, col]
      &averages[1, col] > averages[6, col]
      &averages[1, col] > averages[7, col]
      &averages[1, col] > averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[9, col] > averages[2, col]
      &averages[9, col] > averages[3, col]
      &averages[9, col] > averages[4, col]
      &averages[9, col] > averages[5, col]
      &averages[9, col] > averages[6, col]
      &averages[9, col] > averages[7, col]
      &averages[9, col] > averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
} 

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

#for max
results1 <- list()
results2 <- list()
for (col in colnames(averages)) {
  if (averages[1, col] < averages[2, col]
      &averages[1, col] < averages[3, col]
      &averages[1, col] < averages[4, col]
      &averages[1, col] < averages[5, col]
      &averages[1, col] < averages[6, col]
      &averages[1, col] < averages[7, col]
      &averages[1, col] < averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results1[[col]] <- result
}

for (col in colnames(averages)) {
  if (averages[9, col] < averages[2, col]
      &averages[9, col] < averages[3, col]
      &averages[9, col] < averages[4, col]
      &averages[9, col] < averages[5, col]
      &averages[9, col] < averages[6, col]
      &averages[9, col] < averages[7, col]
      &averages[9, col] < averages[8, col]
  ) {
    result <- (paste("T for column", col))
  } else {
    result <- (paste ("F for column", col))
  }
  results2[[col]] <- result
}

for (col in colnames(averages)) {
  if (results1[[col]] == paste("T for column", col)
      &results2[[col]] == paste("T for column", col)){
    print(paste("T for column", col))
  } 
}

plots <- lapply(names(by.type)[-ncol(by.type)], function(col) {
  ggplot(by.type, aes(x = time, y = .data[[col]])) +
    geom_point() +
    labs(title = col)+
    xlab("")+
    ylab("")+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
    scale_x_continuous(n.breaks = 10)+
    stat_poly_eq(parse=T, aes(label = ..eq.label..), 
                 formula=y ~ x + I(x^2), size = 2.5, color = "red")
  
})
gridExtra::grid.arrange(grobs = plots[1:36])


##### dot plot for all sig genes####
library(readxl)
library(ggplot2)
library(org.Hs.eg.db)

#excel sheet created based on results from trend clasification
all.p <- read_xlsx("time.series.summary.celltype.xlsx", sheet = "by.protein.class.update")
all.p$Classification <- factor(all.p$Classification, levels = c("immune", "transcription", "protein degradation",
                                                                "ribosome/translation", "transferase", "other"))

all.p.plot <- ggplot(all.p, aes(x = Cell.Type, y = Gene, color = X.log.p.value, size = R2))+
  geom_point()+
  facet_grid(Classification ~ Pattern, space = "free", scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12),
        strip.text.y = element_text(angle = 0, size = 12),
        strip.text.x = element_text(size= 12)
  )+
  xlab("Cell Type")+
  scale_size_continuous(range = c(1,4))+
  ggtitle("Trend")+
  labs(color = "-log(p-value)")+ labs(size = "R^2")


ggsave('all.p.plot.update.pdf', all.p.plot, device = "pdf", width = 20, height = 40)

all.p.plot <- ggplot(all.p, aes(x = Cell.Type, fill = Pattern))+
  geom_bar()+
  scale_fill_manual(values = c("darkred","red", "pink","darkblue","blue","lightblue","darkgreen","orange"))+
  labs(x = "Cell Type", y = "Count")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12),
        strip.text.y = element_text(angle = 0, size = 12),
        strip.text.x = element_text(size= 12))

ggsave("cell.type.deg.count.pdf", all.p.plot, device = "pdf", width = 10, height = 40)


# gather control data to plot expression trend with controls #####

covid.subset.controls <- subset(covid, Stage == "Ctrl")
DefaultAssay(covid.subset.controls) <- "RNA"
covid.subset.controls[["integrated"]] <- NULL
covid.subset.controls <- SplitObject(covid.subset.controls, split.by = "batch")

# Calculate average expression for each control
covid.pb.controls1<- AverageExpression(covid.subset.controls$`Ctrl-1`, group.by = "cell_type",
                                       return.seurat = TRUE, slot = "data")
covid.pb.controls1<- covid.pb.controls1@assays$RNA$data
covid.pb.controls2<- AverageExpression(covid.subset.controls$`Ctrl-2`, group.by = "cell_type",
                                       return.seurat = TRUE, slot = "data")
covid.pb.controls2<- covid.pb.controls2@assays$RNA$data
covid.pb.controls3<- AverageExpression(covid.subset.controls$`Ctrl-3`, group.by = "cell_type",
                                       return.seurat = TRUE, slot = "data")
covid.pb.controls3<- covid.pb.controls3@assays$RNA$data


##### plot expression trend with controls for non covid annotated genes in string network ####
#although labeled covid genes, these are actually the genes not annotated in KEGG
library(ggsignif)
covid.genes <- c("EIF3H", "OASL", "RSAD2", "ISG20", "IFI16", "IFIT1",
                 "CDK9", "CX3CR1", "TNFRSF13C", "LGALS9", "IFIT3", "IRF7", "XAF1", "IFI27")
exp.data <- read.csv("all.cells.masigpro.design.2.csv")

#for each gene, cell type color was defined only for cell types who had >= 1 sig gene
#did not automate this process but can if we want to explore specific sets of genes
# for EIF3H ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE"#,
  #"Cylcing T cells" = "#CC6677",
  #"Cytotoxic CD8 T cells" = "#DDCC77"#,
  # "Naive B cells" = "#117733",
  #"Naive T cells" = "#332288",
  #"NKs" = "#AA4499",
  # "Plasma"= "#44AA99",
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

#for EIF3H
a <-ggplot(by.type, aes(x = time, y = EIF3H, 
                                      color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$EIF3H, covid.pb.controls.sig.genes$EIF3H)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = EIF3H, color= cell.type))+
  geom_point(size=2)+
  ylab("EIF3H Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$EIF3H,covid.pb.controls.sig.genes$EIF3H)),labels = scales::number_format(accuracy = 0.01))

EIF3H <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                     padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$EIF3H,
       subset(b$data, cell.type =="Activated CD4 T cells")$EIF3H,
       alternative = "two.sided", paired = F, var.equal = T)

# for OASL ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE"#,
                 #"Cylcing T cells" = "#CC6677",
                 #"Cytotoxic CD8 T cells" = "#DDCC77"#,
                 # "Naive B cells" = "#117733",
                 "Naive T cells" = "#332288"#,
                 #"NKs" = "#AA4499",
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Naive T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Naive T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Naive T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Naive T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = OASL, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$OASL, covid.pb.controls.sig.genes$OASL)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = OASL, color= cell.type))+
  geom_point(size=2)+
  ylab("OASL Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$OASL,covid.pb.controls.sig.genes$OASL)),labels = scales::number_format(accuracy = 0.01))

OASL <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))


# t test
t.test(subset(a$data, time == "3" & cell.type =="Naive T cells")$OASL,
       subset(b$data, cell.type =="Naive T cells")$OASL,
       alternative = "two.sided", paired = F, var.equal = T)


# for RSAD2 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE",
  #"Cylcing T cells" = "#CC6677",
  "Cytotoxic CD8 T cells" = "#DDCC77",
   "Naive B cells" = "#117733",
  #"Naive T cells" = "#332288"#,
  "NKs" = "#AA4499"#,
  # "Plasma"= "#44AA99",
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = RSAD2, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$RSAD2, covid.pb.controls.sig.genes$RSAD2)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = RSAD2, color= cell.type))+
  geom_point(size=2)+
  ylab("RSAD2 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$RSAD2,covid.pb.controls.sig.genes$RSAD2)),labels = scales::number_format(accuracy = 0.01))

RSAD2 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                     padding = unit(0, "lines"))


# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$RSAD2,
       subset(b$data, cell.type =="Activated CD4 T cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$RSAD2,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Naive B cells")$RSAD2,
       subset(b$data, cell.type =="Naive B cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="NKs")$RSAD2,
       subset(b$data, cell.type =="NKs")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)

# for ISG20 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE",
                 #"Cylcing T cells" = "#CC6677",
                 "Cytotoxic CD8 T cells" = "#DDCC77",
                 #"Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 "NKs" = "#AA4499"#,
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "NKs"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "NKs"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "NKs"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "NKs"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = ISG20, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$ISG20, covid.pb.controls.sig.genes$ISG20)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = ISG20, color= cell.type))+
  geom_point(size=2)+
  ylab("ISG20 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$ISG20,covid.pb.controls.sig.genes$ISG20)),labels = scales::number_format(accuracy = 0.01))

ISG20 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))


# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$ISG20,
       subset(b$data, cell.type =="Activated CD4 T cells")$ISG20,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$ISG20,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$ISG20,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="NKs")$ISG20,
       subset(b$data, cell.type =="NKs")$ISG20,
       alternative = "two.sided", paired = F, var.equal = T)

# for IFI16 ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE",
                 #"Cylcing T cells" = "#CC6677",
                 "Cytotoxic CD8 T cells" = "#DDCC77"#,
                 #"Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 #"NKs" = "#AA4499"#,
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IFI16, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFI16, covid.pb.controls.sig.genes$IFI16)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IFI16, color= cell.type))+
  geom_point(size=2)+
  ylab("IFI16 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFI16,covid.pb.controls.sig.genes$IFI16)),labels = scales::number_format(accuracy = 0.01))

IFI16 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))


# t test
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$IFI16,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$IFI16,
       alternative = "two.sided", paired = F, var.equal = T)


# for IFIT1 ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE",
  #"Cylcing T cells" = "#CC6677",
  "Cytotoxic CD8 T cells" = "#DDCC77"#,
  #"Naive B cells" = "#117733",
  #"Naive T cells" = "#332288"#,
  #"NKs" = "#AA4499"#,
  # "Plasma"= "#44AA99",
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IFI16, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFI16, covid.pb.controls.sig.genes$IFI16)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IFI16, color= cell.type))+
  geom_point(size=2)+
  ylab("IFI16 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFI16,covid.pb.controls.sig.genes$IFI16)),labels = scales::number_format(accuracy = 0.01))

IFI16 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))


# t test
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$IFI16,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$IFI16,
       alternative = "two.sided", paired = F, var.equal = T)


# for IFIT1 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE",
                 #"Cylcing T cells" = "#CC6677",
                 "Cytotoxic CD8 T cells" = "#DDCC77",
                 "Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 "NKs" = "#AA4499"#,
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IFIT1, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFIT1, covid.pb.controls.sig.genes$IFIT1)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IFIT1, color= cell.type))+
  geom_point(size=2)+
  ylab("IFIT1 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFIT1,covid.pb.controls.sig.genes$IFIT1)),labels = scales::number_format(accuracy = 0.01))

IFIT1 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$IFIT1,
       subset(b$data, cell.type =="Activated CD4 T cells")$IFIT1,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$IFIT1,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Naive B cells")$IFIT1,
       subset(b$data, cell.type =="Naive B cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="NKs")$IFIT1,
       subset(b$data, cell.type =="NKs")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)


# for CDK9 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE"#,
                 #"Cylcing T cells" = "#CC6677",
                 #"Cytotoxic CD8 T cells" = "#DDCC77",
                 #"Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 #"NKs" = "#AA4499"#,
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = CDK9, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$CDK9, covid.pb.controls.sig.genes$CDK9)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = CDK9, color= cell.type))+
  geom_point(size=2)+
  ylab("CDK9 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$CDK9,covid.pb.controls.sig.genes$CDK9)),labels = scales::number_format(accuracy = 0.01))

CDK9 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$CDK9,
       subset(b$data, cell.type =="Activated CD4 T cells")$CDK9,
       alternative = "two.sided", paired = F, var.equal = T)

# for IFIT1 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE",
                 #"Cylcing T cells" = "#CC6677",
                 "Cytotoxic CD8 T cells" = "#DDCC77",
                 "Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 "NKs" = "#AA4499"#,
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells",
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "NKs"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IFIT1, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFIT1, covid.pb.controls.sig.genes$IFIT1)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IFIT1, color= cell.type))+
  geom_point(size=2)+
  ylab("IFIT1 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFIT1,covid.pb.controls.sig.genes$IFIT1)),labels = scales::number_format(accuracy = 0.01))

IFIT1 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$IFIT1,
       subset(b$data, cell.type =="Activated CD4 T cells")$IFIT1,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$IFIT1,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Naive B cells")$IFIT1,
       subset(b$data, cell.type =="Naive B cells")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="NKs")$IFIT1,
       subset(b$data, cell.type =="NKs")$RSAD2,
       alternative = "two.sided", paired = F, var.equal = T)

# for CDK9 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE"#,
  #"Cylcing T cells" = "#CC6677",
  #"Cytotoxic CD8 T cells" = "#DDCC77"#,
  #"Naive B cells" = "#117733",
  #"Naive T cells" = "#332288"#,
  #"NKs" = "#AA4499"#,
  # "Plasma"= "#44AA99",
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = CDK9, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$CDK9, covid.pb.controls.sig.genes$CDK9)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = CDK9, color= cell.type))+
  geom_point(size=2)+
  ylab("CDK9 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$CDK9,covid.pb.controls.sig.genes$CDK9)),labels = scales::number_format(accuracy = 0.01))

CDK9 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                       padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$CDK9,
       subset(b$data, cell.type =="Activated CD4 T cells")$CDK9,
       alternative = "two.sided", paired = F, var.equal = T)


# for CX3CR1 ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE"#,
                 #"Cylcing T cells" = "#CC6677",
                 "Cytotoxic CD8 T cells" = "#DDCC77"#,
                 #"Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 #"NKs" = "#AA4499"#,
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = CX3CR1, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$CX3CR1, covid.pb.controls.sig.genes$CX3CR1)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = CX3CR1, color= cell.type))+
  geom_point(size=2)+
  ylab("CX3CR1 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$CX3CR1,covid.pb.controls.sig.genes$CX3CR1)),labels = scales::number_format(accuracy = 0.01))

CX3CR1 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                     padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$CX3CR1,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$CX3CR1,
       alternative = "two.sided", paired = F, var.equal = T)

# for TNFRSF13C ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE"#,
  #"Cylcing T cells" = "#CC6677",
  #"Cytotoxic CD8 T cells" = "#DDCC77"#,
  #"Naive B cells" = "#117733",
  #"Naive T cells" = "#332288"#,
  #"NKs" = "#AA4499"#,
  # "Plasma"= "#44AA99",
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = TNFRSF13C, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$TNFRSF13C, covid.pb.controls.sig.genes$TNFRSF13C)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = TNFRSF13C, color= cell.type))+
  geom_point(size=2)+
  ylab("TNFRSF13C Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$TNFRSF13C,covid.pb.controls.sig.genes$TNFRSF13C)),labels = scales::number_format(accuracy = 0.01))

TNFRSF13C <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                       padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$TNFRSF13C,
       subset(b$data, cell.type =="Activated CD4 T cells")$TNFRSF13C,
       alternative = "two.sided", paired = F, var.equal = T)


# for LGALS9 ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE"#,
                 #"Cylcing T cells" = "#CC6677",
                 #"Cytotoxic CD8 T cells" = "#DDCC77"#,
                 #"Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288"#,
                 "NKs" = "#AA4499"#,
                 #"Memory B cells" = "#6699EF",
                 # "Plasma"= "#44AA99",
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "NKs"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "NKs"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "NKs"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "NKs"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = LGALS9, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$LGALS9, covid.pb.controls.sig.genes$LGALS9)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = LGALS9, color= cell.type))+
  geom_point(size=2)+
  ylab("LGALS9 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$LGALS9,covid.pb.controls.sig.genes$LGALS9)),labels = scales::number_format(accuracy = 0.01))

LGALS9 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                          padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="NKs")$LGALS9,
       subset(b$data, cell.type =="NKs")$LGALS9,
       alternative = "two.sided", paired = F, var.equal = T)

# for IFIT3 ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE"#,
  #"Cylcing T cells" = "#CC6677",
  "Cytotoxic CD8 T cells" = "#DDCC77",
  "Naive B cells" = "#117733",
  "Naive T cells" = "#332288",
  "NKs" = "#AA4499",
  "Memory B cells" = "#6699EF"#,
  #"Plasma"= "#44AA99"#,
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "Naive T cells",
  "NKs",
  "Memory B cells"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "Naive T cells",
  "NKs",
  "Memory B cells"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "Naive T cells",
  "NKs",
  "Memory B cells"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Cytotoxic CD8 T cells",
  "Naive B cells",
  "Naive T cells",
  "NKs",
  "Memory B cells"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IFIT3, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFIT3, covid.pb.controls.sig.genes$IFIT3)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IFIT3, color= cell.type))+
  geom_point(size=2)+
  ylab("IFIT3 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFIT3,covid.pb.controls.sig.genes$IFIT3)),labels = scales::number_format(accuracy = 0.01))

IFIT3 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                       padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="NKs")$IFIT3,
       subset(b$data, cell.type =="NKs")$IFIT3,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Cytotoxic CD8 T cells")$IFIT3,
       subset(b$data, cell.type =="Cytotoxic CD8 T cells")$IFIT3,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Naive B cells")$IFIT3,
       subset(b$data, cell.type =="Naive B cells")$IFIT3,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Naive T cells")$IFIT3,
       subset(b$data, cell.type =="Naive T cells")$IFIT3,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Memory B cells")$IFIT3,
       subset(b$data, cell.type =="Memory B cells")$IFIT3,
       alternative = "two.sided", paired = F, var.equal = T)

# for IRF7 ####

# define cell type colors 
cell_colors <- c("Activated CD4 T cells" = "#88CCEE",
  #"Cylcing T cells" = "#CC6677",
  #"Cytotoxic CD8 T cells" = "#DDCC77",
  #"Naive B cells" = "#117733",
  #"Naive T cells" = "#332288",
  "NKs" = "#AA4499"#,
  #"Memory B cells" = "#6699EF"#,
  #"Plasma"= "#44AA99"#,
  # "XCL+ NKs" = "#999933",
  # "DCs" = "#882255",
  # "Megakaryocytes" = "#661100",
  # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Activated CD4 T cells" ,
  "NKs"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells" ,
  "NKs"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells" ,
  "NKs"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Activated CD4 T cells" ,
  "NKs"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IRF7, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IRF7, covid.pb.controls.sig.genes$IRF7)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IRF7, color= cell.type))+
  geom_point(size=2)+
  ylab("IRF7 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IRF7,covid.pb.controls.sig.genes$IRF7)),labels = scales::number_format(accuracy = 0.01))

IRF7 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "3" & cell.type =="NKs")$IRF7,
       subset(b$data, cell.type =="NKs")$IRF7,
       alternative = "two.sided", paired = F, var.equal = T)
t.test(subset(a$data, time == "3" & cell.type =="Activated CD4 T cells")$IRF7,
       subset(b$data, cell.type =="Activated CD4 T cells")$IRF7,
       alternative = "two.sided", paired = F, var.equal = T)

# for IFI27 ####

# define cell type colors 
cell_colors <- c(#"Activated CD4 T cells" = "#88CCEE",
                 #"Cylcing T cells" = "#CC6677",
                 #"Cytotoxic CD8 T cells" = "#DDCC77",
                 #"Naive B cells" = "#117733",
                 #"Naive T cells" = "#332288",
                 #"NKs" = "#AA4499"#,
                 #"Memory B cells" = "#6699EF"#,
                 "Plasma"= "#44AA99"#,
                 # "XCL+ NKs" = "#999933",
                 # "DCs" = "#882255",
                 # "Megakaryocytes" = "#661100",
                 # "Stem cells" = "#6699CC"
)

by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$cell.type <- exp.data$cell.type
by.type$time <- exp.data$time
by.type <- by.type[by.type$cell.type %in% c(
  "Plasma"
), ]

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[16] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[covid.pb.controls1.sig.genes$cell.type %in% c(
  
  "Plasma"
), ]

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[16] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[covid.pb.controls2.sig.genes$cell.type %in% c(
  
  "Plasma"
), ]

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,16] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[16] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[covid.pb.controls3.sig.genes$cell.type %in% c(
  
  "Plasma"
), ]

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = IFI27, 
                        color= cell.type))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFI27, covid.pb.controls.sig.genes$IFI27)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = IFI27, color= cell.type))+
  geom_point(size=2)+
  ylab("IFI27 Expression") + xlab("")+labs(fill="Cell Type")+
  scale_color_manual(values = cell_colors)+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$IFI27,covid.pb.controls.sig.genes$IFI27)),labels = scales::number_format(accuracy = 0.01))

IFI27 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                     padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "7" & cell.type =="Plasma")$IFI27,
       subset(b$data, cell.type =="Plasma")$IFI27,
       alternative = "two.sided", paired = F, var.equal = T)

# for ISG15 ####

# define cell type colors 
covid.genes <- c("ISG15", "ISG20")
by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$time <- exp.data$time

covid.pb.controls1.sig.genes <- covid.pb.controls1[rownames(covid.pb.controls1) %in% covid.genes,]
covid.pb.controls1.sig.genes <- t(covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes <- covid.pb.controls1.sig.genes[,order(colnames(covid.pb.controls1.sig.genes))]
covid.pb.controls1.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls1.sig.genes), covid.pb.controls1.sig.genes)
covid.pb.controls1.sig.genes[,3] <- rep("Baseline", times = 15)
colnames(covid.pb.controls1.sig.genes)[3] = "time"
rownames(covid.pb.controls1.sig.genes) <- paste("Ctrl", 1:15, sep = ".")

covid.pb.controls2.sig.genes <- covid.pb.controls2[rownames(covid.pb.controls2) %in% covid.genes,]
covid.pb.controls2.sig.genes <- t(covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes <- covid.pb.controls2.sig.genes[,order(colnames(covid.pb.controls2.sig.genes))]
covid.pb.controls2.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls2.sig.genes), covid.pb.controls2.sig.genes)
covid.pb.controls2.sig.genes[,3] <- rep("Baseline", times = 15)
colnames(covid.pb.controls2.sig.genes)[3] = "time"
rownames(covid.pb.controls2.sig.genes) <- paste("Ctrl", 1:15, sep = ".")

covid.pb.controls3.sig.genes <- covid.pb.controls3[rownames(covid.pb.controls3) %in% covid.genes,]
covid.pb.controls3.sig.genes <- t(covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes <- covid.pb.controls3.sig.genes[,order(colnames(covid.pb.controls3.sig.genes))]
covid.pb.controls3.sig.genes <- data.frame(cell.type = rownames(covid.pb.controls3.sig.genes), covid.pb.controls3.sig.genes)
covid.pb.controls3.sig.genes[,3] <- rep("Baseline", times = 15)
colnames(covid.pb.controls3.sig.genes)[3] = "time"
rownames(covid.pb.controls3.sig.genes) <- paste("Ctrl", 1:15, sep = ".")

covid.pb.controls.sig.genes <- rbind(covid.pb.controls1.sig.genes,
                                     covid.pb.controls2.sig.genes, 
                                     covid.pb.controls3.sig.genes)

a <-ggplot(by.type, aes(x = time, y = ISG15))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$ISG15, covid.pb.controls.sig.genes$ISG15)),labels = scales::number_format(accuracy = 0.01))

b <-ggplot(covid.pb.controls.sig.genes, aes(x = time, y = ISG15))+
  geom_point(size=2)+
  ylab("ISG15 Expression") + xlab("")+
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$ISG15,covid.pb.controls.sig.genes$ISG15)),labels = scales::number_format(accuracy = 0.01))

ISG15 <- grid.arrange(b, a, widths = c(0.4, 2), nrow= 1, ncol = 2, 
                      padding = unit(0, "lines"))

# t test
t.test(subset(a$data, time == "7")$ISG15,
       subset(b$data)$ISG15,
       alternative = "two.sided", paired = F, var.equal = T)

#standalone legend
# Define the colors and labels for the legend
cell_colors <- c(
  "Activated CD4 T cells" = "#88CCEE",
  # "Cycling T cells" = "#CC6677",
  "Cytotoxic CD8 T cells" = "#DDCC77",
  "Naive B cells" = "#117733",
  "Naive T cells" = "#332288",
  "NKs" = "#AA4499",
  "Plasma"= "#44AA99",
  # "DCs" = "#882255",
  #"Megakaryocytes" = "#661100",
  #"Stem cells" = "#6699CC",
  "Memory B cells" = "#6699EF"
)

# Create an empty plot with a dummy layer
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))

# Add the legend
legend("center", legend = names(cell_colors), fill = cell_colors, title = "Cell Types",
       ncol = 1, box.lwd = 0)

#plot by disease trend 
covid.genes <- c("EIF3H", "OASL", "RSAD2", "ISG20", "IFI16", "IFIT1",
                 "CDK9", "CX3CR1", "TNFRSF13C", "LGALS9", "IFIT3", "IRF7", "XAF1", "IFI27", "ISG15")


# for ISG15 ####
by.type <-as.data.frame(t(covid.masigpro[covid.genes,]))
by.type <- by.type[,order(colnames(by.type))]
by.type$time <- exp.data$time

a <-ggplot(by.type, aes(x = time, y = ISG15))+
  geom_point(size=2)+
  ylab("") + xlab("Days Since Symptom Onset")+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se= F)+
  scale_x_continuous(n.breaks = 10)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,-1.75), "mm"))+
  scale_y_continuous(limits = c(-.1,max(by.type$ISG15, covid.pb.controls.sig.genes$ISG15)),labels = scales::number_format(accuracy = 0.01))

#correlate disease severity with sig genes ####
covid.subset.time <- readRDS("~/Dropbox (EinsteinMed)/Result_from_KevinO/metacell covid paper/code/r.objects/covid.subset.time.RDS")

time_points <- c("3", "7", "9", "10", "13", "15", "16", "22", "25", "28")

setwd("~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test")
# assign all cells a seacell assignment

#find which time points correspond to each sample
covid.subset.time[["sample"]] <- substr(covid.subset.time@meta.data$batch, 1, 5)
table(covid.subset.time@meta.data$sample, covid.subset.time@meta.data$time)

#subset seacell matrix by genes
genes <- c("EIF3H", "OASL", "RSAD2", "ISG20", "IFI16", "IFIT1",
           "CDK9", "CX3CR1", "TNFRSF13C", "LGALS9", "IFIT3", "IRF7", "XAF1", "IFI27", "ISG15")

covid.masigpro.t <- covid.masigpro[rownames(covid.masigpro) %in% genes,]
covid.masigpro.t <- as.data.frame(t(covid.masigpro.t))

#populate covid.masigpro.t with relevent metadata
exp.data <- read.csv("all.cells.masigpro.design.2.csv", row.names = 1)
cell.types <- unique(factor(exp.data$cell.type))

#severe patient is the only sample comprising points 3, 9 and 15
exp.data <- exp.data %>%
  mutate(severity = ifelse(time %in% c(3, 9, 15), "severe", "mild"))

covid.masigpro.t <- covid.masigpro.t[order(rownames(covid.masigpro.t)), , drop = FALSE]
exp.data <- exp.data[order(rownames(exp.data)), , drop = FALSE]
covid.masigpro.t$time <- as.numeric(exp.data$time)
covid.masigpro.t$disease.state <- exp.data$severity

plot_list <- list()
for (gene in genes) {
  current_plot <- ggplot(covid.masigpro.t, aes(x = time, y = .data[[gene]], color = disease.state, group = disease.state)) +
    geom_point() +
    NoLegend() +
    ylab(paste(gene, "Expression", sep = " ")) + xlab("Days Since Symptom Onset")+
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, size = 1)+
    scale_y_continuous(limits = c(-.1,max(covid.masigpro.t[[gene]])),labels = scales::number_format(accuracy = 0.01))
  plot_list[[gene]] <- current_plot
}

plot_list[1:15]


#legend
group <- c(
"mild"= "#FF5054",
"severe" = "#00B6B9")

# Create an empty plot with a dummy layer
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))

# Add the legend
legend("center", legend = names(group), fill = group, title = "Disease Group",
       ncol = 2, box.lwd = 0)

# cell type proportion through time ####
cell_colors <- c(
  "Activated CD4 T cells" = "#88CCEE",
  # "Cycling T cells" = "#CC6677",
  "Cytotoxic CD8 T cells" = "#DDCC77",
  "Naive B cells" = "#117733",
  "Naive T cells" = "#332288",
  "NKs" = "#AA4499",
  "Plasma"= "#44AA99",
  # "DCs" = "#882255",
  #"Megakaryocytes" = "#661100",
  #"Stem cells" = "#6699CC",
  "Memory B cells" = "#6699EF"
)

exp.data <- read.csv("all.cells.masigpro.design.2.csv")

exp.data <- exp.data[exp.data$cell.type %in% c("Activated CD4 T cells", "Cytotoxic CD8 T cells", "Naive B cells",
                                                    "Naive T cells" ,"NKs","Plasma","Memory B cells"), ]

props <- data.frame(prop.table(table(exp.data$cell.type, exp.data$time), margin = 2))
colnames(props) <- c("Cell.Type", "Time", "Proportion")


ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cell_colors)

ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
  geom_smooth(method = "loess", se = F)+
  geom_point()+
  scale_color_manual(values = cell_colors)

ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
  geom_smooth(method = "lm",formula = y ~ poly(x, 2), se = F)+
  geom_point()+
  scale_color_manual(values = cell_colors)

ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
  geom_smooth(method = "lm",formula = y ~ x, se = F)+
  geom_point()+
  scale_color_manual(values = cell_colors)

grid.arrange(ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
               geom_line()+
               geom_point()+
               scale_color_manual(values = cell_colors),
             ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
               geom_smooth(method = "loess", se = F)+
               geom_point()+
               scale_color_manual(values = cell_colors),
             ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
               geom_smooth(method = "lm",formula = y ~ poly(x, 2), se = F)+
               geom_point()+
               scale_color_manual(values = cell_colors),
             ggplot(data = props, aes(x = Time, group = Cell.Type, y = Proportion, color = Cell.Type))+
               geom_smooth(method = "lm",formula = y ~ x, se = F)+
               geom_point()+
               scale_color_manual(values = cell_colors))
