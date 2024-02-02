library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(maSigPro)
library(MASS)
library(dplyr)

#load all cells
covid.subset.time <- readRDS("~/Dropbox (EinsteinMed)/Result_from_KevinO/metacell covid paper/code/r.objects/covid.subset.time.RDS")
covid.subset.time@meta.data[["repl"]] <- covid.subset.time@meta.data$time
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==3] <- 1
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==7] <- 2
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==9] <- 3
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==10] <- 4
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==13] <- 5
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==15] <- 6
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==16] <- 7
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==22] <- 8
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==25] <- 9
covid.subset.time@meta.data$repl[covid.subset.time@meta.data$time==28] <- 10
covid.subset.time@meta.data[,10] <- as.numeric(covid.subset.time@meta.data[,10])
covid.subset.time@meta.data[,11] <- as.numeric(covid.subset.time@meta.data[,11])
meta <- covid.subset.time@meta.data[,c(10,11)]
meta[["group"]] <- "1"
design <- make.design.matrix(meta, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

#masigpro all cells
covid.subset.time@assays$RNA$data <- as.matrix(covid.subset.time@assays$RNA$data)
new.matrix <- covid.subset.time@assays$RNA$data[rowSums(covid.subset.time@assays$RNA$data) != 0,]
NBp <- p.vector(new.matrix, design, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

see.genes(get$sig.genes, k=8)

# masigpro by cell type
cell.types <- levels(covid.subset.time@meta.data$cell_type)

get.list <- list()
for ( type in cell.types ){

design <- meta[rownames(meta) %in% rownames(covid.subset.time@meta.data[covid.subset.time@meta.data$cell_type==type,]),] 

matrix <- new.matrix[,colnames(new.matrix) %in% rownames(covid.subset.time@meta.data[covid.subset.time@meta.data$cell_type==type,])]

design <- make.design.matrix(design, degree = 2, time.col = 1,
                                       repl.col = 2, group.cols = 3)

NBp <- p.vector(matrix, design, counts=FALSE, 
                min.obs = 6)

NBt <- T.fit(NBp)

current_result <- get.siggenes(NBt, vars = "all", rsq = 0.5)

get.list[[type]] <- current_result

}

get.list

