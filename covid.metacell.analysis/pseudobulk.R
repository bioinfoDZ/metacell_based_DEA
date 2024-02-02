#+eval=FALSE

library(Seurat)
library(SeuratData)
library(SeuratDisk)

#####upload data and edit####
setwd("~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test")
#dataset obtained from https://ftp.cngb.org/pub/CNSA/data3/CNP0001102/Single_Cell/CSE0000013/
covid <- readRDS(file = "~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test/covid.RDS")
View(covid@meta.data)
# all cell types, only covid and flu 

covid.subset <- subset(covid, subset = Stage != "Ctrl")
covid.subset <- subset(covid.subset, subset = Stage != "IAV-Stage1")
covid.subset <- subset(covid.subset, subset = Stage != "IAV-Stage2")
covid.subset <- subset(covid.subset, subset = Stage != "IAV-Stage3")

#replace d# w/ days since sx onset 
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
View(covid.subset@meta.data)
covid.subset.time <- covid.subset

#default counts, get rid of integrated idk what it is 
DefaultAssay(covid.subset.time) <- "RNA"
covid.subset.time[["integrated"]] <- NULL

##### for all cells #####
#pseudobulk each time point 
covid.pb <- AverageExpression(covid.subset.time, group.by = "batch",
                              return.seurat = T, slot = "data")

#as.matrix of "data", which is normalized data
covid.masigpro <- as.matrix(covid.pb@assays$RNA$data)
View(covid.masigpro)
#get colnames to create masigpro design matrix in excel
colnames(covid.masigpro)

library(maSigPro)

#load exp design data 
exp.data <- read.csv("covid.masigpro.design.pb.csv")
rownames(exp.data) <- exp.data[,1]
exp.data <- exp.data[,-1]
exp.data <- exp.data[,-3]
exp.data <- exp.data[,-3]
View(exp.data)

#define regression model
design <- make.design.matrix(exp.data, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3) #quadratic in this case to observe fluctuation through 3 time points
View(design)

#masigpro for rna seq 
library(MASS)
NBp <- p.vector(covid.masigpro, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

see.genes(get$sig.genes, k=9)

##### now for cell types  #####
#load exp design data 
exp.data <- read.csv("covid.masigpro.design.pb.csv")
rownames(exp.data) <- exp.data[,1]
exp.data <- exp.data[,-1]
exp.data <- exp.data[,-3]
exp.data <- exp.data[,-3]
View(exp.data)

design <- make.design.matrix(exp.data, degree = 2, time.col = 1,
                             repl.col = 2) 

cell.types <- levels(covid.subset.time@meta.data$cell_type)


covid.subset.list <- list()

for (type in cell.types) {
  covid.subset.type <- subset(covid.subset.time, subset = cell_type == type)
  covid.subset.type.pb <- AverageExpression(covid.subset.type, group.by = "batch",
                                            return.seurat = FALSE, slot = "data")
  
  covid.subset.list[[type]] <- covid.subset.type.pb
}


#find sig genes for each cell type 
#1
NBp <- p.vector(covid.subset.list$`Cytotoxic CD8 T cells`$RNA, design, counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get.cyto.t.cells <-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.cyto.t.cells$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.cyto.t.cells$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-2], function(col) {
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


#2
NBp <- p.vector(covid.subset.list$`Naive T cells`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.naive.t.cells<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.naive.t.cells$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.naive.t.cells$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-2], function(col) {
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

#3
NBp <- p.vector(covid.subset.list$NKs$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.NKs<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.NKs$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.NKs$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-2], function(col) {
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

#4
NBp <- p.vector(covid.subset.list$MAIT$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.MAIT<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary
genes <- get.MAIT$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.MAIT$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-2], function(col) {
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

#5
NBp <- p.vector(covid.subset.list$`Activated CD4 T cells`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.activated.t.cell<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.activated.t.cell$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.activated.t.cell$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-2], function(col) {
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

#6
NBp <- p.vector(covid.subset.list$`Naive B cells`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.naive.b<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.naive.b$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.naive.b$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-3], function(col) {
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

#7
NBp <- p.vector(covid.subset.list$Plasma$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.plasma<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.plasma$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.plasma$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-4], function(col) {
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

#8

NBp <- p.vector(covid.subset.list$`Memory B cells`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.mem.b<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.mem.b$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.mem.b$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-4], function(col) {
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

#9
NBp <- p.vector(covid.subset.list$`XCL+ NKs`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.xcl.nks<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.xcl.nks$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.xcl.nks$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-3], function(col) {
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

#10
NBp <- p.vector(covid.subset.list$`Cycling T cells`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.cycling.t<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.cycling.t$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.cycling.t$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-9], function(col) {
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

#11
exp.data2 <- exp.data[-3,]
exp.data2 <- exp.data.2[-4,]
design2<- make.design.matrix(exp.data2,degree = 2, time.col = 1,
                             repl.col = 2) 

NBp <- p.vector(covid.subset.list$`Monocytes`$RNA, design2,counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get.monocytes<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

genes <- get.monocytes$summary

#get average expression of sig genes for all metacells that of this cell type
by.type <-as.data.frame(t(get.monocytes$sig.genes$sig.profiles[genes,]))
View(by.type)
#add time column that has metacells labeled by time for naive t cells
by.type$time <- exp.data$Time
#find average expression of metacells by time point
averages <- by.type %>% group_by(time) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = T))
plots <- lapply(names(by.type)[-9], function(col) {
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

#some samples not present in this matrix

#12
NBp <- p.vector(covid.subset.list$DCs$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

#some samples not present in this matrix

#13
NBp <- p.vector(covid.subset.list$`Cycling Plasma`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

see.genes(get$sig.genes, k =8)

#14
NBp <- p.vector(covid.subset.list$`Stem cells`$RNA, design, counts=FALSE, 
                min.obs = 6) 
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

#15
NBp <- p.vector(covid.subset.list$Megakaryocytes$RNA, design, counts=FALSE, 
                min.obs = 6) NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary























