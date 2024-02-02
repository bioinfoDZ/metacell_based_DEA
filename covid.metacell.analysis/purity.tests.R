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

# Merge the data matrices
merged.data <- cbind(covid.pb.d3, covid.pb.d7, covid.pb.d9, covid.pb.d10,
                     covid.pb.d13, covid.pb.d15, covid.pb.d16, covid.pb.d22,
                     covid.pb.d25, covid.pb.d28)
merged.data
covid.masigpro <- merged.data
rm(merged.data)

#prep masigpro by cell type, mixed
cell.types <- c("Cytotoxic CD8 T cells", "Naive T cells", "NKs",
                "Activated CD4 T cells", "Naive B cells", "Plasma",
                "Memory B cells", "XCL+ NKs", "MAIT", "Cycling T cells")
# Create  empty lists to store data frames
exp.data.list <- list()
covid.masigpro.list <- list()
design.list <- list()

# Loop through cell types
for (type in cell.types) {
  # Read the CSV file
  exp.data <- read.csv("all.cells.masigpro.design.2.purity.csv")
  
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
  exp.data.type <- exp.data.type[, -c(1, 4)]
  
  # Subset covid.masigpro based on column names
  covid.masigpro.type <- covid.masigpro[, rownames(exp.data.type)]
  
  # Store the data frames in the lists
  covid.masigpro.list[[type]] <- covid.masigpro.type
  
  exp.data.list[[type]] <- exp.data.type
  
}

for (type in cell.types) {
  #define regression model
  exp.data.type <- exp.data.list[[type]]
  design <- make.design.matrix(exp.data.type, degree = 2, time.col = 1,
                               repl.col = 2, group.cols = 3)
  design.list[[type]] <- design
}

get.list.mixed <- list()
for (type in cell.types) {
NBp <- p.vector(covid.masigpro.list[[type]], design.list[[type]], counts=FALSE, 
                min.obs = 6)
NBt <- T.fit(NBp)
get <- get.siggenes(NBt, vars="all", rsq = 0.5)
get.list.mixed[[type]] <- get$summary
}

num.metacells.mixed <- list()
for (type in cell.types) {
  num.metacells.mixed[[type]] <- ncol(covid.masigpro.list[[type]])
} 

#prep masigpro by cell type, pure
cell.types <- c("Cytotoxic CD8 T cells", "Naive T cells", "NKs",
                "Activated CD4 T cells", "Naive B cells", "Plasma",
                "Memory B cells", "XCL+ NKs", "MAIT", "Cycling T cells")
# Create  empty lists to store data frames
exp.data.list <- list()
covid.masigpro.list <- list()
design.list <- list()

# Loop through cell types
for (type in cell.types) {
  # Read the CSV file
  exp.data <- read.csv("all.cells.masigpro.design.2.purity.csv")
  
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
  exp.data.type <- exp.data.type[exp.data.type$purity>0.9,] 
  exp.data.type <- exp.data.type[, -c(1, 4)]
  
  # Subset covid.masigpro based on column names
  covid.masigpro.type <- covid.masigpro[, rownames(exp.data.type)]
  
  # Store the data frames in the lists
  covid.masigpro.list[[type]] <- covid.masigpro.type
  
  exp.data.list[[type]] <- exp.data.type
  
}

for (type in cell.types) {
  #define regression model
  exp.data.type <- exp.data.list[[type]]
  design <- make.design.matrix(exp.data.type, degree = 2, time.col = 1,
                               repl.col = 2, group.cols = 3)
  design.list[[type]] <- design
}

get.list.pure <- list()
for (type in cell.types) {
  NBp <- p.vector(covid.masigpro.list[[type]], design.list[[type]], counts=FALSE, 
                  min.obs = 6)
  NBt <- T.fit(NBp)
  get <- get.siggenes(NBt, vars="all", rsq = 0.5)
  get.list.pure[[type]] <- get$summary
} 

#find number of metacells per cell type
num.metacells.pure <- list()
for (type in cell.types) {
  num.metacells.pure[[type]] <- ncol(covid.masigpro.list[[type]])
} 


#make venn diagrams
library(ggVennDiagram)
venn.diagrams <- list()
for (type in cell.types) {
  x <- list(get.list.pure[[type]], get.list.mixed[[type]])

  plot <- ggVennDiagram(x,
                        label_size = 3,
                        set_size = 2,
                        category.names = c("pure", "mixed"))+
    ggtitle(paste(type))+
    theme(legend.position = "none")
  
  venn.diagrams[[type]] <- plot
  
  
}

x <- list(unique(unlist(get.list.pure)), unique(unlist(get.list.mixed)))

plot <- ggVennDiagram(x,
                      label_size = 3,
                      set_size = 2,
                      category.names = c("pure", "mixed"))+
  ggtitle("Total Unique DEGs")+
  theme(legend.position = "none")

venn.diagrams[[11]] <- plot


library(cowplot)
pure <- list()
mixed <- list()
for (type in cell.types){
  pure[[type]] <- length(get.list.pure[[type]])
  mixed[[type]] <- length(get.list.mixed[[type]])
}

t_test_result <- t.test(unlist(pure), unlist(mixed), alternative = "two.sided", paired = FALSE, var.equal = TRUE)

# Create bar plots
bar_data <- data.frame(
  cell_type = rep(cell.types, times = 2),
  value = c(unlist(pure), unlist(mixed)),
  DEG_type = rep(c("Pure", "Mixed"), each = length(cell.types))
)

add.rows <- data.frame(cell_type=c("Total Unique DEGs", 'Total Unique DEGs'),
                       value=c(length(unique(unlist(get.list.mixed))), length(unique(unlist(get.list.pure)))),
                       DEG_type=c('Mixed', 'Pure'),
                       stringsAsFactors=FALSE)

bar_data <- rbind(bar_data, add.rows)
bar_data$cell_type <- factor(bar_data$cell_type, levels = c(cell.types, "Total Unique DEGs"))

bar1 <- ggplot(bar_data, aes(x = cell_type, y = value, fill = DEG_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs( y = "# DEGs") +
  scale_fill_manual(values = c("Pure" = "blue", "Mixed" = "red")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none")

#number of metacells per type 
num.metacells <- data.frame(
  cell_type = rep(cell.types, times = 2),
  value = c(unlist(num.metacells.pure), unlist(num.metacells.mixed)),
  DEG_type = rep(c("Pure", "Mixed"), each = length(cell.types))
)

add.rows.num <- data.frame(cell_type=c("Total # Metacells", 'Total # Metacells'),
                       value=c(sum(unlist(num.metacells.pure)), sum(unlist(num.metacells.mixed))),
                       DEG_type=c('Pure', 'Mixed'),
                       stringsAsFactors=FALSE)

bar_data.num <- rbind(num.metacells, add.rows.num)

bar_data.num$cell_type <- factor(bar_data.num$cell_type, levels = c(cell.types, "Total # Metacells"))

bar2 <- ggplot(bar_data.num, aes(x = cell_type, y = value, fill = DEG_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs( y = "# Metacells") +
  scale_fill_manual(values = c("Pure" = "blue", "Mixed" = "red")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())

bar <- bar1+bar2

# Combine bar plots and venn diagrams
grid <- plot_grid(
  venn.diagrams[[1]],
  venn.diagrams[[2]],
  venn.diagrams[[3]],
  venn.diagrams[[4]],
  venn.diagrams[[5]],
  venn.diagrams[[6]],
  venn.diagrams[[7]],
  venn.diagrams[[8]],
  venn.diagrams[[9]],
  venn.diagrams[[10]],
  venn.diagrams[[11]],
  bar, ncol = 3
)

print(grid)

#correlation plot
colnames(bar_data) <- c("cell_type", "num.DEGs", "DEG_type")
colnames(bar_data.num) <- c("cell_type", "num.Metacells", "DEG_Type")
combined.bar.data <- cbind(bar_data, bar_data.num)
combined.bar.data <- combined.bar.data[,-4]
combined.bar.data <- combined.bar.data[,-5]
combined.bar.data <- combined.bar.data[-c(21:22),]

ggplot(combined.bar.data, aes(x = num.Metacells, y = num.DEGs, color = DEG_type))+
  geom_point()
