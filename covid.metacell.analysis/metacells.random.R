#####upload data and edit####
setwd("~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test")
#dataset obtained from https://ftp.cngb.org/pub/CNSA/data3/CNP0001102/Single_Cell/CSE0000013/
covid <- readRDS(file = "~/Documents/Einstein/Zheng Lab/Heart Dev/timecourse.test/covid.RDS")
View(covid@meta.data)

# all cell types, only covid and flu 
covid.subset.controls <- subset(covid, subset = Stage != "Ctrl")
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
#View(covid.subset@meta.data)
covid.subset.time <- covid.subset

covid.subset.time[["time"]] <- covid.subset.time$batch
covid.subset.time@meta.data$time <- sub(".*-D", "", covid.subset.time@meta.data$time)
View(covid.subset.time@meta.data)


cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes") 

# random metacells per time point below
#####d3####
covid.subset.time.d3 <- subset(covid.subset.time, subset = time == "3")
DefaultAssay(covid.subset.time.d3) <- "RNA"
covid.subset.time.d3[["integrated"]] <- NULL
summary(covid.subset.time.d3$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 20 cells 
result_list <- list()
for (n in 1:3){
random_columns <- sample(1:ncol(covid.subset.time.d3@assays$RNA@data), 20)
selected_columns <- covid.subset.time.d3@assays$RNA@data[, random_columns]
row_averages <- rowMeans(selected_columns)
new_df <- data.frame(mc = row_averages)
colnames(new_df) <- paste("mc.d3.", n, sep = "")
result_list[[n]] <- new_df
}
d3.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d3.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d3.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d3, subset = cell_type == i))
  covid.subset.time.d3.list[[i]] <- get(subset_name)
}

covid.subset.time.d3.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d3.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d3.new.", i, sep = "")
    covid.subset.time.d3.new.list[[i]] <- covid.subset.time.d3.list[[i]]
  }
}

for (i in names(covid.subset.time.d3.new.list)) {
  d3.subset <- covid.subset.time.d3.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d3.subset@assays$RNA@data), 20)
    selected_columns <- d3.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d3.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d3.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(d3.metacells)
View(`d3.Cytotoxic CD8 T cells`)

#####d7####
covid.subset.time.d7 <- subset(covid.subset.time, subset = time == "7")
DefaultAssay(covid.subset.time.d7) <- "RNA"
covid.subset.time.d7[["integrated"]] <- NULL
summary(covid.subset.time.d7$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d7@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d7@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d7.", n, sep = "")
  result_list[[n]] <- new_df
}
d7.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d7.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d7.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d7, subset = cell_type == i))
  covid.subset.time.d7.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d7.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d7.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d7.new.", i, sep = "")
    covid.subset.time.d7.new.list[[i]] <- covid.subset.time.d7.list[[i]]
  }
}

for (i in names(covid.subset.time.d7.new.list)) {
  d7.subset <- covid.subset.time.d7.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d7.subset@assays$RNA@data), 20)
    selected_columns <- d7.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d7.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d7.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d7.metacells))
View(head(`d7.Cytotoxic CD8 T cells`))

#####d9####
covid.subset.time.d9 <- subset(covid.subset.time, subset = time == "9")
DefaultAssay(covid.subset.time.d9) <- "RNA"
covid.subset.time.d9[["integrated"]] <- NULL
summary(covid.subset.time.d9$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d9@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d9@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d9.", n, sep = "")
  result_list[[n]] <- new_df
}
d9.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d9.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d9.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d9, subset = cell_type == i))
  covid.subset.time.d9.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d9.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d9.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d9.new.", i, sep = "")
    covid.subset.time.d9.new.list[[i]] <- covid.subset.time.d9.list[[i]]
  }
}

for (i in names(covid.subset.time.d9.new.list)) {
  d9.subset <- covid.subset.time.d9.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d9.subset@assays$RNA@data), 20)
    selected_columns <- d9.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d9.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d9.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d9.metacells))
View(head(`d9.XCL+ NKs`))


#####d10####
covid.subset.time.d10 <- subset(covid.subset.time, subset = time == "10")
DefaultAssay(covid.subset.time.d10) <- "RNA"
covid.subset.time.d10[["integrated"]] <- NULL
summary(covid.subset.time.d10$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d10@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d10@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d10.", n, sep = "")
  result_list[[n]] <- new_df
}
d10.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d10.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d10.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d10, subset = cell_type == i))
  covid.subset.time.d10.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d10.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d10.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d10.new.", i, sep = "")
    covid.subset.time.d10.new.list[[i]] <- covid.subset.time.d10.list[[i]]
  }
}

for (i in names(covid.subset.time.d10.new.list)) {
  d10.subset <- covid.subset.time.d10.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d10.subset@assays$RNA@data), 20)
    selected_columns <- d10.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d10.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d10.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d10.metacells))
View(head(`d10.Cytotoxic CD8 T cells`))


#####d13####
covid.subset.time.d13 <- subset(covid.subset.time, subset = time == "13")
DefaultAssay(covid.subset.time.d13) <- "RNA"
covid.subset.time.d13[["integrated"]] <- NULL
summary(covid.subset.time.d13$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d13@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d13@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d13.", n, sep = "")
  result_list[[n]] <- new_df
}
d13.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d13.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d13.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d13, subset = cell_type == i))
  covid.subset.time.d13.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d13.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d13.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d13.new.", i, sep = "")
    covid.subset.time.d13.new.list[[i]] <- covid.subset.time.d13.list[[i]]
  }
}

for (i in names(covid.subset.time.d13.new.list)) {
  d13.subset <- covid.subset.time.d13.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d13.subset@assays$RNA@data), 20)
    selected_columns <- d13.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d13.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d13.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d13.metacells))
View(head(`d13.Cytotoxic CD8 T cells`))


#####d16####
covid.subset.time.d16 <- subset(covid.subset.time, subset = time == "16")
DefaultAssay(covid.subset.time.d16) <- "RNA"
covid.subset.time.d16[["integrated"]] <- NULL
summary(covid.subset.time.d16$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d16@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d16@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d16.", n, sep = "")
  result_list[[n]] <- new_df
}
d16.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d16.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d16.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d16, subset = cell_type == i))
  covid.subset.time.d16.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d16.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d16.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d16.new.", i, sep = "")
    covid.subset.time.d16.new.list[[i]] <- covid.subset.time.d16.list[[i]]
  }
}

for (i in names(covid.subset.time.d16.new.list)) {
  d16.subset <- covid.subset.time.d16.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d16.subset@assays$RNA@data), 20)
    selected_columns <- d16.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d16.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d16.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d16.metacells))
View(head(`d16.Cytotoxic CD8 T cells`))


#####d15####
covid.subset.time.d15 <- subset(covid.subset.time, subset = time == "15")
DefaultAssay(covid.subset.time.d15) <- "RNA"
covid.subset.time.d15[["integrated"]] <- NULL
summary(covid.subset.time.d15$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d15@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d15@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d15.", n, sep = "")
  result_list[[n]] <- new_df
}
d15.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d15.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d15.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d15, subset = cell_type == i))
  covid.subset.time.d15.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d15.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d15.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d15.new.", i, sep = "")
    covid.subset.time.d15.new.list[[i]] <- covid.subset.time.d15.list[[i]]
  }
}

for (i in names(covid.subset.time.d15.new.list)) {
  d15.subset <- covid.subset.time.d15.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d15.subset@assays$RNA@data), 20)
    selected_columns <- d15.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d15.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d15.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d15.metacells))
View(head(`d15.Cytotoxic CD8 T cells`))



#####d22####
covid.subset.time.d22 <- subset(covid.subset.time, subset = time == "22")
DefaultAssay(covid.subset.time.d22) <- "RNA"
covid.subset.time.d22[["integrated"]] <- NULL
summary(covid.subset.time.d22$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d22@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d22@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d22.", n, sep = "")
  result_list[[n]] <- new_df
}
d22.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d22.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d22.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d22, subset = cell_type == i))
  covid.subset.time.d22.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d22.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d22.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d22.new.", i, sep = "")
    covid.subset.time.d22.new.list[[i]] <- covid.subset.time.d22.list[[i]]
  }
}

for (i in names(covid.subset.time.d22.new.list)) {
  d22.subset <- covid.subset.time.d22.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d22.subset@assays$RNA@data), 20)
    selected_columns <- d22.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d22.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d22.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d22.metacells))
View(head(`d22.Cytotoxic CD8 T cells`))


#####d25####
covid.subset.time.d25 <- subset(covid.subset.time, subset = time == "25")
DefaultAssay(covid.subset.time.d25) <- "RNA"
covid.subset.time.d25[["integrated"]] <- NULL
summary(covid.subset.time.d25$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d25@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d25@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d25.", n, sep = "")
  result_list[[n]] <- new_df
}
d25.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d25.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d25.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d25, subset = cell_type == i))
  covid.subset.time.d25.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d25.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d25.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d25.new.", i, sep = "")
    covid.subset.time.d25.new.list[[i]] <- covid.subset.time.d25.list[[i]]
  }
}

for (i in names(covid.subset.time.d25.new.list)) {
  d25.subset <- covid.subset.time.d25.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d25.subset@assays$RNA@data), 20)
    selected_columns <- d25.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d25.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d25.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d25.metacells))
View(head(`d25.Cytotoxic CD8 T cells`))


#####d28####
covid.subset.time.d28 <- subset(covid.subset.time, subset = time == "28")
DefaultAssay(covid.subset.time.d28) <- "RNA"
covid.subset.time.d28[["integrated"]] <- NULL
summary(covid.subset.time.d28$cell_type)
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells",
                "Monocytes","DCs","Cycling Plasma", "Stem cells","Megakaryocytes")
# random subset the average of 75 cells 
result_list <- list()
for (n in 1:3){
  random_columns <- sample(1:ncol(covid.subset.time.d28@assays$RNA@data), 20)
  selected_columns <- covid.subset.time.d28@assays$RNA@data[, random_columns]
  row_averages <- rowMeans(selected_columns)
  new_df <- data.frame(mc = row_averages)
  colnames(new_df) <- paste("mc.d28.", n, sep = "")
  result_list[[n]] <- new_df
}
d28.metacells <- do.call(cbind, result_list)

#subset by cell types and do same thing
covid.subset.time.d28.list <- list()
for (i in cell.types) {
  subset_name <- paste("covid.subset.time.d28.", i, sep = "")
  assign(subset_name, subset(covid.subset.time.d28, subset = cell_type == i))
  covid.subset.time.d28.list[[i]] <- get(subset_name)
}
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

covid.subset.time.d28.new.list <- list()
for (i in cell.types) {
  if (ncol(covid.subset.time.d28.list[[i]]@assays$RNA@data) >= 20) {
    #new_df_name <- i #paste("covid.subset.time.d28.new.", i, sep = "")
    covid.subset.time.d28.new.list[[i]] <- covid.subset.time.d28.list[[i]]
  }
}

for (i in names(covid.subset.time.d28.new.list)) {
  d28.subset <- covid.subset.time.d28.new.list[[i]]
  result_list <- list()
  
  for (n in 1:3) {
    random_columns <- sample(1:ncol(d28.subset@assays$RNA@data), 20)
    selected_columns <- d28.subset@assays$RNA@data[, random_columns]
    row_averages <- rowMeans(selected_columns)
    new_df <- data.frame(mc = row_averages)
    colnames(new_df) <- paste("mc.d28.", n,sep = "")
    result_list[[n]] <- new_df
  }
  
  # Dynamically assign the result to a variable name
  variable_name <- paste("d28.",i, sep = "")
  assign(variable_name, do.call(cbind, result_list))
}

View(head(d28.metacells))
View(head(`d28.Cytotoxic CD8 T cells`))


#created design matrix for each cell type in excel, ie combined.metacell.design.xlsx 
#for each cell type, metacells were output to combined.metacells.csv
#this was purely so i could copy the metacells names into the actual design matrix I was creating for each cell type
#no longer necessary since the design matrix now has all the metacell names and times etc
####merge metacells for all cell types then masigpro  #####
# Get the names of data frames ending with ".metacells"
# metacells_names <- grep("\\.metacells$", ls(), value = TRUE)
# 
# # Subset the data frames
# metacells_list <- lapply(metacells_names, get)
# 
# # Combine the data frames using cbind
# combined_metacells <- do.call(cbind, metacells_list)
# mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 1)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

#remove unnecessary variables
rm(list = ls(pattern = "^covid.subset.time.d"))
rm(list = ls(pattern = "^d3.metacells"))

####merge metacells for cd8 t cells then masigpro  #####

metacells_names <- grep("\\.Cytotoxic CD8 T cells$", ls(), value = TRUE)

# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 2)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.cd8t.rando <-get.siggenes(NBt, vars="all", rsq = 0.5)
get.cd8t.rando$summary

####merge metacells for Naive T cells then masigpro  #####

metacells_names <- grep("\\.Naive T cells$", ls(), value = TRUE)

# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 3)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.naive.t.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.naive.t.rando$summary


####merge metacells for NKs then masigpro  #####

metacells_names <- grep("\\.NKs$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 4)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.nks.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.nks.rando$summary


####merge metacells for cd4 t cells then masigpro  #####

metacells_names <- grep("\\.Activated CD4 T cells$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 5)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.cd4.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.cd4.rando$summary

####merge metacells for naive b then masigpro  #####

metacells_names <- grep("\\.Naive B cells$", ls(), value = TRUE)

# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 6)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.naive.b.rando <-get.siggenes(NBt, vars="all", rsq = 0.5)
get.naive.b.rando$summary

####merge metacells for plasma cells then masigpro  #####

metacells_names <- grep("\\.Plasma$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 7)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.plasma.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.plasma.rando$summary

####merge metacells for memory B then masigpro  #####

metacells_names <- grep("\\.Memory B cells$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 8)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)


NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.mem.b.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.mem.b.rando$summary

####merge metacells for XCL NKs then masigpro  #####

metacells_names <- grep("\\XCL\\+ NKs$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 9)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

combined_metacells <- combined_metacells[,-c(32:33)]
NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.xcl.nks.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.xcl.nks.rando$summary
see.genes <- see.genes(get$sig.genes, k= 5)

####merge metacells for MAIT then masigpro  #####

metacells_names <- grep("\\MAIT$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 10)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

combined_metacells <- combined_metacells[,-41]
NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.mait.rando<-get.siggenes(NBt, vars="all", rsq = 0.5)
get.mait.rando$summary
     
####merge metacells for cycling t  then masigpro  #####

metacells_names <- grep("\\Cycling T cells$", ls(), value = TRUE)
metacells_names <- metacells_names[2:9]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 11)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

combined_metacells <- combined_metacells[,-41]
NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get.cycling.t.rando <-get.siggenes(NBt, vars="all", rsq = 0.5)
get.cycling.t.rando$summary

####merge metacells for monocytes then masigpro  #####

metacells_names <- grep("\\Monocytes$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 12)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

combined_metacells <- combined_metacells[,-41]
NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

####merge metacells for stem cells then masigpro  #####

metacells_names <- grep("\\Stem Cells$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 12)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

combined_metacells <- combined_metacells[,-41]
NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary

####merge metacells for dcs plasma then masigpro  #####

metacells_names <- grep("\\DCs$", ls(), value = TRUE)
metacells_names <- metacells_names[1:10]
# Retrieve the data frames
metacells_list <- mget(metacells_names)

# Combine the data frames using cbind
combined_metacells <- do.call(cbind, metacells_list)
mean(unlist(combined_metacells))
# write.csv(combined_metacells, "combined.metacells.csv")

#masigpro

#define regression model
combined.metacells.design <- read.xlsx("combined.metacells.design.xlsx", sheetIndex = 12)
rownames(combined.metacells.design) <- combined.metacells.design$NA.
combined.metacells.design <- combined.metacells.design[,-1]
design <- make.design.matrix(combined.metacells.design, degree = 2, time.col = 1,
                             repl.col = 2, group.cols = 3)

combined_metacells <- combined_metacells[,-41]
NBp <- p.vector(combined_metacells, design, counts=F, 
                min.obs = 6) #min obs should be degree + 1 x groups + 1 (3x2 for us)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all", rsq = 0.5)
get$summary


#####after running masigpro for each...####
#I copied the gene names from get."name of cell type"$summary to time.series.summary.celltype excel sheet

#####find average # metacells per time point by cell type for seacells #####
avg.list <- list()

for (i in cell.types) {
  type <- design.list[[i]]
  length <- length(table(type$edesign$time))
  avg <- sum(table(type$edesign$time)) / length
  avg.list[[i]] <- avg
}
