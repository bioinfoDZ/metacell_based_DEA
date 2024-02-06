#this code depends on having combined sMetacells, pure sMetacells, rMetacells, and pseudobulked matrices from other scripts

#### determine sd for each cell type regardless of time, plot
cell.types <- c("Cytotoxic CD8 T cells","Naive T cells", "NKs","MAIT", "Activated CD4 T cells",
                "Naive B cells","Plasma","Memory B cells","XCL+ NKs","Cycling T cells")

#for all genes#####
sd.all.list <- list()

for (i in cell.types) {
  subset.sd <- subset(covid.subset.time, subset = cell_type == i)
  subset.sd <- subset.sd@assays$RNA@data
  subset.sd <- apply(subset.sd, 1, sd)
  #subset.sd <- mean(subset.sd)
  sd.all.list[[paste(i, ".sd.all", sep = "")]] <- subset.sd
}

sd.random.list <- list()

for (i in cell.types) {
  if (i == "XCL+ NKs") {
    metacells_names <- grep("\\XCL\\+ NKs$", ls(), value = TRUE)
  } else {
    metacells_names <- grep(paste("\\.", i, "$", sep = ""), ls(),
                            value = TRUE)
  }
  
  metacells_list <- mget(metacells_names)
  combined_metacells <- do.call(cbind, metacells_list)
  combined_metacells <- apply(combined_metacells, 1, sd)
  #combined_metacells <- mean(combined_metacells)
  sd.random.list[[paste(i, ".sd.random", sep = "")]] <- combined_metacells
}

#for seacells 
sd.seacell.list <- list()

for (i in cell.types){
  sd.seacell <- apply(covid.masigpro.list[[i]], 1, sd)
  #sd.seacell <- mean(sd.seacell)
  sd.seacell.list[[paste(i, ".sd.seacell", sep = "")]] <- sd.seacell
}

#for seacells with high purity
sd.seacell.pure.list <- list()

for (i in cell.types){
  sd.seacell <- apply(covid.masigpro.list.pure[[i]], 1, sd)
  #sd.seacell <- mean(sd.seacell)
  sd.seacell.pure.list[[paste(i, ".sd.seacell", sep = "")]] <- sd.seacell
}

#for pb cells
sd.pb.list <- list()

for (i in cell.types) {
  sd.pb <- apply(covid.subset.list[[i]]$RNA, 1, sd)
  #sd.pb <- mean(sd.pb)
  sd.pb.list[[paste(i, ".sd.pb", sep = "")]] <- sd.pb
}

#violin plot for all sds ####
sd.all.unlist <- as.data.frame(sd.all.list)
sd.pb.unlist <- as.data.frame(sd.pb.list)
sd.seacell.unlist <- as.data.frame(sd.seacell.list)
sd.seacell.pure.unlist <- as.data.frame(sd.seacell.pure.list)
sd.random.unlist <- as.data.frame(sd.random.list)
sd.all.unlist$genes <- rownames(sd.all.unlist)
sd.pb.unlist$genes <- rownames(sd.pb.unlist)
sd.seacell.unlist$genes <- rownames(sd.seacell.unlist)
sd.seacell.pure.unlist$genes <- rownames(sd.seacell.pure.unlist)
sd.random.unlist$genes <- rownames(sd.random.unlist)

library(tidyr)

# Convert data frames to long format
sd.all_long <- gather(sd.all.unlist, key = "source", value = "value", -genes)
sd.pb_long <- gather(sd.pb.unlist, key = "source", value = "value", -genes)
sd.seacell_long <- gather(sd.seacell.unlist, key = "source", value = "value", -genes)
sd.seacell.pure_long <- gather(sd.seacell.pure.unlist, key = "source", value = "value", -genes)
sd.random_long <- gather(sd.random.unlist, key = "source", value = "value", -genes)

# Combine the long-format data frames
stacked_data <- rbind(sd.pb_long, sd.seacell_long, sd.random_long,sd.seacell.pure_long, sd.all_long)
write.csv(stacked_data, "stacked.sds.csv")
stacked_data <- read.csv("stacked.sds.csv")

#plot
stacked_data <- stacked_data[stacked_data$Method != "all",]
ggplot(data = stacked_data, aes(x = Cell.Type, y = value, 
                                color = Method))+
  scale_y_continuous(limits =c(0,1))+
  geom_boxplot(width = 0.5, position = position_dodge(0.75),
               outlier.shape = NA)+
  xlab("Cell Type")+
  ylab("SD")+
  scale_color_discrete(labels = c("Pseudobulk", 
                                  "rMetacells", "sMetacells",
                                  "sMetacells Pure"))+
  scale_x_discrete(labels = c("Activated CD4 T Cells",
                                "Cycling T Cells",
                                "Cytotoxic CD8 T Cells",
                                "MAIT",
                                "Memory B Cells",
                                "Naive B Cells",
                                "Naive T cells",
                                'NKs',
                                "PLasma",
                                "XCL+ NKs"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#random ####
for (i in 1:10) {
  cell_type_name <- sub("\\.sd\\.all$", "", names(combined.sds$sd.all.list[i]))
  
  #barplot all cell types 
  barplot(combined.sds)
  # Create the first plot with smaller dots
  plot(combined.sds$sd.pb.list[[i]], combined.sds$sd.random.list[[i]],
       xlim = c(0, 0.5), ylim = c(0, 0.5), xlab = "",
       ylab = "", main = cell_type_name,
       cex = 0.5)
  
  points(combined.sds$sd.all.list[[i]], combined.sds$sd.seacell.list[[i]], 
         col = "red", cex = 0.5)
  
  # Create the second plot with larger dots
  plot(combined.sds$sd.all.list[[i]], combined.sds$sd.random.list[[i]],
       xlim = c(0, 2), ylim = c(0, 2), xlab = "",
       ylab = "", main = cell_type_name,
       cex = 0.5)
  
  points(combined.sds$sd.all.list[[i]], combined.sds$sd.seacell.list[[i]], 
         col = "red", cex = 0.5)
}

for (i in 6:10) {
  cell_type_name <- sub("\\.sd\\.all$", "", names(combined.sds$sd.all.list[i]))
  
  # Create the first plot with smaller dots
  plot(combined.sds$sd.all.list[[i]], combined.sds$sd.random.list[[i]],
       xlim = c(0, 0.5), ylim = c(0, 0.5), xlab = "",
       ylab = "", main = cell_type_name,
       cex = 0.5)
  
  points(combined.sds$sd.all.list[[i]], combined.sds$sd.seacell.list[[i]], 
         col = "red", cex = 0.5)
  
  # Create the second plot with larger dots
  plot(combined.sds$sd.all.list[[i]], combined.sds$sd.random.list[[i]],
       xlim = c(0, 2), ylim = c(0, 2), xlab = "",
       ylab = "", main = cell_type_name,
       cex = 0.5)
  
  points(combined.sds$sd.all.list[[i]], combined.sds$sd.seacell.list[[i]], 
         col = "red", cex = 0.5)
}

# Reset the layout to default
par(mfrow = c(1, 1))


# mean sd across genes ####
sd.all.list <- list()

for (i in cell.types) {
  subset.sd <- subset(covid.subset.time, subset = cell_type == i)
  subset.sd <- subset.sd@assays$RNA@data
  subset.sd <- apply(subset.sd, 1, sd)
  subset.sd <- mean(subset.sd)
  sd.all.list[[paste(i, ".sd.all", sep = "")]] <- subset.sd
}

sd.random.list <- list()

for (i in cell.types) {
  if (i == "XCL+ NKs") {
    metacells_names <- grep("\\XCL\\+ NKs$", ls(), value = TRUE)
  } else {
    metacells_names <- grep(paste("\\.", i, "$", sep = ""), ls(),
                            value = TRUE)
  }
  
  metacells_list <- mget(metacells_names)
  combined_metacells <- do.call(cbind, metacells_list)
  combined_metacells <- apply(combined_metacells, 1, sd)
  combined_metacells <- mean(combined_metacells)
  sd.random.list[[paste(i, ".sd.random", sep = "")]] <- combined_metacells
}

#for seacells 
sd.seacell.list <- list()

for (i in cell.types){
  sd.seacell <- apply(covid.masigpro.list[[i]], 1, sd)
  sd.seacell <- mean(sd.seacell)
  sd.seacell.list[[paste(i, ".sd.seacell", sep = "")]] <- sd.seacell
}

combined.sds <- as.data.frame(cbind(sd.all.list, sd.random.list, sd.seacell.list))





#compare pb to metacell to all cells sd####
# covid.masigpro from all.covid.pt.pb.R, n = 16 (5 patients, 10 time points)
# all cells, n = 25775
# random metacells, n = 30
# seacells metacells, n = 360


#barplot by cell type of mSD

#correlation for all genes by type ######

for (i in 1:10) {
  corr.rando <- cor(combined.sds$sd.all.list[[i]],
                    combined.sds$sd.random.list[[i]],
                    method = "spearman")
  
  corr.seacell <- cor(combined.sds$sd.all.list[[i]],
                      combined.sds$sd.seacell.list[[i]],
                      method = "spearman")
  
  print(paste(i, "random correlation coefficient:", corr.rando))
  print(paste(i, "seacells correlation coefficient:", corr.seacell))
}

