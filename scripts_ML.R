# Code overview: combinations of 10 algorithms for variable selection and prognostic model construction
# 


# Set working directory
work.path <- "/Users/Administrator/Documents/Aljar_machine_learn/101_example/"; setwd(work.path) 

# Set other paths
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

# Create these directories if they do not exist
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# install.packages("randomForestSRC")
# install.packages("snowfall")

# Load the required R packages
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

# Load the scripts for model training and evaluation
source(file.path(code.path, "ML.R"))

## Training Cohort ---------------------------------------------------------
# The training expression matrix has genes (genes of interest) in rows and samples in columns
# (gene identifiers should be in the same format as those in the test set, such as SYMBOL or ENSEMBL)
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_expr <- Train_expr[rowSums(Train_expr > 0) > ncol(Train_expr) * 0.1, ] # Remove genes with excessive zero expression to avoid errors during model fitting
# The training survival data frame has samples in rows and outcome information in columns
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_surv <- Train_surv[Train_surv$OS.time > 0, c("OS", "OS.time")] # Extract samples with OS.time > 0
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
# The test expression matrix has genes (genes of interest) in rows and samples in columns
# (gene identifiers should be in the same format as those in the training set, such as SYMBOL or ENSEMBL)
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# The test survival data frame has samples in rows and outcome information in columns
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_surv <- Test_surv[Test_surv$OS.time > 0, c("Coho","OS", "OS.time")] # Extract samples with OS.time > 0
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# Extract common genes
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # Input expression matrix for the model: rows are samples and columns are genes
Test_expr <- t(Test_expr[comgene,]) # Input expression matrix for the model: rows are samples and columns are genes

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# Record the models to be run here, in the format:
# Algorithm1[parameters] + Algorithm2[parameters]
# Currently, only StepCox and RunEnet support algorithm parameter input
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
head(methods)

## Train the model --------------------------------------------------------
model <- list()
set.seed(seed = 123)
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # Name of the algorithm used in this round
  method <- strsplit(method, "\\+")[[1]] # Algorithm name for each step
  
  Variable = colnames(Train_expr) # Variables ultimately used to build the model
  for (i in 1:length(method)){
    if (i < length(method)){
      selected.var <- RunML(method = method[i], # Machine learning method
                            Train_expr = Train_expr, # Variables with potential predictive value in the training set
                            Train_surv = Train_surv, # Training survival data
                            mode = "Variable",       # Running mode: Variable (feature selection) or Model (model fitting)
                            timeVar = "OS.time", statusVar = "OS") # Survival variables used for training; must exist in Train_surv
      if (length(selected.var) > 5) Variable <- intersect(Variable, selected.var)
    } else {
      model[[method_name]] <- RunML(method = method[i],
                                    Train_expr = Train_expr[, Variable],
                                    Train_surv = Train_surv,
                                    mode = "Model",
                                    timeVar = "OS.time", statusVar = "OS")
    }
  }
}
saveRDS(model, file.path(res.path, "model.rds"))

## Evaluate the model -----------------------------------------------------

# Read the list of models that have already been generated
model <- readRDS(file.path(res.path, "model.rds"))
summary(Train_expr)
summary(Test_expr)

Train_expr <- scale(Train_expr)
Test_expr <- scale(Test_expr)

# Calculate the C-index for each model
#is.finite(Test_expr)
#Test_expr[!is.finite(Test_expr)] <- NA
Cindexlist <- list()
for (method in methods){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # Prognostic model
                                  Test_expr = Test_expr, # Prognostic variables in the test set; should include all variables used in training, otherwise errors may occur
                                  Test_surv = Test_surv, # Test survival data; should contain the required variables, otherwise errors may occur
                                  Train_expr = Train_expr, # Provide training expression data if the training set also needs to be evaluated; otherwise set to NULL
                                  Train_surv = Train_surv, # Provide training survival data if the training set also needs to be evaluated; otherwise set to NULL
                                  Train_name = "TCGA", # Label for the training set if it also needs to be evaluated; otherwise it will be treated as "Training"
                                  cohortVar = "Coho", # Important: variable used to specify the cohort; this column must exist and be correctly assigned [default is "Cohort"], otherwise errors may occur
                                  timeVar = "OS.time", # Survival time used for evaluation; must exist in Test_surv
                                  statusVar = "OS") # Survival status used for evaluation; must exist in Test_surv
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "Cindex_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "Cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- apply(Cindex_mat, 1, mean)           # Calculate the average C-index of each algorithm across all cohorts
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # Sort algorithms by C-index from high to low
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # Reorder the C-index matrix

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # Keep three decimal places
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)

#CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # Set cohort colors

CohortCol <- c("steelblue", "firebrick","green") # You may replace these colors with your preferred ones

names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)

cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              # col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # Yellow-green color scheme
              col = c("#4195C1", "#FFFFFF", "#CB5746"), # Red-blue color scheme
              rect_gp = gpar(col = "black", lwd = 1), # Set borders to black
              cluster_columns = FALSE, cluster_rows = FALSE, # No clustering performed; not meaningful here
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # Add text to each grid cell
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

#pdf(file.path(fig.path, "Cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 3, height = cellheight * nrow(Cindex_mat) * 0.45)
pdf(file.path(fig.path, "Cindex2.pdf"), width = (cellwidth * ncol(Cindex_mat) + 3) * 2, height = cellheight * nrow(Cindex_mat) * 0.45)

draw(hm)
invisible(dev.off())

# Select the optimal genes
for (method in methods){
  
  if (method_name == "StepCox[both]+SuperPC") 
    data1=data.frame(Variable)
    write.table(data1, file = "selected_genes.txt", row.names = FALSE, col.names = FALSE)
  
  
}
