## Developed and Maintained by Sanmati Dugad

rm(list=ls())    #clear the environment in R. Remove all objects (variables, functions, etc.) from the current workspace.
Sys.setenv(TZ='EDT')

ReQ_Modules = c("tidyr","tidyverse", "data.table", "reshape2", "colorspace", "readxl")

for (each in 1:length(ReQ_Modules)) {
  if(ReQ_Modules[each] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_Modules[each])
  }
}

#invisible(lapply(ReQ_Modules, library, character.only = TRUE))

suppressWarnings(library(tidyr))
suppressWarnings(library(tidyverse))
suppressWarnings(library(data.table))
suppressWarnings(library(reshape2))
suppressWarnings(library(colorspace))
suppressWarnings(library(readxl))

args=commandArgs(trailingOnly = T)


scaled_data = NULL    # scaled output
gene_states = NULL    # gene states output
z_statistics = NULL   # z_stats output for a particular sample. 
gene_panel = NULL     # gene panel xls file
panel_name = NULL     # sheet name in the xls file

normal_count = NULL   # number of normal samples
tumor_sample = NULL   # column name of the tumor sample
#title = NULL          # Title to Provide to your Boxplot. It is kind of automated now. Can be removed 
state_type = NULL    ## HIGH, VERY HIGH, ALL, HIGH and VERY HIGH, VERY LOW, LOW, VERY LOW and LOW


if ("-i" %in% args && "-gs" %in% args && "-gp" %in% args && "-n" %in% args && "-t" %in% args && "-p" %in% args && "-st" %in% args && "-gpn" %in% args ) {    
  scaled_data_index = which(args == "-i")    #Getting index of the "--input" option
  gene_states_index = which(args == "-gs")
  z_statistics_index = which(args == "-z")
  gene_panel_index = which(args == "-gp")
  normal_count_index = which(args == "-n")
  tumor_sample_index = which(args == "-t")
  state_type_index = which(args == "st")
  panel_name_index = which(args == "-gpn")
  
  scaled_data <-  args[scaled_data_index + 1] #"union_normalized.txt"
  gene_states <- args[gene_states_index + 1] #"ID_SYMBOL.csv"
  gene_panel <- args[gene_panel_index + 1]
  normal_count = as.numeric(args[normal_count_index + 1])
  tumor_sample = args[tumor_sample_index + 1]
  state_type = args[state_type_index + 1]
  panel_name = args[panel_name_index + 1]
}


## Assigning Variables
setwd("~/3. Boxplot_Analysis/X1107H4D/Input_Files/")

scaled_data = "or4r-scaled.csv"
gene_states = "or4r-gene_states.csv"
z_statistics = "X1107H4D-z_stats.csv"

## Research_Gene_Panel
gene_panel = "~/Gene Panel List and Report Table.xlsx"
panel_name = "Compliment"
normal_count = 24
tumor_sample = "X1107H4D"
state_type = "ALL"    # other options "VERY HIGH, HIGH, VERY HIGH and HIGH, 


title = paste(tumor_sample, ":", panel_name, "Gene Panel -", state_type, "genes")

print(paste("your scaled data file name is:", scaled_data))
print(paste("Your gene states file name is:", gene_states))
print(paste("Your gene panel file name is:", gene_panel))

## reading files
scaled_input = read.csv(scaled_data, row.names = 1)    # scaled data
states_input = read.csv(gene_states, row.names = 1)    # states data
z_input = read.csv(z_statistics, row.names = 1)
gene_list_data = read_xlsx(gene_panel, sheet = panel_name)    # will gene list be always provided as an excel sheet? 

## checking 
genes_list = as.vector(na.omit(unique(trimws(gene_list_data$Gene_ID))))  ## to keep Gene Column name same throughout and without spaces. (SYMBOL)
print(paste("number of genes present in panel are", length(genes_list)))

## Eliminating Gene If Any of them are False Positive or False Negative
True_Positives = z_input$SYMBOL[z_input$Potential_FP_or_FN == 0]
genes_list = genes_list[genes_list %in% True_Positives]



## selecting genes from states data present in the gene panel
states_input.1 = states_input[trimws(states_input$SYMBOL) %in% genes_list, c("SYMBOL",  tumor_sample), drop = FALSE]
states_input.1 = states_input.1[!duplicated(states_input.1$SYMBOL), ]
a = setdiff(states_input.1$SYMBOL, genes_list)
b = setdiff(genes_list, states_input.1$SYMBOL)
print(paste("Genes not selected: ", a, b, FALSE))

## five types of condition to seperate genes
aberant_states = c("VERY HIGH", "HIGH", "NORMAL","LOW", "VERY LOW")
aberant_states_HIGH = c( "HIGH")
aberant_states_VERY_HIGH = c("VERY HIGH")
aberant_states_LOW = c("LOW")
aberant_states_VERY_LOW = c("VERY LOW")
aberant_states_NORMAL = c("NORMAL")


## Categorizing genes based on their gene state for a particular tumor sample
result_genes = c()
result_genes_VERY_HIGH = c()
result_genes_HIGH = c()
result_genes_LOW = c()
result_genes_VERY_LOW = c()
result_genes_NORMAL = c()


for (i in aberant_states) {
  result_genes <- c(result_genes, states_input.1$SYMBOL[states_input.1[, tumor_sample] == i])
}

for (i in aberant_states_VERY_HIGH) {
  result_genes_VERY_HIGH <- sort(c(result_genes_VERY_HIGH, states_input.1$SYMBOL[states_input.1[, tumor_sample] == i]))
}

for (i in aberant_states_HIGH) {
  result_genes_HIGH <- sort(c(result_genes_HIGH, states_input.1$SYMBOL[states_input.1[, tumor_sample] == i]))
}

for (i in aberant_states_LOW) {
  result_genes_LOW <- sort(c(result_genes_LOW, states_input.1$SYMBOL[states_input.1[, tumor_sample] == i]))
}

for (i in aberant_states_VERY_LOW) {
  result_genes_VERY_LOW <- sort(c(result_genes_VERY_LOW, states_input.1$SYMBOL[states_input.1[, tumor_sample] == i]))
}

for (i in aberant_states_NORMAL) {
  result_genes_NORMAL <- sort(c(result_genes_NORMAL, states_input.1$SYMBOL[states_input.1[, tumor_sample] == i]))
}

## selecting genes in scaled data
if (state_type == "HIGH"){
  scaled_input.1 = scaled_input[scaled_input$SYMBOL %in% result_genes_HIGH, ]
}else if (state_type == "VERY HIGH"){
  scaled_input.1 = scaled_input[scaled_input$SYMBOL %in% result_genes_VERY_HIGH, ]
}else if (state_type == "ALL"){
  scaled_input.2 = scaled_input[scaled_input$SYMBOL %in% result_genes, ]
  scaled_input.1 = scaled_input.2[match(result_genes, scaled_input.2$SYMBOL), ]
}else if (state_type == "VERY HIGH and HIGH"){
  scaled_input.2 = scaled_input[scaled_input$SYMBOL %in% c(result_genes_VERY_HIGH, result_genes_HIGH), ]
  scaled_input.1 = scaled_input.2[match(c(result_genes_VERY_HIGH, result_genes_HIGH), scaled_input.2$SYMBOL), ]
}else if (state_type == "HIGH and VERY HIGH"){
  scaled_input.2 = scaled_input[scaled_input$SYMBOL %in% c(result_genes_VERY_HIGH, result_genes_HIGH), ]
  scaled_input.1 = scaled_input.2[match(c(result_genes_VERY_HIGH, result_genes_HIGH), scaled_input.2$SYMBOL), ]
}else if (state_type == "VERY LOW"){
  scaled_input.1 = scaled_input[scaled_input$SYMBOL %in% result_genes_VERY_LOW, ]
}else if (state_type == "LOW"){
  scaled_input.1 = scaled_input[scaled_input$SYMBOL %in% result_genes_LOW, ]
}else if (state_type == "VERY LOW and LOW"){
  scaled_input.2 = scaled_input[scaled_input$SYMBOL %in% c(result_genes_VERY_LOW, result_genes_LOW), ]
  scaled_input.1 = scaled_input.2[match(c(result_genes_LOW, result_genes_VERY_LOW), scaled_input.2$SYMBOL), ]
}else if (state_type == "NORMAL"){
  scaled_input.2 = scaled_input[scaled_input$SYMBOL %in% c(result_genes_NORMAL), ]
  scaled_input.1 = scaled_input.2[match(c(result_genes_NORMAL), scaled_input.2$SYMBOL), ]
}

scaled_input.1 = scaled_input.1[!duplicated(scaled_input.1$SYMBOL), ]
rownames(scaled_input.1) = scaled_input.1$SYMBOL
An = t(scaled_input.1[,2:(normal_count+1)])     #normal samples
An.1 = rownames_to_column(data.frame(An), var = "SAMPLE")

At = t(scaled_input.1[,c(tumor_sample), drop = FALSE])    #selected tumor sample
At.1 = rownames_to_column(data.frame(At), var = "SAMPLE")

#In boxplot its always (value column ~ x-axis column)
A_long <- reshape2::melt(An.1,  id.vars = "SAMPLE", variable.name = "SYMBOL", value.name = "VALUE")
head(A_long)

At_long <-  reshape2::melt(At.1, id.vars = "SAMPLE", variable.name = "SYMBOL", value.name = "VALUE")
head(At_long)

png(paste("Boxplot -", tumor_sample, "-", panel_name, "-", state_type, ".png"), width = 1500, height = 800)
## Generating Boxplots
boxplot(VALUE~SYMBOL, data = A_long, col = "gray70", #col = c("red", "chocolate", "orange", "khaki", "chartreuse", "chartreuse3", "chartreuse4", "aquamarine", "aquamarine4", "cyan", "darkorchid1", "darkorchid4", "magenta") ,
        whiskcol = "black" , staplecol = "black", boxcol = "black", outcol = "chocolate1", outcex = 1.5, outpch = 21 , outbg = "yellow",
        ylim = c(0,25), xaxt = "n",
        xlab = substitute(paste(bold("Gene Name"))), ylab = substitute(paste(bold("log2 Normalized Gene Expression Value"))),
        main = title,
        cex.lab = 1, cex.main = 1.5, cex.axis = 1.1, lwd = 1.5, lty = 1, boxlwd = 1.5)
abline(v = seq(1,40, 1), lty = "solid", col = "gray")
abline(h = seq(0,40, 5), lty = "solid", col = "gray")

boxplot(VALUE~SYMBOL, data = A_long, col = "gray70", #col = c("red", "chocolate", "orange", "khaki", "chartreuse", "chartreuse3", "chartreuse4", "aquamarine", "aquamarine4", "cyan", "darkorchid1", "darkorchid4", "magenta") ,
        whiskcol = "black" , staplecol = "black", boxcol = "black", outcol = "chocolate1", outcex = 1.5, outpch = 21 , outbg = "yellow",
        xaxt = "n",
        xlab = substitute(paste(bold("Gene Name"))), ylab = substitute(paste(bold("log2 Normalized Gene Expression Value"))),
        cex.lab = 1, cex.main = 1.5, cex.axis = 1.1, lwd = 1.5, lty = 1, boxlwd = 1.5, add = T)
legend("topleft", c("Normal_Boxplot","Normal_Outlier", tumor_sample), col = c("lightgray","yellow", "brown1"), pch = c(15, 15, 19), cex = 1.0)
legend("topright", c("Gene Names on X-axis:","VERY HIGH", "HIGH genes","NORMAL" ,"LOW", "VERY LOW"), col = c("black","red","orange", "black","steelblue1", "blue"), pch = c(NA, 15, 15,15, 15, 15), cex = 1.0)

stripchart(VALUE~SYMBOL, data = At_long , col = "brown4", add = T, pch = 1, cex = 1.7, vertical = T)
stripchart(VALUE~SYMBOL, data = At_long , col = "brown1", add = T, pch = 19, cex = 1.5, vertical = T)
#title(sub = "VERY HIGH - red, HIGH - orange, LOW - light blue, VERY LOW - blue")

if (state_type == "HIGH"){
  col = "orange"
}else if (state_type == "VERY HIGH"){
  col = "red"
}else if (state_type == "ALL"){
  col = c(rep("red", length(result_genes_VERY_HIGH)), rep("orange", length(result_genes_HIGH)), rep("black", length(result_genes_NORMAL)),rep("steelblue1", length(result_genes_LOW)), rep("blue", length(result_genes_VERY_LOW)))
}else if (state_type == "VERY HIGH and HIGH"){
  col = c(rep("red", length(result_genes_VERY_HIGH)), rep("orange", length(result_genes_HIGH)))
}else if (state_type == "HIGH and VERY HIGH"){
  col = c(rep("red", length(result_genes_VERY_HIGH)), rep("orange", length(result_genes_HIGH)))
}else if (state_type == "VERY LOW"){
  col = "darkblue"
}else if (state_type == "LOW"){
  col = "steelblue1"
}else if (state_type == "VERY LOW and LOW"){
  col = c(rep("steelblue1", length(result_genes_LOW)), rep("darkblue", length(result_genes_VERY_LOW)))
}


#axis(1, 1:ncol(An), colnames(An), cex.axis =1.0,font = 1,xaxt = "n")
mtext(side = 1, text = colnames(An), at= 1:ncol(An), col = col, line = 0.1, font = 2, las = 2, cex = 0.85)

#mtext(side = 1, text = colnames(An), at= 1:ncol(An), col = col, line = 0.5, font = 2)

# mtext(side = 1, text = "State of Gene is Color Coordinated as: ", col = "black", line = 4, font = 1, at = 1.0)
# mtext(side = 1, text = "VERY HIGH", col = "red", line = 4, font = 2, at = 3.0)
# mtext(side = 1, text = ", ", col = "black", line = 4, font = 2, at = 3.5)
# mtext(side = 1, text = "HIGH", col = "orange", line = 4, font = 2, at = 4.0)
# mtext(side = 1, text = ",", col = "black", line = 4, font = 2, at = 4.5)
# mtext(side = 1, text = "LOW", col = "steelblue1", line = 4, font = 2, at = 5.0)
# mtext(side = 1, text = ",", col = "black", line = 4, font = 2, at = 5.5)
# mtext(side = 1, text = "VERY LOW", col = "darkblue", line = 4, font = 2, at = 6.5)

dev.off()

