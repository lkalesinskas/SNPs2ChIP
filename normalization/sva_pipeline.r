#!/bin/R

# sva.r
# Description: Performs Surrogate Variable Analysis (SVA) and validation
# Code author: Craig Smail

## Load required libraries
library(data.table)
library(sva)
library(ggplot2)
library(ggfortify)
library(stringr)
library(dplyr)

## Set working directory
setwd("/users/csmail/bmi212/")

## Read data
dat <- fread("ourMatrix.csv")

## Read rownames in dat
row_names <- fread("row_names.csv", header=FALSE)
row_names <- unlist(strsplit(as.character(unlist(row_names)), '[_]'))[c(F,F,F,T)]
row_names <- gsub(".bed", "", row_names)

## Load metadata
metadata <- fread("exp.txt") %>% filter(V1 %in% row_names)
colnames(metadata) <- gsub(" ", "", colnames(metadata))
stopifnot(metadata$V1==row_names) # make sure sample numbers match in data and metadata

## Split metadata column
metadata_split <- lapply(metadata$Metadata, function(x) {
			split_col <- data.frame(strsplit(strsplit(gsub("?Ref=", "", x), "\\|\\|")[[1]], "="))
			colnames(split_col) <- tolower(gsub(" ", "", as.character(unlist(split_col[1, ]))))
			split_col <- split_col[2, ]
			return(split_col)	
		})

metadata_combine <- data.frame(metadata$V1, stringsAsFactors=F)
colnames(metadata_combine) <- "sample_id"

for (i in 1:length(metadata_split)) {
	cols <- unique(colnames(metadata_split[[i]]))
	for (j in cols) {
		metadata_combine[i, j] <- gsub(" ", "", as.character(unlist(metadata_split[[i]][, j])))
	}
}

## Sparsity of metadata
metadata_sparse <- data.frame(colnames(metadata_combine))
colnames(metadata_sparse) <- "col_name"
metadata_sparse$percent_missing <- apply(metadata_combine, 2, function(x) length(which(is.na(x)))/length(x)*100)
metadata_sparse <- metadata_sparse[order(metadata_sparse$percent_missing, decreasing=T), ]

## Remove genes with zero variance
dat_remove_zero_var <- dat[, apply(dat, 2, function(x) !var(x) == 0), with=FALSE] 

## Unit variance and center for cols with var != 0
dat_scale <- apply(dat_remove_zero_var, 2, function(x) scale(x, center=TRUE, scale=TRUE))

## SVA
mod <- model.matrix(~1, data=as.data.frame(dat_scale))
sva_fit <- sva(t(dat_scale), mod, method="two-step")

## Regress out SVs
modsv <- cbind(mod, sva_fit$sv)
fitsv <- lm.fit(modsv, dat_scale)

## Center and scale
fitsv_scale <- apply(fitsv$residuals, 2, function(x) scale(x, center=TRUE, scale=T))
rownames(fitsv_scale) <- row_names

## Write corrected data
write.table(fitsv_scale, "corrected_data_21may18.txt", sep="\t", row.names=F, col.names=F) # write corrected matrix
write.table(as.data.frame(sva_fit$sv), "sig_svs_21may18.txt", sep="\t", row.names=F, col.names=F) # write SVs

## <----- TODO: find metadata correlated with PCs...

## Plot PCA
pca_object <- prcomp(fitsv$residuals, center=TRUE, scale=TRUE)
combined_nonnum <- cbind(fitsv$residuals, as.data.frame(metadata)) # non-numeric holder

p1 <- autoplot(pca_object, scale=0, data=combined_nonnum, colour="Antigenclass", alpha=0.9) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) #+ xlim(-200, 210) + ylim(-140, 100) # optionally manually set axes length
ggsave("gene_sva.pdf", p1, width=5, height=4)

## Plot PCA (no correction - for comparison)
pca_object_none <- prcomp(dat_scale, center=TRUE, scale=TRUE)

p2 <- autoplot(pca_object_none, scale=0, data=combined_nonnum, colour="Antigenclass", alpha=0.9) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)  # + xlim(-150, 100) + ylim(-125, 100) # optionally manually set axes length
ggsave("gene_no_sva.pdf", p2, width=5, height=4)

