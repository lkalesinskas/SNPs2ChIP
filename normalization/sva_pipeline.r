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
matrix_type <- "score" # choose score or binary matrix
dat <- fread(paste0("data/", matrix_type, "/final_1/matrix.csv"), header=F)
colnames(dat) <- gsub(" ", "", apply(fread(paste0("data/", matrix_type, "/final_1/columns.csv"), header=F), 1, function(i) paste0(i[1], "_", i[4])))
row.names <- gsub(".bed", "", unlist(strsplit(as.character(unlist(fread(paste0("data/", matrix_type, "/final_1/rows.csv"), header=F))), "_"))[c(F,F,F,F,F,F,T)])

for (i in 2:22) {
	message(paste("Chr:", i))
	tmp <- fread(paste0("data/", matrix_type, "/final_", i, "/matrix.csv"), header=F)
	colnames(tmp) <- gsub(" ", "", apply(fread(paste0("data/", matrix_type, "/final_", i, "/columns.csv"), header=F), 1, function(i) paste0(i[1], "_", i[4])))
	tmp_row_names <- gsub(".bed", "", unlist(strsplit(as.character(unlist(fread(paste0("data/", matrix_type, "/final_", i, "/rows.csv"), header=F))), "_"))[c(F,F,F,F,F,F,T)])
	stopifnot(tmp_row_names==row.names)
	dat <- cbind(dat, tmp)
}

rownames(dat) <- row.names

## Load metadata
metadata <- fread("exp.txt") %>% filter(V1 %in% row.names)
colnames(metadata) <- gsub(" ", "", colnames(metadata))
stopifnot(metadata$V1==row.names) # make sure sample numbers match in data and metadata

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
		metadata_combine[i, j] <- tolower(gsub(" ", "", as.character(unlist(metadata_split[[i]][, j]))))
	}
}

## Combine population and ethnicity columns
metadata_combine$pop_name <- tolower(ifelse(is.na(metadata_combine$ethnicity), metadata_combine$populationname, metadata_combine$ethnicity))

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
dat_scale_t <- t(dat_scale)
sva_fit <- sva(dat_scale_t, mod, method="two-step")

## Find SVs significantly associated with covariate of interest
meta_index <- which(!is.na(metadata_combine[, "antibody"])) # index of samples with metadata for the covariate specified
sv_lm <- apply(sva_fit$sv, 2, function(x) summary(lm(x[meta_index] ~ factor(metadata_combine[meta_index, "antibody"])))[[4]][-1,4])
sig_sv_lm <- which(apply(sv_lm, 2, function(x) any(x <= 1e-80)))

## Regress out SVs
modsv <- cbind(mod, sva_fit$sv[, -sig_sv_lm])
fitsv <- lm.fit(modsv, dat_scale)

## Center and scale
fitsv_scale <- apply(fitsv$residuals, 2, function(x) scale(x, center=TRUE, scale=T))
rownames(fitsv_scale) <- row.names

## Write corrected data
write.table(fitsv_scale, "corrected_data_genome_score_31may18.txt", sep="\t", row.names=F, col.names=T) # write corrected matrix
write.table(as.data.frame(sva_fit$sv), "sig_svs_21may18.txt", sep="\t", row.names=F, col.names=F) # write SVs

## Subset data and plot PCs
meta_choice <- "cellline"
meta_index <- which(!is.na(metadata_combine[, meta_choice])) # index of samples with metadata for the covariate specified

## Plot PCA
pca_object <- prcomp(fitsv$residuals[meta_index, ], center=TRUE, scale=TRUE)
combined_nonnum <- cbind(fitsv$residuals[meta_index, ], as.data.frame(metadata_combine[meta_index, ])) # non-numeric holder

p1 <- autoplot(pca_object, scale=0, data=combined_nonnum, colour="cellline") + #, alpha="pop_name") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	guides(colour=FALSE) +
	guides(alpha=FALSE) +
	#scale_colour_manual(values=c("#F8766D", "gray50", "#00BFC4")) +
	#scale_alpha_manual(values=c(1, 0.1, 1)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) #+ xlim(-200, 210) + ylim(-140, 100) # optionally manually set axes length
ggsave("gene_sva.pdf", p1, width=4, height=3)

## Plot PCA (no correction - for comparison)
pca_object_none <- prcomp(dat_scale[meta_index, ], center=F, scale=F)

p2 <- autoplot(pca_object_none, scale=0, data=combined_nonnum, colour="cellline") + #, alpha="pop_name") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	guides(colour=FALSE) +
	guides(alpha=FALSE) +
	#scale_colour_manual(values=c("#F8766D", "gray50", "#00BFC4")) +
	#scale_alpha_manual(values=c(1, 0.1, 1)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)  # + xlim(-150, 100) + ylim(-125, 100) # optionally manually set axes length
ggsave("gene_no_sva.pdf", p2, width=4, height=3)


## Plot all points, gray for points with no metadata value
pca_object <- prcomp(fitsv$residuals, center=TRUE, scale=TRUE)
pca_object_none <- prcomp(dat_scale, center=F, scale=F)
combined_nonnum <- cbind(fitsv$residuals, as.data.frame(metadata_combine)) # non-numeric holder

# Find PCs most correlated with covariate of interest
variable_of_interest <- "pop_name"
combined_nonnum[is.na(combined_nonnum[, variable_of_interest]), variable_of_interest] <- "unknown"
meta_index <- which(metadata_combine[, variable_of_interest] != "unknown") # index of samples with metadata for the covariate specified
pca_filter <- pca_object_none$x[meta_index, ]
sig_pc <- apply(pca_filter, 2, function(x) summary(lm(x ~ factor(metadata_combine[meta_index, variable_of_interest])))[[4]][-1,4])
sig_pick <- apply(sig_pc, 2, function(x) sum(x <= 1e-30))
sig_pick[order(sig_pick, decreasing=T)][1:10]
which(sig_pc < 1e-3)

# Plot
pca_object_current <- pca_object

p2 <- ggplot(pca_object_current$x) +
	geom_point(aes(x=PC1, y=PC2), colour="gray20", alpha=0.3) +
	geom_point(data=pca_object_current$x[-which(combined_nonnum[, variable_of_interest]=="unknown"), ], aes(x=PC1, y=PC2, 
		colour=combined_nonnum[-which(combined_nonnum[, variable_of_interest]=="unknown"), variable_of_interest])) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	guides(colour=FALSE) +
	guides(alpha=FALSE) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) + xlim(-600, 600) + ylim(-400, 400) # optionally manually set axes length
ggsave("plot_pc_all_points.pdf", p2, width=4, height=3)


## Plot metadata sparsity
metadata_sparse_filter <- metadata_sparse[-nrow(metadata_sparse), ]
metadata_sparse_filter$rownum <- seq(1,nrow(metadata_sparse_filter), 1)

p3 <- ggplot(metadata_sparse_filter) +
	geom_col(aes(x=rownum, y=percent_missing)) +
	scale_x_continuous(breaks=seq(1,nrow(metadata_sparse_filter),1), labels=as.character(unlist(metadata_sparse_filter$col_name))) +
	labs(x="Metadata Name", y="Percent Missing (%)") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text.x=element_text(angle=45, hjust=1),
		axis.text=element_text(size=8),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) 
ggsave("metdata_sparse.pdf", p3, width=6, height=3.5)

## Plot SV sig associated with antibody

summary(lm(sva_fit$sv ~ factor(metadata_combine$pop_name))) 
sv_correlation <- data.frame(sv=sva_fit$sv[, 4], meta=factor(metadata_combine$pop_name))
sv_correlation <- subset(sv_correlation, !is.na(meta))
#sv_correlation <- subset(sv_correlation, meta != "none" & meta != "none(input)")

p4 <- ggplot(sv_correlation, aes(x=meta, y=sv)) +
	geom_jitter(alpha=0.1) +
	geom_boxplot(aes(fill=meta), alpha=0.7) +
	labs(x="Ancestry", y="") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text.x=element_text(angle=45, hjust=1),
		axis.text=element_text(size=9),
		axis.title=element_text(size=10),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(fill=FALSE)
ggsave("sv_correlation.pdf", p4, width=3, height=2.5)




