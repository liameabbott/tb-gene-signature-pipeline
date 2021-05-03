
# sys args
# 1 - data_dir
# 2 - group 1 category ('HC', 'LTBI', 'ATB', 'OD')
# 3 - group 0 category
# 4 - GSM sample inclusion pattern
# 5 - GSE ID 
# 6 - platform type (microarray or rnaseq)
#
# args if platform_type == 'microarray':
# 7 - micrarrary platform 
#
# args if platform_type == 'rnaseq':
# 7 - RNA-seq counts dataset path
# 8 - number of samples in dataset
# 9 - TMM normalize
# 10 - gene symbol column name
# 11 - column index of first sample
# 12 - column index of last sample

# Example invocation:
#
# #!/bin/bash
# data_dir="/efs/liam/tb-gene-signature-data/roger-data"
# group_0="HC"
# group_1="ATB"
# gsms="X101X00X10X1XX01X01X1X1X1X11X01000X0X1X0XX"
# gse_id="GSE19439"
# microarray_platform="GPL6947"
# platform_type="microarray"
#
# Rscript --vanilla differential-expression-analysis.R \
# $data_dir $group_0 $group_1 $gsms $gse_id $microarray_platform
#
print(.libPaths())
args = commandArgs(trailingOnly=TRUE)
data_dir <- args[1]
groups <- args[2]
gsms <- args[3]
GSE_ID <- args[4]
platform_type <- args[5]
if (platform_type == 'microarray') {
    microarray_platform <- args[6]
} else if (platform_type == 'rnaseq') {
    counts_dataset_path = args[6]
    num_samples = as.integer(args[7])
    TMM_norm = args[8] %in% c('False', 'false', 'FALSE', 1, '1')
    gene_symbol_col = args[9]
    first_sample_column_i = as.integer(args[10])
    last_sample_column_i = as.integer(args[11])
}

groups_split = strsplit(groups, '-')
group_0 = groups_split[[1]][2]
group_1 = groups_split[[1]][1]

library(Biobase)
library(GEOquery)
library(limma)
library(repr)
library(edgeR)
library(preprocessCore)

#define operator
`%notin%` <- Negate(`%in%`)

reduce_expression_matrix_rows <- function(gset, microarray_platform) {
    
    #filter rows of exp matrix, keep only one probe per gene (symbol) with the highest expression sum across samples
    
    #useful commands for de-bugging
    #print(gset)
    #print(varLabels(featureData(gset)))

    #store expression matrix & get the sum for each row/probe from expression matrix
    exp_matrix <- exprs(gset)
    probe_sums <- rowSums(exp_matrix)
    probe_IDs <- names(probe_sums)

    #iterate through all rows/probes and store non-redundant gene symbols (first check to see that string is not empty) store in vector

    #initialize list that will hold information for each row/probe: probe ID (index), list[vector(gene symbols), float(row/probe exp sum)]
    probe_info <- list()

    if ((microarray_platform == 'GPL6947') | (microarray_platform == 'GPL4133') | (microarray_platform == 'GPL10558') | (microarray_platform == 'GPL570') | (microarray_platform == 'GPL11532') | (microarray_platform == 'GPL6102') | (microarray_platform == 'GPL6480') | (microarray_platform == 'GPL6883')) {

        gene_symbol_labels <- fData(gset)$Gene.symbol #get the gene symbols mapped to each row/probe

        i <- 1 #indexes rows/probes
        for (gene_symbol_i in gene_symbol_labels) {

            #get info for row
            probe_i_ID = probe_IDs[i]
            probe_i_sum = probe_sums[i]
            probe_i_gene_symbols = vector()

            #multiple genes
            if ((gene_symbol_i != '') & grepl('///', gene_symbol_i, fixed = TRUE)) {

                for (gene_symbol_j in unlist(strsplit(gene_symbol_i, '///', fixed = TRUE))) {
                    probe_i_gene_symbols <- c(probe_i_gene_symbols, trimws(gene_symbol_j))
                }
            #one gene
            } else if (gene_symbol_i != '') {

                probe_i_gene_symbols <- c(probe_i_gene_symbols, trimws(gene_symbol_i))
            }
            probe_info[[probe_i_ID]] <- list(probe_i_gene_symbols, probe_i_sum)
            i = i + 1
        }
    } else if (microarray_platform == 'GPL5175') {

        gene_symbol_labels <- fData(gset)$gene_assignment #get the gene symbols mapped to each row/probe

        i <- 1 #indexes rows/probes
        for (gene_assignment_i in gene_symbol_labels) {

            #get info for row
            probe_i_ID = probe_IDs[i]
            probe_i_sum = probe_sums[i]
            probe_i_gene_symbols = vector()

            #multiple genes
            if ((gene_assignment_i != '') & grepl('///', gene_assignment_i, fixed = TRUE)) {

                for (gene_assignment_j in unlist(strsplit(gene_assignment_i , '///' , fixed = TRUE))) {
                    gene_symbol_i <- unlist(strsplit(gene_assignment_j, '//'))[2]
                    probe_i_gene_symbols <- c(probe_i_gene_symbols, trimws(gene_symbol_i))
                }
            #one gene
            } else if ((gene_assignment_i != '') & grepl('//', gene_assignment_i, fixed = TRUE)) {

                gene_symbol_i <- unlist(strsplit(gene_assignment_i, '//'))[2]
                probe_i_gene_symbols <- c(probe_i_gene_symbols, trimws(gene_symbol_i))
            }   
            probe_info[[probe_i_ID]] <- list(probe_i_gene_symbols, probe_i_sum)
            i = i + 1
        }
    } else if (microarray_platform == 'GPL16951') {

        gene_symbol_labels <- fData(gset)$Gene_symbol #get the gene symbols mapped to each row/probe

        i <- 1 #indexes rows/probes
        for (gene_symbol_i in gene_symbol_labels) {
            
            if (is.na(gene_symbol_i)) { gene_symbol_i = '' }

            #get info for row
            probe_i_ID = probe_IDs[i]
            probe_i_sum = probe_sums[i]
            probe_i_gene_symbols = vector()

            #multiple genes
            if ((gene_symbol_i != 'previous version conserved probe') & (gene_symbol_i != '') & grepl('|', gene_symbol_i, fixed = TRUE)) {

                for (gene_symbol_j in unlist(strsplit(gene_symbol_i, '|', fixed = TRUE))) {
                    probe_i_gene_symbols <- c(probe_i_gene_symbols, trimws(gene_symbol_j))
                }
            #one gene
            } else if ((gene_symbol_i != 'previous version conserved probe') & (gene_symbol_i != '')) {

                probe_i_gene_symbols <- c(probe_i_gene_symbols, trimws(gene_symbol_i))
            }   
            probe_info[[probe_i_ID]] <- list(probe_i_gene_symbols, probe_i_sum)
            i = i + 1
        }
    }

    #iterate through each probe - gene symbol list - gene exp sum, keep track of the probe with the largest sum (across samples) for each gene symbol
    #> each gene symbol will get mapped to only the row/probe with the highest expression sum across samples
    probe_exp_sum_by_gene = list()
    for (probe_i in names(probe_info)) {

        gene_symbol_vec_probe_i = probe_info[probe_i][[1]][[1]] #holds the vector of gene symbols
        gene_symbol_exp_sum = unname(probe_info[probe_i][[1]][[2]]) #holds the probe exp sum (total across samples)

        #iterate through gene symbols that correspond to this probe (if there are any)
        if (length(gene_symbol_vec_probe_i) > 0) { #there's at least 1 gene symbol
            for (gene_symbol_i in gene_symbol_vec_probe_i) { #iterate through gene symbols

                #if gene symbol not present in list, then store gene symbol along with (probe exp sum) & (probe ID)
                if (gene_symbol_i %notin% names(probe_exp_sum_by_gene)) {
                    probe_exp_sum_by_gene[[gene_symbol_i]] <- list(probe_i, gene_symbol_exp_sum)

                #if gene symbol is present in list, then replace gene symbol along with (probe exp sum) & (probe ID) IF (probe exp sum) > than the (probe exp sum) currently in list
                } else if ((gene_symbol_i %in% names(probe_exp_sum_by_gene)) & (gene_symbol_exp_sum > probe_exp_sum_by_gene[[gene_symbol_i]][[2]])) {
                    probe_exp_sum_by_gene[[gene_symbol_i]] <- list(probe_i, gene_symbol_exp_sum)
                }
            }
        }
    }

    #create a vector of probe IDs we want to keep
    probe_IDs_to_keep <- vector()
    for (gene_symbol_i in names(probe_exp_sum_by_gene)) {
        probe_IDs_to_keep <- c(probe_IDs_to_keep, probe_exp_sum_by_gene[[gene_symbol_i]][[1]])
    }
    #some probe IDs may be duplicated since some probes map to multiple genes
    probe_IDs_to_keep <- unique(probe_IDs_to_keep)

    #create a boolean vector to filter Gene Expression matrix (TRUE = keep row/probe , FALSE = drop probe)
    exp_probe_bool_filter = vector()
    for (probe_id in row.names(exp_matrix)) {
        if (probe_id %in% probe_IDs_to_keep) {
            exp_probe_bool_filter <- c(exp_probe_bool_filter, TRUE)
        } else {
            exp_probe_bool_filter <- c(exp_probe_bool_filter, FALSE)
        }
    }

    #subset to rows that have the highest expression sums for each gene
    sel <- which(exp_probe_bool_filter == TRUE)
    gset <- gset[sel, ]
    return(gset)
}

microarray_diff_exp_analysis_for_samples <- function(GSE_ID, microarray_platform, gsms, normalized_expression_path, differential_expression_path, qnorm) {
    
    # Differential expression analysis with limma

    # load series and platform data from GEO
    gset <- getGEO(GSE_ID, GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep(microarray_platform, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # make proper column names to match toptable 
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    # Reduce rows to one probe per gene (keep probes with the highest expression sum across samples)
    gset <- reduce_expression_matrix_rows(gset, microarray_platform)
    
    # group names for all samples
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # eliminate samples marked as "X"
    sel <- which(sml != "X")
    sml <- sml[sel]
    gset <- gset[ ,sel]

    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
              (qx[6]-qx[1] > 50 && qx[2] > 0) ||
              (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
      exprs(gset) <- log2(ex) }

    # perform quantile normalization
    if (qnorm == TRUE) {
        col_names <- colnames(exprs(gset))
        row_names <- rownames(exprs(gset))
        exprs(gset) <- preprocessCore::normalize.quantiles(exprs(gset))
        colnames(exprs(gset)) <- col_names
        rownames(exprs(gset)) <- row_names }
    
    ############################################################
    #store normalized gene expression data
    ############################################################
    norm_exp_matrix <- exprs(gset) 
    
    #get the gene symbols mapped to each row/probe
    if ((microarray_platform == 'GPL6947') | (microarray_platform == 'GPL4133') | (microarray_platform == 'GPL10558') | (microarray_platform == 'GPL570') | (microarray_platform == 'GPL11532') | (microarray_platform == 'GPL6102') | (microarray_platform == 'GPL6480') | (microarray_platform == 'GPL6883')) {
        gene_symbol_labels <- fData(gset)$Gene.symbol
    } else if (microarray_platform == 'GPL5175') {
        gene_symbol_labels <- fData(gset)$gene_assignment
    } else if (microarray_platform == 'GPL16951') {
        gene_symbol_labels <- fData(gset)$Gene_symbol
    }
    
    #label each sample according to the group
    group_type_labels <- list()
    for (sample_i_group in sml) {
        if (sample_i_group == '0') {
            group_type_labels <- append(group_type_labels, group_0) 

        } else if (sample_i_group == '1') {
            group_type_labels <- append(group_type_labels, group_1)
        }   
    }

    #re-label rows according gene symbols & columns according to sample group
    rownames(norm_exp_matrix) <- gene_symbol_labels
    colnames(norm_exp_matrix) <- group_type_labels

    #transpose rows & columns
    norm_exp_matrix <- t(norm_exp_matrix)

    write.csv(norm_exp_matrix, normalized_expression_path)
    ############################################################
    
    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort="none" , n=Inf)

    write.csv(tT, differential_expression_path)
}

rnaseq_diff_exp_analysis_for_samples <- function(GSE_ID, gsms, counts_dataset_path, num_samples, TMM_norm, first_sample_column_i, last_sample_column_i, gene_symbol_col) {
    
    # Differential expression analysis with limma

    #read in counts matrix
    countsMatrix_filepath = file.path(counts_dataset_path)
    if (endsWith(countsMatrix_filepath, '.csv') | endsWith(countsMatrix_filepath, '.csv.gz')) {
        countsMatrix_full <- read.csv(file = countsMatrix_filepath, header=TRUE) #retrieve full counts matrix
    } else {
        countsMatrix_full <- read.table(file = countsMatrix_filepath, header=TRUE) #retrieve full counts matrix
    }

    annotation <- countsMatrix_full[[gene_symbol_col]]
    countsMatrix <- countsMatrix_full[first_sample_column_i:last_sample_column_i] #select columns that correspond to samples

    #convert counts matrix to a dataframe
    countsMatrix <- as.data.frame(countsMatrix)

    #construct design matrix
    group <- vector()
    for (sample_i in strsplit(gsms, "")[[1]]) {
        if (sample_i == '0') {
            group <- c(group, 1)
        }
        else if (sample_i == '1') {
            group <- c(group, 2)
        }
    }
    group <- factor(group)
    design <- model.matrix(~ group + 0)

    # group names for all samples
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # eliminate samples marked as "X"
    sel <- which(sml != "X")
    sml <- sml[sel]
    countsMatrix <- countsMatrix[ ,sel]

    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names

    # rename columns of design matrix to match group names
    colnames(design) <- levels(as.factor(sml)) 

    # Put the counts, groups, and annotation into a 'DGEList' object 
    x <- DGEList(counts=countsMatrix, group=sml, genes=annotation)

    # TMM Normalization (use if counts matrix not already normalized)
    if (TMM_norm == TRUE) {
        x <- calcNormFactors(x)
    }

    # Filter out genes that don't have at least 1 count-per-million in at least number of samples / 2
    isexpr <- rowSums(cpm(x) > 1) >= num_samples / 2
    x <- x[isexpr,]
    
    ############################################################
    #store normalized gene expression data
    ############################################################
    norm_exp_matrix <- x[['counts']]

    #get the gene symbols mapped to each row/probe
    gene_symbol_labels <- x[['genes']]
    gene_symbol_labels <- gene_symbol_labels[,"genes"]

    #label each sample according to the group
    group_type_labels <- list()
    for (sample_i_group in x$samples$group) {
        if (sample_i_group == 'G0') {
            group_type_labels <- append(group_type_labels, group_0) 

        } else if (sample_i_group == 'G1') {
            group_type_labels <- append(group_type_labels, group_1)
        }   
    }

    #re-label rows according gene symbols & columns according to sample group
    rownames(norm_exp_matrix) <- gene_symbol_labels
    colnames(norm_exp_matrix) <- group_type_labels

    #transpose rows & columns
    norm_exp_matrix <- t(norm_exp_matrix)

    write.csv(norm_exp_matrix, normalized_expression_path)
    ############################################################

    # Run voom
    voomResults <- voom(x,design)

    #Run limma
    fit <- lmFit(voomResults, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort="none" , n=Inf)

    #save results
    write.csv(tT, differential_expression_path)
}

normalized_expression_dir = file.path(
    data_dir, 'normalized-expression-data', groups)
dir.create(normalized_expression_dir, recursive=TRUE)

normalized_expression_path = file.path(
    normalized_expression_dir,
    paste0(GSE_ID, '_', platform_type, '_Exp_EachGene.csv'))

differential_expression_dir = file.path(
    data_dir, paste0(platform_type, '-differential-gene-expression'), groups)
dir.create(differential_expression_dir, recursive=TRUE)

differential_expression_path = file.path(
    differential_expression_dir, paste0(GSE_ID, '.csv'))

if (platform_type == 'microarray') {

    microarray_diff_exp_analysis_for_samples(
        GSE_ID, microarray_platform, gsms,
        normalized_expression_path,
        differential_expression_path, qnorm=TRUE)
    
} else if (platform_type == 'rnaseq') {

    rnaseq_diff_exp_analysis_for_samples(
        GSE_ID, gsms, counts_dataset_path, num_samples, TMM_norm,
        first_sample_column_i, last_sample_column_i, gene_symbol_col)
    
}

cat(
    normalized_expression_path,
    differential_expression_path
)
