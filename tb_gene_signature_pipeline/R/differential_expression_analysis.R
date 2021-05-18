#!/shared/software/anaconda3/envs/tb-gene-signature/bin/Rscript

quietly = function(x) {
    return(suppressWarnings(suppressMessages(x)))
}

quietly(library(GEOquery))
quietly(library(rjson))
quietly(library(dplyr))
quietly(library(tidyr))
quietly(library(tibble))
quietly(library(stringr))
quietly(library(limma))
quietly(library(here))
quietly(library(edgeR))
quietly(library(yaml))

args = commandArgs(trailingOnly=TRUE)

if (args[1] == 'true') {
    log_transform = TRUE
} else if (args[1] == 'false') {
    log_transform = FALSE
}

config = read_yaml(file=here('config.yml'))
datasets = fromJSON(file=here('data', 'datasets.json'))
comparisons = fromJSON(file=here('data', 'comparisons.json'))
platforms = fromJSON(file=here('data', 'platforms.json'))
        
gene_id_symbol_map = read.table(
    here('data', 'gene_id_symbol_map.tsv'), sep='\t', header=TRUE)
srr_gsm_id_map = read.table(
    here('data', 'srr_gsm_id_map.tsv'), sep='\t', header=TRUE) %>% 
    column_to_rownames('srr_id')

data_dir = config$data_directory
rnaseq_dir = file.path(data_dir, 'rnaseq-datasets')

is_log_normalized = function(gset) {

    exprs = exprs(gset)
    
    # non-missing, variant columns
    cols = apply(exprs, 2, function(x) sd(x, na.rm=TRUE) > 0)
    
    # test using first column
    test_vec = exprs[,cols][,1]
                 
    # add offset to make all values positive, if necessary
    min_value = min(test_vec, na.rm=TRUE)
    if (min_value <= 0) {
        test_vec = test_vec - min_value + 1
    }
                 
    real = qqnorm(test_vec, plot.it=FALSE)
    exp = qqnorm(exp(test_vec), plot.it=FALSE)
    log = qqnorm(log(test_vec), plot.it=FALSE)
                 
    exp$x[which(is.infinite(exp$x))] = NA
    exp$y[which(is.infinite(exp$y))] = NA

    cor_real = cor(real$x, real$y, use='pairwise.complete.obs')
    cor_exp = cor(exp$x, exp$y, use='pairwise.complete.obs')
    cor_log = cor(log$x, log$y, use='pairwise.complete.obs')  
    
    if (all(is.na(cor_exp))) {
        score = cor_log - cor_real
    } else {
        score = log2(abs(cor_real**2 - cor_log**2)/abs(cor_real**2 - cor_exp**2))
    }
        
    return(score <= 0)

}
                 
quantile_normalize = function(p, exprs, sample_phenotypes) {
    samples = names(Filter(function(x) { x == p }, sample_phenotypes))
    exprs_norm = preprocessCore::normalize.quantiles(
        exprs[,samples]) %>% data.frame
    rownames(exprs_norm) = rownames(exprs)
    colnames(exprs_norm) = samples
    return(exprs_norm)
}
                 
code_phenotype = function(gsm_id, phenotypes, control, case) {
    if (phenotypes[[gsm_id]] == control) { return(0) }
    else if (phenotypes[[gsm_id]] == case) { return(1) }
    else { return(NA) }
}
                 
normalize_expression_matrix = function(exprs, dataset) {
    
    sample_phenotypes = as.factor(unlist(dataset$phenotypes))
    phenotype_levels = levels(sample_phenotypes)

    if (dataset$platform == 'rnaseq') {
        design_matrix = model.matrix(~0 + sample_phenotypes)
        rownames(design_matrix) = names(sample_phenotypes)
        colnames(design_matrix) = phenotype_levels
        
        dge_list = DGEList(
            counts=exprs, group=sample_phenotypes, genes=rownames(exprs))
        norm_factors = calcNormFactors(dge_list)
        exprs = voom(dge_list, design_matrix)$E
    }
    
    else if (startsWith(dataset$platform, 'GPL')) {
        # quantile normalize within phenotype groups
        norm_dfs = lapply(phenotype_levels, quantile_normalize, exprs, sample_phenotypes)
        exprs = as.matrix(bind_cols(norm_dfs))
        #genes = rownames(exprs)
        #exprs = preprocessCore::normalize.quantiles(exprs)
        #rownames(exprs) = genes
        #colnames(exprs) = names(sample_phenotypes)
    }

    return (exprs)
}
                 
fit_differential_expression_model = function(comparison, exprs, dataset, data_dir) {

    coded_phenotypes = sapply(
        colnames(exprs), code_phenotype, dataset$phenotypes,
        comparison$control, comparison$case)

    coded_phenotypes = as.factor(coded_phenotypes[!is.na(coded_phenotypes)])

    if (length(unique(coded_phenotypes)) == 2) { 
 
        design_matrix = model.matrix(~0 + coded_phenotypes)
        rownames(design_matrix) = names(coded_phenotypes)[!is.na(coded_phenotypes)]
        colnames(design_matrix) = c('control', 'case')
        exprs = exprs[,rownames(design_matrix)]
        
        if (dataset$platform == 'rnaseq') {
            
            # normalize data
            dge_list = DGEList(
                counts=exprs, group=coded_phenotypes, genes=rownames(exprs))
            norm_factors = calcNormFactors(dge_list)
            exprs = voom(dge_list, design_matrix)
           
        } 

        limma_fit = lmFit(exprs, design_matrix)
        contrast_matrix = makeContrasts(case - control, levels=design_matrix)
        contrast_fit = contrasts.fit(limma_fit, contrast_matrix)
        contrast_fit_ebayes = eBayes(contrast_fit, proportion=0.01)
        
        results_df = topTable(
            contrast_fit_ebayes, genelist=rownames(exprs),
            adjust='fdr', sort='none', n=Inf)

        results_df = results_df %>%
            rename(gene_symbol='ID',
                   log_fc=logFC,
                   adj_p_val=adj.P.Val) %>%
            mutate(dataset=dataset$gse_id,
                   control=comparison$control,
                   case=comparison$case) %>%
            select(dataset,
                   control,
                   case,
                   gene_symbol,
                   log_fc,
                   adj_p_val) 
        
        # write out results to disk
        results_path = file.path(
            data_dir, 'differential-expression-results',
            paste(dataset$gse_id, paste0(comparison$control, '_vs_', comparison$case),
                  'diff_expr_results', 'tsv', sep='.'))
        write.table(results_df, results_path, sep='\t', row.names=FALSE, quote=FALSE)
        
        return(results_df)
    }
}
                 
differential_expression_analysis = function(dataset, comparisons, platforms, data_dir, qnorm=TRUE) {
    
    # process microarray datasets
    if (startsWith(dataset$platform, 'GPL')) {

        # fetch dataset from GEO
        gset = quietly(getGEO(dataset$gse_id, AnnotGPL=TRUE))
        if (length(gset) > 1) {
            idx = grep(dataset$platform, attr(gset, 'names'))
        } else {
            idx = 1
        }
        gset = gset[[idx]]
        
        # log normalize expression matrix values if not already
        if (log_transform || !is_log_normalized(gset)) {
            
            # values may still be negative due to background subtraction algo, add offset
            min_val = min(exprs(gset), na.rm=TRUE)
            if (min_val < 0) {
                exprs(gset) = exprs(gset) - min_val + 1
            }
            
            # log2 normalize
            exprs(gset) = log2(exprs(gset))
        }
        
        # format feature annotation names
        fvarLabels(gset) = make.names(fvarLabels(gset))
        
        # name of gene symbol column varies by microarray platform
        gene_col_name = platforms[[dataset$platform]]$gene_col_name

        # create expression dataframe
        exprs_df = fData(gset) %>%
            rename(gene_symbol = all_of(gene_col_name)) %>%
            select(gene_symbol) %>%
            rownames_to_column('probe_id') %>%
        
            # join gene symbol to expression matrix
            inner_join(
                exprs(gset) %>%
                    data.frame %>%
                    rownames_to_column('probe_id'),
                by='probe_id') %>%
            select(-probe_id) %>%
        
            # compute total expression per probe ID across samples
            mutate(probe_sum = rowSums(across(where(is.numeric)))) %>%
        
            # explode rows with multiple gene symbols
            separate_rows(gene_symbol, sep='///') %>%
            mutate(gene_symbol = str_trim(gene_symbol)) %>%
        
            # drop missing gene symbols
            mutate(gene_symbol = na_if(gene_symbol, '')) %>%
            mutate(gene_symbol = na_if(gene_symbol, '---')) %>%
            mutate(gene_symbol = na_if(gene_symbol, 'previous version conserved probe')) %>%
            drop_na(gene_symbol) %>%
        
            # parse out gene symbol
            extract(gene_symbol, into=c('gene_symbol', 'gene_symbol1'),
                    regex='^[\\s]*([[:graph:]]+)[\\s//\\s]*([[:graph:]]+)*') %>%
            mutate(gene_symbol1 = na_if(gene_symbol1, '')) %>%
            mutate(gene_symbol = ifelse(is.na(gene_symbol1), gene_symbol, gene_symbol1)) %>%
            select(-gene_symbol1) %>%
        
            # explode rows with multiple gene symbols (again)
            separate_rows(gene_symbol, sep='\\|') %>%
            mutate(gene_symbol = str_trim(gene_symbol)) %>%
        
            # keep one row per gene symbol (row with highest expression across samples)
            group_by(gene_symbol) %>%
            slice_max(n=1, order_by=probe_sum, with_ties=FALSE) %>%
            ungroup() %>%
            select(-probe_sum) %>%
            column_to_rownames('gene_symbol')
        
        # convert to matrix
        exprs = as.matrix(exprs_df) 
        
        # remove null phenotype samples
        is_null = sapply(dataset$phenotypes, is.null)
        keep_samples = intersect(
            names(is_null)[!is_null], colnames(exprs))
        dataset$phenotypes = dataset$phenotypes[keep_samples]
        exprs = exprs[,keep_samples]
        
        # set of unique phenotypes in the dataset
        #phenotypes = unique(dataset$phenotypes)
        
        # quantile normalize within phenotype groups
        #norm_dfs = lapply(phenotypes, quantile_normalize, exprs_df, dataset$phenotypes)
        #exprs = as.matrix(bind_cols(norm_dfs))
    }
    
    # process rnaseq datasets
    else if (dataset$platform == 'rnaseq') {
        
        filename = paste(
            dataset$gse_id, 'salmon', 'genes', 'read_count',
            'tximport', 'tsv', sep='.')
        filepath = file.path(rnaseq_dir, filename)
        
        df_exprs = read.table(filepath, sep='\t', header=TRUE)
        df_exprs = inner_join(gene_id_symbol_map, df_exprs, by='gene_id') %>%
            select(-gene_id)
        colnames(df_exprs)[2:ncol(df_exprs)] = sapply(
            colnames(df_exprs)[2:ncol(df_exprs)], function(x) srr_gsm_id_map[[x, 'gsm_id']])
        
        # for each gene symbol, keep only row with highest total expression across samples
        df_exprs = df_exprs %>%
            group_by(gene_symbol) %>%
            arrange(desc(rowSums(across(where(is.numeric))))) %>%
            filter(row_number() == 1) %>%
            column_to_rownames('gene_symbol')
            
        # convert to matrix
        exprs = as.matrix(df_exprs) 
            
        # remove null phenotype samples
        is_null = sapply(dataset$phenotypes, is.null)
        keep_samples = intersect(
            names(is_null)[!is_null], colnames(exprs))
        dataset$phenotypes = dataset$phenotypes[keep_samples]
        exprs = exprs[,keep_samples]
            
        # filter out genes that don't have at least 1 count-per-million in at least half the samples
        exprs = exprs[rowSums(cpm(exprs)) >= (ncol(exprs)/2),]
    }
    
    # write out preprocessed expression matrices
    exprs_out = exprs %>% data.frame %>% rownames_to_column('gene_symbol')
    exprs_path = file.path(
        data_dir, 'transformed-expression-matrices',
        paste(dataset$gse_id, dataset$platform, 'transformed_expr_matrix', 'tsv', sep='.'))
    write.table(exprs_out, exprs_path, sep='\t', row.names=FALSE, quote=FALSE)

    # normalize expression data
    exprs_normalized = normalize_expression_matrix(exprs, dataset)
            
    # write out normalized expression matrices
    exprs_norm_out = exprs_normalized %>% data.frame %>% rownames_to_column('gene_symbol')
    exprs_norm_path = file.path(
        data_dir, 'normalized-expression-matrices',
        paste(dataset$gse_id, dataset$platform, 'normalized_expr_matrix', 'tsv', sep='.'))
    write.table(exprs_norm_out, exprs_norm_path, sep='\t', row.names=FALSE, quote=FALSE)
        
    gse_id = dataset$gse_id
    platform = dataset$platform
            
    # combine preprocessed and normalized values in long dataframe to return
    expr_return = inner_join(
        exprs_out %>%
            pivot_longer(cols=!gene_symbol, names_to='gsm_id', values_to='transformed'),   
        exprs_norm_out %>%
            pivot_longer(cols=!gene_symbol, names_to='gsm_id', values_to='normalized'),
        by=c('gene_symbol', 'gsm_id')) %>%
            
        # add sample phenotypes as column
        inner_join(
            data.frame(phenotype=unlist(dataset$phenotypes)) %>%
                rownames_to_column('gsm_id'),
            by='gsm_id') %>%
            
        # order samples by phenotype level
        arrange(phenotype) %>%
        mutate(gsm_id=factor(gsm_id, levels=unique(gsm_id))) %>%
            
        # add dataset and platform as columns
        mutate(dataset=gse_id,
               platform=platform)
            
    # run differential expression analysis
    dataset_results_list = lapply(
        comparisons, fit_differential_expression_model,
        exprs, dataset, file.path(data_dir))
    dataset_results_df = bind_rows(dataset_results_list)
    rownames(dataset_results_df) = NULL
            
    return_obj = list(
        expr_df=expr_return,
        results_df=dataset_results_df)

    return(return_obj)
}

all_datasets_output = lapply(
    datasets, differential_expression_analysis,
    comparisons, platforms, data_dir)
            
all_exprs_df = bind_rows(
    lapply(all_datasets_output, function(x) x$expr_df)) %>%
    select(dataset, platform, gene_symbol, gsm_id, phenotype, transformed, normalized) %>%
    arrange(dataset, gene_symbol, gsm_id)
            
all_results_df = bind_rows(
    lapply(all_datasets_output, function(x) x$results_df)) %>%
    select(dataset, control, case, gene_symbol, log_fc, adj_p_val) %>%
    arrange(dataset, control, case, gene_symbol)

exprs_file = file.path(data_dir, 'differential_expression_values.tsv')
results_file = file.path(data_dir, 'differential_expression_results.tsv')
           
write.table(all_exprs_df, exprs_file, sep='\t', row.names=FALSE, quote=FALSE)       
write.table(all_results_df, results_file, sep='\t', row.names=FALSE, quote=FALSE)
