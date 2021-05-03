
import os
import pickle
import subprocess
from subprocess import check_output, CalledProcessError
import numpy as np
import pandas as pd
#from collections import Counter
#import fastcluster
#from scipy.cluster.hierarchy import dendrogram
#import scipy.cluster.hierarchy as hier
import networkx as nx
#import networkx.algorithms.community as nxcom
#import community
#from scipy.cluster.hierarchy import dendrogram
#import scipy.cluster.hierarchy as hier
#import gseapy as gp
import numpy as np
from pathlib import Path
import yaml
import json

with (Path.cwd().resolve().parent / 'config.yml').open('r') as f:
    config = yaml.safe_load(f)
    
with (Path.cwd().resolve().parent / 'data' / 'datasets.json').open('r') as f:
    datasets = json.load(f)
    
with (Path.cwd().resolve().parent / 'data' / 'platforms.json').open('r') as f:
    platforms = json.load(f)
    
with (Path.cwd().resolve().parent / 'data' / 'comparisons.json').open('r') as f:
    comparisons = json.load(f)
    
this_dir = Path(__file__).resolve().parent

data_dir = Path(config['data_directory'])
data_dir.mkdir(parents=True, exist_ok=True)

rnaseq_dataset_dir = Path(config['rnaseq_dataset_directory'])

# notebook 1
def run_differential_expression_analysis():
    """
    Run marginal, per-gene differential expression analyses.
    
    :func:`run_differential_expression_analysis` runs the differential analysis script
    `gene-signature-pipeline-ml-pipeline/R/differential_expression_analysis.R`.
    
    An independent case/control differential expression analysis is run for
    each dataset defined in `data/datasets.json`, and for each group comparison
    defined in `data/comparisons.json`.
    
    Returns
    -------
    :class:`pd.DataFrame`
        A dataframe with the schema
    
    Examples
    --------
    
    """
    r_script = this_dir / 'R' / 'differential_expression_analysis.R'

    results_file = Path(
        check_output(r_script).decode('utf8'))

    df = pd.read_table(results_file, sep='\t')
    results_file.unlink()
    
    return df

# notebook 2
def merge_differential_expression_results(
    differential_expression_df, pval_thresh=0.05,
        log_fc_thresh=np.log2(1.5)):
    """
    """
    df = differential_expression_df.copy()
    
    df = df.loc[
        (df['log_fc'] >= log_fc_thresh) & (df['adj_p_val'] <= pval_thresh)]
    
    df = (df
        .groupby(['control', 'case'])
        .apply(
            lambda x: x.pivot(
                index='gene_symbol', columns='dataset', values='log_fc'))
        .fillna(0.)
        .reset_index())
    
    return df

# notebook 3
def run_edge_weight_distribution_analysis(data_dir, group, qval_thresh):
    """
    """
    diff_exp_data_dir = os.path.join(data_dir, 'pooled-differential-gene-expression')
    
    def edge_weight_distribution_from_pooled_logFC(logFC_table_filename):

        # read in CSV file with significant logFC changes for genes between comparison X samples collected from various datasets
        #diff_gene_exp_df = pd.read_csv(diff_exp_data_dir + logFC_table_filename).set_index('Unnamed: 0')
        diff_gene_exp_df = pd.read.csv(os.path.join(diff_exp_data_dir, logFC_table_filename)).set_index('Unnamed: 0')
        diff_gene_exp_df.rename_axis('gene' , inplace = True)

        # Construct simplified matrix of logFC direction from DataFrame with significant logFC changes across all analyses by converting values:
        # +1 if logFC > 0
        # 0 if logFC = 0
        # -1 if logFC < 0

        # store copy of array from dataframe with sig. logFC values (rows = genes, columns = GSE ID)
        direc_diff_gene_exp_matrix = diff_gene_exp_df.copy().values 

        # replace values in logFC matrix
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix > 0.0] = 1
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix == 0.0] = 0
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix < 0.0] = -1

        # convert to lower memory int8 datatype
        direc_diff_gene_exp_matrix = direc_diff_gene_exp_matrix.astype('int8')

        # compute the dot product between every two pairs of gene vectors (will calculate the edges weights for our network)
        # multiply direction logFC matrix by its transpose to get the dot products between all pairs of rows
        network_edge_weight_matrix = direc_diff_gene_exp_matrix.dot(direc_diff_gene_exp_matrix.T)

        # the row/column annotation (genes) can be copied from the logFC differential gene expression DataFrame
        network_edge_weight_matrix_labels = pd.Series(list(diff_gene_exp_df.index) , index = range(0 , len(diff_gene_exp_df.index)))

        #number of rows / columns
        num_genes = np.shape(network_edge_weight_matrix)[0]

        # retrieve the distribution of the Edge Weights by returning the upper triangular part of the matrix
        edge_weight_array = network_edge_weight_matrix[np.triu_indices(num_genes, k = 0)]

        #convert array to a Counter dict to save space (keys: edge weight values, values: count of edge weights in edge weight distribution)
        edge_weight_distr_counter_dict = Counter(list(edge_weight_array))

        return edge_weight_distr_counter_dict

    qval_thresh = 0.001
    edge_weight_distr_counter_dict = edge_weight_distribution_from_pooled_logFC(
        f'{group}_signif_logFC_qval_thresh_{qval_thresh}.csv')
    #plot_distribution_of_edge_weights(axes[0,0], edge_weight_distr_counter_dict, f'q-val threshold = {qval_thresh}')

    return

# notebook 4
def construct_edge_weight_null_distribution(data_dir, group, qval_thresh, n_iter=25):
    """
    """
    diff_exp_data_dir = os.path.join(data_dir, 'pooled-differential-gene-expression')
    
    def edge_weight_distribution_from_pooled_logFC(logFC_table_filename):

        # read in CSV file with significant logFC changes for genes between comparison X samples collected from various datasets
        diff_gene_exp_df = pd.read_csv(
            os.path.join(diff_exp_data_dir, logFC_table_filename)).set_index('Unnamed: 0')
        diff_gene_exp_df.rename_axis('gene' , inplace = True)

        # Construct simplified matrix of logFC direction from DataFrame with significant logFC changes across all analyses by converting values:
        # +1 if logFC > 0
        #  0 if logFC = 0
        # -1 if logFC < 0

        # store copy of array from dataframe with sig. logFC values (rows = genes, columns = GSE ID)
        direc_diff_gene_exp_matrix = diff_gene_exp_df.copy().values 

        # replace values in logFC matrix
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix > 0.0] = 1
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix == 0.0] = 0
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix < 0.0] = -1

        # convert to lower memory int8 datatype
        direc_diff_gene_exp_matrix = direc_diff_gene_exp_matrix.astype('int8')

        # compute the dot product between every two pairs of gene vectors (will calculate the edges weights for our network)
        # multiply direction logFC matrix by its transpose to get the dot products between all pairs of rows
        network_edge_weight_matrix = direc_diff_gene_exp_matrix.dot(direc_diff_gene_exp_matrix.T)

        # the row/column annotation (genes) can be copied from the logFC differential gene expression DataFrame
        network_edge_weight_matrix_labels = pd.Series(list(diff_gene_exp_df.index) , index = range(0 , len(diff_gene_exp_df.index)))

        #number of rows / columns
        num_genes = np.shape(network_edge_weight_matrix)[0]

        # retrieve the distribution of the Edge Weights by returning the upper triangular part of the matrix
        edge_weight_array = network_edge_weight_matrix[np.triu_indices(num_genes, k = 0)]

        #convert array to a Counter dict to save space (keys: edge weight values, values: count of edge weights in edge weight distribution)
        edge_weight_distr_counter_dict = Counter(list(edge_weight_array))

        return edge_weight_distr_counter_dict
    
    def shuffle_rows_within_each_column(array):
        """
        This function takes in a 2D-array and shuffles the rows of the array seperately for each column.
        """
        #get number of rows & columns for input array
        nrows, ncols = array.shape

        #get the column indices for each element (will keep the same)
        cols = np.indices((nrows, ncols))[1]

        #permute the row indices for each column
        rows = np.array([np.random.permutation(nrows) for _ in range(ncols)]).T

        #re-arrange elements in each column according to the chosen row indices for that column
        shuffled_array = array[rows, cols]

        return shuffled_array
    
    def edge_weight_null_distribution_from_pooled_logFC(logFC_table_filename, N):
    
        # read in CSV file with significant logFC changes for genes between comparison X samples collected from various datasets
        diff_gene_exp_df = pd.read_csv(
            os.path.join(diff_exp_data_dir, logFC_table_filename)).set_index('Unnamed: 0')
        diff_gene_exp_df.rename_axis('gene' , inplace = True)

        # Construct simplified matrix of logFC direction from DataFrame with significant logFC changes across all analyses by converting values:
        # +1 if logFC > 0
        # 0 if logFC = 0
        # -1 if logFC < 0

        # store copy of array from dataframe with sig. logFC values (rows = genes, columns = GSE ID)
        direc_diff_gene_exp_matrix = diff_gene_exp_df.copy().values 

        # replace values in logFC matrix
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix > 0.0] = 1
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix == 0.0] = 0
        direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix < 0.0] = -1

        # convert to lower memory int8 datatype
        direc_diff_gene_exp_matrix = direc_diff_gene_exp_matrix.astype('int8')

        # construct null distribution of shuffling the columns within each row, extracting edge weights, then adding to bucket of edge weights from many iterations
        edge_weight_null_distr_counter_dict = Counter([]) #initialize empty Counter dict

        # (shuffle columns within each row > get edge weights) x N > add to bucket of edge weights (null distribution)
        for iter_i in range(1, N+1):

            # permute the rows within each column once
            direc_diff_gene_exp_matrix_shuffled = shuffle_rows_within_each_column(direc_diff_gene_exp_matrix)

            # compute the dot product between every two pairs of gene vectors (will calculate the edges weights for our network)
            # multiply direction logFC matrix by its transpose to get the dot products between all pairs of rows
            network_edge_weight_matrix = direc_diff_gene_exp_matrix_shuffled.dot(direc_diff_gene_exp_matrix_shuffled.T)

            # the row/column annotation (genes) can be copied from the logFC differential gene expression DataFrame
            network_edge_weight_matrix_labels = pd.Series(list(diff_gene_exp_df.index) , index = range(0 , len(diff_gene_exp_df.index)))

            #number of rows / columns
            num_genes = np.shape(network_edge_weight_matrix)[0]

            # retrieve the distribution of the Edge Weights by returning the upper triangular part of the matrix
            edge_weight_array = network_edge_weight_matrix[np.triu_indices(num_genes, k = 0)]

            #convert array to a Counter dict to save space (keys: edge weight values, values: count of edge weights in edge weight distribution)
            edge_weight_distr_counter_dict = Counter(list(edge_weight_array))

            # append to counter with distribution of edge weights
            edge_weight_null_distr_counter_dict = edge_weight_null_distr_counter_dict + edge_weight_distr_counter_dict

            if iter_i % 1 == 0:
                print(f'finished loop {iter_i}')

        return edge_weight_null_distr_counter_dict
    
    null_distr_data_dir = os.path.join(
        data_dir, 'null-distributions', f'qval-thresh-{qval_thresh}')
    if not os.path.exists(null_distr_data_dir):
        os.makedirs(null_distr_data_dir)
        
    logFC_table_filename = f'{group}-signif-logFC-qval-thresh-{qval_thresh}.csv'

    #get the null distribution of edge weights (shuffle columns within each row, repeat process and add to bucket N times)
    edge_weight_null_distr_counter_dict = edge_weight_null_distribution_from_pooled_logFC(logFC_table_filename, n_iter)

    #output the null distribution as a pickled object for later analysis
    edge_weights_path = os.path.join(null_distr_data_dir, f'{group}.pickle')
    with open(edge_weights_path, 'wb') as handle:
        pickle.dump(
            edge_weight_null_distr_counter_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return edge_weights_path

# notebook 5
def construct_network(data_dir, group, qval_thresh):
    """
    """
    num_datasets_per_comparison = {
        'atb-ltbi': 16,
        'atb-hc' : 15,
        'atb-od': 10,
        'ltbi-hc': 9
    }
    
    diff_exp_data_dir = os.path.join(data_dir, 'pooled-differential-gene-expression')
    
    pickled_objects_dir = os.path.join(data_dir, 'pickled-files')
    if not os.path.exists(pickled_objects_dir):
        os.makedirs(pickled_objects_dir)

    #for degree series
    degree_series_path = os.path.join(
        pickled_objects_dir, 'network-files', 'degree-series')
    if not os.path.exists(degree_series_path):
        os.makedirs(degree_series_path)

    #for weighted degree series
    weighted_degree_path = os.path.join(
        pickled_objects_dir, 'network-files', 'weighted-degree-series')
    if not os.path.exists(weighted_degree_path):
        os.makedirs(weighted_degree_path)                

    #for eigenvector centrality series
    eigenvector_path = os.path.join(
        pickled_objects_dir, 'network-files', 'eigenvector-centrality-series')
    if not os.path.exists(eigenvector_path):
        os.makedirs(eigenvector_path)

    #for networks
    networks_path = os.path.join(
        pickled_objects_dir, 'network-files', 'networks')
    if not os.path.exists(networks_path):
        os.makedirs(networks_path)

    #for mean logFC series
    mean_logFC_path = os.path.join(
        pickled_objects_dir, 'network-files', 'mean-logFC-network-nodes-series')
    if not os.path.exists(mean_logFC_path):
        os.makedirs(mean_logFC_path)

    # read in CSV file with significant logFC changes for genes between comparison X samples collected from various datasets
    diff_gene_exp_df = pd.read_csv(os.path.join(
        diff_exp_data_dir, f'{group}-signif-logFC-qval-thresh-{qval_thresh}.csv')).set_index('Unnamed: 0')
    diff_gene_exp_df.rename_axis('gene' , inplace = True)

    # Construct simplified matrix of logFC direction from DataFrame with significant logFC changes across all analyses by converting values:
    # +1 if logFC > 0
    #  0 if logFC = 0
    # -1 if logFC < 0

    # store copy of array from dataframe with sig. logFC values (rows = genes, columns = GSE ID)
    direc_diff_gene_exp_matrix = diff_gene_exp_df.copy().values 

    # replace values in logFC matrix
    direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix > 0.0] = 1
    direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix == 0.0] = 0
    direc_diff_gene_exp_matrix[direc_diff_gene_exp_matrix < 0.0] = -1

    # convert to lower memory int8 datatype
    direc_diff_gene_exp_matrix = direc_diff_gene_exp_matrix.astype('int8')

    # compute the dot product between every two pairs of gene vectors (will calculate the edges weights for our network)
    # Compute $M \cdot M^T \text{ to get } i \cdot j \text{ , } \forall \text{ pairs of rows } i, j \text{ in matrix } M$
    # multiply direction logFC matrix by its transpose to get the dot products between all pairs of rows
    network_edge_weight_matrix = direc_diff_gene_exp_matrix.dot(direc_diff_gene_exp_matrix.T)

    # the row/column annotation (genes) can be copied from the logFC differential gene expression DataFrame
    network_edge_weight_matrix_labels = pd.Series(
        list(diff_gene_exp_df.index) , index = range(0 , len(diff_gene_exp_df.index)))

    # DISTRIBUTION OF EDGE WEIGHTS

    #number of rows / columns
    num_genes = np.shape(network_edge_weight_matrix)[0]

    # retrieve the distribution of the Edge Weights by returning the upper triangular part of the matrix
    edge_weight_array = network_edge_weight_matrix[np.triu_indices(num_genes, k = 0)]

    #convert array to a Counter dict to save space (keys: edge weight values, values: count of edge weights in edge weight distribution)
    edge_weight_distr_counter_dict = Counter(list(edge_weight_array))
    
    #Return the upper triangular part of the matrix with elements in lower part ZEROED out
    upper_tri_network_edge_weight_matrix = np.triu(network_edge_weight_matrix, k = 0)

    #Return a boolean for elements in the upper triangular part of the matrix for elments that are <= -3 OR >= 3
    upper_tri_network_edge_weight_matrix_bool = abs(upper_tri_network_edge_weight_matrix) >= 3

    #get the indices for the elements in the upper triangle where elements (edge weights) <= -3 OR >= 3
    node_i_indices = upper_tri_network_edge_weight_matrix_bool.nonzero()[0]
    node_j_indices = upper_tri_network_edge_weight_matrix_bool.nonzero()[1]

    #get the normalization factor (number of datasets used to construct edge weights = maximum possible weight)
    edge_weight_norm_factor = float(num_datasets_per_comparison[group])

    #Create list of edges for NetworkX graph by iterating through numpy (adjancency) matrix (with edge weights) + node labels (rows/columns of matrix) & storing edges with weights <= -3 OR >= 3
    G_edge_list = [ ( network_edge_weight_matrix_labels[node_i], network_edge_weight_matrix_labels[node_j], (float(abs(network_edge_weight_matrix[node_i, node_j]))/edge_weight_norm_factor) ) for node_i, node_j in zip(node_i_indices, node_j_indices)]
    
    G = nx.Graph()
    G.add_weighted_edges_from(G_edge_list)
    
    nx.write_gpickle(G, os.path.join(networks_path, f'{group}.pkl'))
                     
    node_list = [node_deg[0] for node_deg in list(G.degree())]
    degree_list = [node_deg[1] for node_deg in list(G.degree())]
    degree_series = pd.Series(degree_list , index = node_list).sort_values(ascending = False)
    
    node_list = [node_deg[0] for node_deg in list(G.degree(weight = 'weight'))]
    degree_list = [node_deg[1] for node_deg in list(G.degree(weight = 'weight'))]
    weighted_degree_series = pd.Series(degree_list , index = node_list).sort_values(ascending = False)
    
    eigenvector_centrality_series = pd.Series(
        nx.eigenvector_centrality(G, weight = 'weight', max_iter=10000)).sort_values(ascending = False) #takes edge weight into account
                     
    degree_series.to_pickle(
        os.path.join(degree_series_path, f'{group}.pkl'))
    weighted_degree_series.to_pickle(
        os.path.join(weighted_degree_path, f'{group}.pkl'))
    eigenvector_centrality_series.to_pickle(
        os.path.join(eigenvector_path, f'{group}.pkl'))
                     
    genes_in_network = list(G.nodes()) #get all nodes in the network
    mean_logFC_series = pd.Series(index = genes_in_network) #create series

    for gene_i in mean_logFC_series.index:

        #use code below to calculate average logFC for gene across ALL datasets/studies
        gene_i_mean_logFC = diff_gene_exp_df.loc[gene_i , :].mean()

        #append to series
        mean_logFC_series[gene_i] = gene_i_mean_logFC

    mean_logFC_series.sort_values(inplace = True, ascending = False)
                     
    mean_logFC_series.to_pickle(os.path.join(mean_logFC_path, f'{group}.pkl'))
                                
    return {
        'network-files': networks_path,
        'degree-series': degree_series_path,
        'weighted-degree-series': weighted_degree_path,
        'eigenvector-centrality-series': eigenvector_path,
        'mean-logFC-network_nodes-series': mean_logFC_path
    }

# notebook 6
def compare_across_networks(data_dir):
    """
    """
    #directoy for Pickled Objects
    pickled_objects_dir = os.path.join(data_dir, 'pickled-files')

    #specify directory for Gene Set Files
    Gene_Set_dir = os.path.join(data_dir, 'gene-set-files', 'csv-files')

    #specify directory for Data
    Data_files_dir = data_dir
    
    top_nodes_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists')                                   
    if not os.path.exists(top_nodes_path): 
        os.makedirs(top_nodes_path)
                                       
    G_ATB_HC = nx.read_gpickle(
        os.path.join(pickled_objects_dir, 'network-files', 'networks', 'atb-hc.pkl'))
    G_ATB_HC = nx.read_gpickle(
        os.path.join(pickled_objects_dir, 'network-files', 'networks', 'atb-ltbi.pkl'))
    G_ATB_HC = nx.read_gpickle(
        os.path.join(pickled_objects_dir, 'network-files', 'networks', 'atb-od.pkl'))
                                       
    weighted_deg_ATB_HC_series = pd.read_pickle(
        os.path.join(pickled_objects_dir, 'network-files', 'weighted-degree-series', 'atb-hc.pkl'))
    weighted_deg_ATB_LTBI_series = pd.read_pickle(
        os.path.join(pickled_objects_dir, 'network-files', 'weighted-degree-series', 'atb-ltbi.pkl'))
    weighted_deg_ATB_OD_series = pd.read_pickle(
        os.path.join(pickled_objects_dir, 'network-files', 'weighted-degree-series', 'atb-od.pkl'))
                                       
    mean_logFC_ATB_HC_series = pd.read_pickle(
        os.path.join(pickled_objects_dir, 'network-files', 'mean-logFC-network-nodes-series', 'atb-hc.pkl'))
    mean_logFC_ATB_LTBI_series = pd.read_pickle(
        os.path.join(pickled_objects_dir, 'network-files', 'mean-logFC-network-nodes-series', 'atb-ltbi.pkl'))
    mean_logFC_ATB_OD_series = pd.read_pickle(
        os.path.join(pickled_objects_dir, 'network-files', 'mean-logFC-network-nodes-series', 'atb-od.pkl'))
    
    top_N_nodes = 100
    
    ATB_HC_df = pd.DataFrame(index = mean_logFC_ATB_HC_series.index)
    ATB_HC_df['mean_log2FC'] = mean_logFC_ATB_HC_series.values
    ATB_HC_df['weighted_degree'] = weighted_deg_ATB_HC_series[mean_logFC_ATB_HC_series.index].values
    
    ATB_LTBI_df = pd.DataFrame(index = mean_logFC_ATB_LTBI_series.index)
    ATB_LTBI_df['mean_log2FC'] = mean_logFC_ATB_LTBI_series.values
    ATB_LTBI_df['weighted_degree'] = weighted_deg_ATB_LTBI_series[mean_logFC_ATB_LTBI_series.index].values
        
    ATB_OD_df = pd.DataFrame(index = mean_logFC_ATB_OD_series.index)
    ATB_OD_df['mean_log2FC'] = mean_logFC_ATB_OD_series.values
    ATB_OD_df['weighted_degree'] = weighted_deg_ATB_OD_series[mean_logFC_ATB_OD_series.index].values
        
    top_N_nodes_ATB_HC = set(
        ATB_HC_df.sort_values(by = 'weighted_degree', ascending = False).head(n=top_N_nodes).index)
    top_N_nodes_ATB_LTBI = set(
        ATB_LTBI_df.sort_values(by = 'weighted_degree', ascending = False).head(n=top_N_nodes).index)
    top_N_nodes_ATB_OD = set(
        ATB_OD_df.sort_values(by = 'weighted_degree', ascending = False).head(n=top_N_nodes).index)
    
    ATB_HC_only = top_N_nodes_ATB_HC - top_N_nodes_ATB_LTBI - top_N_nodes_ATB_OD

    #pickle list of nodes for downstream analysis
    atb_hc_only_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-hc-only.pkl')
    with open(atb_hc_only_path, 'wb') as f:
        pickle.dump(list(ATB_HC_only), f)
        
    ATB_LTBI_only = top_N_nodes_ATB_LTBI - top_N_nodes_ATB_HC - top_N_nodes_ATB_OD

    #pickle list of nodes for downstream analysis
    atb_ltbi_only_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-ltbi-only.pkl')
    with open(atb_ltbi_only_path, 'wb') as f:
        pickle.dump(list(ATB_LTBI_only), f)
        
    ATB_OD_only = top_N_nodes_ATB_OD - top_N_nodes_ATB_HC - top_N_nodes_ATB_LTBI

    #pickle list of nodes for downstream analysis
    atb_od_only_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-od-only.pkl')
    with open(atb_od_only_path, 'wb') as f:
        pickle.dump(list(ATB_OD_only), f)
        
    ATB_HC_and_ATB_LTBI = top_N_nodes_ATB_HC.intersection(top_N_nodes_ATB_LTBI) - top_N_nodes_ATB_OD

    #pickle list of nodes for downstream analysis
    atb_hc_atb_ltbi_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-hc-and-atb-ltbi.pkl')
    with open(atb_hc_atb_ltbi_path, 'wb') as f:
        pickle.dump(list(ATB_HC_and_ATB_LTBI), f)
        
    ATB_HC_and_ATB_OD = top_N_nodes_ATB_HC.intersection(top_N_nodes_ATB_OD) - top_N_nodes_ATB_LTBI

    #pickle list of nodes for downstream analysis
    atb_hc_atb_od_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-hc-and-atb-od.pkl')
    with open(atb_hc_atb_od_path, 'wb') as f:
        pickle.dump(list(ATB_HC_and_ATB_OD), f)
        
    ATB_LTBI_and_ATB_OD = top_N_nodes_ATB_LTBI.intersection(top_N_nodes_ATB_OD) - top_N_nodes_ATB_HC

    #pickle list of nodes for downstream analysis
    atb_ltbi_atb_od_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-ltbi-and-atb-od.pkl')
    with open(atb_ltbi_atb_od_path, 'wb') as f:
        pickle.dump(list(ATB_LTBI_and_ATB_OD), f)
        
    ATB_HC_and_ATB_LTBI_and_ATB_OD = top_N_nodes_ATB_HC.intersection(top_N_nodes_ATB_LTBI.intersection(top_N_nodes_ATB_OD))

    #pickle list of nodes for downstream analysis
    atb_hc_atb_ltbi_atb_od_path = os.path.join(
        pickled_objects_dir, 'network-files', 'top-weighted-node-lists', 'atb-hc-and-atb-ltbi-and-atb-od.pkl')
    with open(atb_hc_atb_ltbi_atb_od_path, 'wb') as f:
        pickle.dump(list(ATB_HC_and_ATB_LTBI_and_ATB_OD), f)
        
    return {
        'atb-hc-only': ATB_HC_only,
        'atb-ltbi-only': ATB_LTBI_only,
        'atb-od': ATB_OD_only,
        'atb-hc-and-atb-ltbi': ATB_HC_and_ATB_LTBI,
        'atb-hc-and-atb-od': ATB_HC_and_ATB_OD,
        'atb-ltbi-and-atb-od': ATB_LTBI_and_ATB_OD,
        'atb-hc-and-atb-ltbi-and-atb-od': ATB_HC_and_ATB_LTBI_and_ATB_OD
    }                            
        
# notebook 7
def generate_gene_set_heatmaps(data_dir):
    """
    """
    pass

# notebook 8
def run_gene_enrichment_analysis(data_dir, group):
    """
    """
    
    
    
    return
