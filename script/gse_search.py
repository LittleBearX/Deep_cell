import os
import pickle
import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr
from scipy.stats import spearmanr


def read_index(gse_index):
    """read pre-built gse_index"""

    with open(gse_index, 'rb') as f:
        gse_index_dict = pickle.load(f)

    return gse_index_dict


def search(project, gse_index_dict, organism, gene_info, gene_ident, count_by='mean', n_threshold=5, spearman=True):
    """find similar clusters among database"""

    data_log1p = project.data_log1p.T
    data_log1p.fillna(0, inplace=True)

    if project.all_marker_genes is not None:
        genes_selected = project.all_marker_genes.gene.unique()
    else:
        genes_selected = data_log1p.columns

    if gene_ident == 'geneSymbol':
        # transfer geneSymbol to ensemblID
        with open(gene_info, encoding='utf8') as f:
            gene_data = pd.read_csv(f, sep='\t')
        columns = gene_data.columns
        assert 'Gene_stable_ID' in columns
        assert 'Gene_name' in columns
        gene_data = gene_data.loc[gene_data['Organism'] == organism]
        genes_dict = dict(zip(gene_data['Gene_name'], gene_data['Gene_stable_ID']))
        genes_selected = [i for i in genes_selected if i in genes_dict]
        gene_id = [genes_dict[i] for i in genes_selected]
        data_log1p = data_log1p[genes_selected]
        data_log1p.columns = gene_id
    else:
        data_log1p = data_log1p[genes_selected]

    data_log1p['cluster'] = project.cells_ident
    target_clusters = None
    if count_by == 'mean':
        target_clusters = data_log1p.groupby('cluster').mean()
    elif count_by == 'median':
        target_clusters = data_log1p.groupby('cluster').median()
    else:
        assert count_by in {'mean', 'median'}

    genes = target_clusters.columns

    all_result_list = []
    for gse in tqdm(gse_index_dict[organism][count_by]):
        join_genes = genes & gse.columns
        if len(join_genes) > n_threshold:
            join_target_clusters = target_clusters[join_genes]
            join_search_clusters = gse[join_genes]
            for i in join_target_clusters.index:
                for j in join_search_clusters.index:
                    pearson_result = pearsonr(join_target_clusters.loc[i], join_search_clusters.loc[j])
                    result_list = [i, j, pearson_result[0], pearson_result[1]]
                    if spearman:
                        spearman_result = spearmanr(join_target_clusters.loc[i], join_search_clusters.loc[j])
                        result_list.extend([spearman_result.correlation, spearman_result.pvalue])
                    all_result_list.append(result_list)

    columns = ['cluster', 'GSE_cluster', 'pearson_corr', 'pearson_pvalue']
    if spearman:
        columns.extend(['spearman_corr', 'spearman_pvalue'])
    all_result = pd.DataFrame(all_result_list, columns=columns)
    if spearman:
        all_result.sort_values('spearman_corr', ascending=False, inplace=True)
    else:
        all_result.sort_values('pearson_corr', ascending=False, inplace=True)

    search_result_file = os.path.join(project.path, 'search_result.csv')
    with open(search_result_file, 'w', encoding='utf8') as f:
        all_result.to_csv(f, index=False)


def main(project, gse_index, organism, gene_info, gene_ident, search_by, spearman):
    gse_index_dict = read_index(gse_index)
    search(project, gse_index_dict, organism, gene_info, gene_ident, count_by=search_by, spearman=spearman)
