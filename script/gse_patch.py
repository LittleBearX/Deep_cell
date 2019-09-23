import os
import pandas as pd


def main(gse_data, gene_info, gene_ident):
    gene_info = pd.read_csv(gene_info, '\t')
    if gene_ident == 'geneSymbol':
        gene_info_dict = dict(zip(gene_info.Gene_name, gene_info.Gene_stable_ID))
    else:
        gene_info_dict = dict(zip(gene_info.Gene_stable_ID, gene_info.Gene_name))

    for gse in gse_data:
        gse_dir = gse
        go_file = os.path.join(gse_dir, 'go.csv')
        coexpression_file = os.path.join(gse_dir, 'coexpression_go.csv')
        if os.path.isfile(go_file):
            with open(go_file) as f:
                data1 = pd.read_csv(f, '\t')
                data1_hits = data1['hits']
                data1_hits_name = ['|'.join([gene_info_dict[gene] for gene in gene_list.split('|')])
                                   for gene_list in data1['hits']]
            if gene_ident == 'geneSymbol':
                data1_hits, data1_hits_name = data1_hits_name, data1_hits

            data1['hits'] = data1_hits
            data1['hits_name'] = data1_hits_name
            data1 = data1[['cluster', 'go', 'type', 'name', 'definition', 'p_value', 'odds_ratio', 'n_exp_hit',
                           'n_hit', 'n_cluster', 'n_go', 'n_universe', 'hits', 'hits_name']]
            with open(go_file, 'w') as f:
                data1.to_csv(f, '\t', index=False)

        if os.path.isfile(coexpression_file):
            with open(coexpression_file) as f:
                data2 = pd.read_csv(f, '\t')

            data2_ncoding = data2['ncoding']
            data2_ncoding_name = [gene_info_dict[gene] for gene in data2['ncoding']]
            data2_co_genes = data2['co_genes']
            data2_co_genes_name = ['|'.join([gene_info_dict[gene] for gene in gene_list.split('|')])
                                   for gene_list in data2['co_genes']]
            data2_hits = data2['hits']
            data2_hits_name = ['|'.join([gene_info_dict[gene] for gene in gene_list.split('|')])
                               for gene_list in data2['hits']]
            if gene_ident == 'geneSymbol':
                data2_ncoding, data2_ncoding_name = data2_ncoding_name, data2_ncoding
                data2_co_genes, data2_co_genes_name = data2_co_genes_name, data2_co_genes
                data2_hits, data2_hits_name = data2_hits_name, data2_hits

            data2['ncoding'] = data2_ncoding
            data2['ncoding_name'] = data2_ncoding_name
            data2['co_genes'] = data2_co_genes
            data2['co_genes_name'] = data2_co_genes_name
            data2['hits'] = data2_hits
            data2['hits_name'] = data2_hits_name
            data2 = data2[['cluster', 'ncoding', 'ncoding_name', 'co_genes', 'co_genes_name', 'go', 'type',
                           'name', 'definition', 'p_value',  'odds_ratio', 'n_exp_hit', 'n_hit', 'n_cluster',
                           'n_go', 'n_universe', 'hits', 'hits_name']]
            with open(coexpression_file, 'w') as f:
                data2.to_csv(f, '\t', index=False)
