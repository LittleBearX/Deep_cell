import os
import pandas as pd
from scipy import stats


def get_cell_marker_dict(cell_marker, gene_ident='ensemblID', sep='\t'):
    """read cell_marker, return organism marker gene of every cell type"""

    print('Reading cell_marker:', cell_marker)
    with open(cell_marker, encoding='utf8') as f:
        cell_marker_data = pd.read_csv(f, sep=sep)

    columns = cell_marker_data.columns
    for col_name in ['speciesType', 'tissueType', 'cellName', 'ensemblID', 'geneSymbol']:
        if col_name not in columns:
            raise AttributeError('cell_marker must have column: ' + col_name)

    cell_marker_dict = {}
    for organism, data in cell_marker_data.groupby('speciesType'):
        organism_cell_dict = {}
        for _, cell in data.iterrows():
            cell_type = (cell['tissueType'], cell['cellName'])
            if gene_ident == 'ensemblID':
                gene_set = set(cell['ensemblID'].strip().split(','))
            elif gene_ident == 'geneSymbol':
                gene_set = set(cell['geneSymbol'].strip().split(','))
            else:
                raise AttributeError('gene_ident must be in {geneSymbol, ensemblID}')
            organism_cell_dict[cell_type] = gene_set
        cell_marker_dict[organism] = organism_cell_dict

    print('get cell type marker info from these organisms: \n%s' % set(cell_marker_dict.keys()))
    if len(cell_marker_dict) == 0:
        exit()

    # count all marker gene of all organisms
    sum_dict = {}
    for organism in cell_marker_dict:
        organism_all_gene_set = set()
        for cell_type, gene_set in cell_marker_dict[organism].items():
            organism_all_gene_set.update(gene_set)
        sum_dict[organism] = organism_all_gene_set
    cell_marker_dict['all'] = sum_dict

    return cell_marker_dict


def gse_marker_handle(gse_data, organism, cell_marker_dict, odds_ratio_threshold=2,
                      p_value_threshold=0.01, method='greater'):
    """Fisher exact text for every cluster"""

    assert method in {'two-sided', 'less', 'greater'}
    all_gse_data = gse_data
    for count, gse in enumerate(all_gse_data, 1):
        marker_genes_file = os.path.join(gse, 'marker_genes.csv')
        if os.path.isdir(gse) and not os.path.isfile(marker_genes_file):
            text = f'Missing: {marker_genes_file}!'
            print(text)
        else:
            if organism not in cell_marker_dict:
                text = f'{gse}: Did not find marker genes.txt of {organism} in cell_marker!'
                print(text)
                continue

            text = f'Handling: {gse} {organism} ({count}/{len(all_gse_data)})'
            print(text)
            with open(marker_genes_file, 'r', encoding='utf8') as f:
                marker_genes_data = pd.read_csv(f, sep=',')

            item_list = []
            all_marker = cell_marker_dict['all'][organism]  # all marker
            for cluster, data in marker_genes_data.groupby('cluster'):
                cluster_marker = set(data['gene']) & all_marker  # marker in one cluster
                n_all_marker = len(all_marker)
                n_cluster_marker = len(cluster_marker)
                if n_cluster_marker == 0:
                    continue
                cluster_marker_prop = n_cluster_marker / n_all_marker  # proportion of cluster marker in all marker
                for cell_type, cell_type_marker in cell_marker_dict[organism].items():
                    n_cell_type_marker = len(cell_type_marker)  # marker in one cell type
                    # expected hit in random condition
                    n_expected_hit = cluster_marker_prop * n_cell_type_marker
                    hit_genes = cluster_marker & cell_type_marker
                    n_hit = len(hit_genes)
                    odds_ratio = n_hit / n_expected_hit
                    if odds_ratio > odds_ratio_threshold:
                        n_non_hit_cell_type_marker = n_cell_type_marker - n_hit
                        n_non_hit_cluster_marker = n_cluster_marker - n_hit
                        n_other_marker = n_all_marker - n_hit - n_non_hit_cell_type_marker - n_non_hit_cluster_marker
                        table = [[n_other_marker, n_non_hit_cell_type_marker], [n_non_hit_cluster_marker, n_hit]]
                        p_value = stats.fisher_exact(table, method)[1]
                        if p_value < p_value_threshold:
                            item = [cluster, organism, cell_type[0], cell_type[1], n_all_marker, n_cluster_marker,
                                    n_cell_type_marker, n_hit, n_expected_hit, odds_ratio, p_value, '|'.join(hit_genes)]
                            item_list.append(item)
            if item_list:
                item_data = pd.DataFrame(item_list)
                columns = ['cluster', 'organism', 'tissueType', 'cellName', 'n_all_marker', 'n_cluster_marker',
                           'n_cell_type_marker', 'n_hit', 'n_expected_hit', 'odds_ratio', 'p_value', 'hits']
                item_data.columns = columns
                item_data.sort_values(by=['cluster', 'p_value'], inplace=True)

                cells_type_file = os.path.join(gse, 'cells_type.csv')
                with open(cells_type_file, 'w', encoding='utf8') as f:
                    item_data.to_csv(f, index=False)
                text = f'Finished: {gse}'
                print(text)
            else:
                text = f'Not cluster can be marked to cell type: {gse}!'
                print(text)


def main(cell_marker, gse_data, organism, gene_ident):
    cell_marker_dict = get_cell_marker_dict(cell_marker, gene_ident, sep='\t')
    gse_marker_handle(gse_data, organism, cell_marker_dict)
