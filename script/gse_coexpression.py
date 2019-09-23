import os
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from rpy2.robjects import r
import matplotlib.pyplot as plt
from rpy2.robjects import StrVector
from rpy2.robjects import NA_Logical
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr


warnings.filterwarnings('ignore')


def get_organism_genes_dict(gene_info, gene_ident, sep='\t'):
    """get organism gene info from GENE_info"""

    print('Reading GENE_info: ', gene_info)
    with open(gene_info, encoding='utf8') as f:
        gene_data = pd.read_csv(f, sep=sep)

    columns = gene_data.columns
    assert 'Gene_stable_ID' in columns
    assert 'Gene_name' in columns
    assert 'Organism' in columns
    assert 'Coding' in columns
    organism_genes_dict = {}
    for organism, data in gene_data.groupby('Organism'):
        if gene_ident == 'ensemblID':
            foo = {'coding': set(data['Gene_stable_ID'][data['Coding'] == 'coding']),
                   'ncoding': set(data['Gene_stable_ID'][data['Coding'] == 'ncoding'])}
        elif gene_ident == 'geneSymbol':
            foo = {'coding': set(data['Gene_name'][data['Coding'] == 'coding']),
                   'ncoding': set(data['Gene_name'][data['Coding'] == 'ncoding'])}
        else:
            raise AttributeError("gene_ident must be in {'geneSymbol', 'ensemblID'}")
        organism_genes_dict[organism] = foo
    if len(organism_genes_dict) == 0:
        raise AttributeError('Not info in GENE_info: ', gene_info)

    return organism_genes_dict


def get_go_info_dict(go_info, sep='\t'):
    """get go info dict"""

    print('Reading GO_info: ', go_info)
    with open(go_info, encoding='utf8') as f:
        go_info_data = pd.read_csv(f, sep=sep)
    columns = go_info_data.columns
    assert 'go' in columns
    assert 'old_go' in columns
    assert 'name' in columns
    assert 'definition' in columns
    go, old_go, name, definition = \
        go_info_data['go'], go_info_data['old_go'], go_info_data['name'], go_info_data['definition']
    go_info_dict = {b: [a, c, d] for a, b, c, d in zip(go, old_go, name, definition)}
    print('Reading %d go information' % len(go_info_dict))
    if len(go_info_dict) == 0:
        exit()

    return go_info_dict


def import_r_packages(gse_annotation):
    """import R packages, R functions and annotation"""

    r_packages = {'AnnotationDbi', 'GOstats'}
    r_functions = {'new', 'unique', 'unlist', 'mget', 'hyperGTest'}

    print('import R packages: %s' % r_packages)
    for r_package in r_packages:
        try:
            importr(r_package)
            print('√\t%s' % r_package)
        except RRuntimeError:
            print('×\t%s' % r_package)
            exit()

    print('import R function: %s' % r_functions)
    r_functions_dict = {}
    for r_function in r_functions:
        try:
            r_functions_dict[r_function] = r(r_function)
            print('√\t%s' % r_function)
        except RRuntimeError:
            print('×\t%s' % r_function)
            exit()

    print('import annotation: %s' % gse_annotation)
    try:
        importr(gse_annotation)
        print('√\t%s' % gse_annotation)
    except RRuntimeError:
        print('×\t%s' % gse_annotation)
        exit()

    return r_functions_dict


def co_expression_go(gse_data, organism_genes_dict, r_functions_dict, go_info_dict, organism, gse_annotation,
                     gene_ident, spearman=True, samples=1000, coef_threshold=0.5, odds_ratio_threshold=1,
                     p_value_threshold=0.01, method='greater', conditional=False, direction='over'):
    """Inference of go function of non-coding genes from similar coding genes"""

    new = r_functions_dict['new']
    mget = r_functions_dict['mget']
    unique = r_functions_dict['unique']
    unlist = r_functions_dict['unlist']
    hyper_test = r_functions_dict['hyperGTest']

    all_gse_data = gse_data
    for count, gse in enumerate(all_gse_data, 1):
        print('-----------------------')
        gse_file = os.path.join(gse, 'matrix.csv')
        marker_genes_file = os.path.join(gse, 'marker_genes.csv')
        co_file = os.path.join(gse, 'coexpression_go.csv')
        co_hist = os.path.join(gse, 'coexpression_hist.png')
        if os.path.isdir(gse):
            if not os.path.isfile(gse_file):
                text = f'Missing: {gse_file}!'
                print(text)
                continue
            if not os.path.isfile(marker_genes_file):
                text = f'Missing: {marker_genes_file}!'
                print(text)
                continue

            file_size = '%.3fM' % (os.path.getsize(gse_file)/(10**6))
            text = f'Handling: {gse} {organism} {file_size} ({count}/{len(all_gse_data)})'
            print(text)

            with open(marker_genes_file, 'r', encoding='utf8') as f:
                marker_genes_data = pd.read_csv(f, sep=',')
                marker_genes_data['cluster'] = marker_genes_data['cluster'].astype('int').astype('str')

            with open(gse_file, 'r', encoding='utf8') as f:
                matrix_data = pd.read_csv(f, index_col=0)
                matrix_data.fillna(0, inplace=True)
                n_genes, n_cells = matrix_data.shape
                cells_count = (matrix_data > 0.00001).sum(axis=0)
                genes_count = (matrix_data > 0.00001).sum(axis=1)
                min_cells_count = cells_count.sort_values()[int(n_cells*0.05)]
                min_genes_count = int(n_cells * 0.1)
                matrix_data = matrix_data.loc[:, cells_count > min_cells_count]
                data_norm = (matrix_data * 10000 / matrix_data.sum(axis=0))
                data_log1p = np.log1p(data_norm)

            coding = [gene for gene in data_log1p.index if gene in organism_genes_dict[organism]['coding']]
            ncoding = [gene for gene in data_log1p.index if gene in organism_genes_dict[organism]['ncoding']]
            coding_marker = [gene for gene in marker_genes_data['gene'].unique() if gene in coding]
            ncoding_marker = [gene for gene in marker_genes_data['gene'].unique() if gene in ncoding]
            selected_marker = [gene for gene in marker_genes_data['gene'].unique()
                               if (gene in coding_marker or gene in ncoding_marker)]
            sub_data = data_log1p.loc[selected_marker].T
            if not coding_marker:
                print('Lack of coding genes in marker genes!')
                continue
            elif not ncoding_marker:
                print('Lack of non coding genes in marker genes!')
                continue

            target_genes = set(genes_count[genes_count > min_genes_count].index)
            foo_coding = [i for i in coding if i in target_genes]
            foo_ncoding = [i for i in ncoding if i in target_genes]
            samples = min([samples, len(foo_coding), len(foo_ncoding)])
            coding_genes = list(np.random.choice(foo_coding, samples, replace=False))
            ncoding_genes = list(np.random.choice(foo_ncoding, samples, replace=False))
            selected_genes = coding_genes + ncoding_genes
            sub_data2 = data_log1p.loc[selected_genes].T

            if spearman:
                corr_matrix = stats.spearmanr(sub_data)[0]
                corr_matrix2 = stats.spearmanr(sub_data2)[0]
            else:
                corr_matrix = np.corrcoef(sub_data.T)
                corr_matrix2 = np.corrcoef(sub_data2.T)

            corr_matrix = pd.DataFrame(corr_matrix, index=selected_marker, columns=selected_marker)
            corr_matrix = corr_matrix.loc[coding_marker, ncoding_marker]
            corr_matrix2 = pd.DataFrame(corr_matrix2, index=selected_genes, columns=selected_genes)
            corr_matrix2 = corr_matrix2.loc[coding_genes, ncoding_genes]

            coef_of_same_cluster = []
            coef_of_diff_cluster = []
            for cluster, genes in marker_genes_data.groupby('cluster')['gene']:
                cluster_coding_marker = {gene for gene in genes if gene in coding_marker}
                cluster_ncoding_marker = {gene for gene in genes if gene in ncoding_marker}
                non_cluster_coding_marker = set(coding_marker) - cluster_coding_marker
                coef_of_cluster = corr_matrix.loc[cluster_coding_marker, cluster_ncoding_marker]
                coef_of_non_cluster = corr_matrix.loc[non_cluster_coding_marker, cluster_ncoding_marker]
                coef_of_cluster = np.array(coef_of_cluster).flatten()
                coef_of_non_cluster = np.array(coef_of_non_cluster).flatten()
                coef_of_same_cluster.extend(coef_of_cluster)
                coef_of_diff_cluster.extend(coef_of_non_cluster)
            coef_of_random = list(np.array(corr_matrix2).flatten())
            plt.figure(figsize=(8, 8))
            sns.distplot(coef_of_same_cluster, hist=False, color='g', kde_kws={'shade': True}, label='same')
            sns.distplot(coef_of_diff_cluster, hist=False, color='r', kde_kws={'shade': True}, label='diff')
            sns.distplot(coef_of_random, hist=False, color='y', kde_kws={'shade': True}, label='random')
            plt.yticks([])
            plt.savefig(co_hist, dpi=300)
            plt.close('all')

            try:
                if gene_ident == 'geneSymbol':
                    annotation_entrez = r(gse_annotation.replace('.db', 'SYMBOL2EG'))
                    annotation_id = r(gse_annotation.replace('.db', 'SYMBOL'))
                else:
                    annotation_entrez = r(gse_annotation.replace('.db', 'ENSEMBL2EG'))
                    annotation_id = r(gse_annotation.replace('.db', 'ENSEMBL'))
                universe_genes = coding_marker
                r_universe_genes = StrVector(universe_genes)
                r_universe_entrez = unique(unlist(mget(r_universe_genes, annotation_entrez, ifnotfound=NA_Logical)))

                all_go_item_dict = {}
                for go_type in ['BP', 'CC', 'MF']:
                    params = new("GOHyperGParams", geneIds=r_universe_entrez, universeGeneIds=r_universe_entrez,
                                 annotation=gse_annotation, ontology=go_type, pvalueCutoff=p_value_threshold,
                                 conditional=conditional, testDirection=direction)
                    hyper_result = hyper_test(params)

                    foo_data = hyper_result.do_slot('goDag').do_slot('nodeData').do_slot('data').items()
                    go_item_dict = {i[0]: {j[0] for j in mget(i[1][0], annotation_id)} for i in foo_data}
                    all_go_item_dict[go_type] = go_item_dict
                    print('universe %s: %s' % (go_type, len(go_item_dict)))
            except RRuntimeError as e:
                print('Error：%s' % e)
                continue

            item_list = []
            for go_type in ['BP', 'CC', 'MF']:
                go_universe_genes = {gene for go in all_go_item_dict[go_type] for gene in all_go_item_dict[go_type][go]}
                go_universe_genes = go_universe_genes & set(universe_genes)
                n_universe_genes = len(go_universe_genes)
                for ncoding_gene in ncoding_marker:
                    clusters = marker_genes_data[marker_genes_data['gene'] == ncoding_gene]['cluster'].tolist()
                    gene_item_list = []
                    foo = corr_matrix.loc[go_universe_genes, ncoding_gene]
                    foo = foo[foo > coef_threshold]
                    co_genes = set(foo.index)
                    n_co_genes = len(co_genes)
                    if n_co_genes == 0:
                        continue
                    co_genes_prop = n_co_genes / n_universe_genes
                    for go, go_genes in all_go_item_dict[go_type].items():
                        if go in go_info_dict:
                            _, name, definition = go_info_dict[go]
                        else:
                            continue
                        n_go_genes = len(go_genes)
                        n_expected_hit = co_genes_prop * n_go_genes
                        hit_genes = co_genes & go_genes
                        n_hit = len(hit_genes)
                        odds_ratio = n_hit / n_expected_hit
                        if odds_ratio > odds_ratio_threshold:
                            n_non_hit_go_genes = n_go_genes - n_hit
                            n_non_hit_co_genes = n_co_genes - n_hit
                            n_other_genes = n_universe_genes - n_hit - n_non_hit_go_genes - n_non_hit_co_genes
                            table = [[n_other_genes, n_non_hit_go_genes], [n_non_hit_co_genes, n_hit]]
                            p_value = stats.fisher_exact(table, method)[1]
                            if p_value < p_value_threshold:
                                item = ['|'.join(clusters), ncoding_gene, '|'.join(co_genes), go, go_type, name,
                                        definition, p_value, odds_ratio, n_expected_hit, n_hit, n_co_genes,
                                        n_go_genes, n_universe_genes, '|'.join(hit_genes)]
                                gene_item_list.append(item)
                    item_list.extend(gene_item_list)

            # save_data
            if item_list:
                co_data = pd.DataFrame(item_list)
                columns = ['cluster', 'ncoding', 'co_genes', 'go', 'type', 'name', 'definition', 'p_value',
                           'odds_ratio', 'n_exp_hit', 'n_hit', 'n_cluster', 'n_go', 'n_universe', 'hits']
                co_data.columns = columns
                co_data.sort_values(['cluster', 'ncoding', 'type', 'p_value'], inplace=True)
                with open(co_file, 'w', encoding='utf8') as f:
                    co_data.to_csv(f, '\t', index=False)
                text = f'Task:{count}/{len(all_gse_data)}\tsuccess:{gse}' \
                       f'\tnon coding count:{len(set(co_data["ncoding"]))}\tgo items:{len(item_list)}'
                print(text)
            else:
                text = f'Task:{count}/{len(all_gse_data)}\tfail:{gse}' \
                       f'\tnon coding count:{len(set(co_data["ncoding"]))}\tgo items:{len(item_list)}'
                print(text)


def main(gse_data, gene_info, organism, gse_annotation, gene_ident, go_info, spearman):
    organism_genes_dict = get_organism_genes_dict(gene_info, gene_ident)
    go_info_dict = get_go_info_dict(go_info)
    r_functions_dict = import_r_packages(gse_annotation)
    co_expression_go(gse_data, organism_genes_dict, r_functions_dict, go_info_dict, organism, gse_annotation,
                     gene_ident, spearman=spearman, coef_threshold=0.5, odds_ratio_threshold=1, p_value_threshold=0.01,
                     method='greater')
