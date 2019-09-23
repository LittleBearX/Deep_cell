import os
import warnings
import pandas as pd
from scipy import stats
from rpy2.robjects import r
from rpy2.robjects import StrVector
from rpy2.robjects import NA_Logical
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr

warnings.filterwarnings('ignore')


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


def go_annotation(gse_data, gene_ident, gse_annotation, r_functions_dict, go_info_dict, odds_ratio_threshold=1,
                  p_value_threshold=0.01, method='greater', conditional=False, direction='over'):
    """use GO_stats"""

    assert direction in {'over', 'under'}
    new = r_functions_dict['new']
    mget = r_functions_dict['mget']
    unique = r_functions_dict['unique']
    unlist = r_functions_dict['unlist']
    hyper_test = r_functions_dict['hyperGTest']

    all_gse_data = gse_data
    for count, gse in enumerate(all_gse_data, 1):
        print('-----------------------')
        marker_genes_file = os.path.join(gse, 'marker_genes.csv')
        go_file = os.path.join(gse, 'go.csv')

        if os.path.isdir(gse) and not os.path.isfile(marker_genes_file):
            text = f'Missing: {marker_genes_file}!'
            print(text)
        else:
            text = f'Handling: {gse} {gse_annotation} ({count}/{len(all_gse_data)})'
            print(text)

            with open(marker_genes_file, 'r', encoding='utf8') as f:
                marker_genes_data = pd.read_csv(f, sep=',')
                universe_genes = marker_genes_data['gene'].unique()

            try:
                if gene_ident == 'geneSymbol':
                    annotation_entrez = r(gse_annotation.replace('.db', 'SYMBOL2EG'))
                    annotation_id = r(gse_annotation.replace('.db', 'SYMBOL'))
                else:
                    annotation_entrez = r(gse_annotation.replace('.db', 'ENSEMBL2EG'))
                    annotation_id = r(gse_annotation.replace('.db', 'ENSEMBL'))
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
                for cluster, cluster_genes in marker_genes_data.groupby('cluster')['gene']:
                    cluster_item_list = []
                    cluster_genes = set(cluster_genes) & go_universe_genes
                    n_cluster_genes = len(cluster_genes)
                    if n_cluster_genes == 0:
                        continue
                    cluster_genes_prop = n_cluster_genes / n_universe_genes
                    for go, go_genes in all_go_item_dict[go_type].items():
                        if go in go_info_dict:
                            _, name, definition = go_info_dict[go]
                        else:
                            continue
                        n_go_genes = len(go_genes)
                        n_expected_hit = cluster_genes_prop * n_go_genes
                        hit_genes = cluster_genes & go_genes
                        n_hit = len(hit_genes)
                        odds_ratio = n_hit / n_expected_hit
                        if odds_ratio > odds_ratio_threshold:
                            n_non_hit_go_genes = n_go_genes - n_hit
                            n_non_hit_cluster_genes = n_cluster_genes - n_hit
                            n_other_genes = n_universe_genes - n_hit - n_non_hit_go_genes - n_non_hit_cluster_genes
                            table = [[n_other_genes, n_non_hit_go_genes], [n_non_hit_cluster_genes, n_hit]]
                            p_value = stats.fisher_exact(table, method)[1]
                            if p_value < p_value_threshold:
                                item = [cluster, go, go_type, name, definition, p_value, odds_ratio, n_expected_hit,
                                        n_hit, n_cluster_genes, n_go_genes, n_universe_genes, '|'.join(hit_genes)]
                                cluster_item_list.append(item)
                    item_list.extend(cluster_item_list)
                    print('cluster: %s %s: %s' % (cluster, go_type, len(cluster_item_list)))

            # save_data
            if item_list:
                go_data = pd.DataFrame(item_list)
                columns = ['cluster', 'go', 'type', 'name', 'definition', 'p_value', 'odds_ratio', 'n_exp_hit',
                           'n_hit', 'n_cluster', 'n_go', 'n_universe', 'hits']
                go_data.columns = columns
                go_data.sort_values(['cluster', 'type', 'p_value'], inplace=True)
                with open(go_file, 'w', encoding='utf8') as f:
                    go_data.to_csv(f, '\t', index=False)
                text = f'Finished: {gse} ({count}/{len(all_gse_data)})'
                print(text)
            else:
                text = f'No go: {gse} ({count}/{len(all_gse_data)})'
                print(text)


def main(gse_data, gene_ident, gse_annotation, go_info):
    go_info_dict = get_go_info_dict(go_info)
    r_functions_dict = import_r_packages(gse_annotation)
    go_annotation(gse_data, gene_ident, gse_annotation, r_functions_dict, go_info_dict)
