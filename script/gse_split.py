import os
import shutil
import pandas as pd


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


def make_dir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass


def split_gse(matrix, organism, organism_genes_dict, output):
    """split matrix into coding and non-coding by gene info"""

    coding_data = os.path.join(output, 'CODING_data')
    ncoding_data = os.path.join(output, 'NCODING_data')
    all_data = os.path.join(output, 'ALL_data')
    make_dir(coding_data)
    make_dir(ncoding_data)
    make_dir(all_data)

    organism = organism.replace(' ', '_')

    if organism not in organism_genes_dict:
        text = f'organism must be in {set(organism_genes_dict.keys())}'
        raise AttributeError(text)

    file_size = '%.3fM' % (os.path.getsize(matrix) / (10 ** 6))
    text = f'Handling: {matrix}\t{organism}\t{file_size}'
    print(text)
    with open(matrix) as f:
        matrix_data = pd.read_csv(f, index_col=0)

    coding = organism_genes_dict[organism]['coding']
    ncoding = organism_genes_dict[organism]['ncoding']
    coding_genes = matrix_data.index[matrix_data.index.isin(coding)]
    ncoding_genes = matrix_data.index[matrix_data.index.isin(ncoding)]

    if len(coding_genes):
        print('Coding genes.txt: %d' % len(coding_genes))
        coding_file = os.path.join(coding_data, 'matrix.csv')
        with open(coding_file, 'w') as f:
            foo = matrix_data.loc[coding_genes, :]
            foo.to_csv(f, sep=',')
    else:
        text = f"Didn't find coding genes.txt!"
        print(text)

    if len(ncoding_genes):
        print('Non-coding genes.txt: %d' % len(ncoding_genes))
        ncoding_file = os.path.join(ncoding_data, 'matrix.csv')
        with open(ncoding_file, 'w') as f:
            foo = matrix_data.loc[ncoding_genes, :]
            foo.to_csv(f, sep=',')
    else:
        text = f"Didn't find non-coding genes.txt!"
        print(text)

    shutil.copy(matrix, all_data)


def main(gene_info, matrix, organism, output, gene_ident):
    organism_genes_dict = get_organism_genes_dict(gene_info, gene_ident)
    split_gse(matrix, organism, organism_genes_dict, output)
