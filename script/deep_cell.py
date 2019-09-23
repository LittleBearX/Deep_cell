import os
import sys
import getopt
import warnings

path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path)
warnings.filterwarnings('ignore')


from gse_worker import Project
from gse_go import main as gse_go_main
from gse_split import main as gse_split_main
from gse_search import main as gse_search_main
from gse_patch import main as gse_patch_main
from gse_marker import main as gse_marker_main
from gse_coexpression import main as gse_coexpression_main


annotation_dict = {'Homo_sapiens': 'org.Hs.eg.db', 'Mus_musculus': 'org.Mm.eg.db',
                   'Drosophila_melanogaster': 'org.Dm.eg.db', 'Danio_rerio': 'org.Dr.eg.db',
                   'Rattus_norvegicus': 'org.Rn.eg.db'}


def main(gene_info, go_info, matrix, organism, output, gse_index, cell_marker, gene_ident, do_filter=True,
         scale_factor=10000, y_low_cutoff=0.5, scale_max=10., only_hvg=False, n_pca='auto', do_corr=True, n_knn=30,
         do_fast=0, do_snn=True, resolution='auto', n_iter_tsne=1000, three_tsne=True, only_pos=True, min_pvalue=0.01,
         top_diff=10, sort_by='mean.diff', save_data=True, save_pdf=True, search=True, search_by='mean', spearman=True,
         go=True, coexpression=True, cell_type=True):
    print('\n=======================')
    print('Splitting info')
    gse_split_main(gene_info, matrix, organism, output, gene_ident)
    print('Finished!')

    all_file = ['ALL_data', 'CODING_data', 'NCODING_data']
    print('\n=======================')
    print('Using algorithm in data')
    for i in all_file:
        print('-----------------------')
        file = os.path.join(output, i, 'matrix.csv')
        if os.path.isfile(file):
            print('Handling file: %s' % file)
            project = Project(file)
            project.check()
            if do_filter:
                project.filter_data(min_genes_percent=0.001, min_cells_percent=0.005, max_cells_percent=0.995)
            else:
                project.filter_data(min_genes_percent=0, min_cells_percent=0, max_cells_percent=1)
            project.normalize_data(scale_factor=scale_factor)
            project.find_variable_genes(bin_nums=20, x_low_cutoff=0.0125, x_high_cutoff=3, y_low_cutoff=y_low_cutoff)
            project.scale_data(scale_max=scale_max)
            project.run_pca(only_hvg=only_hvg, n=n_pca, corr=do_corr)
            project.find_clusters(k=n_knn, do_fast=do_fast, snn=do_snn, resolution=resolution, modularity=1,
                                  algorithm=1, n_start=100, n_iter=10)
            project.run_tsne(n_iter=n_iter_tsne, three=three_tsne)
            project.find_all_marker_genes(min_pct=0.25, min_diff_pct=0.0, log_threshold=0.25, only_pos=only_pos,
                                          min_pvalue=min_pvalue)
            project.run_heat_map(top=top_diff, sort_by=sort_by)
            if save_data:
                project.save_data()
            if save_pdf:
                project.save_picture(pdf=True)

            if search:
                gse_search_main(project, gse_index, organism, gene_info, gene_ident, search_by, spearman)

        else:
            print('Warning: not such file: %s' % file)
    print('Finished!')

    if go:
        print('\n=======================')
        print('GO analysis of DE genes for clusters')
        try:
            gse_annotation = annotation_dict[organism.replace(' ', '_')]
            gse_go_main([os.path.join(output, i) for i in all_file], gene_ident, gse_annotation, go_info)
            print('Finished!')
        except KeyError:
            print('No annotation for %s' % organism)

    if coexpression:
        print('\n=======================')
        print('Inference of go function of non-coding genes from similar coding genes')
        try:
            gse_annotation = annotation_dict[organism.replace(' ', '_')]
            gse_coexpression_main([os.path.join(output, i) for i in all_file],
                                  gene_info, organism, gse_annotation, gene_ident, go_info, spearman)
            print('Finished!')
        except KeyError:
            print('No annotation for %s' % organism)

    try:
        gse_patch_main([os.path.join(output, i) for i in all_file], gene_info, gene_ident)
    except Exception as e:
        print(e)

    if cell_type:
        print('\n=======================')
        print('Finding cell type in cluster')
        gse_marker_main(cell_marker, [os.path.join(output, i) for i in all_file], organism, gene_ident)
        print('Finished!')


if __name__ == '__main__':
    param_dict = {'geneInfo': 'File contains the coding genes and non-coding genes.',
                  'goInfo': 'File contains go information.',
                  'matrix': 'File of gene expression matrix, rows for samples and columns for genes, CSV format.',
                  'organism': 'Organism.',
                  'output': 'Directory to put the output data.',
                  'gseIndex': 'Pre-built index for database.',
                  'cellMarker': 'File contains marker genes for cell type.',
                  'geneIdent': 'Genes with ensemblID or geneSymbol in matrix file.',
                  'do_filter': 'If true, low expression samples and genes will be filtered.',
                  'scale_factor': 'Scale the data.',
                  'scale_max': 'Clip data after scaling.',
                  'y_low_cutoff': 'Threshold of high variable genes(0.25-1 recommended).',
                  'only_hvg': 'If true, only high variable genes will be taken into account.',
                  'n_pca': 'Number of principal component.',
                  'do_corr': 'If true, Pearson correlation coefficient matrix will be calculated.',
                  'n_knn': 'N samples with shared nearest neighbor.',
                  'do_fast': 'do_fast.',
                  'do_snn': 'If true, shared nearest neighbor matrix will be calculated.',
                  'resolution': 'Higher resolution for more clusters(0.6-2 recommended).',
                  'n_iter_tsne': 'Number of iterations for t-SNE.',
                  'three_tsne': 'If true, t-SNE of three dimensions will be calculated.',
                  'only_pos': 'If true, only up-regulated genes will be taken into account.',
                  'min_pvalue': 'Threshold of adjusted p-value for differential expression genes.',
                  'top_diff': 'Choose the top n of the most different expression genes to draw heatmap.',
                  'sort_by': 'The most different genes is according to: '
                             '\n\t\t1. p_value: compute by rank sum test. '
                             '\n\t\t2. count.diff: percentage of gene expression difference.'
                             '\n\t\t3. mean.diff: gene expression level difference.',
                  'save_data': 'If true, save data.',
                  'save_pdf': 'If true, save pdf.',
                  'search': 'If true, search similar clusters among database.',
                  'search_by': 'Use mean or median expression for cluster while searching.',
                  'spearman': 'If true, calculate spearman correlation coefficient while searching and co-expression.',
                  'go': 'If true, go analysis.',
                  'coExpression': 'If true, infer go function of non-coding genes from co-expression coding genes.',
                  'cell_type': 'If true, calculate cell type.'}

    if '--help' in sys.argv[1:]:
        print('Parameters list:')
        for param in param_dict:
            print(f'\t- {param}: {param_dict[param]}')
        exit()

    options, _ = getopt.getopt(sys.argv[1:], '', [param + '=' for param in param_dict])
    options_dict = dict(options)
    miss_key = {'--' + param for param in param_dict} - options_dict.keys()
    if miss_key:
        raise AttributeError('missing parameter: %s' % str(miss_key))

    GENE_INFO = options_dict['--geneInfo']
    GO_INFO = options_dict['--goInfo']
    MATRIX = options_dict['--matrix']
    ORGANISM = options_dict['--organism']
    OUTPUT = options_dict['--output']
    GSE_INDEX = options_dict['--gseIndex']
    CELL_MARKER = options_dict['--cellMarker']
    GENE_IDENT = options_dict['--geneIdent']
    DO_FILTER = True if options_dict['--do_filter'] == 'True' else False
    SCALE_FACTOR = int(options_dict['--scale_factor'])
    Y_LOW_CUTOFF = float(options_dict['--y_low_cutoff'])
    SCALE_MAX = float(options_dict['--scale_max'])
    ONLY_HVG = True if options_dict['--only_hvg'] == 'True' else False
    N_PCA = 'auto' if options_dict['--n_pca'] == 'auto' else int(options_dict['--n_pca'])
    DO_CORR = True if options_dict['--do_corr'] == 'True' else False
    N_KNN = int(options_dict['--n_knn'])
    DO_FAST = int(options_dict['--do_fast'])
    DO_SNN = True if options_dict['--do_snn'] == 'True' else False
    RESOLUTION = 'auto' if options_dict['--resolution'] == 'auto' else float(options_dict['--resolution'])
    N_ITER_TSNE = int(options_dict['--n_iter_tsne'])
    THREE_TSNE = True if options_dict['--three_tsne'] == 'True' else False
    ONLY_POS = True if options_dict['--only_pos'] == 'True' else False
    MIN_PVALUE = float(options_dict['--min_pvalue'])
    TOP_DIFF = int(options_dict['--top_diff'])
    SORT_BY = options_dict['--sort_by']
    SAVE_DATA = True if options_dict['--save_data'] == 'True' else False
    SAVE_PDF = True if options_dict['--save_pdf'] == 'True' else False
    SEARCH = True if options_dict['--search'] == 'True' else False
    SEARCH_BY = options_dict['--search_by']
    SPEARMAN = True if options_dict['--spearman'] == 'True' else False
    GO = True if options_dict['--go'] == 'True' else False
    CO_EXPRESSION = True if options_dict['--coExpression'] == 'True' else False
    CELL_TYPE = True if options_dict['--cell_type'] == 'True' else False

    main(GENE_INFO, GO_INFO, MATRIX, ORGANISM, OUTPUT, GSE_INDEX, CELL_MARKER, GENE_IDENT, do_filter=DO_FILTER,
         scale_factor=SCALE_FACTOR, y_low_cutoff=Y_LOW_CUTOFF, scale_max=SCALE_MAX, only_hvg=ONLY_HVG, n_pca=N_PCA,
         do_corr=DO_CORR, n_knn=N_KNN, do_fast=DO_FAST, do_snn=DO_SNN, resolution=RESOLUTION, n_iter_tsne=N_ITER_TSNE,
         three_tsne=THREE_TSNE, only_pos=ONLY_POS, min_pvalue=MIN_PVALUE, top_diff=TOP_DIFF, sort_by=SORT_BY,
         save_data=SAVE_DATA, save_pdf=SAVE_PDF, search=SEARCH, search_by=SEARCH_BY, spearman=SPEARMAN, go=GO,
         coexpression=CO_EXPRESSION, cell_type=CELL_TYPE)
