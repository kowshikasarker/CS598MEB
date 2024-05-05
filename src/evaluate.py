# compute unweighted RF distance for all methods on the celegans data

import dendropy
from dendropy.calculate import treecompare
import pandas as pd
import os

mu_list = [0.03, 0.25, 0.1, 0.2, 0.05, 0.3, 0.15]
pd_list = [0, 1]
run_list = list(range(1, 11))
criterion_list = ['s1', 'd1', 'd2', 'd3']
threshold_list = [str(x/100) for x in range(0, 101, 10)]

true_tree = '/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/ground_truth_cell_lineage.nwk'

filename_map = {
    'Startle-NNI': ['Startle_363_', '_bin_tree.newick'],
    'Cassiopeia-Greedy': ['Cas_363_', '_bin_tree.newick'],
    'LinRace-IST': ['LinRace_363_', '_bin_tree.newick'],
    'LinTIMaT': ['LinTIMaT_363_', '_bin_tree.newick'],
    'DCLEAR-kmer-NJ': ['DCLEAR_363_', '_bin_tree.newick'],
}

result_dict = {
    'Mutation_Rate': [],
    'Dropout': [],
    'Run': [],
    'Startle-NNI': [],
    'Cassiopeia-Greedy': [],
    'LinRace-IST': [],
    'LinTIMaT': [],
    'DCLEAR-kmer-NJ': [],
}


for criterion in criterion_list:
    for th in threshold_list:
        result_dict['CS598_' + criterion + '_th=' + th] = []

print(result_dict)


def calc_rf_dist(tree1_path, tree2_path):
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get_from_path(tree1_path, schema="newick", preserve_underscores=True, taxon_namespace=tns)
    tree2 = dendropy.Tree.get_from_path(tree2_path, schema="newick", preserve_underscores=True, taxon_namespace=tns)

    rf_dist = treecompare.symmetric_difference(tree1, tree2)
    normalizer = len(tree1.internal_nodes()) + len(tree2.internal_nodes())
    rf_dist = round(rf_dist / normalizer, 2)
    return rf_dist

for mu in mu_list:
    for dp in pd_list:
        for run in run_list:
            barcode_params = 'mu_' + str(mu) + '_pd_' + str(dp) + '_Nchar_9_run_' + str(run)
            print(barcode_params)
            if(barcode_params == 'mu_0.03_pd_1_Nchar_9_run_5'):
                continue
            outer_dir = '/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/' + barcode_params

            result_dict['Mutation_Rate'].append(mu)
            result_dict['Dropout'].append(dp)
            result_dict['Run'].append(run)
            
            for method, ext_pair in filename_map.items():
                pred_tree = outer_dir + '/' + ext_pair[0] + barcode_params + ext_pair[1]
                print(method, pred_tree)
                result_dict[method].append(calc_rf_dist(pred_tree, true_tree))
            for criterion in criterion_list:
                for threshold in threshold_list:
                    inner_dir = outer_dir + '/refined/Startle/' + criterion + '/threshold_' + threshold
                    pred_tree = inner_dir + '/refined.nwk'
                    method = 'CS598_' + criterion + '_th=' + threshold
                    result_dict[method].append(calc_rf_dist(pred_tree, true_tree))
result_df = pd.DataFrame.from_dict(result_dict)
result_df.to_csv('/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/RF-Distance.tsv', sep='\t', index=False)


