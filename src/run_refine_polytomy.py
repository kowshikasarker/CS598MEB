import os
from multiprocessing import Pool

mu_list = [0.03, 0.25, 0.1, 0.2, 0.05, 0.3, 0.15]
pd_list = [0, 1]
run_list = list(range(1, 11))
similarity_list = ['s1', 's2', 's3']
threshold_list = [x/100 for x in range(0, 101, 10)]

cell_column = 'cell_id'
cluster_column = 'ClusterIdent'
lineage_path = '/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/sl_elegans.nwk'
script = '/shared/nas/data/m1/ksarker2/CS598MEB/refine_polytomy.py'

base_common_command = 'python3 ' + script + ' -c1 ' + cell_column + ' -c2 ' + cluster_column + ' -l ' + lineage_path
command_list = []

for mu in mu_list:
    for pd in pd_list:
        for run in run_list:
            barcode_params = 'mu_' + str(mu) + '_pd_' + str(pd) + '_Nchar_9_run_' + str(run)
            if (barcode_params == 'mu_0.03_pd_1_Nchar_9_run_5'):
                continue
            in_dir = '/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/' + barcode_params
            base_out_dir = in_dir + '/refined/Startle'

            in_tree = in_dir + '/Startle_363_' + barcode_params + '_bin_tree.newick'
            cell_state = in_dir + '/profile_363_' + barcode_params + '.csv'

            common_command =  base_common_command + ' -i ' + in_tree + ' -e ' + cell_state
            
            for similarity in similarity_list:
                for threshold in threshold_list:
                    out_dir = base_out_dir + '/' + similarity + '/threshold_' + str(threshold)
                    command = common_command + ' -o ' + out_dir + ' -s ' + similarity + ' -th ' + str(threshold)
                    command_list.append(command)

print(len(command_list), 'commands')

command_list = [command_list[x:x+40] for x in range(0, len(command_list), 40)]
print(len(command_list))

def run(command):
    print('Started', command)
    exit_code = os.system(command)
    if(exit_code != 0):
        print('ERRROR CODE:', exit_code, command)
    else:
        print('Completed', command)
    print()

for sublist in command_list:
    pool = Pool()
    pool.map(run, sublist)
    pool.close()
    pool.join()
