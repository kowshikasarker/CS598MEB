import os
from multiprocessing import Pool, cpu_count

mu_list = [0.03, 0.25, 0.1, 0.2, 0.05, 0.3, 0.15]
#mu_list = [0.03, 0.15, 0.3]
pd_list = [0, 1]
run_list = list(range(1, 11))
#run_list = list(range(2, 11, 2))
distance_list = ['jensen_shannon', 'cluster_overlap', 'cluster_path']
threshold_list = [str(x/100) for x in range(10, 101, 10)]
#threshold_list = [0.3, 0.6, 0.9]

cell_column = 'cell_id'
cluster_column = 'ClusterIdent'
lineage_path = '/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/sl_elegans.nwk'
script = '/shared/nas/data/m1/ksarker2/CS598MEB/resolve_states.py'

base_common_command = 'python3 ' + script + ' -c1 ' + cell_column + ' -c2 ' + cluster_column + ' -l ' + lineage_path

command_list = []

for mu in mu_list:
    for pd in pd_list:
        for run in run_list:
            barcode_params = 'mu_' + str(mu) + '_pd_' + str(pd) + '_Nchar_9_run_' + str(run)
            base_dir = '/shared/nas/data/m1/ksarker2/CS598MEB/Data/Experiments/celegans/' + barcode_params
            cell_state = base_dir + '/profile_363_' + barcode_params + '.csv'
            common_command =  base_common_command + ' -s ' + cell_state
            
            for distance in distance_list:
                for threshold in threshold_list:
                    working_dir = base_dir + '/refined/Startle/' + distance + '/threshold_' + str(threshold)
                    #print('working_dir', working_dir)
                    in_tree = working_dir + '/refined.nwk'
                    #print('in_tree', in_tree)
                    command = common_command + ' -i ' + in_tree + ' -o ' + working_dir
                    command_list.append(command)

print(len(command_list), 'commands')

#cpu_count = cpu_count()
cpu_count = 40

command_list = [command_list[x:x+cpu_count] for x in range(0, len(command_list), cpu_count)]
print(command_list)

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