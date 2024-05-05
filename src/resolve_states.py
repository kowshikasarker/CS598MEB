import argparse
from pathlib import Path
import os
import dendropy
import pandas as pd
import json
import time

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_tree", type=str,
                        help="path to the input tree in newick format",
                        required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="path to the output folder",
                        required=True, default=None)
    parser.add_argument("-s", "--cell_state", type=str,
                        help="path to cell to cluster mapping",
                        required=True, default=None)
    parser.add_argument("-c1", "--cell_column", type=str,
                        help="name of the cell column in the cell to cluster mapping",
                        required=True, default=None)
    parser.add_argument("-c2", "--cluster_column", type=str,
                        help="name of the cluster column in the cell to cluster mapping",
                        required=True, default=None)
    parser.add_argument("-l", "--lineage_path", type=str,
                        help="path to the cluster lineage in newick format",
                        required=False, default=None)
    args = parser.parse_args()
    return args

def assign_taxon_name_to_internal_nodes(tree):
    #print(tree.__str__())
    serial = 1
    for node in tree.internal_nodes():
        node.taxon = dendropy.Taxon('Int_' + str(serial))
        serial += 1

    #print(tree.__str__())
    return tree

'''def calculate_state_change_cost(cluster_lineage, cost_map=None):
    if(cost_map is None):
        cost_map = {}
    
    node = cluster_lineage.seed_node
    if(node.is_leaf()):
        node_name = int(node.taxon.label)
    else:
        node_name = int(node.label)
    #print('node', node, node.label)
    cost_map[(node_name, node_name)] = 0
    for child in node.child_nodes():
        #print('child', child, child.label)
        if(child.is_leaf()):
            child_name = int(child.taxon.label)
        else:
            child_name = int(child.label)
        new_cost_map = {}
        for state_pair, cost in cost_map.items():
            if(state_pair[1] == node_name):
                new_cost_map[(int(state_pair[0]), int(child_name))] = cost + 1
                new_cost_map[(int(child_name), int(state_pair[0]))] = cost + 1
        new_cost_map[(int(node_name), int(child_name))] = 1
        new_cost_map[(int(child_name), int(node_name))] = 1
        cost_map.update(new_cost_map)
        
        if not child.is_leaf():
            #print(child)
            #print(cluster_lineage)
            subtree = dendropy.Tree(seed_node=child.extract_subtree(suppress_unifurcations=False))
            cost_map = calculate_state_change_cost(subtree, cost_map=cost_map)
        else:
            cost_map[(child_name, child_name)] = 0
         
    #print('cost_map', cost_map) 
    return cost_map'''

def cluster_path_len_neg_exp_similarity_matrix(cluster_lineage, similarity=None, node_names=[], first_call=True): # similarity = exp(-(path_length))
    if(similarity is None):
        similarity = {}

    sim_of_parent_child = math.exp(-1) # path_length = 1
        
    node = cluster_lineage.seed_node
    
    if(node.is_leaf()):
        node_name = int(node.taxon.label)
    else:
        node_name = int(node.label)
    node_names.append(node_name)   
    
    similarity[(node_name, node_name)] = 1 # exp(-0) = 1
    
    for child in node.child_nodes():
        if(child.is_leaf()):
            child_name = int(child.taxon.label)
        else:
            child_name = int(child.label)
            
        new_similarity = {}
        
        for state_pair, sim in similarity.items():
            if(state_pair[1] == node_name):
                new_similarity[(state_pair[0], child_name)] = sim * sim_of_parent_child
                new_similarity[(child_name, state_pair[0])] = 0 # descendant -> ancestor similarity = 0 = exp(-inf)
        
        new_similarity[(node_name, child_name)] = sim_of_parent_child # path length of parent -> child = 1
        new_similarity[(child_name, node_name)] = 0 # similarity (child -> parent) = 0 = exp(-inf)
        
        similarity.update(new_similarity)
        
        if not child.is_leaf():
            subtree = dendropy.Tree(seed_node=child.extract_subtree(suppress_unifurcations=False))
            similarity = cluster_path_len_neg_exp_similarity_matrix(subtree, similarity=similarity, node_names=node_names, first_call=False)
        else:
            similarity[(child_name, child_name)] = 1

    #print('cost_map', cost_map) 
    print('nodes on different root-to-leaf paths')
    if (first_call):
        all_node_pairs = list(product(node_names, repeat=2))
        for node_pair in all_node_pairs:
            if not (node_pair in similarity.keys()):
                print(node_pair)
                mrca = cluster_lineage.mrca(node_pair).label
                similarity[node_pair] = 0
                similarity[reversed(node_pair)] = 0
         
    #print('similarity', similarity) 
    return similarity


def sankoff_backward_pass(tree, pointers, assigned_cell_states):
    root_name = tree.seed_node.taxon.label
    print('sankoff_backward_pass', 'root_name', root_name)
    root_state = assigned_cell_states[root_name] - 1 # -1 for 0-based index
    '''if(tree.seed_node.__getattribute__('state') != assigned_cell_states[root_name]):
        #print('Something wrong', tree.seed_node.__getattribute__('state'), assigned_cell_states[root_name])'''
    #print()
    root_pointer = pointers[root_name][root_state] # dict
    print('root_pointer', root_pointer)
    for child in tree.seed_node.child_nodes():
        
        child_name = child.taxon.label
        print('child_name', child_name)
        assigned_cell_states[child_name] = root_pointer[child_name] + 1
        #child.label = root_pointer[child_name] + 1
        #child.__setattr__('state', root_pointer[child_name] + 1) 
        if child.is_internal():
            subtree = dendropy.Tree(seed_node=child.extract_subtree())
            subtree, assigned_cell_states = sankoff_backward_pass(subtree, pointers, assigned_cell_states)
            child.set_child_nodes(subtree.seed_node.child_nodes())
    #print('tree', tree.__str__())
    return tree, assigned_cell_states

def sankoff_forward_pass(tree, cost_map, dp_memory, last_cluster, pointers):
    '''print()
    print('sankoff_forward_pass')
    print('tree', tree)
    print('cost_map', cost_map)
    print('dp_memory', dp_memory)
    print('pointers', pointers)
    print('root', tree.seed_node.taxon.label)'''
    
    for child in tree.seed_node.child_nodes():
        #print(child.taxon.label)
        if not (child.taxon.label in dp_memory.keys()):
            subtree = dendropy.Tree(seed_node=child.extract_subtree())
            dp_memory, pointers = sankoff_forward_pass(subtree, cost_map, dp_memory, last_cluster, pointers)
    #print('dp_memory', dp_memory)
    best_cost = [float('inf')] * last_cluster # best cost per state assignment # state = i
    pointer = [-1] * last_cluster # backtrack pt for each state assignment best cost
    
    for i in range(last_cluster):
        parent_cost = 0
        pointer_i = {}
        #print('i', i)
        for child in tree.seed_node.child_nodes():
            #print('child', child.taxon.label)
            best_child_cost = float('inf')
            child_dp_memory = dp_memory[child.taxon.label]
            #print('child_dp_memory', child_dp_memory)
            for j in range(last_cluster):
                child_cost = child_dp_memory[j]
                if((i, j) in cost_map.keys()):
                    child_cost = child_cost + cost_map[(i, j)]
                else:
                    child_cost = float('inf')
                #print('j', j, 'child_cost', child_cost)
                if(child_cost <= best_child_cost):
                    best_child_cost = child_cost
                    pointer_i[child.taxon.label] = j
                #print('best_child_cost', best_child_cost, pointer_i)
            parent_cost += best_child_cost
        best_cost[i] = parent_cost
        pointer[i] = pointer_i
    #print('best_cost', best_cost)
    #print('pointer[i]', pointer)
    pointers[tree.seed_node.taxon.label] = pointer
    dp_memory[tree.seed_node.taxon.label] = best_cost
    #print()
    return dp_memory, pointers
    
def sankoff_algorithm(in_tree, cell_to_cluster, cluster_lineage):
    in_tree = assign_taxon_name_to_internal_nodes(in_tree)
    in_tree.write(path="internal.nwk", schema="newick")
    #cost_map = calculate_state_change_cost(cluster_lineage)
    similarity = 
    print('cost_map', cost_map) 
    max_cost = max(cost_map.values())
    cost_map = {cluster_pair: cost / max_cost for cluster_pair, cost in cost_map.items()}
    print('cost_map', cost_map) 
    with open('state_change_cost_map.tsv', 'w') as file:
        file.write('State1' + '\t' + 'State2' + '\t' + 'Cost' + '\n')
        for state_pair, cost in cost_map.items():
            file.write(str(state_pair[0]) + '\t' + str(state_pair[1]) + '\t' + str(cost) + '\n')
            file.flush()
    last_cluster = max(cell_to_cluster.values())

    dp_memory = {}
    pointers = {}

    for node in in_tree.leaf_nodes():
        node_name = node.taxon.label
        best_node_cost = [float('inf')] * last_cluster
        best_node_cost[cell_to_cluster[node_name]-1] = 0 # -1 for 0-indexed list
        dp_memory[node_name] = best_node_cost
            
    #print('dp_memory', dp_memory)
    

    dp_memory, pointers = sankoff_forward_pass(in_tree, cost_map, dp_memory, last_cluster, pointers)
    

    assigned_cell_states = {}
    
    root_name = in_tree.seed_node.taxon.label
    root_memeory = dp_memory[root_name]
    best_cost = min(root_memeory)
    #print('Tree Best Cost', best_cost)
    best_cluster = root_memeory.index(best_cost) + 1 # +1 for 0-indexed list
    assigned_cell_states[root_name] = best_cluster
    #in_tree.seed_node.__setattr__('state', best_cluster) 

    
    new_dp_memory = {}
    for key, value in dp_memory.items():
        new_value = []
        for i in range(len(value)):
            if(value[i] == float('inf')):
                new_value.append(-1)
            else:
                new_value.append(value[i])
        #print(value, new_value)
        new_dp_memory[key] = new_value
    #print('dp_memory', dp_memory)
    with open("dp_memory.json", "w") as file:
        json.dump(new_dp_memory, file, indent=4)

    with open("pointers.json", "w") as file:
        json.dump(pointers, file, indent=4)

    resolved_tree, assigned_cell_states = sankoff_backward_pass(in_tree, pointers, assigned_cell_states)
    #print('resolved_tree', type(resolved_tree), resolved_tree.__str__())
    resolved_tree.write(path="resolved.nwk", schema="newick")
    
    
    
    with open("assigned_internal_states.json", "w") as file:
        json.dump(assigned_cell_states, file, indent=4)

    
def resolve_states(in_tree_path, cell_state_path, cell_column, cluster_column, lineage_path):
    in_tree = dendropy.Tree.get(path=in_tree_path, schema="newick", preserve_underscores=True)
    '''print('Input Tree')
    print(in_tree.__str__())
    print()'''

    state_df = pd.read_csv(cell_state_path, sep='\t', usecols=[cell_column, cluster_column])
    cell_to_cluster = dict(zip(state_df[cell_column], state_df[cluster_column]))
    '''print('cell_to_cluster')
    print(cell_to_cluster)
    print()'''

    cluster_lineage = dendropy.Tree.get(path=lineage_path, schema="newick", preserve_underscores=True)
    '''print('cluster_lineage')
    print(cluster_lineage)
    print()'''

    start_time = time.time()
    sankoff_algorithm(in_tree, cell_to_cluster, cluster_lineage)
    print('Runtime', time.time() - start_time)

    
def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    os.chdir(args.out_dir)

    resolve_states(args.in_tree, args.cell_state, args.cell_column, args.cluster_column, args.lineage_path)

if __name__ == "__main__":
    main(parse_args())