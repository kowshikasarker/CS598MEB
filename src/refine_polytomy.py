import argparse
import dendropy
import pandas as pd
import os
import sys
from pathlib import Path
import time
from scipy.spatial.distance import jensenshannon
import math
from itertools import product

cost_matrix = None

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_tree", type=str,
                        help="path to the input tree in newick format",
                        required=True, default=None)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="path to the output directory",
                        required=True, default=None)
    parser.add_argument("-e", "--expr_state", type=str,
                        help="path to cell to cluster mapping in tsv format (column for cell = 'cell id', column for cluster id = 'ClusterIdent')",
                        required=True, default=None)
    parser.add_argument("-l", "--lineage_path", type=str,
                        help="path to the cluster lineage in newick format",
                        required=False, default=None)
    parser.add_argument("-c", "--criterion", type=str,
                        help="name of the criterion to compute cost between subtree pairs",
                        required=True, default=None)
    parser.add_argument("-t", "--threshold", type=float,
                        help="criterion threshold to refine polytomy",
                        required=True, default=None)
    
    args = parser.parse_args()
    return args

def calulcate_cluster_distribution(leafset, cell_to_cluster):
    max_cluster_no = max(cell_to_cluster.values())
    counts = [0 for i in range(max_cluster_no)]
    
    for leaf in leafset:
        cluster = cell_to_cluster[leaf]
        counts[cluster-1] += 1
    
    total = sum(counts)
    if(total != len(leafset)):
        raise Exception('Something went wrong')
    
    distribution = []
    for i in range(max_cluster_no):
        distribution.append(counts[i] / len(leafset))
        
    return distribution

def jensen_shannon_cost(leafset1, leafset2, cell_to_cluster):
    p = calulcate_cluster_distribution(leafset1, cell_to_cluster)
    q = calulcate_cluster_distribution(leafset2, cell_to_cluster)
    return jensenshannon(p, q, base=2)

def cluster_overlap_cost(leafset1, leafset2, cell_to_cluster):
    p = calulcate_cluster_distribution(leafset1, cell_to_cluster)
    q = calulcate_cluster_distribution(leafset2, cell_to_cluster)
    similarity = 0
    for i in range(len(p)):
        similarity += min(p[i], q[i])
    cost = 1 - similarity
    return cost

def state_change_cost(leafset1, leafset2, cell_to_cluster):
    leaf_pairs = list(product(leafset1, leafset2))
    total_cost = 0
    for leaf_pair in leaf_pairs:
        cluster_1 = cell_to_cluster[leaf_pair[0]]
        cluster_2 = cell_to_cluster[leaf_pair[1]]
        if (cluster_1 != cluster_2):
            total_cost += 1
    avg_cost = total_cost / len(leaf_pairs)
    return avg_cost

def state_distance_matrix(cluster_lineage, distance=None):
    if(distance is None):
        distance = {}
        
    node = cluster_lineage.seed_node
    
    if(node.is_leaf()):
        node_name = int(node.taxon.label)
    else:
        node_name = int(node.label)
    
    distance[(node_name, node_name)] = 0
    
    for child in node.child_nodes():
        if(child.is_leaf()):
            child_name = int(child.taxon.label)
        else:
            child_name = int(child.label)
            
        new_distance = {}
        
        for state_pair, dist in distance.items():
            if(state_pair[1] == node_name):
                new_distance[(state_pair[0], child_name)] = dist + 1
        
        new_distance[(node_name, child_name)] = 1 
        
        distance.update(new_distance)
        
        if child.is_leaf():
            distance[(child_name, child_name)] = 0
        else:
            subtree = dendropy.Tree(seed_node=child.extract_subtree(suppress_unifurcations=False))
            distance.update(state_distance_matrix(subtree, distance=distance))

    return distance

def scaled_sigmoid(x):
  return 2 * (1 / (1 + math.exp(-x)) - 0.5)

def lineage_linear_cost(leafset1, leafset2, cell_to_cluster, cluster_lineage, inf_value=2.0):
    global cost_matrix
    if cost_matrix is None:
        cost_matrix = {}
        distance = state_distance_matrix(cluster_lineage)
        max_dist = max(distance.values()) 
        states = list(range(1, max(cell_to_cluster.values())+1))
        state_pairs = list(product(states, repeat=2))
        for state_pair in state_pairs:
            if not state_pair in distance.keys():
                cost_matrix[state_pair] = inf_value
            else:
                cost_matrix[state_pair] = distance[state_pair] / max_dist
    print('cost_matrix', cost_matrix)
    
    leaf_pairs = list(product(leafset1, leafset2))
    total_cost = 0
    for leaf_pair in leaf_pairs:
        cluster_pair = (cell_to_cluster[leaf_pair[0]], cell_to_cluster[leaf_pair[1]]) 
        total_cost += cost_matrix[cluster_pair]
    avg_cost = total_cost / len(leaf_pairs)
    return avg_cost


def lineage_sigmoid_cost(leafset1, leafset2, cell_to_cluster, cluster_lineage):
    global cost_matrix
    
    if cost_matrix is None:
        cost_matrix = {}
        distance = state_distance_matrix(cluster_lineage)
        states = list(range(1, max(cell_to_cluster.values())+1))
        state_pairs = list(product(states, repeat=2))
        
        for state_pair in state_pairs:
            if not state_pair in distance.keys():
                cost_matrix[state_pair] = 1
            else:
                cost_matrix[state_pair] = scaled_sigmoid(distance[state_pair])
    
    leaf_pairs = list(product(leafset1, leafset2))
    total_cost = 0
    
    for leaf_pair in leaf_pairs:
        cluster_pair = (cell_to_cluster[leaf_pair[0]], cell_to_cluster[leaf_pair[1]]) 
        total_cost += cost_matrix[cluster_pair]
    avg_cost = total_cost / len(leaf_pairs)
    return avg_cost


def cost_between_leafset_pairs(leafset1, leafset2, criterion, cell_to_cluster, cluster_lineage):
    # c1 = 'jensen_shannon_cost'
    # c2 = 'cluster_overlap_cost'
    # c3 = 'state_change_cost'
    # c4 = 'lineage_linear_cost'
    # c5 = 'lineage_sigmoid_cost'
    if(criterion == 'c1'):
        return jensen_shannon_cost(leafset1, leafset2, cell_to_cluster)
    if(criterion == 'c2'):
        return cluster_overlap_cost(leafset1, leafset2, cell_to_cluster)
    if(criterion == 'c3'):
        return state_change_cost(leafset1, leafset2, cell_to_cluster)
    if(criterion == 'c4'):
        return lineage_linear_cost(leafset1, leafset2, cell_to_cluster, cluster_lineage)
    if(criterion == 'c5'):
        return lineage_sigmoid_cost(leafset1, leafset2, cell_to_cluster, cluster_lineage)
    else:
        raise Exception('Unrecognized criterion', criterion)
    
def recursively_refine_polytomy(tree, cell_to_cluster, cluster_lineage, criterion, threshold, criterion_dict=None, save_criterion=True, refined_dict=None, save_refined=True):
    if (tree.seed_node.is_leaf()):
        return None

    if criterion_dict is None:
        criterion_dict = {
            'Max_Subtree_Size': [],
            'Min_Subtree_Size': [],
            'Max_Subtree': [],
            'Min_Subtree': [],
            criterion: []
        }

    if refined_dict is None:
        refined_dict = {
            'Unrefined_Degree': [],
            'Refined_Degree': [],
            'Refined_Topology': []
        }

    child_nodes = tree.seed_node.child_nodes()
    for i in range(len(child_nodes)):
        child = child_nodes[i]
        child.label = ''
        subtree = dendropy.Tree(seed_node=child.extract_subtree())
        refined = recursively_refine_polytomy(subtree, cell_to_cluster, cluster_lineage, criterion, threshold, criterion_dict=criterion_dict, save_criterion=False, refined_dict=refined_dict, save_refined=False)
        if not child.is_leaf():
            child.set_child_nodes(refined.seed_node.child_nodes())

    unrefined_degree = len(child_nodes)
    while(True):
        flag = False # if any refinement happened

        children = tree.seed_node.child_nodes()
        
        if(len(children) < 3):
            #print('<3 children.')
            break
            
        leafsets = []
        topologies = []
        
        for i in range(len(children)):
            topologies.append(dendropy.Tree(seed_node=children[i].extract_subtree()).__str__())
            leafset = [leaf.taxon.label for leaf in children[i].leaf_nodes()]
            leafsets.append(leafset)
        #print('leafsets', leafsets)
            
        pairwise_cost = {}
        for i in range(len(children)):
            for j in range(i+1, len(children)):
                cost = cost_between_leafset_pairs(leafsets[i], leafsets[j], criterion, cell_to_cluster, cluster_lineage)
                pairwise_cost[(i, j)] = cost
                if(len(leafsets[i]) > len(leafsets[j])):
                    criterion_dict['Max_Subtree_Size'].append(len(leafsets[i]))
                    criterion_dict['Max_Subtree'].append(topologies[i])
                    criterion_dict['Min_Subtree_Size'].append(len(leafsets[j]))
                    criterion_dict['Min_Subtree'].append(topologies[j])
                else:
                    criterion_dict['Max_Subtree_Size'].append(len(leafsets[j]))
                    criterion_dict['Max_Subtree'].append(topologies[j])
                    criterion_dict['Min_Subtree_Size'].append(len(leafsets[i]))
                    criterion_dict['Min_Subtree'].append(topologies[i])                
                criterion_dict[criterion].append(cost)
        
        while(True):
            if not pairwise_cost:
                break
            closest_pair = min(pairwise_cost, key=pairwise_cost.get)
            min_cost = pairwise_cost[closest_pair]
            if(min_cost <= threshold):
                flag = True
                new_node = dendropy.Node()
                new_node.parent_node = tree.seed_node
                children[closest_pair[0]].parent_node = new_node
                children[closest_pair[1]].parent_node = new_node
                new_node.label = '(' + children[closest_pair[0]].label + ',' + children[closest_pair[1]].label + ')'
                children[closest_pair[0]].label = None
                children[closest_pair[1]].label = None
                valid_pairs = [pair for pair in pairwise_cost.keys() if not ((closest_pair[0] in pair) or (closest_pair[1] in pair))]  
                pairwise_cost = {key: value for key, value in pairwise_cost.items() if key in valid_pairs}
            else:
                break
        if not flag:
            break
     
    child_nodes = tree.seed_node.child_nodes()
    refined_degree = len(child_nodes)
    if(unrefined_degree > refined_degree):
        refined_topology = '('
        for child in child_nodes:
            refined_topology += child.label + ','
            child.label = None       
        refined_topology = refined_topology[:-1] + ')' # -1 to remove the last comma
    else:
        refined_topology = '-'
    refined_dict['Unrefined_Degree'].append(unrefined_degree)
    refined_dict['Refined_Degree'].append(refined_degree)
    refined_dict['Refined_Topology'].append(refined_topology)
    if save_criterion:
        criterion_df = pd.DataFrame.from_dict(criterion_dict)
        criterion_df.to_csv(criterion + '.tsv', sep='\t', index=False)
    if save_refined:
        refined_df = pd.DataFrame.from_dict(refined_dict)
        refined_df.to_csv('refined.tsv', sep='\t', index=False)

    return tree

def refine_polytomy(tree_path, state_path, cluster_lineage_path, criterion, threshold):    
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    
    state_df = pd.read_csv(state_path, sep='\t', usecols=['cell_id', 'ClusterIdent'])
    cell_to_cluster = dict(zip(state_df['cell_id'], state_df['ClusterIdent']))
    
    cluster_lineage = dendropy.Tree.get(path=cluster_lineage_path, schema="newick", preserve_underscores=True)
    
    refined_tree = recursively_refine_polytomy(tree, cell_to_cluster, cluster_lineage, criterion, threshold)
    refined_tree_nwk = refined_tree.__str__().replace("'", "") + ";"
    with open("refined.nwk", "w") as outfile:
        outfile.write(refined_tree_nwk)
        outfile.flush()
    #refined_tree.write(path="refined.nwk", schema="newick")

     
def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    os.chdir(args.out_dir)
    
    log_file = open('log.txt', 'w')
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = log_file
    sys.stderr = log_file

    print(args)

    start = time.process_time()
    
    refine_polytomy(args.in_tree, args.expr_state, args.lineage_path, args.criterion, args.threshold)

    print('Runtime', time.process_time() - start)

    log_file.flush()
    log_file.close()
    sys.stderr = original_stderr
    sys.stdout = original_stdout
    
if __name__ == "__main__":
    main(parse_args())
