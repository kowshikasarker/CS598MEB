import argparse
import dendropy
import pandas as pd
from scipy.stats import entropy
import os
import sys
from pathlib import Path
import time
from scipy.spatial.distance import jensenshannon
import itertools
import math
from itertools import product

state_change_matrix = None

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
    parser.add_argument("-c", "--criterion", type=str,
                        help="name of the criterion to compute distance/similarity between subtree pairs",
                        required=True, default=None)
    parser.add_argument("-th", "--criterion_threshold", type=float,
                        help="distance threshold to refine polytomy",
                        required=True, default=None)
    
    args = parser.parse_args()
    return args

def calulcate_cluster_distribution(leafset, cell_to_cluster):
    max_cluster_no = max(cell_to_cluster.values())
    counts = [0 for i in range(max_cluster_no)]
    
    for leaf in leafset:
        #print(type(node.taxon.label))
        cluster = cell_to_cluster[leaf]
        counts[cluster-1] += 1
    
    total = sum(counts)
    if(total != len(leafset)):
        print('Something went wrong')
    
    distribution = []
    for i in range(max_cluster_no):
        distribution.append(counts[i] / len(leafset))
        
    #print('counts', counts)
    print('distribution', distribution)
    return distribution

def jensen_shannon_distance(leafset1, leafset2, cell_to_cluster):
    print('leafset1', leafset1)
    print('leafset2', leafset2)
    p = calulcate_cluster_distribution(leafset1, cell_to_cluster)
    q = calulcate_cluster_distribution(leafset2, cell_to_cluster)
    return jensenshannon(p, q, base=2)

def cluster_overlap_distance(leafset1, leafset2, cell_to_cluster):
    p = calulcate_cluster_distribution(leafset1, cell_to_cluster)
    q = calulcate_cluster_distribution(leafset2, cell_to_cluster)
    similarity = 0
    for i in range(len(p)):
        similarity += min(p[i], q[i])
    distance = 1 - similarity
    return distance


def cluster_path_len_distance_matrix(cluster_lineage, distance=None, node_names=[], first_call=True):
    if(distance is None):
        distance = {}
    
    node = cluster_lineage.seed_node
    if(node.is_leaf()):
        node_name = int(node.taxon.label)
    else:
        node_name = int(node.label)
    node_names.append(node_name)
    
    #print('node', node, node.label)
    distance[(node_name, node_name)] = 0
    for child in node.child_nodes():
        #print('child', child, child.label)
        if(child.is_leaf()):
            child_name = int(child.taxon.label)
        else:
            child_name = int(child.label)
        new_distance = {}
        for state_pair, cost in distance.items():
            if(state_pair[1] == node_name):
                new_distance[(int(state_pair[0]), int(child_name))] = cost + 1
                new_distance[(int(child_name), int(state_pair[0]))] = cost + 1
        new_distance[(int(node_name), int(child_name))] = 1
        new_distance[(int(child_name), int(node_name))] = 1
        distance.update(new_distance)
        
        if not child.is_leaf():
            #print(child)
            #print(cluster_lineage)
            subtree = dendropy.Tree(seed_node=child.extract_subtree(suppress_unifurcations=False))
            distance = cluster_path_len_distance_matrix(subtree, distance=distance, node_names=node_names, first_call=False)
        else:
            distance[(child_name, child_name)] = 0
         
    #print('distance', distance) 
    print('nodes on different root-to-leaf paths')
    if (first_call):
        all_node_pairs = list(product(node_names, repeat=2))
        for node_pair in all_node_pairs:
            if not (node_pair in distance.keys()):
                print(node_pair)
                mrca = cluster_lineage.mrca(node_pair).label
                distance[node_pair] = distance[(mrca, node_pair[0])] + distance[(mrca, node_pair[1])]
                distance[reversed(node_pair)] = distance[node_pair]
                    
    return distance


def cluster_path_len_distance(leafset1, leafset2, cell_to_cluster, cluster_lineage):
    global state_change_matrix
    if(state_change_matrix is None):
        dist_matrix = cluster_path_len_distance_matrix(cluster_lineage)
        max_dist = max(dist_matrix.values())
        state_change_matrix = {cluster_pair: dist / max_dist for cluster_pair, dist in dist_matrix.items()}
        print('state_change_matrix', state_change_matrix)
    
    leaf_pairs = list(itertools.product(leafset1, leafset2))
    total_dist = 0
    for leaf_pair in leaf_pairs:
        cluster_pair = (cell_to_cluster[leaf_pair[0]], cell_to_cluster[leaf_pair[1]]) 
        total_dist += state_change_matrix[cluster_pair]
    avg_dist = total_dist / len(leaf_pairs)
    print('total_dist', total_dist, 'avg_dist', avg_dist)
    return avg_dist

def cell_path_distance(leafset1, leafset2, cell_lineage):
    pass
    
def cell_pseudotemporal_distance(leafset1, leafset2, cell_pseudotime):
    pass

def cluster_pseudotemporal_distance(leafset1, leafset2, cell_pseudotime):
    pass

def distance_between_leafset_pairs(leafset1, leafset2, distance, cell_to_cluster, cluster_lineage):
    # d1 = 'jensen_shannon'
    # d2 = 'cluster_overlap'
    # d3 = 'cluster_path'
    if(distance == 'd1'):
        return jensen_shannon_distance(leafset1, leafset2, cell_to_cluster)
    if(distance == 'd2'):
        return cluster_overlap_distance(leafset1, leafset2, cell_to_cluster)
    if(distance == 'd3'):
        return cluster_path_len_distance(leafset1, leafset2, cell_to_cluster, cluster_lineage)

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


def cluster_path_len_neg_exp_similarity(leafset1, leafset2, cell_to_cluster, cluster_lineage):
    global state_change_matrix
    if state_change_matrix is None:
        state_change_matrix = cluster_path_len_neg_exp_similarity_matrix(cluster_lineage)
    print('cluster_path_len_neg_exp_similarity', state_change_matrix)
    

    leaf_pairs = list(itertools.product(leafset1, leafset2))
    total_sim = 0
    for leaf_pair in leaf_pairs:
        cluster_pair = (cell_to_cluster[leaf_pair[0]], cell_to_cluster[leaf_pair[1]]) 
        total_sim += state_change_matrix[cluster_pair]
    avg_sim = total_sim / len(leaf_pairs)
    print('total_sim', total_sim, 'avg_sim', avg_sim)
    return avg_sim

def similarity_between_leafset_pairs(leafset1, leafset2, similarity, cell_to_cluster, cluster_lineage):
    # s1 = 'cluster_similarity'
    if(similarity == 's1'):
        return cluster_path_len_neg_exp_similarity(leafset1, leafset2, cell_to_cluster, cluster_lineage)
    
def recursively_refine_polytomy_using_distance(tree, cell_to_cluster, cluster_lineage, distance, threshold, distance_dict=None, save_distance=True):
    #print('Function: recursively_refine_polytomy_using_distance', tree.__str__())
    if(tree.seed_node.is_leaf() == 1):
        #print('Leaf node.')
        return None

    if distance_dict is None:
        distance_dict = {
            'Max_Subtree_Size': [],
            'Min_Subtree_Size': [],
            'Max_Subtree': [],
            'Min_Subtree': [],
            'Distance': []
        }

    for child in tree.seed_node.child_nodes():
        subtree = dendropy.Tree(seed_node=child.extract_subtree())
        refined = recursively_refine_polytomy_using_distance(subtree, cell_to_cluster, cluster_lineage, distance, threshold, distance_dict=distance_dict, save_distance=False)
        if not refined is None:
            child.set_child_nodes(refined.seed_node.child_nodes())
    

    print('tree.seed_node.child_nodes()', tree.seed_node.child_nodes())
    print('tree', tree.__str__())
    

    while(True):
        flag = False # if any refinement happened

        children = tree.seed_node.child_nodes()
        print('children', children)
        
        if(len(children) < 3):
            print('<3 children.')
            break
            
        leafsets = []
        topologies = []
        
        for i in range(len(children)):
            topologies.append(dendropy.Tree(seed_node=children[i].extract_subtree()).__str__())
            leafset = [leaf.taxon.label for leaf in children[i].leaf_nodes()]
            leafsets.append(leafset)
        print('leafsets', leafsets)
            
        pairwise_distance = {}
        for i in range(len(children)):
            for j in range(i+1, len(children)):
                dist = distance_between_leafset_pairs(leafsets[i], leafsets[j], distance, cell_to_cluster, cluster_lineage)
                pairwise_distance[(i, j)] = dist
                if(len(leafsets[i]) > len(leafsets[j])):
                    distance_dict['Max_Subtree_Size'].append(len(leafsets[i]))
                    distance_dict['Max_Subtree'].append(topologies[i])
                    distance_dict['Min_Subtree_Size'].append(len(leafsets[j]))
                    distance_dict['Min_Subtree'].append(topologies[j])
                else:
                    distance_dict['Max_Subtree_Size'].append(len(leafsets[j]))
                    distance_dict['Max_Subtree'].append(topologies[j])
                    distance_dict['Min_Subtree_Size'].append(len(leafsets[i]))
                    distance_dict['Min_Subtree'].append(topologies[i])                
                distance_dict['Distance'].append(dist)
                print(distance_dict)
        print('pairwise_distance', pairwise_distance)
        
        while(True):
            if not pairwise_distance:
                break
            closest_pair = min(pairwise_distance, key=pairwise_distance.get)
            print('closest_pair', closest_pair, children[closest_pair[0]], children[closest_pair[1]])
            min_dist = pairwise_distance[closest_pair]
            print('min_dist', min_dist)
            if(min_dist <= threshold):
                flag = True
                new_node = dendropy.Node()
                new_node.parent_node = tree.seed_node
                children[closest_pair[0]].parent_node = new_node
                children[closest_pair[1]].parent_node = new_node
                valid_pairs = [pair for pair in pairwise_distance.keys() if not ((closest_pair[0] in pair) or (closest_pair[1] in pair))]  
                pairwise_distance = {key: value for key, value in pairwise_distance.items() if key in valid_pairs}
            else:
                break
        print('One level refinement complete')
        print('tree', tree.__str__())
        if not flag:
            break
    if save_distance:
        distance_df = pd.DataFrame.from_dict(distance_dict)
        distance_df.to_csv(distance + '.tsv', sep='\t', index=False)

    return tree

def recursively_refine_polytomy_using_similarity(tree, cell_to_cluster, cluster_lineage, similarity, threshold, similarity_dict=None, save_similarity=True):
    #print('Function: recursively_refine_polytomy_using_similarity', tree.__str__())
    if (tree.seed_node.is_leaf()):
        return None

    if similarity_dict is None:
        similarity_dict = {
            'Max_Subtree_Size': [],
            'Min_Subtree_Size': [],
            'Max_Subtree': [],
            'Min_Subtree': [],
            'Similarity': []
        }

    for child in tree.seed_node.child_nodes():
        subtree = dendropy.Tree(seed_node=child.extract_subtree())
        refined = recursively_refine_polytomy_using_similarity(subtree, cell_to_cluster, cluster_lineage, similarity, threshold, similarity_dict=similarity_dict, save_similarity=False)
        if not child.is_leaf():
            child.set_child_nodes(refined.seed_node.child_nodes())
    

    print('tree.seed_node.child_nodes()', tree.seed_node.child_nodes())
    print('tree', tree.__str__())
    

    while(True):
        flag = False # if any refinement happened

        children = tree.seed_node.child_nodes()
        print('children', children)
        
        if(len(children) < 3):
            print('<3 children.')
            break
            
        leafsets = []
        topologies = []
        
        for i in range(len(children)):
            topologies.append(dendropy.Tree(seed_node=children[i].extract_subtree()).__str__())
            leafset = [leaf.taxon.label for leaf in children[i].leaf_nodes()]
            leafsets.append(leafset)
        print('leafsets', leafsets)
            
        pairwise_similarity = {}
        for i in range(len(children)):
            for j in range(i+1, len(children)):
                sim = similarity_between_leafset_pairs(leafsets[i], leafsets[j], similarity, cell_to_cluster, cluster_lineage)
                pairwise_similarity[(i, j)] = sim
                if(len(leafsets[i]) > len(leafsets[j])):
                    similarity_dict['Max_Subtree_Size'].append(len(leafsets[i]))
                    similarity_dict['Max_Subtree'].append(topologies[i])
                    similarity_dict['Min_Subtree_Size'].append(len(leafsets[j]))
                    similarity_dict['Min_Subtree'].append(topologies[j])
                else:
                    similarity_dict['Max_Subtree_Size'].append(len(leafsets[j]))
                    similarity_dict['Max_Subtree'].append(topologies[j])
                    similarity_dict['Min_Subtree_Size'].append(len(leafsets[i]))
                    similarity_dict['Min_Subtree'].append(topologies[i])                
                similarity_dict['Similarity'].append(sim)
                print(similarity_dict)
        print('pairwise_similarity', pairwise_similarity)
        
        while(True):
            if not pairwise_similarity:
                break
            closest_pair = max(pairwise_similarity, key=pairwise_similarity.get)
            print('closest_pair', closest_pair, children[closest_pair[0]], children[closest_pair[1]])
            max_sim = pairwise_similarity[closest_pair]
            print('max_sim', max_sim)
            if(max_sim >= threshold):
                flag = True
                new_node = dendropy.Node()
                new_node.parent_node = tree.seed_node
                children[closest_pair[0]].parent_node = new_node
                children[closest_pair[1]].parent_node = new_node
                valid_pairs = [pair for pair in pairwise_similarity.keys() if not ((closest_pair[0] in pair) or (closest_pair[1] in pair))]  
                pairwise_similarity = {key: value for key, value in pairwise_similarity.items() if key in valid_pairs}
            else:
                break
        print('One level refinement complete')
        print('tree', tree.__str__())
        if not flag:
            break
    if save_similarity:
        similarity_df = pd.DataFrame.from_dict(similarity_dict)
        similarity_df.to_csv(similarity + '.tsv', sep='\t', index=False)

    return tree

def recursively_refine_polytomy(tree, cell_to_cluster, cluster_lineage, criterion, threshold):
    if(criterion.startswith('s')):
        return recursively_refine_polytomy_using_similarity(tree, cell_to_cluster, cluster_lineage, criterion, threshold)
    if(criterion.startswith('d')):
        return recursively_refine_polytomy_using_distance(tree, cell_to_cluster, cluster_lineage, criterion, threshold)
    print('Criterion', criterion, 'unrecognized')
    return None

def refine_polytomy(tree_path, state_path, cell_col, cluster_col, cluster_lineage_path, criterion, threshold):
    print('Function: refine_polytomy')
    
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    
    print('Input Tree')
    print(tree.__str__())
    print()

    state_df = pd.read_csv(state_path, sep='\t', usecols=[cell_col, cluster_col])
    cell_to_cluster = dict(zip(state_df[cell_col], state_df[cluster_col]))
    print('cell_to_cluster')
    print(cell_to_cluster)
    print()
    
    cluster_lineage = dendropy.Tree.get(path=cluster_lineage_path, schema="newick", preserve_underscores=True)
    print('cluster_lineage')
    print(cluster_lineage)
    print()

    refined_tree = recursively_refine_polytomy(tree, cell_to_cluster, cluster_lineage, criterion, threshold)
    refined_tree.write(path="refined.nwk", schema="newick")
    print('refined_tree')
    print(refined_tree.__str__())
    print()

    print('criterion', criterion)
    print('threshold', threshold)
     
def main(args):
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    os.chdir(args.out_dir)
    
    log_file = open('log.txt', 'w')
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = log_file
    sys.stderr = log_file

    start = time.process_time()
    
    refine_polytomy(args.in_tree, args.cell_state, args.cell_column, args.cluster_column, args.lineage_path, args.criterion, args.criterion_threshold)

    print('Runtime', time.process_time() - start)

    log_file.flush()
    log_file.close()
    sys.stderr = original_stderr
    sys.stdout = original_stdout
    
if __name__ == "__main__":
    main(parse_args())