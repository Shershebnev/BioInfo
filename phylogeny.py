from random import choice as choice
from copy import deepcopy as deepcopy



def distance_between_leaves(file_input, file_output):
    a = open(file_input, 'r')
    strings = a.readlines()
    a.close()
    n = strings[0].rstrip('\n')
    adj_list = strings[1:]
    adj_list = [item.rstrip('\n') for item in adj_list]
    
    edges = {}
    weights = {}
    
    for item in adj_list:
        try:
            edges[int(item.split('->')[0])].append(int(item.split('->')[1].split(':')[0]))
        except KeyError:
            edges[int(item.split('->')[0])] = [int(item.split('->')[1].split(':')[0])]
        weights[item.split(':')[0]] = int(item.split(':')[1])
        
    leaves = tree_leaves(edges)
    
    matrix = [['0'] * len(leaves) for k in range(len(leaves))]
    
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):            
            tree_copy = deepcopy(edges)
            path = pathfinder(tree_copy, leaves[i], leaves[j])
            weight = 0
            for k in range(len(path) - 1):
                weight += weights[str(path[k]) + '->' + str(path[k + 1])]
            matrix[i][j] = str(weight)
            matrix[j][i] = matrix[i][j]
            
    a = open(file_output, 'a+')
    for i in range(len(matrix)):
        a.write('\t'.join(matrix[i]))
        a.write('\n')
    return 'Success'


def tree_leaves(tree):
    in_nodes = []
    out_nodes = []
    for key, value in tree.items():
        in_nodes.append(key)
        out_nodes.extend(list(value))
    leaves_list = []
    for node in in_nodes:
        if out_nodes.count(node) == 1:
            leaves_list.append(node)
    return leaves_list


def pathfinder(tree, start_edge, end_edge):
    path = []
    stack = [start_edge]
    next_edge = choice(tree[start_edge])
    tree[start_edge].pop(tree[start_edge].index(next_edge))
    tree[next_edge].pop(tree[next_edge].index(start_edge))
    path.append(start_edge)
    if next_edge == end_edge:
        path.append(end_edge)
        return path
    while next_edge != end_edge:
        if not tree[next_edge]:
            stack.pop(-1)
            next_edge = path.pop(-1)
        else:
            stack.append(next_edge)
            path.append(next_edge)
            old_edge = next_edge
            next_edge = choice(tree[next_edge])            
            tree[old_edge].pop(tree[old_edge].index(next_edge))
            tree[next_edge].pop(tree[next_edge].index(old_edge))
    path.append(next_edge)
    return path


def limb_length_from_file(file):
    a = open(file, 'r')
    strings = a.readlines()
    strings = [string.rstrip('\n') for string in strings]
    n = strings.pop(0)
    j = int(strings.pop(0))
    matrix = []
    for item in strings:
        matrix.append([int(weight) for weight in item.split(' ')])
    length = float('inf')
    for i in range(len(matrix)):
        if i != j:
            for k in range(len(matrix)):
                if k != j and i != k:
                    limb = (matrix[i][j] + matrix[j][k] - matrix[i][k]) / 2
                    if limb < length:
                        length = limb
    return int(length)


def limb_length(j, matrix):
    length = float('inf')
    min_i = 0
    min_k = 0
    for i in range(len(matrix)):
        if i != j:
            for k in range(len(matrix)):
                if k != j and i != k:
                    limb = (matrix[i][j] + matrix[j][k] - matrix[i][k]) / 2
                    if limb < length:
                        length = limb
                        min_i = i
                        min_k = k
    return int(length), min_i, min_k



def additive_phylogeny(n, matrix):
    if n == 2:
        return matrix[1][1]
    length_limb = limb_length(n, matrix)
    for j in range(n - 1):
        matrix[j][n] = matrix[j][n] - length_limb[0]
        matrix[n][j] = matrix[j][n]
    i = length_limb[1]
    k = length_limb[2]
    x = matrix[i][n]
    a = matrix.pop(n)
    for row in range(len(matrix)):
        a = matrix[row].pop(n)
    T = additive_phylogeny(n - 1, matrix)
    
    
