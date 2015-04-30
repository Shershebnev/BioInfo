from random import choice as choice
from copy import deepcopy as deepcopy

matrix_size = 0


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


##def adjacency(first_node, second_node):
##    return str(first_node) + '->' + str(second_node)
##
##
##def additive_phylogeny(matrix):
##    global matrix_size, current_node_number
##    if matrix_size == 0:
##        matrix_size = len(matrix)
##        current_node_number = matrix_size + 1
##    n = len(matrix)
##    if n == 2:
##        return {'0->1' : matrix[0][1]}, ['0->1']
##    limb = limb_length(n, matrix)
##    for i in range(len(matrix) - 1):
##        matrix[i][n] -= limb[0]
##        matrix[n][i] = matrix[i][n]
##    ink = [0]*3
##    for i in range(len(matrix) - 1):
##        for n in range(len(matrix) - 1):
##            for k in range(len(matrix) - 1):
##                if matrix[i][k] == (matrix[i][n] + matrix[n][k]) and \
##                   matrix[i][k] * matrix[i][n] * matrix[n][k] != 0:
##                    ink[0] = i
##                    ink[1] = n
##                    ink[2] = k
##    x = matrix[ink[0]][ink[1]]
##    matrix_copy = deepcopy(matrix)
##    matrix_copy.pop(ink[1])
##    for item in matrix_copy:
##        item.pop(ink[1])
##    T = additive_phylogeny(matrix_copy)
##    if len(matrix) == matrix_size:
##        v = matrix_size - 1
##    else:
##        v = current_node_number
##        current_node_number += 1
##    if len(T[0]) == 1:
##        adj = T[1][0]
##        adj_nodes = adj.split('->')
##        weight = list(T[0].items())[0]
##        T[0].clear()
##        T[0][adjacency(adj_nodes[0], v)] = weight - (matrix[0][-2] - matrix[0][-1])
##        T[0][adjacency(adj_nodes[1], v)] = matrix[0][-2] - matrix[0][-1]
##        T[1].clear()
##        T[1].extend(adjacency(adj_nodes[0], v), adjacency(adj_nodes[1], v))
##        if last_nodes[0] < last_nodes[1]:
##            last_node = last_nodes[1]
##        else:
##            last_node = last_nodes[0]
##        T[0][adjacency(last_node, current_node_number)] = x
##        T[1].append(adjacency(last_node, current_node_number))
##    else if len(matrix) != matrix_size:
##        last_added = T[1][-1]
##        last_nodes = last_added.split('->')
##        if last_nodes[0] < last_nodes[1]:
##            last_node = last_nodes[1]
##        else:
##            last_node = last_nodes[0]
##        T[0][adjacency(last_node, v)] = matrix[0][-2] - matrix[0][-1]
##        for key in T[0].keys():
##            if last_node in key:
##                T[0][key] -= T[0][adjacency(last_node, v)]
##        T[1].append(adjacency(last_node, v))
##        


def upgma(matrix):
    clusters = [i for i in range(len(matrix))]  
