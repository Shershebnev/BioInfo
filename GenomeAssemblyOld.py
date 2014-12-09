def StringComposition(file_input, file_output, k):

    '''
    Create a list of all kmers in the string in lexicographic order
    '''
    
    genome_file = open(file_input, 'r')
    genome = genome_file.readline()
    genome = genome.rstrip('\n')
    
    strings = open(file_output, 'w')

    kmers = []
    
    for i in range(len(genome) - k + 1):
        kmers.append(genome[i : i + k])

    kmers.sort()

    for item in kmers:
         c = strings.write(item + '\n')

    genome_file.close()
    strings.close()

    return 'Mission accomplished!'


def OverlapGraph(file_input, file_output):

    '''
    Generate a list of overlaped kmers. Overlapping occur if suffix of on kmer == prefix
    of another kmer. Prefix = kmer[ : len(kmer) - 1], suffix = kmer[ 1 : ]
    Output in file in form kmer -> kmer
    '''
    
    pattern_list = []

    a = open(file_input, 'r')
    pattern = a.readline()
    pattern = pattern.rstrip('\n')
    
    while pattern != '':
        pattern_list.append(pattern)
        pattern = a.readline()
        pattern = pattern.rstrip('\n')

    prefix = {}
    suffix = {}
    overlap = {}
    overlap_list = []

    # generate dictionaries of prefixes and suffixes
    for item in pattern_list:
        prefix[item[ : len(item) - 1]] = item[-1]
        suffix[item[1 : ]] = item[0]

    # look for overlapping. If found, add to overlap dict with keys being kmers with suffix
    # and values being kmers with prefix
    for item in suffix:
        if item in prefix:
            overlap[suffix[item] + item] = item + prefix[item]

    # for sorting
    for key in overlap:
        overlap_list.append(key)

    overlap_list.sort()
    
    b = open(file_output, 'w')
    
    for item in overlap_list:
        c = b.write(item + ' -> ' + overlap[item] + '\n')

    a.close()
    b.close()
    return 'Mission accomplished!'


def DeBruijnGraph(string, file_output, k):
    
    bruijn = {}
    bruijn_list = []
    
    for i in range(len(string) - k + 1):
        kmer = string[i : i + k - 1]
        if not kmer in bruijn:
            bruijn[kmer] = [string[i + 1 : i + k]]
        else:
            bruijn[kmer].append(string[i + 1 : i + k])
    
    for item in bruijn:
        bruijn[item].sort()
        bruijn_list.append(item)

    bruijn_list.sort()

    a = open(file_output, 'w')
    
    for item in bruijn_list:
        string = ''
        for value in bruijn[item]:
            string += value + ','
        string = string.rstrip(',')
        c = a.write(item + ' -> ' + string + '\n')

    a.close()
    return 'Mission accomplished!'


def DeBruijnGraphFromKmer(file_input, file_output):
    
    pattern_list = []

    a = open(file_input, 'r')
    pattern = a.readline()
    pattern = pattern.rstrip('\n')
    
    while pattern != '':
        pattern_list.append(pattern)
        pattern = a.readline()
        pattern = pattern.rstrip('\n')
    
    k = len(pattern_list[0])
    
    pref_suf = {}

    for item in pattern_list:
        if not item[:k - 1] in pref_suf:
            pref_suf[item[:k - 1]] = [item[1:]]
        else:
            pref_suf[item[:k - 1]].append(item[1:])

    prefix = []

    for key in pref_suf:
        prefix.append(key)

    prefix.sort()

    b = open(file_output, 'w')
    
    for value in prefix:
        string = ''
        for item in pref_suf[value]:
            string += item + ','
        string = string.rstrip(',')
        c = b.write(value + ' -> ' + string + '\n')

    a.close()
    b.close()
    return 'Mission accomplished!'


##def EulerianCycle(file):
##
##    a = open(file, 'r')
##    b = a.readline()
##    b = b.rstrip('\n')
##    dict_of_paths = {}
##    while b:
##        first_parse = b.split()
##        dict_of_paths[int(first_parse[0])] = []
##        
##        second_parse = first_parse[-1].split(',')
##        for item in second_parse:
##            dict_of_paths[int(first_parse[0])].append(int(item))
##
##        b = a.readline()
##        b = b.rstrip('\n')
##
##    a.close()
##    
##    from random import choice as choice
##    from copy import deepcopy as deepcopy    
##
##    keys = dict_of_paths.keys()
##    nodes = [key for key in keys]
##
##    node = choice(nodes)
##    path = [node]
##    dict_of_paths_copy = deepcopy(dict_of_paths)
##    '''
##    if len(dict_of_paths_copy[node]) > 1:
##        c = dict_of_paths_copy[node].pop(index(node))
##    else:
##        c = dict_of_paths_copy.pop(node)
##    '''
##    
##    while dict_of_paths_copy:
##        if node in dict_of_paths_copy and len(dict_of_paths_copy[node]) > 1:
##            next_node = choice(dict_of_paths_copy[node])
##            path += '->' + str(next_node)
##            c = dict_of_paths_copy[node].pop(dict_of_paths_copy[node].index(next_node))
##            node = next_node
##            print(path + ' node is ' + str(node))
##        elif node in dict_of_paths_copy and len(dict_of_paths_copy[node]) == 1:
##            next_node = dict_of_paths_copy[node][0]
##            c = dict_of_paths_copy.pop(node)
##            path += '->' + str(next_node)
##            node = next_node
##            print(path + ' node is ' + str(node))
##        elif not node in dict_of_paths_copy:
##            dict_of_paths_copy = deepcopy(dict_of_paths)
##            node = choice(nodes)
##            path = str(node)
##            print(path + ' node is ' + str(node))
##    
##    return path
##
##
##def find_path(graph, start, end, path=[]):
##        path = path + [start]
##        if start == end:
##            return path
##        if not start in graph:
##            return None
##        for node in graph[start]:
##            if node not in path:
##                newpath = find_path(graph, node, end, path)
##                if newpath: return newpath
##        return None
##
##def find_all_paths(graph, start, end, path=[]):
##        path = path + [start]
##        if start == end:
##            return [path]
##        if not start in graph:
##            return []
##        paths = []
##        for node in graph[start]:
##            if node not in path:
##                newpaths = find_all_paths(graph, node, end, path)
##                for newpath in newpaths:
##                    paths.append(newpath)
##        return paths
##
##
##
##def EulerianCycle2(file):
##
##    a = open(file, 'r')
##    b = a.readline()
##    b = b.rstrip('\n')
##    dict_of_paths = {}
##    while b:
##        first_parse = b.split()
##        dict_of_paths[int(first_parse[0])] = []
##        
##        second_parse = first_parse[-1].split(',')
##        for item in second_parse:
##            dict_of_paths[int(first_parse[0])].append(int(item))
##
##        b = a.readline()
##        b = b.rstrip('\n')
##
##    a.close()
##    
##    from random import choice as choice
##    from copy import deepcopy as deepcopy
##
##    forks = {}
##    list_of_forks = []
##    steps = 0
##
##    keys = dict_of_paths.keys()
##    nodes = [key for key in keys]
##
##    node = choice(nodes)
##    path = [node]
##    dict_of_paths_copy = deepcopy(dict_of_paths)
##
##    while dict_of_paths_copy:
##        if node in dict_of_paths_copy and len(dict_of_paths_copy[node]) == 1:
##            steps += 1
##            next_node = dict_of_paths_copy[node][0]
##            c = dict_of_paths_copy.pop(node)
##            node = next_node
##            path.append(node)
##            print('path is ',path, 'fork is ', forks)
##        elif node in dict_of_paths_copy and len(dict_of_paths_copy[node]) > 1:
##            next_node = choice(dict_of_paths_copy[node])
##            if not node in forks:
##                forks[node] = [[next_node], [steps], deepcopy(dict_of_paths_copy)]
##                list_of_forks.append(node)
##            else:
##                forks[node][0].append(node)
##                forks[node][1].append(steps)
##                list_of_forks.append(node)
##            path.append(next_node)
##            c = dict_of_paths_copy[node].pop(dict_of_paths_copy[node].index(next_node))
##            node = next_node
##            steps = 0
##            print('path is ',path, 'fork is ', forks)
##        elif not node in dict_of_paths_copy:
##            dict_of_paths_copy = forks[list_of_forks[-2]][2]
##            steps = 0
##            next_node = choice(dict_of_paths_copy[list_of_forks[-2]])
##            while not next_node in forks[list_of_forks[-2]][0]:
##                next_node = choice(dict_of_paths_copy[list_of_forks[-2]])
##            node = next_node
##            for i in range(forks[list_of_forks[-2]][1][-1]):
##                path.pop()
##            print('path is ',path, 'fork is ', forks)
##    return path
##
def GraphForInt(file):

    a = open(file, 'r')
    b = a.readline()
    b = b.rstrip('\n')
    dict_of_paths = {}
    while b:
        first_parse = b.split()
        dict_of_paths[int(first_parse[0])] = []
        
        second_parse = first_parse[-1].split(',')
        for item in second_parse:
            dict_of_paths[int(first_parse[0])].append(int(item))

        b = a.readline()
        b = b.rstrip('\n')

    a.close()
    return dict_of_paths

def Graph(file):

    a = open(file, 'r')
    b = a.readline()
    b = b.rstrip('\n')
    dict_of_paths = {}
    while b:
        first_parse = b.split()
        dict_of_paths[first_parse[0]] = []
        
        second_parse = first_parse[-1].split(',')
        for item in second_parse:
            dict_of_paths[first_parse[0]].append(item)

        b = a.readline()
        b = b.rstrip('\n')

    a.close()
    return dict_of_paths

def EulerianCycleForInt(file):

    a = open(file, 'r')
    b = a.readline()
    b = b.rstrip('\n')
    dict_of_paths = {}
    while b:
        first_parse = b.split()
        dict_of_paths[int(first_parse[0])] = []
        
        second_parse = first_parse[-1].split(',')
        for item in second_parse:
            dict_of_paths[int(first_parse[0])].append(int(item))

        b = a.readline()
        b = b.rstrip('\n')

    a.close()
    
    from random import choice as choice

    keys = dict_of_paths.keys()
    nodes = [key for key in keys]

    path = []
    
    stack = []
    stack.append(choice(nodes))

    while stack:
        V = len(dict_of_paths[stack[-1]])
        if V == 0:
            path.append(stack[-1])
            stack.pop()
        else:
            node = choice(dict_of_paths[stack[-1]])
            c = dict_of_paths[stack[-1]].pop(dict_of_paths[stack[-1]].index(node))
            stack.append(node)

    path.reverse()

    path_string = ''
    for item in path:
        path_string += str(item) + '->'
        
    return path_string.rstrip('->')

def EulerianCycle(file):

    a = open(file, 'r')
    b = a.readline()
    b = b.rstrip('\n')
    dict_of_paths = {}
    while b:
        first_parse = b.split()
        dict_of_paths[first_parse[0]] = []
        
        second_parse = first_parse[-1].split(',')
        for item in second_parse:
            dict_of_paths[first_parse[0]].append(item)

        b = a.readline()
        b = b.rstrip('\n')

    a.close()
    
    from random import choice as choice

    keys = dict_of_paths.keys()
    nodes = [key for key in keys]

    path = []
    
    stack = []
    stack.append(choice(nodes))

    while stack:
        V = len(dict_of_paths[stack[-1]])
        if V == 0:
            path.append(stack[-1])
            stack.pop()
        else:
            node = choice(dict_of_paths[stack[-1]])
            c = dict_of_paths[stack[-1]].pop(dict_of_paths[stack[-1]].index(node))
            stack.append(node)

    path.reverse()

    path_string = ''
    for item in path:
        path_string += str(item) + '->'
        
    return path_string.rstrip('->')

def EulerianCycleFromDict(dict_of_paths):

    from random import choice as choice

    keys = dict_of_paths.keys()
    nodes = [key for key in keys]

    path = []
    
    stack = []
    stack.append(choice(nodes))

    while stack:
        V = len(dict_of_paths[stack[-1]])
        if V == 0:
            path.append(stack[-1])
            stack.pop()
        else:
            node = choice(dict_of_paths[stack[-1]])
            c = dict_of_paths[stack[-1]].pop(dict_of_paths[stack[-1]].index(node))
            stack.append(node)

    path.reverse()

    path_string = ''
    for item in path:
        path_string += str(item) + '->'
        
    return path_string.rstrip('->')

def EulerianPathFromInt(file):

    a = open(file, 'r')
    b = a.readline()
    b = b.rstrip('\n')
    dict_of_paths = {}
    while b:
        first_parse = b.split()
        dict_of_paths[int(first_parse[0])] = []
        
        second_parse = first_parse[-1].split(',')
        for item in second_parse:
            dict_of_paths[int(first_parse[0])].append(int(item))

        b = a.readline()
        b = b.rstrip('\n')

    a.close()

    c = dict_of_paths.values()
    
    values = []
    for item in c:
        values.extend(item)
    
            
    all_keys = []

    for value in values:
        if not value in all_keys:
            all_keys.append(value)
    for key in dict_of_paths:
        if not key in all_keys:
            all_keys.append(key)
            
    in_count = {}
    out_count = {}

    for key in all_keys:
        if not key in dict_of_paths:
            out_count[key] = 0
        else:
            out_count[key] = len(dict_of_paths[key])

    for key in all_keys:
        if not key in in_count:
            in_count[key] = values.count(key)

    unbalanced_for_out = []
    unbalanced_for_in = []

    for key in in_count:
        if in_count[key] > out_count[key]:
            unbalanced_for_out.append(key)
        elif in_count[key] < out_count[key]:
            unbalanced_for_in.append(key)
    
    if not unbalanced_for_out[0] in dict_of_paths:
        dict_of_paths[unbalanced_for_out[0]] = [unbalanced_for_in[0]]
    elif unbalanced_for_out[0] in dict_of_paths:
        dict_of_paths[unbalanced_for_out[0]].append(unbalanced_for_in[0])
    
    d = EulerianCycleFromDict(dict_of_paths)

    e = d.split(str(unbalanced_for_out[0]) + '->' + str(unbalanced_for_in[0]))

    f = e[1].split('->')
    if len(e[1]) == 0:
        g = ''
    else:
        g = str(unbalanced_for_in[0])
    g += e[1]+ e[0][len(f[-1]):]
    if len(e[0]) != 0:
        g += str(unbalanced_for_out[0])

    return g


def EulerianPath(file):

    a = open(file, 'r')
    b = a.readline()
    b = b.rstrip('\n')
    dict_of_paths = {}
    while b:
        first_parse = b.split()
        dict_of_paths[first_parse[0]] = []
        
        second_parse = first_parse[-1].split(',')
        for item in second_parse:
            dict_of_paths[first_parse[0]].append(item)

        b = a.readline()
        b = b.rstrip('\n')

    a.close()

    c = dict_of_paths.values()
    
    values = []
    for item in c:
        values.extend(item)
    
            
    all_keys = []

    for value in values:
        if not value in all_keys:
            all_keys.append(value)
    for key in dict_of_paths:
        if not key in all_keys:
            all_keys.append(key)
            
    in_count = {}
    out_count = {}

    for key in all_keys:
        if not key in dict_of_paths:
            out_count[key] = 0
        else:
            out_count[key] = len(dict_of_paths[key])

    for key in all_keys:
        if not key in in_count:
            in_count[key] = values.count(key)

    unbalanced_for_out = []
    unbalanced_for_in = []

    for key in in_count:
        if in_count[key] > out_count[key]:
            unbalanced_for_out.append(key)
        elif in_count[key] < out_count[key]:
            unbalanced_for_in.append(key)
    
    if not unbalanced_for_out[0] in dict_of_paths:
        dict_of_paths[unbalanced_for_out[0]] = [unbalanced_for_in[0]]
    elif unbalanced_for_out[0] in dict_of_paths:
        dict_of_paths[unbalanced_for_out[0]].append(unbalanced_for_in[0])
    
    d = EulerianCycleFromDict(dict_of_paths)

    e = d.split(unbalanced_for_out[0] + '->' + unbalanced_for_in[0])

    f = e[1].split('->')
    if len(e[1]) == 0:
        g = ''
    else:
        g = unbalanced_for_in[0]
    g += e[1]+ e[0][len(f[-1]):]
    if len(e[0]) != 0:
        g += unbalanced_for_out[0]

    h = g.split('->')

    j = h[0]

    h.pop(0)

    for item in h:
        j += item[-1]

    return j
    

def HashTable(k, file_output):

    mutations = {'1' : '0', '0' : '1'}
    
    m = 0

    kstring_list = []

    start_kstring = '0' * k

    end_kstring = '1' * k

    kstring_list.append(start_kstring)

    while not kstring_list[-1] == end_kstring:
        kstring_list_length = len(kstring_list)
        for i in range(kstring_list_length):
            new_kstring = kstring_list[i][:m] + mutations[kstring_list[i][m]] + kstring_list[i][m+1:]
            kstring_list.append(new_kstring)
        m += 1

    a = open(file_output, 'w')
    for item in kstring_list:
        a.write(item + '\n')
    a.close()
    return 'Mission accomplished!'

def UniversalString(k, file_path, file_bruijn, file_output):
    
    HashTable(k, file_path)
    DeBruijnGraphFromKmer(file_path, file_bruijn)
    a = EulerianCycle(file_bruijn)

    b = a.split('->')
    value = []
    for item in b:
        if not item in value:
            value.append(item)
            value.append(item)
    c = value[0]

    value.pop(0)
    for item in value:
        c += item[-1]

    d = open(file_output, 'w')
    d.write(c)
    d.close()
    return 'Mission accomplished!'


def de_bruijn(k, n):
    """
    De Bruijn sequence for alphabet size k 
    and subsequences of length n.
    """
    a = [0] * k * n
    sequence = []
    def db(t, p):
        if t > n:
            if n % p == 0:
                for j in range(1, p + 1):
                    sequence.append(a[j])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)
    db(1, 1)

    string = ''
    for item in sequence:
        string += str(item)

    return string
