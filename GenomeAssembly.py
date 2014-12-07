from random import choice as choice

def Prefix(string):
    return string[:-1]

def Suffix(string):
    return string[1:]

# Generates list of all kmers present in a string, sorted alphabetically
# Rewrites the output file. If needs otherwise, change mode to "a"
def Composition(input_file, output_file):
    '''
    (int, str) -> list
    
    Generates list of all kmers present in a string, sorted alphabetically
    NB: REWRITES THE OUTPUT FILE. 
    '''
    a = open(input_file, 'r')
    k = int(a.readline().rstrip('\n'))
    string = a.readline().rstrip('\n')
    
    kmers = []
    for i in range(len(string) - k + 1):
        kmers.append(string[i : i + k])

    kmers = sorted(kmers)
    
    b = open(output_file, 'w')
    
    for item in kmers:
        b.write(item + '\n')
        
    a.close()
    b.close()
    
    return "Azaza"


# Generates string based on list of kmers present (kmers are ordered as they
# appear in the string)
def GenomePathString(input_file, output_file):
    '''
    (list) -> str
    
    Generates string based on list of kmers present, kmers are ordered as
    they appear in the string

    ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'] -> ACCGAAGCT
    '''
    a = open(input_file, 'r')
    string = a.readline().rstrip('\n')
    kmers = []
    while string:
        kmers.append(string)
        string = a.readline().rstrip('\n')
        
    genome = ''

    # Takes first kmer entirely, for every other kmer just add one last
    # nucleotide to the string
    for item in kmers:
        if kmers.index(item) == 0:
            genome += item
        else:
            genome += item[-1]
            
    b = open(output_file, 'w')
    
    b.write(genome + '\n')
    
    a.close()
    b.close()
    
    return "Azaza"

# Generates an adjacency list, which shows the connections between nodes
def OverlapGraph(file_input, file_output):
    """
    Generates an adjacency list
    """
    
    a = open(file_input, 'r')
    kmer = a.readline().rstrip('\n')
    kmers = []
    while kmer:
        kmers.append(kmer)
        kmer = a.readline().rstrip('\n')

    # Generates list of prefixes and list of suffixes for each kmer
    prefixes = []
    suffixes = []
    for item in kmers:
        prefixes.append(Prefix(item))
        suffixes.append(Suffix(item))
        
    connections = []

    # For each prefix itereates through suffix and if prefix == suffix,
    # updates the list of connections
    for i in range(len(prefixes)):
        for j in range(len(suffixes)):
            if prefixes[i] == suffixes[j]:
                connections.append(kmers[j] + ' -> ' + kmers[i])

    b = open(file_output, 'w')
    
    for item in connections:
        b.write(item + '\n')
        
    a.close()
    b.close()
    
    return "Azaza"


# From a string generates De Brujin graph in the form of adjacency
# list, where each edge is a kmer, node at the beginning of edge is
# Prefix(kmer) and node at the end is Suffix(kmer)
def DeBrujinGraph(file_input, file_output):
    """
    Generates De Bujin graph in the form of adjacency list
    Edges == kmers, nodes == Prefix(kmer) and Suffix(kmer)
    """
    a = open(file_input, 'r')
    k = int(a.readline().rstrip('\n'))
    string = a.readline().rstrip('\n')

    connections = {}

    # Walk through the string, for each kmer if Prefix(kmer) not in dict
    # add it as a key with value == Suffix(kmer). If in dict: append
    # Suffix(kmer) to corresponding list
    for i in range(len(string) - k + 1):
        if not Prefix(string[i : i + k]) in connections:
            connections[Prefix(string[i : i + k])] = [Suffix(string[i : i + k])]
        else:
            connections[Prefix(string[i : i + k])].append(Suffix(string[i : i + k]))

    cons = []

    # Output format is: node -> node,node,node...
    for key, value in connections.items():
        connected = key + ' -> '
        for item in value:
            connected += item + ','
        cons.append(connected.rstrip(','))
        
    b = open(file_output, 'w')
    
    for item in cons:
        b.write(item + '\n')
        
    a.close()
    b.close()
    
    return "Azaza"           

# Generates De Brujin Graph from list of kmers in the form of adjacency list
# General principle is the same as in DeBrujinGraph
def DeBrujinGraphFromKmers(file_input, file_output):
    '''
    Generates De Brujin Graph from list of kmers in the form of adjacency list
    '''

    a = open(file_input, 'r')
    string = a.readline().rstrip('\n')
    kmers = []
    while string:
        kmers.append(string)
        string = a.readline().rstrip('\n')

    connections = {}

    # For each kmer in the list, if Prefix(kmer) not in dict, add it as key
    # with value == Suffix(kmer). Otherwise, appends Suffix(kmer) to
    # corresponding list
    for item in kmers:
        if not Prefix(item) in connections:
            connections[Prefix(item)] = [Suffix(item)]
        else:
            connections[Prefix(item)].append(Suffix(item))

    cons = []

    # Output format is: node -> node,node,node...
    for key, value in connections.items():
        connected = key + ' -> '
        for item in value:
            connected += item + ','
        cons.append(connected.rstrip(','))
        
    b = open(file_output, 'w')
    
    for item in cons:
        b.write(item + '\n')
        
    a.close()
    b.close()
    
    return "Azaza" 


# Checks if given node in adjacency list has unexplored edges
# (if length of corresponding list is 0 -> no unexplored edges)
def CheckV(adjacency_list, node):
    """
    Checks if given node in adjacency list has unexplored edges
    """
    return len(adjacency_list[node])


# Takes adjacency list, retuns Eulerian Cycle from this list
def EulerianCycle(file_input, file_output):
    """
    Generates Eulerian Cycle from adjacency list
    """
    
    a = open(file_input, 'r')
    adj_list = {}
    string = a.readline().rstrip('\n')

    # File formating is: 1 -> 2,3
    # First parse based on " -> " element, from this uses the first
    # element as key. For values takes second element and splits it again
    # at ','. Each element is converted into int
    while string:
        adj_list[int(string.split(' -> ')[0])] = list(map(int, string.split(' -> ')[1].split(',')))
        string = a.readline().rstrip('\n')
    
    path = []

    # Contains elements obtained from random walk on the graph
    stack = [choice(list(adj_list.keys()))]

    # If during random walk ends up at the node which doesn't have any
    # out edges starts "returning back" to the nearest edge with
    # available out edges. At each step back latest node in stack
    # is popped and inserted at index 0 in the path.
    # At nodes, which has at least one out edge, randomly choose one of
    # available edges, which is then popped from adjacency list
    while stack:
        if CheckV(adj_list, stack[-1]) == 0:
            path.insert(0, stack.pop(-1)) # step back
        else:
            stack.append(choice(adj_list[stack[-1]])) #random choice
            # after "moving" one node forward uses previous node to
            # delete current node from the list of available edges
            adj_list[stack[-2]].pop(adj_list[stack[-2]].index(stack[-1]))

    # Output format is node->node->node...
    path_string = ''

    for item in path:
        path_string += str(item) + '->'
        
    b = open(file_output, 'w')
    
    b.write(path_string.rstrip('->'))

    a.close()
    b.close()
    
    return "Azaza"


##>>> in_edge = 0
##>>> out_edge = 0
##>>> for key in a.keys():
##	for value in a.values():
##		if key in value:
##			in_edge += 1
##	out_edge = len(a[key])
##	if in_edge != out_edge:
##		print(key)
##	in_edge = 0
##	out_edge = 0
##
##	
##6
# ^ first node

# for value in values:
# for item in values:
# if not item in keys:
# item -> end node
