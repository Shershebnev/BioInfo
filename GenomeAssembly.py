from random import choice as choice

def Prefix(string):
    return string[:-1]

def Suffix(string):
    return string[1:]

def PrefixPaired(string):
    k = len(string.split('|')[0])
    return (string[ : k - 1] + '|' + string[k + 1 : 2 * k])

def SuffixPaired(string):
    k = len(string.split('|')[0])
    return (string[1 : k] + '|' + string[k + 2 : 2 * k + 1])

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
            
    a.close()
    
    kmers = []
    for i in range(len(string) - k + 1):
        kmers.append(string[i : i + k])

    kmers = sorted(kmers)
    
    b = open(output_file, 'w')
    
    for item in kmers:
        b.write(item + '\n')

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
        
    a.close()
    
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

    b.close()
    
    return "Azaza"


def GenomePathStringFromList(kmers):
    genome = ''

    # Takes first kmer entirely, for every other kmer just add one last
    # nucleotide to the string
    for item in kmers:
        genome += item[-1]
    return genome


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
            
    a.close()

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
        
    a.close()
    
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

    a.close()
        
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
        
    b.close()
    
    return "Azaza" 


# Checks if given node in adjacency list has unexplored edges
# (if length of corresponding list is 0 -> no unexplored edges)
def CheckV(adjacency_list, node):
    """
    Checks if given node in adjacency list has unexplored edges
    """
    return len(adjacency_list[node])


### Takes adjacency list, retuns Eulerian Cycle from this list
##def EulerianCycle(file_input, file_output):
##    """
##    Generates Eulerian Cycle from adjacency list
##    """
##    
##    a = open(file_input, 'r')
##    adj_list = {}
##    string = a.readline().rstrip('\n')
##
##    # File formating is: 1 -> 2,3
##    # First parse based on " -> " element, from this uses the first
##    # element as key. For values takes second element and splits it again
##    # at ','. Each element is converted into int
##    while string:
##        adj_list[int(string.split(' -> ')[0])] = list(map(int, string.split(' -> ')[1].split(',')))
##        string = a.readline().rstrip('\n')
##    
##    a.close()
##    
##    path = []
##
##    # Contains elements obtained from random walk on the graph
##    stack = [choice(list(adj_list.keys()))]
##
##    # If during random walk ends up at the node which doesn't have any
##    # out edges starts "returning back" to the nearest edge with
##    # available out edges. At each step back latest node in stack
##    # is popped and inserted at index 0 in the path.
##    # At nodes, which has at least one out edge, randomly choose one of
##    # available edges, which is then popped from adjacency list
##    while stack:
##        if CheckV(adj_list, stack[-1]) == 0:
##            path.insert(0, stack.pop(-1)) # step back
##        else:
##            stack.append(choice(adj_list[stack[-1]])) #random choice
##            # after "moving" one node forward uses previous node to
##            # delete current node from the list of available edges
##            adj_list[stack[-2]].pop(adj_list[stack[-2]].index(stack[-1]))
##
##    # Output format is node->node->node...
##    path_string = ''
##
##    for item in path:
##        path_string += str(item) + '->'
##        
##    b = open(file_output, 'w')
##    
##    b.write(path_string.rstrip('->'))
##
##    b.close()
##    
##    return "Azaza"


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
        adj_list[string.split(' -> ')[0]] = string.split(' -> ')[1].split(',')
        string = a.readline().rstrip('\n')
    
    a.close()
    
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
        path_string += item + '->'
        
    b = open(file_output, 'w')
    
    b.write(path_string.rstrip('->'))

    b.close()
    
    return "Azaza"


# From adjacency list generates Eulerian Path. Need optimization %)
def EulerianPath(file_input, file_output):
    """
    From adjacency list generates Eulerian Path
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
        
    a.close()

    # Search for first and last edge. Compares ammount of in and out edges
    # for each node present as a key
    for key in adj_list.keys():
        in_edge = 0
        out_edge = 0
        for value in adj_list.values():
            if key in value:
                in_edge += 1
            out_edge = len(adj_list[key])
        if in_edge > out_edge:
            end_edge = key
        elif in_edge < out_edge:
            start_edge = key
            
    # In case node is presented only in values (i.e. has no in edges)
    for value in adj_list.values():
        for item in value:
            if not item in adj_list.keys():
                end_edge = item
                
    if not end_edge and start_edge:
        return EulerianCycle(file_input, file_output)
    
    # Add edge end_edge -> start_edge
    if not end_edge in adj_list:
        adj_list[end_edge] = [start_edge]
    else:
        adj_list[end_edge].append(start_edge)

    # Search eulerian cycle
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
    
    path_string = path_string.rstrip('->')

    # Parse Eulerian cycle around added edge, then reconnect to form
    # Eulerian path
    c = path_string.split(str(end_edge) + '->' + str(start_edge))

    # In case added edge is the first step, split will generate following
    # ['', '->7->8->9->6->3->0->2->1->3->4'], which only lack of start edge
    if not c[0]:
        eu_path = str(start_edge) + c[1]
    # If added edge is the last step, it will only lack the last edge
    # ['6->7->8->9->6->3->0->2->1->3->', '']
    elif not c[1]:
        eu_path = c[0] + str(end_edge)
    # In all other cases: delete the "wrap-around" node (2 in example case)
    # add first and last nodes and merge.
    # Example case: ['2->1->3->', '->7->8->9->6->3->0->2']
    else:
        eu_path = str(start_edge) + c[1].rsplit('->', 1)[0] + '->' + c[0] + str(end_edge)

    b = open(file_output, 'w')
    b.write(eu_path)
    b.close()
    
    return 'Azaza'       

# Same as EulerianPath, but works with strings (kmers), instead of int
def EulerianPathStr(file_input, file_output):
    """
    From adjacency list generates Eulerian Path
    """
    
    a = open(file_input, 'r')
    adj_list = {}
    string = a.readline().rstrip('\n')

    # File formating is: 1 -> 2,3
    # First parse based on " -> " element, from this uses the first
    # element as key. For values takes second element and splits it again
    # at ','. Each element is converted into int
    while string:
        adj_list[string.split(' -> ')[0]] = string.split(' -> ')[1].split(',')
        string = a.readline().rstrip('\n')
        
    a.close()

    # Search for first and last edge. Compares ammount of in and out edges
    # for each node present as a key
    for key in adj_list.keys():
        in_edge = 0
        out_edge = 0
        for value in adj_list.values():
            if key in value:
                in_edge += 1
            out_edge = len(adj_list[key])
        if in_edge > out_edge:
            end_edge = key
        elif in_edge < out_edge:
            start_edge = key
            
    # In case node is presented only in values (i.e. has no in edges)
    for value in adj_list.values():
        for item in value:
            if not item in adj_list.keys():
                end_edge = item
                
    if not end_edge and start_edge:
        return EulerianCycle(file_input, file_output)
    
    # Add edge end_edge -> start_edge
    if not end_edge in adj_list:
        adj_list[end_edge] = [start_edge]
    else:
        adj_list[end_edge].append(start_edge)

    # Search eulerian cycle
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
        path_string += item + '->'
    
    path_string = path_string.rstrip('->')

    # Parse Eulerian cycle around added edge, then reconnect to form
    # Eulerian path
    c = path_string.split(end_edge + '->' + start_edge)

    # In case added edge is the first step, split will generate following
    # ['', '->7->8->9->6->3->0->2->1->3->4'], which only lack of start edge
    if not c[0]:
        eu_path = start_edge + c[1]
    # If added edge is the last step, it will only lack the last edge
    # ['6->7->8->9->6->3->0->2->1->3->', '']
    elif not c[1]:
        eu_path = c[0] + end_edge
    # In all other cases: delete the "wrap-around" node (2 in example case)
    # add first and last nodes and merge.
    # Example case: ['2->1->3->', '->7->8->9->6->3->0->2']
    else:
        eu_path = start_edge + c[1].rsplit('->', 1)[0] + '->' + c[0] + end_edge

    b = open(file_output, 'w')
    b.write(eu_path)
    b.close()
    
    return 'Azaza'


# Reconstruct genome from the list of kmers. NB: there is number on the first
# line of datasets, which is not used in this implementation, so delete it
# before using
def StringReconstruction(file_input, file_output):
    '''
    Reconstruct string from the list of kmers
    '''
    # Generate adjacency list
    DeBrujinGraphFromKmers(file_input, file_output)
    # Generate Eulerian path
    EulerianPathStr(file_output, file_output)

    # From Eulerian path generate a list of kmers in the same order as
    # they appear in the path
    a = open(file_output, 'r')
    path = a.readline().rstrip('\n')
    a.close()
    kmers = path.split('->')

    # From this list reconstruct the string
    return GenomePathStringFromList(kmers)

# Generates k-universal binary circular string
def UniversalStringReconstruction(file_input, file_output):
    '''
    Generates k-universal binary circular string
    '''
    
    a = open(file_input, 'r')
    k = int(a.readline().strip('\n'))
    a.close()
    
    a = open(file_output, 'w')
    
    # Generates all possible binary kmers buy converting decimal integers
    # to binary and adding zeros at the beginning so that length of each
    # number is k
    for i in range(2**k):
        c = a.write(str(0)*(k - len(bin(i)[2:])) + bin(i)[2:] + '\n')
        
    a.close()
    # Generate adjacency list
    DeBrujinGraphFromKmers(file_output, file_output)

    # Generate Eulerian Cycle
    EulerianCycle(file_output, file_output)

    # From Eulerian cycle generate a list of kmers in the same order as
    # they appear in the path
    a = open(file_output, 'r')
    path = a.readline().rstrip('\n')
    a.close()

    kmers = path.split('->')
    kmers.pop(0) # repeat deletion (we need circular uni.string, not linear)
    string = GenomePathStringFromList(kmers)
    a = open(file_output, 'w')
    c = a.write(string)
    a.close()
    return 'Azaza'

# Generate string from path obtained from Paired De Bruijn graph.
# Format of list is ['kmer|kmer', ...].
# Generates two strings, first - from first nucleotides of first kmer,
# last - from first nucleotide of second kmer + all nucleotides from the
# very last second kmer. As an output uses first k + d + 1 nucleotides from
# first string, then the hole 'last' string
def StringFromPairedDeBrujin(path, k, d):
    """
    Generate string from path obtained from Paired De Bruijn graph.
    """
    
    first = ''
    last = ''
    
    for item in path:
        first += item[0]
        last += item[k + 1]

    # Add nucleotides from the very last second kmer
    last += path[-1][k + 2:]
    
    return first[: k + d + 1] + last
    

# Take paired reads in format 'read|read\n', reconstruct a string from them
# Uses the same algorithm as in EulerianPath: find start and end nodes, add
# edge between them, find Eulerian Cycle, delete last 'step' in cycle
def PairedDeBrujin(file_input, file_output):
    '''
    Reconstruct string from paired reads. Format for reads in
    file: read|read\n
    '''
    a = open(file_input, 'r')
    k_and_d = a.readline().rstrip('\n').split(' ')
    paired = a.readline().rstrip('\n')
    
    kmers = []
    while paired:
        kmers.append(paired)
        paired = a.readline().rstrip('\n')
        
    a.close()
    
    k = int(k_and_d[0])
    d = int(k_and_d[1])

    connections = {}

    # Generates connections between Prefixes as keys and suffixes as values
    for kmer in kmers:
        if not PrefixPaired(kmer) in connections:
            connections[PrefixPaired(kmer)] = [SuffixPaired(kmer)]
        else:
            connections[PrefixPaired(kmer)].append(SuffixPaired(kmer))


    # Search for first and last edge. Compares ammount of in and out edges
    # for each node present as a key
    for key in connections.keys():
        in_edge = 0
        out_edge = 0
        for value in connections.values():
            if key in value:
                in_edge += 1
            out_edge = len(connections[key])
        if in_edge > out_edge:
            end_edge = key
        elif in_edge < out_edge:
            start_edge = key
            
    # In case node is presented only in values (i.e. has no in edges)
    for value in connections.values():
        for item in value:
            if not item in connections.keys():
                end_edge = item
                
    # Add edge end_edge -> start_edge
    if not end_edge in connections:
        connections[end_edge] = [start_edge]
    else:
        connections[end_edge].append(start_edge)
        
    path = []

    stack = [start_edge]

    # Same algorithm as in above functions
    while stack:
        if CheckV(connections, stack[-1]) == 0:
            path.insert(0, stack.pop(-1)) # step back
        else:
            stack.append(choice(connections[stack[-1]])) #random choice
            # after "moving" one node forward uses previous node to
            # delete current node from the list of available edges
            connections[stack[-2]].pop(connections[stack[-2]].index(stack[-1]))
            
    c = path.pop(-1)

    # Function to reconstruct string from the path
    string = StringFromPairedDeBrujin(path, k, d)

    a = open(file_output, 'w')
    g = a.write(string)
    a.close()

    return "Azaza"
