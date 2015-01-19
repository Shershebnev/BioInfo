import sys
import re
import random
from DNA import ReverseComplement as Reverse
import time

def Permutation(P, k):
    '''
    (string, int) -> list

    Produces rearrangement, integers in rearranged region change sign
    >>> Permutation('(+1 -4 +3 +5 -2)', 2)
    [1, 2, -5, -3, 4]
    '''
    try:
        index = P.index(k)
    except ValueError:
        index = P.index(-k)
        
    rearrange = P[k - 1 : index + 1]
    rearrange = [x*(-1) for x in rearrange]
    rearrange.reverse()
    
    P_new = P[:k - 1] + rearrange + P[index + 1:]
    return P_new

def PermOutput(P):
    string = '('
    for item in P:
        if item > 0:
            string += '+' + str(item) + ' '
        else:
            string += str(item) + ' '
    string = string.rstrip(' ')
    return string + ')' + '\n'
            

def GreedySorting(file_input, file_output):
    
    a = open(file_input, 'r')
    P = a.readline().rstrip('\n')
    P = P[1:-1]
    P = P.split(' ')
    P = [int(x) for x in P]
    a.close()
    
    perm_length = len(P)

    a = open(file_output, 'a')
    
    for k in range(1, perm_length + 1):
        if abs(P[k-1]) != k:
            P = Permutation(P, k)
            a.write(PermOutput(P))
        if P[k-1] == -k:
            P[k-1] = P[k-1] * (-1)
            a.write(PermOutput(P))
    return 'azaza'
        

def BreakpointCounter(file_input):

    a = open(file_input, 'r')
    P = a.readline().rstrip('\n')
    P = P[1:-1]
    P = P.split(' ')
    P = [int(x) for x in P]
    a.close()
    P.insert(0, 0)
    P.append(len(P))
    
    count = 0
    for i in range(len(P) - 1):
        if P[i + 1] - P[i] != 1:
            count += 1
    return count


# Functions from Charge station on representing 2-breaks
def ChromosomeToCycle(P):
    '''
    (str) -> str

    Transforms chromosome of n block into integers
    Charging station
    
    >>> ChromosomeToCycle('(+1 -2 -3 +4)')
    '(1 2 4 3 6 5 7 8)'
    '''
    P = P[1:-1]
    P = P.split(' ')
    P = [int(x) for x in P]
    node = []
    
    for i in range(1, len(P) + 1):
        j = P[i - 1]
        if j > 0:
            node.append(2 * j - 1)
            node.append(2 * j)
        else:
            node.append(-2 * j)
            node.append(-2 * j - 1)
            
    nodes = '('
    for item in node:
        nodes += str(item) + ' '
    nodes = nodes.rstrip(' ')
    
    return nodes + ')'
        

def CycleToChromosome(nodes):
    '''
    (str) -> str

    Reverse ChromosomeToCycle, converts integers into n blocks

    >>> CycleToChromosome('(1 2 4 3 6 5 7 8)')
    '(+1 -2 -3 +4)'
    '''
    
    nodes = nodes[1:-1]
    nodes = nodes.split(' ')
    nodes = [int(x) for x in nodes]
    chromosome = []
    
    for i in range(0, len(nodes), 2):
        if nodes[i + 1] > nodes[i]:
            if i == 0:
                if nodes[i] == 1:
                    chromosome.append(1)
                else:
                    chromosome.append(nodes[i] // 2 + nodes[i] % 2)
            else:
                chromosome.append(max(nodes[i + 1], nodes[i]) // 2)
        else:
            if i == 0:
                if nodes[i] == 2:
                    chromosome.append(-1)
                else:
                    chromosome.append(-(nodes[i] // 2 + nodes[i] % 2))
            else:
                chromosome.append((-(max(nodes[i + 1], nodes[i])) // 2))

    string = '('
    for item in chromosome:
        if item > 0:
            string += '+' + str(item) + ' '
        else:
            string += str(item) + ' '
    return string.rstrip(' ') + ')'


def ColoredEdges(P):
    '''
    (str) -> str
    
    Converts a set of chromosomes from blocks to integers and returns
    colored edges (edges between synteny blocks)

    >>> ColoredEdges('(+1 -2 -3)(+4 +5 -6)')
    '(2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)'
    '''
    
    chromosomes = P.split(')(')
    for item in chromosomes:
        if item[-1] != ')' and item[0] != '(':
            chromosomes[chromosomes.index(item)] = '(' + item + ')'
        elif item[-1] != ')':
            chromosomes[chromosomes.index(item)] = item + ')'
        elif item[0] != '(':
            chromosomes[chromosomes.index(item)] = '(' + item

    colored_edges = []
    
    for item in chromosomes:
        item = ChromosomeToCycle(item)
        item = item[1:-1].split(' ')
        item.append(item[0])
        
        for i in range(1, len(item), 2):
            string = '(' + str(item[i]) + ', ' + str(item[i + 1]) + ')'
            colored_edges.append(string)
        
    string = ''
    for item in colored_edges:
        string += item + ', '
    return string.rstrip(', ')


def GraphToGenome(graph):
    '''
    Reverse ColoredEdges, from list of colored edges returns list of chromosomes
    After splitting and converting into int: logic behind the algoritm is -
    second integer in each final colored edge is smaller than first one, i.e.
    in (2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8) two final edges, which
    are (5, 1) and (11, 8). At this edges cycle population is stopped, last
    integer is moved to the very first position in the cycle list, which transforms
    the cycle into list of black edges (i.e. (2, 4), (3, 6), (5, 1) - colored edges,
    (1, 2), (4, 3), (6, 5) - black edges)

    >>> GraphToGenome('(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)')
    '(+1 -2 -3)(-4 +5 -6)'
    '''
    
    graph = graph[1:-1].split('), (')
    nodes = []
    for item in graph:
        nodes.append(item.split(', '))
    for item in nodes: # converts each str into int [[2, 4], [3, 8], [7, 5], [6, 1]]
        for i in range(len(item)):
            item[i] = int(item[i])
    cycles = []
    cycle = []
    for i in range(len(nodes)):
        if nodes[i][1] > nodes[i][0]:
            cycle.append(nodes[i][0])
            cycle.append(nodes[i][1])
        else:
            cycle.append(nodes[i][0])
            cycle.insert(0, nodes[i][1])
            cycles.append(cycle)
            cycle = []
    chromosomes = []
    for item in cycles:
        item = [str(x) for x in item]
        chromosomes.append(CycleToChromosome('(' + ' '.join(item) + ')'))
    return ''.join(chromosomes)
    

def TwoBreakOnGenomeGraph(graph, a, b, c, d, output = 'string'):
    '''
    >>> TwoBreak('(2, 4), (3, 8), (7, 5), (6, 1)', 1, 6, 3, 8)
    '(2, 4), (3, 1), (7, 5), (6, 8)'
    '''
    graph = graph[1:-1].split('), (')
    nodes = []
    for item in graph:
        nodes.append(item.split(', '))
    for item in nodes:
        for i in range(len(item)):
            item[i] = int(item[i])
    for item in nodes:
        if a in item:
            i = nodes.index(item)
            if item[0] == a:
                nodes.insert(i, [d, b])
            elif item[1] == a:
                nodes.insert(i, [b, d])
            nodes.pop(i + 1)
        elif c in item:
            i = nodes.index(item)
            if item[0] == c:
                nodes.insert(i, [c, a])
            elif item[1] == c:
                nodes.insert(i, [a, c])
            nodes.pop(i + 1)
    if output == 'string':       
        string = ''
        for item in nodes:
            string += '(' + str(item[0]) + ', ' + str(item[1]) + '), '
        return string.rstrip(', ')
    else:
        return nodes


##def TwoBreakOnGenome(P, a, b, c, d):
##
##    # converts genome into integers, introduces 2-breaks, create list
##    # chromosome_parts consisting of two chromosome parts
##    # [[[2, 4], [3, 1]], [[7, 5], [6, 8]]]
##    edges = ColoredEdges(P)
##    print(edges)
##    new_edges = TwoBreakOnGenomeGraph(edges, a, b, c, d, 'list')
##    print(new_edges)
##    for item in new_edges:
##        if a in item and new_edges.index(item) != len(new_edges):
##            index = new_edges.index(item)
##        elif b in item and new_edges.index(item) != len(new_edges):
##            index = new_edges.index(item)
##    chromosome_parts = [new_edges[:index - 1], new_edges[index - 1:]]
##    print(chromosome_parts)
##
##    # converts each chromosome part into string in format, which is used in
##    # CycleToChromosome, i.e. chromosomes = [('2 4 3 1'), ('7 5 6 8')]
##    chromosomes = []
##    for item in chromosome_parts:
##        string = []
##        for list_item in item:
##            list_item = [str(x) for x in list_item]
##            a = ' '.join(list_item)
##            string.append(a)
##        chromosome = ' '.join(string)
##        chromosomes.append('(' + chromosome + ')')
##    print(chromosomes)
##
##    # converts each cycle into chromosome and appends it to the string for output
##    string = ''
##    for part in chromosomes:
##        string += CycleToChromosome(part)
##    return string

# Logic is following: 2-breaks occur in that manner that new chr.parts can be
# obtained simply by subsetting the list of colored edges, i.e. for
# genome = (+1 -2 -4 + 3) list of colored edges = [2, 4, 3, 8, 7, 5, 6, 1] and
# with breaks at edges 3-8 and 6-1 new cycles are [2, 4, 3, 1] and [8, 7, 5, 6]
# which can be obtained by subsetting based on indexes of breaks coordinates.
# From these two elements the latter is already in form of black edges and
# the former requires the last element of the list to be moved to the 0-th
# position of this list, i.e. [2, 4, 3, 1] -> [1, 2, 4, 3]
def TwoBreakOnGenome(P, a, b, c, d):
    '''
    (str, int, int, int, int) -> str

    Introduces 2-breaks and return new chromosome parts, obtained after break
    >>> TwoBreakOnGenome('(+1 -2 -4 +3)', 1, 6, 3, 8)
    '(+1 -2)(-4 +3)'
    '''
    
    edges = ColoredEdges(P)
    edges = edges.split(' ')
    for item in edges:
        try:
            edges[edges.index(item)] = re.search('\((\d+)\,', item).group(1)
        except AttributeError:
            edges[edges.index(item)] = re.search('(\d+)[),]+', item).group(1)

    # didn't made conversion str -> int on previous step because later at
    # .join step it would require backward transformation
    index_first = min(edges.index(str(a)), edges.index(str(b)),
                      edges.index(str(c)), edges.index(str(d)))
    index_last = max(edges.index(str(a)), edges.index(str(b)),
                      edges.index(str(c)), edges.index(str(d)))

    # checks just in case
    if not {index_first, index_last}.difference({a, c}) or {index_first, index_last}.difference({b, d}):
        next
    else:
        sys.exit('sorry')

    chromosomes = [edges[:index_first + 1] + edges[index_last:], edges[index_first + 1:index_last]]
    chromosomes[0].insert(0, chromosomes[0].pop(-1)) # move last element to 0-th pos
    
    cycles = []
    for item in chromosomes:
        cycles.append('(' + ' '.join(item) + ')')
    string = ''
    for item in cycles:
        string += CycleToChromosome(item)
    return string


def TwoBreakDistance(P, Q):

    colored_p = ColoredEdges(P)
    colored_q = ColoredEdges(Q)

    colored_p = colored_p.split(', ')
    colored_q = colored_q.split(', ')

    for item in colored_p:
        try:
            colored_p[colored_p.index(item)] = re.search('\((\d+)', item).group(1)
        except AttributeError:
            colored_p[colored_p.index(item)] = re.search('(\d+)\)', item).group(1)
    for item in colored_q:
        try:
            colored_q[colored_q.index(item)] = re.search('\((\d+)', item).group(1)
        except AttributeError:
            colored_q[colored_q.index(item)] = re.search('(\d+)\)', item).group(1)

    blocks_count = len(colored_p) // 2
    cycles_count = 0
    
    while colored_p:
        start_edge = colored_p[0]
        current_edge = start_edge
        next_edge = colored_p[1]
        colored_p.pop(0)
        colored_p.pop(0)
        cur_genome = 'p'
        while next_edge != start_edge:
            current_edge = next_edge
            if cur_genome == 'p':
                if colored_q.index(current_edge) % 2 == 1:
                    next_edge = colored_q[colored_q.index(current_edge) - 1]
                    cur_genome = 'q'
                    colored_q.pop(colored_q.index(current_edge))
                    colored_q.pop(colored_q.index(next_edge))
                else:
                    next_edge = colored_q[colored_q.index(current_edge) + 1]
                    cur_genome = 'q'
                    colored_q.pop(colored_q.index(current_edge))
                    colored_q.pop(colored_q.index(next_edge))
            elif cur_genome == 'q':
                if colored_p.index(current_edge) % 2 == 1:
                    next_edge = colored_p[colored_p.index(current_edge) - 1]
                    cur_genome = 'p'
                    colored_p.pop(colored_p.index(current_edge))
                    colored_p.pop(colored_p.index(next_edge))
                else:
                    next_edge = colored_p[colored_p.index(current_edge) + 1]
                    cur_genome = 'p'
                    colored_p.pop(colored_p.index(current_edge))
                    colored_p.pop(colored_p.index(next_edge))
        cycles_count += 1
        if cycles_count % 10 == 0:
            print(cycles_count)
    return blocks_count - cycles_count

##def BreakpointGraph(colored_p, colored_q):
##
##    colored_p = colored_p.split(', ')
##    colored_q = colored_q.split(', ')
##
##    for item in colored_p:
##        try:
##            colored_p[colored_p.index(item)] = re.search('\((\d+)', item).group(1)
##        except AttributeError:
##            colored_p[colored_p.index(item)] = re.search('(\d+)\)', item).group(1)
##    for item in colored_q:
##        try:
##            colored_q[colored_q.index(item)] = re.search('\((\d+)', item).group(1)
##        except AttributeError:
##            colored_q[colored_q.index(item)] = re.search('(\d+)\)', item).group(1)
##
##    breakpoint_graph = []
##    
##    while colored_p:
##        start_edge = colored_p[0]
##        current_edge = start_edge
##        next_edge = colored_p[1]
##        colored_p.pop(0)
##        colored_p.pop(0)
##        cur_genome = 'p'
##        breakpoint_graph.extend([current_edge, next_edge])
##        while next_edge != start_edge:
##            current_edge = next_edge
##            if cur_genome == 'p':
##                if colored_q.index(current_edge) % 2 == 1:
##                    next_edge = colored_q[colored_q.index(current_edge) - 1]
##                    breakpoint_graph.append(next_edge)
##                    cur_genome = 'q'
##                    colored_q.pop(colored_q.index(current_edge))
##                    colored_q.pop(colored_q.index(next_edge))
##                else:
##                    next_edge = colored_q[colored_q.index(current_edge) + 1]
##                    cur_genome = 'q'
##                    breakpoint_graph.append(next_edge)
##                    colored_q.pop(colored_q.index(current_edge))
##                    colored_q.pop(colored_q.index(next_edge))
##            elif cur_genome == 'q':
##                if colored_p.index(current_edge) % 2 == 1:
##                    next_edge = colored_p[colored_p.index(current_edge) - 1]
##                    cur_genome = 'p'
##                    breakpoint_graph.append(next_edge)
##                    colored_p.pop(colored_p.index(current_edge))
##                    colored_p.pop(colored_p.index(next_edge))
##                else:
##                    next_edge = colored_p[colored_p.index(current_edge) + 1]
##                    cur_genome = 'p'
##                    breakpoint_graph.append(next_edge)
##                    colored_p.pop(colored_p.index(current_edge))
##                    colored_p.pop(colored_p.index(next_edge))
##                    
##    return breakpoint_graph                    
##
##
##def NonTrivialCycle(breakpoint_graph):
##    for i in range(0, len(breakpoint_graph) - 2, 2):
##        if breakpoint_graph[i] != breakpoint_graph[i + 2]:
##            return True
##
##def TrivialCycle(breakpoint_graph):
##    '''
##    >>> TrivialCycle(BreakpointGraph('(+1 +2 +3 +4 +5 +6)', '(+1 -3 -6 -5)(+2 -4)'))
##    ['(10, 11)', '(11, 10)']
##    '''
##    trivial_cycles = []
##    for i in range(0, len(breakpoint_graph) - 2, 2):
##        if breakpoint_graph[i] == breakpoint_graph[i + 2]:
##            trivial_cycles.extend(['(' + breakpoint_graph[i] +
##                                  ', ' + breakpoint_graph[i + 1] + ')',
##                                   '(' + breakpoint_graph[i + 1] +
##                                  ', ' + breakpoint_graph[i + 2] + ')'])
##    return trivial_cycles
##
##def RandomEdge(blue_edges):
##    '''
##    >>> RandomEdge(ColoredEdges('(+1 +2 -4 -3)'))
##    '2, 3'
##    '''
##    blue_edges = re.split('\), \(', blue_edges)
##    return blue_edges[random.randrange(0, len(blue_edges))].strip(')').strip('(')
##
##def DeleteTrivial(edges, trivial_cycles):
##    '''
##    >>> DeleteTrivial(ColoredEdges('(+1 +2 +3 +4 +5 +6)'), TrivialCycle(BreakpointGraph('(+1 +2 +3 +4 +5 +6)', '(+1 -3 -6 -5)(+2 -4)')))
##    '(2, 3), (4, 5), (6, 7), (8, 9), (12, 1)']
##    '''
##    for item in trivial_cycles:
##        if item in edges:
##            edges = edges[ : edges.find(item)] + edges[edges.find(item) +
##                                                     len(item) + 2 : ]
##    return edges
##def TwoBreakSorting(P, Q):
##    
####    red_edges = ColoredEdges(P)
####    blue_edges = ColoredEdges(Q)
##    
##    breakpoint_graph = BreakpointGraph(ColoredEdges(P), ColoredEdges(Q))
##
##    red_edges = re.split('\), \(', DeleteTrivial(ColoredEdges(P), TrivialCycle(breakpoint_graph)))
##    red_edges[0] = red_edges[0].strip('(')
##    red_edges[-1] = red_edges[-1].strip(')')
##    blue_edges = DeleteTrivial(ColoredEdges(Q), TrivialCycle(breakpoint_graph))
##
##    while NonTrivialCycle(breakpoint_graph):
##        
##        random_blue_edge = RandomEdge(blue_edges).split(', ')
##        c = random_blue_edge[0]
##        b = random_blue_edge[1]
##        
##        for item in red_edges:
##            print(item, c)
##            if c == item[- len(c) : ]:
##                index_c = red_edges.index(item)
##                break
##        print(index_c)
##        red_edge = red_edges[index_c].split(', ')
##        a = red_edge[0]
##        
##        for item in red_edges:
##            if b == item[- len(b) : ]:
##                index_b = red_edges.index(item)
##                break
##        red_edge = red_edges[index_b].split(', ')
##        d = red_edge[0]
##        red_edges.pop(index_b)
##        red_edges.insert(index_b, d + ', ' + a)
##        red_edges.pop(index_c)
##        red_edges.insert(index_c, c + ', ' + b)
##
##        red_edges = ['(' + item + ')' for item in red_edges]
##        print(red_edges)
##        breakpoint_graph = BreakpointGraph(', '.join(red_edges), blue_edges)
##        red_edges = [item.strip('(').strip(')') for item in red_edges]
##        print(red_edges)
##        

for_hash = {'A' : 1, 'C' : 2, 'T' : 3, 'G' : 4}
def Hash(kmer):
    hash_sum = 0
    for i in range(len(kmer)):
        hash_sum += for_hash[kmer[i]]
    return hash_sum

def SharedKmers(file_input, file_output):
    timer = time.time()
    a = open(file_input, 'r')
    k = int(a.readline().rstrip('\n'))
    string1 = a.readline().rstrip('\n')
    string2 = a.readline().rstrip('\n')

    # creates hash table for each kmer and its reverse complement from string2
    string2_hash = {}
    for i in range(len(string2) - k + 1):
        kmer = string2[i : i + k]
        reverse_kmer = Reverse(kmer)
        hash_sum = Hash(kmer)
        hash_sum_reverse = Hash(reverse_kmer)
        try:
            string2_hash[hash_sum].append([kmer, i])
        except KeyError:
            string2_hash[hash_sum] = [[kmer, i]]
        try:
            string2_hash[hash_sum_reverse].append([reverse_kmer, i])
        except KeyError:
            string2_hash[hash_sum_reverse] = [[reverse_kmer, i]]
    #print(string2_hash)

    a.close()
    print(len(string2_hash))
    coords = []
    for i in range(len(string1) - k + 1):
        if i % 10000 == 0:
            print(i)
            print(time.time() - timer)
            timer = time.time()
        kmer = string1[i : i + k]
        hash_sum = Hash(kmer)
        try:
            potential_kmers = string2_hash[hash_sum]
            for item in potential_kmers:
                if kmer == item[0]:
                    coords.append('(' + str(i) + ', ' + str(item[1]) + ')')
        except KeyError:
            next

    a = open(file_output, 'a')
    for item in coords:
        a.write(item + '\n')
    a.close()
    return 'azaza'
