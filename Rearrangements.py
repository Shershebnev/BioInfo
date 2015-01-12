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


# Functions from Chrage station on representing 2-breaks
def ChromosomeToCycle(P):

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

    nodes = nodes[1:-1]
    nodes = nodes.split(' ')
    nodes = [int(x) for x in nodes]
    chromosome = []
    
    for i in range(0, len(nodes), 2):
        if nodes[i + 1] > nodes[i]:
            if i == 0:
                chromosome.append(1)
            else:
                chromosome.append(max(nodes[i + 1], nodes[i]) // 2)
        else:
            if i == 0:
                chromosome.append(-1)
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
