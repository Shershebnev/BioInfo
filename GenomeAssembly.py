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

    # for each prefix itereates through suffix and if prefix == suffix,
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


def DeBrujinGraph(file_input, file_output):

    a = open(file_input, 'r')
    k = int(a.readline().rstrip('\n'))
    string = a.readline().rstrip('\n')

    connections = {}

    for i in range(len(string) - k + 1):
        if not Prefix(string[i : i + k]) in connections:
            connections[Prefix(string[i : i + k])] = [Suffix(string[i : i + k])]
        else:
            connections[Prefix(string[i : i + k])].append(Suffix(string[i : i + k]))

    cons = []
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

def DeBrujinGraphFromKmers(file_input, file_output):

    a = open(file_input, 'r')
    string = a.readline().rstrip('\n')
    kmers = []
    while string:
        kmers.append(string)
        string = a.readline().rstrip('\n')

    connections = {}

    for item in kmers:
        if not Prefix(item) in connections:
            connections[Prefix(item)] = [Suffix(item)]
        else:
            connections[Prefix(item)].append(Suffix(item))

    cons = []
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
