# list of possible nucleotide changes
mutations = {'A' : ['T', 'G', 'C'], 'T' : ['A', 'G', 'C'], 'G' : ['A', 'T', 'C'], 'C' : ['A', 'T', 'G']}


# Counts number of times pattern appears in the text
def PatternCount(genome, pattern):
    '''
    (str, str) -> int

    >>> PatternCount("GCGCG", "GCG")
    2
    '''

    # Converts genome to upper case, just in case
    if not genome.isupper():
        genome = genome.upper()

    # Same with pattetn
    if not pattern.isupper():
        pattern = pattern.upper()

    # Store length of pattern in variable so that it is not counter every time
    k = len(pattern)

    # Initial count
    count = 0

    # Compares part of the genome to the pattern for strong match, then moves
    # one symbol right and repeat
    for i in range(len(genome) - k + 1):
        if genome[i : k + i] == pattern:
            count += 1
    
    return count

# Same as PatternCount, but reads data from file. In file: two strings,
# first is genome, second is pattern
def PatternCountFile(file):

    # Opens file with attribute read
    string = open(file, 'r')

    # Strip off special character from the end of the genome string
    genome = string.readline().strip("\n")

    # Same with pattern string
    pattern = string.readline().strip("\n")

    # Close connection
    string.close()

    # Converts genome to upper case, just in case
    if not genome.isupper():
        genome = genome.upper()

    # Same with pattetn
    if not pattern.isupper():
        pattern = pattern.upper()

    # Store length of pattern in variable so that it is not counter every time
    k = len(pattern)

    # Initial count
    count = 0

    # Compares part of the genome to the pattern for strong match, then moves
    # one symbol right and repeat
    for i in range(len(genome) - k + 1):
        if genome[i : k + i] == pattern:
            count += 1

    # Asks user if he/she wants the genome and pattern to be printed or just
    # the count
    def UserInput():
        user_input = input("Do you want to print genome and string? Y/N ")
        if user_input.upper() == "Y":
            print("For pattern " + pattern + " " + str(count) +
                  " matches were found in genome " + genome)
        elif user_input.upper() == "N":
            print("Count is " + str(count))
        else:
            print("Sorry, wrong input")
            UserInput()

    UserInput()

# produce a dict, keys == kmer, values == frequency of appearance of this kmer in sequence
def count(text, k):
    '''
    (str, int) -> list of str
    Find n most frequent k-mers in a string.

    >>> count(ATATA, 3)
    ATA
    >>> count(ACGTTGCATGTCGCATGATGCATGAGAGCT, 4)
    CATG GCAT
    '''
    if not text.isupper():
        text = text.upper()
        
    kmers = {}
    for i in range(len(text) - k):
        kmer = ''
        for m in range(i, k+i):
            kmer = kmer + text[m]
        if kmer in kmers:
            kmers[kmer] = kmers[kmer] + 1
        else:
            kmers[kmer] = 1
    return kmers

# inverse dict from previous function, keys == frequency of appearance of kmers in sequence, values == list of kmers of this frequency of appearance
def sort_kmers(dictio):
    '''
    inverse keys and values in dictionary
    '''

    value_to_kmer = {}
    for kmer in dictio:
        ammount = dictio[kmer]

        if not (ammount in value_to_kmer):
            value_to_kmer[ammount] = [kmer]

        else:
            value_to_kmer[ammount].append(kmer)
    return value_to_kmer

#mix of functions above, input - sequence and k (length of kmer)
#output - dict, keys == frequency of appearance of kmers in sequence
#values == list of kmers of this frequency of appearance
def list_of_kmers(text, k):
    '''
    (str, int) -> dict
    Find all k-mers in string and append it to a dictionary. Keys are frequency of appearance of kmer, value == list of kmers
    '''
    if not text.isupper():
        text = text.upper()
        
    kmers = {}
    for i in range(len(text) - k):
        kmer = ''
        for m in range(i, k+i):
            kmer = kmer + text[m]
        if kmer in kmers:
            kmers[kmer] = kmers[kmer] + 1
        else:
            kmers[kmer] = 1

    return sort_kmers(kmers)


# Produces a complementary strand in 3' -> 5' direction, then reverse it to 5' -> 3'
def reverse_complementary_strand(genome):
    '''
    (str) -> str
    Creates a reverse complementary strand in 5' -> 3' direction
    >>> reverse_complementary_strand('ATGC')
    'GCAT'
    '''
    if not genome.isupper():
        genome = genome.upper()
        
    #dict of complementary nucleotides
    comp_nucleotides = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    reverse_comp_str = ''
    #creates reverse complementary strand in 3' -> 5' direction
    for i in range(len(genome)):
        reverse_comp_str = reverse_comp_str + comp_nucleotides[genome[-1 - i]]

    return reverse_comp_str


# Finds position(s) of the pattern given and its complementary pattern in the genome
# Returns number(s) meaning the position in the genome of the first nucleotide in pattern
def position_of_both_patterns(pattern, genome):
    '''
    (str, str) -> str
     Find all occurrences of a pattern in a genome. Pattern may be either as entered
     or complementary, return a string of positions

     >>> position_matches('ATAT', 'GATATATGCATATACTT')
     '1 3 9'
     '''
    if not genome.isupper():
        genome = genome.upper()

    if not pattern.isupper():
        pattern = pattern.upper()
        
    # Creates complementary pattern
    comp_pattern = reverse_complementary_strand(pattern)
    
    # Creates a dict
    pattern_position = {pattern: []}
    
    for i in range(len(genome) - len(pattern)):
        kmer = ''
        
    # Creates a kmer of length (len(pattern))
        for m in range(i, len(pattern) + i):
            kmer = kmer + genome[m]
            
    # Checks whether kmer is pattern or comp_pattern. If True, add
    # position i to the list of positions in dict
        if kmer == pattern or kmer == comp_pattern:
            pattern_position[pattern].append(i)
            
    # Creates a string of positions with space between two position
    #(i.e. '0 3 6 9')
    pattern_positions_string = ''
    for i in range(len(pattern_position[pattern])):
        pattern_positions_string = pattern_positions_string + str(pattern_position[pattern][i]) + ' '

    return pattern_positions_string.strip()


# Finds position(s) of the pattern given
# Returns number(s) meaning the position in the genome of the first nucleotide in pattern
def position_of_pattern(pattern, genome):
    '''
    (str, str) -> str
     Find all occurrences of a pattern in a genome. Return a string of positions

     >>> position_matches('ATAT', 'GATATATGCATATACTT')
     '1 3 9'
     '''
    
    if not genome.isupper():
        genome = genome.upper()

    if not pattern.isupper():
        pattern = pattern.upper()
        
    # Creates a dict
    pattern_position = {pattern: []}
    
    for i in range(len(genome) - len(pattern)):
        kmer = ''
        
    # Creates a kmer of length (len(pattern))
        for m in range(i, len(pattern) + i):
            kmer = kmer + genome[m]
            
    # Checks whether kmer is pattern. If True, add
    # position i to the list of positions in dict
        if kmer == pattern:
            pattern_position[pattern].append(i)
            
    # Creates a string of positions with space between two position
    # (i.e. '0 3 6 9')
    pattern_positions_string = ''
    for i in range(len(pattern_position[pattern])):
        pattern_positions_string = pattern_positions_string + str(pattern_position[pattern][i]) + ' '

    return pattern_positions_string.strip()

import time

# Find repeated k-mers of given frequency in the part of the genome with the length L
# Upd. 17/10/14. Dafuq was I thinking writing this????
def clump_find(genome, k, L, t):
    '''
    (str, int, int, int) -> list
    Find patterns forming clumps in a string. Given integers L and t, a k-mer Pattern forms
    an (L, t)-clump inside a (larger) string Genome if there is an interval of Genome
    of length L in which this k-mer appears at least t times.

    >>> clump_find('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4)
    ['CGACA', 'GAAGA']
    '''
    if not genome.isupper():
        genome = genome.upper()
    
    kmers_list = {t: []}
    
    # Generate a string of len = 50
    for m in range(len(genome) - L):
        genome_part = ''
        for i in range(m, L + m):
            genome_part = genome_part + genome[i]
            
    # Create a dict of k-mers in generated string, keys == frequency, values == k-mers
        kmers_pattern = list_of_kmers(genome_part, k)
        
    # Checks whether k-mers of frequency t is already in kmers_list. If False, append it
    # to kmers_list
        if t in kmers_pattern:
            for x in range(len(kmers_pattern[t])):
                if not kmers_pattern[t][x] in kmers_list[t]:
                    kmers_list[t].append(kmers_pattern[t][x])

    return kmers_list[t]


# Define the difference between the total number of occurences of G and C in Genome
# If C -> skew = skew - 1; if G -> skew = skew + 1
# If total number increases, the cytosine frequency is low meaning we are on the forward
# half-strand. If total number decreases, the guanine frequency is low meaning we are on
# the reverse half_strand
def skew(string):
    '''
    (str) -> str

    Returns a string of skews

    >>> skew('CATGGGCATCGGCCATACGCC')
    '0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2'
    '''

    if not string.isupper():
        string = string.upper()
        
    m = 0
    m_string = '0 '
    for i in range(len(string)):
        if string[i] == 'C':
            m -= 1
        elif string[i] == 'G':
            m += 1
        m_string += str(m) + ' '
    return m_string.strip()


#Counts skew (look above for the description of skew), finds position(-s), where skew is
# minimum, thus locating oriC (the skew achieve a minimum at the position where the
# reverse half-strand ends and the forward half-strand begins, which is exactly the
# location of oriC)
def minimum_skew(string):
    '''
    Looks for oriC based on minimum skew, returns positions of this minimum(s)

    >>> minimum_skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
    '11 24'
    '''

    if not string.isupper():
        string = string.upper()
        
    # sets the initial skew to 0
    m = 0
    
    # sets the least skew to 0
    n = 0
    list_of_positions = []
    
    # Checks each charachter and changes skew accordingly
    
    for i in range(len(string)):
        if string[i] == 'C':
            m -= 1
        elif string[i] == 'G':
            m += 1
            
    # if new skew less than old least skew (m < n), erases a list_of_positions,
    # appends a position of nucleotide, sets value of the least skew n equal to m
        if m < n:
            list_of_positions = []
            list_of_positions.append(i+1)
            n = m
            
    # If new skew is equal to the least skew, appends its position to the list
        elif m == n:
            list_of_positions.append(i+1)
            
    # Creates a string from the list_of_positions
    positions = ''
    for i in range (len(list_of_positions)):
        positions += str(list_of_positions[i]) + ' '
    return positions.strip()


# Finds position(s) in the genome, where pattern appears with at most d mismatches
def approximate_pattern_matching(pattern, genome, d):
    '''
    (str, str, int) -> int
    Find positions in the genome, where pattern appears with at most d mismatches (differences)
    CGGAT and CGGAC have one mismatch (The last nucleotide)

    >>> approximate_pattern_matching('ATTCTGGA',
                                     'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',
                                     3)
    '6 7 26 27'

    '''

    if not genome.isupper():
        genome = genome.upper()

    if not pattern.isupper():
        pattern = pattern.upper()
        
    position = ''
    for i in range(len(genome) - len(pattern)+1):
        mismatch = 0
        
        # creates a pattern from genome with the same length as pattern
        genome_pattern = genome[i : i + len(pattern)]
        
        # compares two patterns, count mismatches
        for m in range(len(pattern)):
            if not pattern[m] == genome_pattern[m]:
                mismatch += 1
                
        # add position i to the resulting string in case ammount of mismatches
        # at this position between two pattern is less than d
        if mismatch <= d:
            position = position + str(i) + ' '
    return position.strip()


def EndKmerGenerator(kmer):
    '''
    Generate final kmer in dict in case all kmers of length k are generated
    
    AAAA -> CCCC
    ATGC -> CCCG
    '''
    if not kmer.isupper():
        kmer = kmer.upper()

    end_kmer = ''
    
    for i in range(len(kmer)):
        end_kmer = end_kmer + mutations[kmer[i]][2]

    return end_kmer

def KmerMutation(kmer, d):
    '''
    Mutate kmer for d generations, so that at the end a dict is obtained,
    consisted of all possible kmers with up to d mismatches from the original kmer

    >>> KmerMutation('AA', 1)
    >>> {'AA': 0, 'GA': 0, 'AT': 0, 'AC': 0, 'CA': 0, 'AG': 0, 'TA': 0}
    '''

    if not kmer.isupper():
        kmer = kmer.upper()

    kmer_list = []
    kmer_dict = {}

    kmer_list.append(kmer)
    kmer_dict[kmer] = 0
    m = 0 # number of generation
    
    k = len(kmer)
    generation = 1 # first generation
    
    while generation <= d: # checks for number of generation
        kmer_list_length = len(kmer_list)
        
        for i in range(kmer_list_length): # i-th element of list
            
            for j in range(m, k): # j-th symbol of i-th element of list
                
                # generates three new kmers with changes at j-th nucleotide
                new_kmer1 = kmer_list[i][:j] + mutations[kmer_list[i][j]][0] + kmer_list[i][j+1:]
                if not new_kmer1 in kmer_dict:
                    kmer_dict[new_kmer1] = 0
                    kmer_list.append(new_kmer1)

                new_kmer2 = kmer_list[i][:j] + mutations[kmer_list[i][j]][1] + kmer_list[i][j+1:]
                if not new_kmer2 in kmer_dict:
                    kmer_dict[new_kmer2] = 0
                    kmer_list.append(new_kmer2)

                new_kmer3 = kmer_list[i][:j] + mutations[kmer_list[i][j]][2] + kmer_list[i][j+1:]
                if not new_kmer3 in kmer_dict:
                    kmer_dict[new_kmer3] = 0
                    kmer_list.append(new_kmer3)
                    
        generation += 1
        m += 1
        
    return kmer_dict
        
def PossibleKmers(genome, k, d):
    '''
    Generates all possible kmers of length k from genome with up to d mismatches
    For that takes each kmer from genome and generate all possible kmers with up to d
    mismatches with the help of KmerMutation function
    '''

    if not genome.isupper():
        genome = genome.upper()
        
    all_possible_kmers = {}

    for i in range(len(genome) - k + 1):
        
        # generates kmer from genome
        a = KmerMutation(genome[i : i + k], d)
        
        # generates all possible kmers with up to d mismatch for kmer from genome
        for item in a:            
            if not item in all_possible_kmers:
                all_possible_kmers[item] = 0
    return all_possible_kmers


def PossibleKmersReverse(genome, k, d):
    '''
    Generates all possible kmers of length k from genome with up to d mismatches
    For that takes each kmer from genome and generate all possible kmers with up to d
    mismatches with the help of KmerMutation function
    '''

    if not genome.isupper():
        genome = genome.upper()
        
    all_possible_kmers = PossibleKmers(genome, k, d)
    all_possible_kmers_reverse = {}

    for item in all_possible_kmers:
        reverse = reverse_complementary_strand(item)
        if not reverse in all_possible_kmers_reverse:
            all_possible_kmers_reverse[reverse] = 0
    return all_possible_kmers_reverse


def KmerHashTable(k):
    
    '''
    Generates the dict of all possible kmers with the length k
    '''
    
    # number of nucleotide exposed to mutation
    m = 0
    
    kmer_list = [] 
    kmer_dict = {}
       
    # kmer to begin with
    start_kmer = 'A' * k
    
    # final kmer, which ends the kmer_list
    end_kmer = 'C' * k
    
    kmer_list.append(start_kmer)
    kmer_dict[start_kmer] = 0

    # while-loop checks for appearance of end_kmer as the last element of kmer_list
    while not kmer_list[-1] == end_kmer:

        # checks for the length of changed kmer_list
        kmer_list_length = len(kmer_list)

        # dissects the kmer changing the nucleotide at position m to three corresponding
        # nucleotides according to mutations dict
        for i in range(kmer_list_length):

            # inserts mutations[item][0] generating new_kmer1
            new_kmer1 = kmer_list[i][:m] + mutations[kmer_list[i][m]][0] + kmer_list[i][m+1:]

            # appending to the list
            kmer_list.append(new_kmer1)

            # adding to kmer_dict if wasn't in already
            if not new_kmer1 in kmer_dict:
                kmer_dict[new_kmer1] = 0

            # inserts mutations[item][1] generating new_kmer2
            new_kmer2 = kmer_list[i][:m] + mutations[kmer_list[i][m]][1] + kmer_list[i][m+1:]

            # appending to the list
            kmer_list.append(new_kmer2)
            
            # adding to kmer_dict if wasn't in already
            if not new_kmer2 in kmer_dict:
                kmer_dict[new_kmer2] = 0
                
            # inserts mutations[item][2] generating new_kmer3
            new_kmer3 = kmer_list[i][:m] + mutations[kmer_list[i][m]][2] + kmer_list[i][m+1:]
            
            # appending to the list
            kmer_list.append(new_kmer3)

            # adding to kmer_dict if wasn't in already
            if not new_kmer3 in kmer_dict:
                kmer_dict[new_kmer3] = 0

        # moving to the next nucleotide position to induce mutation        
        m += 1

    return kmer_dict



def genome_to_dict(genome, k):
    '''
    Converts genome to the dict of all presented kmers of length k
    '''
    if not genome.isupper():
        genome = genome.upper()
        
    genome_dict = {}
    for i in range(len(genome) - k + 1):
        genome_str = genome[i : i + k]
        if not genome_str in genome_dict:
            genome_dict[genome_str] = 0

    return genome_dict


def check_freq_word_with_mismatch(genome, k, d):
    '''
    Generates a dict of all possible kmers with up to d mismatches,
    first checks the first half of pattern of length l for number of mismatches,
    if <= d, then compares the whole kmer from table with pattern from genome
    '''

    if not genome.isupper():
        genome = genome.upper()
        
    table = PossibleKmers(genome, k, d)
    
    for item in table:
        
        # compare first half of the possible kmers and pattern from the genome
        for i in range(len(genome) + k - 1):
            
            # generates first half of pattern
            string = genome[i : (i + k) // 2]
            mismatch = 0

            # compares strings
            for m in range(len(string)):
                if not item[m] == string[m]:
                    mismatch += 1

            # if number of mismatches overlimits, 'for' for this i breaks execution
            if mismatch > d:
                break

            # if there was no overlimit for first halves compares the whole kmer and pattern
            else:
                full_string = genome[i : i + k]
                mismatch_full = 0
                for j in range(len(full_string)):
                    if not item[j] == full_string[j]:
                        mismatch_full += 1
                if mismatch_full <= d:
                    table[item] = table[item] + 1
                    
    # inverse keys and values in table            
    return sort_kmers(table)


def reverse_table(dictionary):
    '''
    creates a dictionary of reverse complement for each kmer in dictionary given
    '''
    
    reverse_table = {}
    for item in dictionary:
        if not reverse_complementary_strand(item) in reverse_table:
            reverse_table[reverse_complementary_strand(item)] = 0
    return reverse_table


def freq_words_misrev(genome, k, d):
    '''
    Generates a dict of all possible kmers with up to d mismatches,
    first checks the first half of pattern of length l for number of mismatches,
    if <= d, then compares the whole kmer from table with pattern from genome
    '''

    if not genome.isupper():
        genome = genome.upper()
        
    kmer_table = PossibleKmers(genome, k, d)
    reverse_kmer_table = reverse_table(kmer_table)
    total_table = {}

    # compares each kmer in kmer_table with kmers from genome
    for item in kmer_table:
        for i in range(len(genome) - k + 1):
            mismatch = 0
            genome_pattern = genome[i : i + k]
            for m in range(k):
                if not item[m] == genome_pattern[m]:
                    mismatch += 1
            if mismatch <= d:
                kmer_table[item] = kmer_table[item] + 1
    
    # compares each kmer in reverse_kmer_table with kmers from genome
    for item in reverse_kmer_table:
        for i in range(len(genome) - k + 1):
            mismatch = 0
            genome_pattern = genome[i : i + k]
            for m in range(k):
                if not item[m] == genome_pattern[m]:
                    mismatch += 1
            if mismatch <= d:
                reverse_kmer_table[item] = reverse_kmer_table[item] + 1
                
    # sum frequencies of kmer and its' reverse complement and add it to total_table
    for item in kmer_table:
        total_table[item + ' and ' + reverse_complementary_strand(item)] = kmer_table[item] + reverse_kmer_table[reverse_complementary_strand(item)]

    # inverse keys and values in table
    return sort_kmers(total_table)


# Converts fasta genomes to the string
def FastaToString(file):
    
    '''
    Converts files in fasta format or txt files in fasta like formatting.
    Output is a dictionary, keys are Taxons (string after '>' in file) and values
    are genomes
    '''
    
    genome_file = open(file, 'r')
    line = genome_file.readline()
    genomes = {}
    raw_genomes = []
    lines = []
    genome = []
    final_lines = []
    raw_genome = ''

    # parse the file, generates list of taxon names in lines
    # and each genome but the last one in raw_genomes
    while line != '':
        if line[0] == '>':
            lines.append(line)
            raw_genomes.append(raw_genome)
            raw_genome = ''
        else:
            raw_genome += line
        line = genome_file.readline()

    # appends last genome to raw_genomes
    raw_genomes.append(raw_genome)

    # deletes first empty entry in raw_genomes. Appears due to first "raw_genome.append"
    # in while-cond
    if raw_genomes[0] == '':
        raw_genomes.pop(0)

    # get rid of '\n' in genomes, generates final list - genome
    for item in raw_genomes:
        final_genome = ''
        for i in range(len(item)):
            if item[i] in 'ATGC':
                final_genome += item[i]
        genome.append(final_genome)

    # fill in a dictionary
    if len(lines) == 0: # if there is only genome sequence in file
        genomes['Unknown Taxon'] = genome[0]
    else:
        for m in range(len(lines)):

            # [1 : -1] for getting rid of '>' and '\n'
            genomes[lines[m][1 : -1]] = genome[m]
    
    return genomes
