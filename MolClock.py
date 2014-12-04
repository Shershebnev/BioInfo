import DNA
from DNA import Nucleotides as Nucleotides
from random import randrange as randrange


# Finds patterns of length k, which appear in each string with up to d mismatches.
# Uses Neighbors and ApproximatePatternCount functions from DNA.py
def MotifEnumeration(DNA_list, k ,d):
    '''
    (list of str, int, int) -> str

    Finds patterns of length k, which appear in each string with up to
    d mismatches
    >>> MotifEnumeration(['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1)
    'ATT TTT ATA GTT'
    '''

    # List for storing kmers
    final_list = []

    # From each string from the list generates kmers, then generates kmer's
    # neighbors, each of them is then checked for occuring with up to
    # d mismatches in each string. If appears in each string -> adds to
    # final_list
    for item in DNA_list:
        for i in range(len(item) - k + 1):
            pattern = item[i : i + k]
            neighbors = DNA.Neighbors(pattern, d)
            for patterns in neighbors:
                add = 1 # this will change to 0 during following for if pattern
                        # didn't appear in the string. This will also break for-loop
                for string in DNA_list:
                    count = DNA.ApproximatePatternCount(string, patterns, d)
                    if count == 0:
                        add = 0
                        break
                if add == 1:
                    final_list.append(patterns)

    # Cut the number of copies of each kmer down to 1        
    final_list = list(set(final_list))

    output = ""
    for item in final_list:
        output += item + ' '

    return output.rstrip(' ')

### Creates matrix in a form of hash (equal to Profile)
##matrix = {'A' : [0]*len(strings[0]), "G" : [0]*len(strings[0]),
##          "T" : [0]*len(strings[0]), "C" : [0]*len(strings[0])}
##
### For each nucleotide in string increase corresponding count in matrix by 0.1
##for item in strings:
##	for i in range(len(item)):
##		matrix[item[i].upper()] += 1
##
### Round because of such - 0.9999999999999999 - values
##for values in matrix.values():
##	for i in range(len(values)):
##		values[i] = round(values[i], 1)
##
##entropy = 0
##
### i - index of element in the list, j - index of list. For each element i of
### list iterates through each list
##for i in range(len(matrix['T'])):
##	for j in range(len(matrix)):            
### Calculates entropy
##for value in matrix.values():
##	for i in range(len(values)):
##		if value[i] != 0:
##			entropy += value[i] * math.log2(value[i])
##
### For counting Score instead of entropy it is necessary to distinguish which of
### them has max value so that it won't added to score. Possible solution:
### in the if-loop append each score to the list, then pop max value and then
### sum up every entry. Though it looks a little bit stupid for me at the moment


# Calculates minimal distance between pattern and each string in the list
def Distance(pattern, dna_list):
    '''
    (str, list) -> int

    Calculates minimal distance between pattern and each string in the list
    '''
    total_distance = 0

    for item in dna_list:
        distance = float('inf')
        for i in range(len(item) - len(pattern) + 1):
            if DNA.HammingDistance(pattern, item[i : i + len(pattern)]) < distance:
                distance = DNA.HammingDistance(pattern, item[i : i + len(pattern)])
        total_distance += distance

    return total_distance


# Generates all possible kmers of length k. NEED MORE EFFICIENT WAY!!!
def Kmers(k):
    '''
    (int) -> list

    Generates all possible kmers of length k   
    '''
    kmer_list = list(Nucleotides)

    # Takes kmer_list and expands it buy adding each nucleotide to each element
    # until length of k is achieved
    while len(kmer_list[0]) < k:
        
        new_list = []

        for item in kmer_list:
            for nucleotide in Nucleotides:
                new_list.append(item + nucleotide)

        kmer_list = list(new_list)
        
    return kmer_list

# Finds median string, i.e. kmer, which has minimal distance between itself
# and collection of strings
def MedianString(dna_list, k):

    minimal_distance = float('inf')

    kmer_list = Kmers(k)
    median = ""
    for kmer in kmer_list:
        distance = Distance(kmer, dna_list)
        if distance <= minimal_distance:
            minimal_distance = distance
            median = kmer

    return median


# Given matrix of probabilities, finds the most probable kmer in string.
# File contains:
# DNA string\n
# k\n
# Probabilities of A as a string of integers, separated by space\n
# Probabilities of C as a string of integers, separated by space\n
# Probabilities of G as a string of integers, separated by space\n
# Probabilities of T as a string of integers, separated by space\n
def MostProbableFile(file):
    '''
    (file) -> str

    Given matrix of probabilities, finds the most probable kmer in string
    
    >>> a = "C:/Users/Александр/Desktop/Downloads/profile_most_1.txt"
    # ^ extra dataset from the course, Correct answer is TGTCGC
    >>> MostProbable(a)
    'TGTCGC'
    '''

    # File parsing
    a = open(file, 'r')
    pattern = a.readline().rstrip('\n')
    k = int(a.readline().rstrip('\n'))
    b = a.readline().rstrip('\n')
    probs = []
    i = 0
    while b:
        probabilities = list(map(float, b.split(' ')))
        probs.append(probabilities)
        b = a.readline().rstrip('\n')
        i += 1
    a.close()
    
    # Variable for storing maximum probability
    max_prob = 0

    # Variable for storing kmer with maximum probability
    kmer = ''

    # Iterates through the string,
    for j in range(len(pattern) - k + 1):
        string = pattern[j : j + k]
        prob = 1

        # calculates probability for each kmer
        for m in range(len(string)):
            prob *= probs[Nucleotides.index(string[m])][m]

        # and compares it to current maximum probability. If new prob is bigger
        # updates max_prob and kmer
        if prob > max_prob:
            max_prob = prob
            kmer = string

    return kmer


# Same as MostProbableFile, but input data is not read from file, but entered
# "by hand"
def MostProbable(probabilities, k, string):
    '''
    Same as MostProbableFile, but input data is not read from file, but entered
    "by hand"
    '''
    
    # Variable for storing maximum probability
    max_prob = 0

    # Variable for storing kmer with maximum probability
    kmer = ''

    # Iterates through the string,
    for j in range(len(string) - k + 1):
        pattern = string[j : j + k]
        prob = 1

        # calculates probability for each kmer
        for m in range(len(pattern)):
            prob *= probabilities[Nucleotides.index(pattern[m])][m]

        # and compares it to current maximum probability. If new prob is bigger
        # updates max_prob and kmer
        if prob > max_prob:
            max_prob = prob
            kmer = pattern

    return kmer


def Probabilities(profile, k, string):

    probs = []

    for j in range(len(string) - k + 1):
        pattern = string[j : j + k]
        prob = 1

        # calculates probability for each kmer
        for m in range(len(pattern)):
            prob *= profile[Nucleotides.index(pattern[m])][m]
        probs.append(prob)
    return probs


# Generates the probability profile based on list of strings given. Output is in
# form of list of length 4 with each element is list of length len(string)
def Profile(string_list):
    '''
    (list) -> list

    Generates the probability profile based on list of strings given

    >>> Profile(['ATGC', "CGGC"])
    [[0.5, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 1.0], [0.0, 0.5, 1.0, 0.0], [0.0, 0.5, 0.0, 0.0]]
    '''

    # List for storing profile
    profile = [[0]*len(string_list[0]),
               [0]*len(string_list[0]),
               [0]*len(string_list[0]),
               [0]*len(string_list[0])]

    # List for storing number of each nucleotide at i-th position in each string
    count = [0] * 4

    for i in range(len(string_list[0])):
        # Counts nucleotides at i-th position of each string
        for string in string_list:
            count[Nucleotides.index(string[i])] += 1

        # Updates profile with probabilities of each nucleotide at i-th position
        for j in range(len(count)):
            profile[j][i] = round(count[j]/sum(count), 3)
            
        count = [0] * 4

    return profile

# Counts the score of given list of DNA strings
def Score(dna_list):
    '''
    (list) -> int

    Counts the score of given list of DNA strings. Score is defined as a sum
    of number of appearances of nucleotides different from the most popular
    nucleotide at each position
    >>> Score(['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC', 'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
    30
    '''
    
    count = {'A' : [0]*len(dna_list[0]), 'C' : [0]*len(dna_list[0]),
              'G' : [0]*len(dna_list[0]), 'T' : [0]*len(dna_list[0])}

    # Counts appearance of each nucleotide at each position
    for item in dna_list:
        for i in range(len(item)):
            count[item[i]][i] += 1

    score = 0

    # Take number for each nucleotide at each position and calculates the score
    # for this position
    for j in range(len(count['A'])):
        a = []
        for key in count.keys():
            a.append(count[key][j])
        score += sum(a) - max(a)
    
    return score


# Looks for a list of motifs best_motifs, that minimizes Score(best_motifs)
# Iterates through kmers in the first string in dna_list, for each builds
# probabilities profile, which is used to select most probable kmer from
# other strings. For each string finds most probable (or, if not found, takes
# the first) kmer, adds it to the list, generates new profile based on this list
def GreedyMotifSearch(dna_list, k, t):
    """
    (list, int, int) -> str

    Looks for a list of motifs best_motifs, that minimizes Score(best_motifs)
    >>> GreedyMotifSearch(['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG'], 3, 5)
    'CAG CAG CAA CAA CAA'
    """

    best_motifs = []
    for item in dna_list:
        best_motifs.append(item[:k])
    best_score = Score(best_motifs)
    
    for i in range(len(dna_list[0]) - k + 1):
        motifs = [dna_list[0][i : i + k]]
        for j in range(1, t):
            profile = Profile(motifs)
            most_probable = MostProbable(profile, k, dna_list[j])
            if not most_probable:
                motifs.append(dna_list[j][:k])
            else:
                motifs.append(most_probable)

        if Score(motifs) < best_score:
            best_score = Score(motifs)
            best_motifs = motifs

    # best_motifs_string = ''
    #for item in best_motifs:
        #best_motifs_string += item + ' '
            
    for item in best_motifs:
        print(item)

    #return best_motifs_string.rstrip(' ')

# Same as Profile, but uses pseudocounts
def ProfilePseudo(string_list):
    '''
    (list) -> list

    Generates the probability profile based on list of strings given
    Uses pseudocounts
    '''

    # List for storing profile
    profile = [[0]*len(string_list[0]),
               [0]*len(string_list[0]),
               [0]*len(string_list[0]),
               [0]*len(string_list[0])]

    # List for storing number of each nucleotide at i-th position in each string
    '''
    The only difference from Profile, default values in the list are now 1
    instead of 0 in Profile
    '''
    count = [1] * 4

    for i in range(len(string_list[0])):
        # Counts nucleotides at i-th position of each string
        for string in string_list:
            count[Nucleotides.index(string[i])] += 1

        # Updates profile with probabilities of each nucleotide at i-th position
        for j in range(len(count)):
            profile[j][i] = round(count[j]/sum(count), 3)
            
        count = [1] * 4

    return profile

# Same as GreedyMotifSearch, but uses ProfilePseudo for building probabilities
# profile
def GreedyMotifSearchPseudo(dna_list, k, t):
    """
    (list, int, int) -> str

    Looks for a list of motifs best_motifs, that minimizes Score(best_motifs)
    Uses pseudocounts
    >>> GreedyMotifSearchPseudo(['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG'], 3, 5)
    'TTC ATC TTC ATC TTC'
    """

    best_motifs = []
    for item in dna_list:
        best_motifs.append(item[:k])
    best_score = Score(best_motifs)
    
    for i in range(len(dna_list[0]) - k + 1):
        motifs = [dna_list[0][i : i + k]]
        for j in range(1, t):
            profile = ProfilePseudo(motifs)
            most_probable = MostProbable(profile, k, dna_list[j])
            if not most_probable:
                motifs.append(dna_list[j][:k])
            else:
                motifs.append(most_probable)

        if Score(motifs) < best_score:
            best_score = Score(motifs)
            best_motifs = motifs

    # best_motifs_string = ''
    #for item in best_motifs:
        #best_motifs_string += item + ' '
            
    for item in best_motifs:
        print(item)

    #return best_motifs_string.rstrip(' ')


def RandomizedMotifSearch(dna_list, k):

    best_motifs = []
    for item in dna_list:
        pos = randrange(len(item) - k + 1)
        best_motifs.append(item[pos : pos + k])

    best_score = Score(best_motifs)
    
    while float('inf'):
        profile = ProfilePseudo(best_motifs)
        motifs = []
        for item in dna_list:
            motifs.append(MostProbable(profile, k, item))
        motifs_score = Score(motifs)
        if motifs_score < best_score:
            best_motifs = motifs
            best_score = motifs_score
        else:
            return best_motifs, best_score
            

def RandomizedRuns(dna_list, k, n):
    best_motifs = []
    best_score = float('inf')
    for i in range(n):
        a = RandomizedMotifSearch(dna_list, k)
        if a[1] < best_score:
            best_motifs = a[0]
            best_score = a[1]
        if i % 50 == 0:
            print(i)
    for item in best_motifs:
        print(item)
    #return best_motifs
