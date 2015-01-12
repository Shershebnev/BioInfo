Complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
Nucleotides = ['A', 'C', 'G', 'T']

##    a[2].append(a[1].pop(a[1].index("atc")))

import time
import sys
##
##def Timing2():
##    start_time = time.time()
##    Neighbors("ATGCATGCAT", 2)
##    print("--- %s seconds ---" % (time.time() - start_time))

### Prints some output. For programms, which take long to process data. Shows
### it's working and didn't hang up =)
### Requires start_time = time.time() at the beginning of the programm
##    
### Prints "Processing" for the first round
##if i == 0:
##    print("Processing...\n")
##
### Prints how much time was spent, what i is, what maximum i is,
### % of completeness, estimated time left
##if i%50 == 0:
##    print("--- %s seconds ---" % round((time.time() - start_time), 2)
##          + " i is %i " % i + ' , max i is ' + str(len(genome) - k) +
##          '. ' + str(round(i/(len(genome) - k)*100, 2)) + "% completed."
##          + "\nEstimated time left: %s" % round((timing *
##                                                 (len(genome) - k) - (time.time()
##                                                                      - start_time)), 2)
##          + " seconds.\n")


# Counts number of times pattern appears in the text
def PatternCount(genome, pattern):
    '''
    (str, str) -> int

    Counts number of times pattern appears in the text

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
    '''
    Same as PatternCount, but reads data from file. In file: two strings, first
    is genome, second is pattern
    '''
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
    def UserInputPC():
        user_input = input("Do you want to print genome and string? Y/N ")
        if user_input.upper() == "Y":
            print("For pattern " + pattern + " " + str(count) +
                  " matches were found in genome " + genome)
        elif user_input.upper() == "N":
            print("Count is " + str(count))
        else:
            print("Sorry, wrong input")
            UserInputPC()

    UserInputPC()
    

# Looks for the most frequent k-mer(-s)
def FrequentWords(genome, k):
    '''
    (str, int) -> list of str
    
    Looks for the most frequent k-mer(-s)
    
    >>> FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
    ['GCAT', 'CATG']
    '''

    if not genome.isupper():
        genome = genome.upper()

    # Initializes dict
    kmer_array = {}

    # Counts k-mers. If k-mer already in dict, it's count (i.e. value)
    # is increased by 1. If not in dict, add pair k-mer (key) : 1 (value)
    for i in range(len(genome) - k + 1):
        if genome[i : k + i] in kmer_array:
            kmer_array[genome[i : k + i]] += 1
        else:
            kmer_array[genome[i : k + i]] = 1

    # Finds the maximum ammount of matches. First converts object of class
    # dict_items to class list, then looks for the maximum value, then converts
    # it to the string
    max_match = max(list(kmer_array.values()))

    # Initialize arrow where keys (i.e. patterns) will be stored
    keys = []

    # Searches for keys, which has values containing max_match value

    # Converts all entries into tuples (key : value)
    for key, value in kmer_array.items():

        # Compares value with max_match. If True, appends corresponding
        # key (i.e. pattern) to the array
        if value == max_match:
            keys.append(key)
            
    return keys

# Produces the reverse complement of a string. NB result is printed in reverse
# direction (i.e. for ATGC the result will be not TACG, but GCAT)
def ReverseComplement(pattern):
    '''
    (str) -> str

    Produces the reverse complement of a string. NB result is printed in reverse
    direction (i.e. for ATGC the result will be not TACG, but GCAT)

    >>> ReverseComplement(AAAACCCGGT)
    ACCGGGTTTT
    '''

    # Just in case
    if not pattern.isupper():
        pattern = pattern.upper()

    # Initializes empty string
    reverse = ""

    # For each letter in pattern gets it's complementary nucleotide from the
    # hash Complement (see the beginning of the file). Starts from the end of
    # the pattern in order to obtain reverse complement.
    for i in range(1, len(pattern) + 1):
        reverse += Complement[pattern[-i]]

    return reverse

# Returns positions [starting] of the genome at which the given pattern
# is located. Result is printed as a string.
def PatternMatching(genome, pattern):
    '''
    (str, str) -> str

    Returns positions [starting] of the genome at which the given pattern is
    located. Result is printed as a string.
    
    >>> PatternMatching("GATATATGCATATACTT", "ATAT")
    1 3 9
    '''

    # Just in case
    if not genome.isupper():
        genome = genome.upper()

    if not pattern.isupper():
        pattern = pattern.upper()

    # Initializes empty string, where all the positions will be kept
    positions = ''

    # For less typing =)
    k = len(pattern)

    # Takes portion of the genome string of length, equal to pattern length,
    # and compares it with the pattern. If they match, position is added to
    # the string
    for i in range(len(genome) - k + 1):
        if genome[i : i + k] == pattern:
            positions += str(i) + ' '

    # Removes last space
    return positions.rstrip()


# Searches for the k-mers, which occurs at least t times in each part of
# the genome of length L
def ClumpFinding(genome, k, L, t):
    '''
    (str, int, int, int) -> str

    Searches for the k-mers, which occurs at least t times in each part of
    the genome of length L

    >>> ClumpFinding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4)
    CGACA GAAGA
    '''
    
    start_time = time.time()
 
    # Just in case
    if not genome.isupper():
        genome = genome.upper()

    # Initializes array for patterns found
    keys = []

    # Modification of FrequentWords function
    for i in range(len(genome) - L + 1):

        # Part of the genome of given length L
        pattern = genome[i : i + L]
        
        # Initializes dict
        kmer_array = {}

        # Counts k-mers. If k-mer already in dict, it's count (i.e. value)
        # is increased by 1. If not in dict, add pair k-mer (key) : 1 (value)
        for j in range(L - k + 1):
            if pattern[j : k + j] in kmer_array:
                kmer_array[pattern[j : k + j]] += 1
            else:
                kmer_array[pattern[j : k + j]] = 1

        # Searches for keys, which has values more or equal to t

        # Converts all entries into tuples (key : value)
        for key, value in kmer_array.items():

            # Compares value with t. If True and not in the array keys,
            # appends corresponding key (i.e. pattern) to the array
            if value >= t and not key in keys:
                keys.append(key)

        if (i % 10000) == 0:
            print("--- %s seconds ---" % (time.time() - start_time) + " i is " + str(i))
        
        
    # Initializes empty string, where all found pattern will be added for
    # result returning
    patterns = ""

    # Joining all the patterns found into one string
    for item in keys:
        patterns += item + " "

    # Deletes the last space from the string
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    return patterns.rstrip()

# REWRITE. SOME OVERFLOW I GUESS. TIME FOR EACH 10K WINDOW GROWS EXPONENTIALLY
##def BetterClumpFinding(genome, k, L, t):
##
##    start_time = time.time()
##    
##    if not genome.isupper():
##        genome = genome.upper()
##
##    keys = []
##    
##    kmer_array = {}
##
##    start_pattern = genome[ : L]
##
##    for j in range(L - k + 1):
##        if start_pattern[j : k + j] in kmer_array:
##            kmer_array[start_pattern[j : k + j]] += 1
##        else:
##            kmer_array[start_pattern[j : k + j]] = 1
##
##    for key, value in kmer_array.items():
##        if value >= t and not key in keys:
##            keys.append(key)
##
##    for i in range(len(genome) - L):
##        
##        new_pattern = start_pattern[1 : ] + genome[L + i]
##
##        if not kmer_array:
##            kmer_array[start_pattern[ : k]] -= 1
##
##        if not new_pattern[-k : ] in kmer_array:
##            kmer_array[new_pattern[-k : ]] = 1
##        else:
##            kmer_array[new_pattern[-k : ]] += 1
##
##        for key, value in kmer_array.items():
##            if value >= t and not key in keys:
##                keys.append(key)
##
##        start_pattern = new_pattern
##
##        if (i % 10000) == 0:
##            print("--- %s seconds ---" % (time.time() - start_time) + " i is " + str(i))
##    
##    return keys


# Faster algorithm then ClumpFinding. For the first window similar to
# ClumpFinding, but for following windows it only decreases frequency of
# first kmer of previous window by 1, increasing frequency of last kmer of new
# window by 1 and then check if new frequency of latter is equal to or
# above threshold (t).
def BetterClumpFinding(genome, k, L, t):
    '''
    Faster algorithm then ClumpFinding
    '''
    
    # Look at the end of the script for explanation
    #start_time = time.time()

    # Again, just in case
    if not genome.isupper():
        genome = genome.upper()

    # Initializing array
    kmers = []

    # Initializing hash table
    kmer_hash = {}

    # First window of length L
    pattern = genome[ : L]

    # Same logic as in ClumpFinding function, only for 1 window
    # Checks each kmer
    for i in range(L - k + 1):

        # If exists in hash table, increases its frequency to 1
        if pattern[i : i + k] in kmer_hash:
            kmer_hash[pattern[i : i + k]] += 1

        # Otherwise adds kmer with frequency equal to 1
        else:
            kmer_hash[pattern[i : i + k]] = 1

    # Looks for kmers having values equal to or above the threshold t
    for key, value in kmer_hash.items():

        # If finds and it is not in the array, append it to the array
        if value >= t and not key in kmers:
            kmers.append(key)

    # Iterating through the windows, increasing frequency for new kmer by 1
    # and decreasing for the gone kmer by 1. Also checks the newly counted
    # kmer whether its frequency is above threshold. If so and it is not in
    # the array, adds it there
    for j in range(L, len(genome)):

        # Generats new pattern
        new_pattern = pattern[1 : ] + genome[j]

        # Decreases the frequency of the gone kmer
        kmer_hash[pattern[ : k]] -= 1

        # Checks whether new kmer is in hash table. If so, increases frequency
        # by 1, if not - adds it with frequency equal to 1
        if not new_pattern[-k : ] in kmer_hash:
            kmer_hash[new_pattern[-k : ]] = 1
        else:
            kmer_hash[new_pattern[-k : ]] += 1

        # Checks the threshold and adds kmer if equal or above and not in the array
        if kmer_hash[new_pattern[-k : ]] >= t and not new_pattern[-k : ] in kmers:
            kmers.append(new_pattern[-k : ])

        # Assign pattern used in this round of iteration to a variable
        # for storing previous pattern
        pattern = new_pattern

# Tool for checking time spent on the algorithm. Prints how many time has
# passed since the start of the algorithm is seconds for every 10000 iteration.
# Module "time" should be imported for this to work (see the beginning of the
# file. and the start_time variable should be initialized (see the beginning of
# this script)
##        if (j % 10000) == 0:
##            print("--- %s seconds ---" % (time.time() - start_time) +
##                  " j is " + str(j))
    return kmers


# Creates a string of skews of the genome given. For each C it decreases by one
# and for each G it increases by 1. For A and T no changes for skew
def Skew(string):
    '''
    (str) -> str

    Creates a string of skews of the genome given. For each C it decreases by one
    and for each G it increases by 1. For A and T no changes for skew
    
    >>> Skew("CATGGGCATCGGCCATACGCC")
    0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
    '''
    # Initializes string for storing skew. Starts from 0
    skew_line = "0 "

    # Variable for storing current skew
    current_skew = 0

    # Iterates throw the genome, increases, decreases or doesn't change the skew
    # and add a new skew to the string
    for i in range(len(string)):
        
        if string[i] == "C":
            current_skew -= 1
        elif string[i] == "G":
            current_skew += 1

        skew_line += str(current_skew) + " "

    # Deletes the last space
    return skew_line.rstrip()


# Finds the positions at which skew reaches its minimum
def MinimumSkew(string):
    '''
    (str) -> str

    Finds the positions at which skew reaches its minimum
    
    >>> MinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
    11 24
    '''

    # Initializes variables for keeping current skew, minimal skew, list of
    # positions and string for results print out
    current_skew = 0
    
    minimal_skew = 0

    position = []
    
    position_string = ""

    # For each nucleotide first check current skew vs minimal skew and updates
    # list of positions according to the comparison results, then shifts one
    # nucleotide right and update current skew
    for i in range(len(string)):

        if current_skew < minimal_skew:
            minimal_skew = current_skew
            position.clear()
            position.append(i)
        elif current_skew == minimal_skew:
            position.append(i)
            
        if string[i] == "C":
            current_skew -= 1
        elif string[i] == "G":
            current_skew += 1

    # Sticks together positions from the list
    for item in position:
        position_string += str(item) + ' '

    # Deletes trailing space
    return position_string.rstrip()


# Counts Hamming Distance, i.e. number of positions, at which given strings
# differs from each other. Strings should be of equal size
def HammingDistance(string1, string2):
    '''
    (str, str) -> int

    Counts Hamming Distance, i.e. number of positions, at which given strings
    differs from each other. Strings should be of equal size
    
    >>> HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC")
    3
    '''

    if not len(string1) == len(string2):
        print("Sorry, but strings should be of equal length :(")
        sys.exit()

    # Initializes variable for number of differences
    distance = 0

    # Compares each position between two strings. If differ, increases distance
    # by 1
    for i in range(len(string1)):
        if not string1[i] == string2[i]:
            distance += 1

    return distance


# Looks for positions at which substring from genome of length len(pattern)
# differs from pattern at d or less positions, i.e. number of mismatches <= d
def ApproximatePatternMatching(genome, pattern, d):
    '''
    (str, str, int) -> str

    Looks for positions at which substring from genome of length len(pattern)
    differs from pattern at d or less positions, i.e. number of mismatches <= d
    
    >>> ApproxPatternMatching('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',
                              'ATTCTGGA', 3)
    6 7 26 27
    '''

    # Initializes string, which will print the resulting list of positions
    positions = ""

    # Iterates through genome generating substrings genome_part, then applies
    # HammingDistance function to substring and pattern. If result of this
    # comparison is less than or equal to d, adds position to the positions
    # string.
    for i in range(len(genome) - len(pattern) + 1):

        # Zeroes the distance
        distance = 0

        # Generate new substring
        genome_part = genome[i : i + len(pattern)]

        # Compares counted distance to the threshold, adds to the string, if True
        if HammingDistance(genome_part, pattern) <= d:
            positions += str(i) + ' '

    # Removes space from right side of the string
    return positions.rstrip()


# Counts how many times pattern appears in genome with up to d mismatch
def ApproximatePatternCount(genome, pattern, d):
    '''
    (str, str, int) -> int

    Counts how many times pattern appears in genome with up to d mismatch
    
    >>> ApproximatePatternCount('TTTAGAGCCTTCAGAGG', 'GAGG', 2)
    4
    '''
    
    # Initializes variable to keep the count
    count = 0

    # Iterates through genome, comparing substring to the pattern. If difference
    # is less than or equal to d, increases count by 1
    for i in range(len(genome) - len(pattern) + 1):

        # Zeroes the distance
        distance = 0

        # Generates new substring
        genome_part = genome[i : i + len(pattern)]

        # Calculates the distance and compares it to the threshold d. Increases
        # count if distance is less than or equal to d
        if HammingDistance(genome_part, pattern) <= d:
            count += 1

    return count

# Generates all the neighbors of pattern, i.e. all kmers with Hamming Distance
# less than or equal to d
def Neighbors(pattern, d):
    '''
    (str, int) -> list

    Generates all the neighbors of pattern, i.e. all kmers with Hamming Distance
    less than or equal to d
    
    >>> Neighbors("ACG", 1)
    ['ACA', 'ACT', 'AAG', 'ATG', 'AGG', 'ACG', 'TCG', 'GCG', 'CCG', 'ACC']
    '''

    # Without any mismatches
    if d == 0:
        return pattern

    # Obvious, I guess. 
    if len(pattern) == 1:
        return Nucleotides

    # Initializes array for storing neighbor kmers
    neighborhood = []

    # Recursively generate suffixes of pattern
    suffix_neighbors = Neighbors(pattern[1:], d)

    # Compares Hamming Distance between each item from suffix_neighbors and
    # suffix(pattern)
    for item in suffix_neighbors:

        # If distance is less than d, add each nucleotide to the item and append
        # it to the neighborhood array, if it is not already there
        if HammingDistance(pattern[1:], item) < d:
            for nucleotide in Nucleotides:
                if not (nucleotide + item) in neighborhood:
                    neighborhood.append(nucleotide + item)

        # If distance is equal to d, take first nucleotide from the pattern, add
        # item to it and append to the neighborhood array, if it is not already there
        else:
            if not (pattern[0] + item) in neighborhood:
                neighborhood.append(pattern[0] + item)

    return neighborhood

# Does the same as FWwMSorting (see below), but takes too long (> 5 minutes).
# Based on genome generates dict of all possbile kmers with up to d mismatch,
# then counts ocurence of each one in the genome (with up to d mismatch)
def FrequentWordsWithMismatch(genome, k, d):
    '''
    Looks for the most frequent kmers of length k with up to d mismatch.

    Too long, faster version is FWwMSorting
    '''
    kmer_array = {}
    for i in range(len(genome) - k + 1):
        kmers = Neighbors(genome[i : i + k], d)
        for kmer in kmers:
            if not kmer in kmer_array:
                count = ApproximatePatternCount(genome, kmer, d)
                kmer_array[kmer] = count

    # Finds the maximum ammount of matches. First converts object of class
    # dict_items to class list, then looks for the maximum value, then converts
    # it to the string
    max_match = max(list(kmer_array.values()))

    # Initialize arrow where keys (i.e. patterns) will be stored
    keys = []

    # Searches for keys, which has values containing max_match value

    # Converts all entries into tuples (key : value)
    for key, value in kmer_array.items():

        # Compares value with max_match. If True, appends corresponding
        # key (i.e. pattern) to the array
        if value == max_match:
            keys.append(key)
            
    return keys


# Looks for the most frequent kmers of length k with up to d mismatch. Uses
# sorting
def FWwMSorting(genome, k, d):
    '''
    (str, int, int) -> list or str (depends on user input)

    Looks for the most frequent kmers of length k with up to d mismatch. Uses
    sorting.
    
    >>> FWwMSorting("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
    GATG ATGC ATGT
    >>> FWwMSorting("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
    ['GATG', 'ATGC', 'ATGT']
    '''

    start_time = time.time()
    
    # Initializes list for storing kmers
    kmer_list = []

    # For each part of genome of length k generates all possible kmers with
    # up to d mismatch, then add them to the list of kmers kmer_list
    for i in range(len(genome) - k + 1):
        kmers = Neighbors(genome[i : i + k], d)
        for kmer in kmers:
            kmer_list.append(kmer)
        
    # Sorts all kmers in the list
    kmer_list = sorted(kmer_list)

    # Initializes dict for storing kmers with counts as values
    kmer_array = {}

    # Counts kmers in kmer_list. If item i and item i + 1 are equal increases
    # count for this kmer in kmer_array by 1
    for i in range(len(kmer_list) - 1):
        if kmer_list[i] == kmer_list[i+1]:
            if not kmer_list[i] in kmer_array:
                kmer_array[kmer_list[i]] = 1
            else:
                kmer_array[kmer_list[i]] += 1

    # Finds the maximum ammount of matches. First converts object of class
    # dict_items to class list, then looks for the maximum value, then converts
    # it to the string
    max_match = max(list(kmer_array.values()))

    # Initialize arrow where keys (i.e. patterns) will be stored
    keys = []

    # Searches for keys, which has values containing max_match value

    # Converts all entries into tuples (key : value)
    for key, value in kmer_array.items():

        # Compares value with max_match. If True, appends corresponding
        # key (i.e. pattern) to the array
        if value == max_match:
            keys.append(key)

    def UserInputFWwM():
        user_input = input("In what form do you want the result, list(L) or string(S)? ")
        if user_input.lower() == "l":
            return 1
        elif user_input.lower() == "s":
            return 2
        else:
            b = UserInputFWwM()
            return b

    a = UserInputFWwM()

    kmer_string = ''
    
    if a == 1:
        return keys
    elif a == 2:
        for item in keys:
            kmer_string += item + ' '
        return kmer_string.rstrip()

# Does the same as FWwMSorting, but also accounts for reverse complements.
# Finds the pair kmer - ReverseComplement(kmer), which maximizes the count of
# matches (with up to d mismatch)
def FWwMAndReverse(genome, k, d):
    '''
    (str, int, int) -> list or str (depends on user input)

    Does the same as FWwMSorting, but also accounts for reverse complements
    
    >>> FWwMAndReverse("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
    ['ACAT ATGT']
    
    >>> FWwMAndReverse("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
    'ACAT ATGT'
    '''

# For time stamps
    start_time = time.time()

    # This part is same as in FWwMSorting function
    kmer_list = []

    for i in range(len(genome) - k + 1):
            
        kmers = Neighbors(genome[i : i + k], d)
        for kmer in kmers:
            kmer_list.append(kmer)
            
        """
        Uncomment if time stamps needed
        """
        # Prints "Processing" for the first round
        if i == 0:
            timing = time.time() - start_time
            print("Processing...\n")

        # Prints how much time was spent, what i is, what maximum i is,
        # % of completeness, estimated time left
        if i%50 == 0:
            print("--- %s seconds ---" % round((time.time() - start_time), 2)
                  + " i is %i " % i + ' , max i is ' + str(len(genome) - k) +
                  '. ' + str(round(i/(len(genome) - k)*100, 2)) + "% completed."
                  + "\nEstimated time left: %s" % round(timing *(len(genome)-k)
                                                        - (time.time() -
                                                         start_time), 2)
                  + " seconds.\n")

    kmer_list = sorted(kmer_list)

    kmer_hash = {}
    
    for i in range(len(kmer_list) - 1):
        if kmer_list[i] == kmer_list[i+1]:
            if not kmer_list[i] in kmer_hash:
                kmer_hash[kmer_list[i]] = 1
            else:
                kmer_hash[kmer_list[i]] += 1

    # Converts all the keys in kmer_hash into a list, which consists of all the
    # kmers, which matched at least once
    a = list(kmer_hash.keys())

    # Initializes the hash where keys will be in a form of string, containing
    # kmer and its reverse complement separated by a space, and values will be
    # the sum of scores for corresponding kmer and its reverse
    kmer_reverse_hash = {}

    # For each kmer in a checks, whether the string of format
    # "kmer ReverseComplement(kmer)" (or vice versa) is already in hash-table.
    # If so, updates the count for this entry.
    # Second if and elif checks whether both kmer and its reverse were present
    # in the first hash kmer_hash (which contains all kmers matched). If both
    # present - sum their scores, otherwise add only score for current kmer
    """
    Part with second if and elif seems ебанутой какой-то to me currently,
    perhaps it can be rewritten %) Or maybe I just need more sleep. Or both...
    """
    for kmer in a:
        if not (kmer + ' ' + ReverseComplement(kmer) in kmer_reverse_hash or
                ReverseComplement(kmer) + ' ' + kmer in kmer_reverse_hash):
            if kmer in kmer_hash and ReverseComplement(kmer) in kmer_hash:
                kmer_reverse_hash[kmer + ' ' + ReverseComplement(kmer)] = kmer_hash[kmer] + kmer_hash[ReverseComplement(kmer)]
            elif kmer in kmer_hash and not ReverseComplement(kmer) in kmer_hash:
                kmer_reverse_hash[kmer + ' ' + ReverseComplement(kmer)] = kmer_hash[kmer]
    
    # Finds the maximum ammount of matches. First converts object of class
    # dict_items to class list, then looks for the maximum value, then converts
    # it to the string
    max_match = max(list(kmer_reverse_hash.values()))

    # Initialize arrow where keys (i.e. patterns) will be stored
    keys = []

    # Searches for keys, which has values containing max_match value

    # Converts all entries into tuples (key : value)
    for key, value in kmer_reverse_hash.items():

        # Compares value with max_match. If True, appends corresponding
        # key (i.e. pattern) to the array
        if value == max_match:
            keys.append(key)

    # Asks user about output format, list or string
    def UserInputFWwMR():
        user_input = input("In what form do you want the result, list(L) or string(S)? ")
        if user_input.lower() == "l":
            return 1
        elif user_input.lower() == "s":
            return 2
        else:
            b = UserInputFWwMR()
            return b

    a = UserInputFWwMR()
    
    kmer_string = ''
    
    if a == 1:
        return keys
    elif a == 2:
        for item in keys:
            kmer_string += item + ', '
        return kmer_string.rstrip(', ')
