# Hash of aminoacids, keys are codons, values are one-letter abreviations for
# AA or * for stop codons
Aminoacids = {'UAC': 'Y', 'AAC': 'N', 'GCA': 'A', 'UCA': 'S', 'UGG': 'W', 'CAA': 'Q', 'UCG': 'S', 'GUU': 'V', 'CUU': 'L', 'CGA': 'R', 'GUC': 'V', 'GGC': 'G', 'ACA': 'T', 'AUC': 'I', 'UUG': 'L', 'CAG': 'Q', 'AGA': 'R', 'GCC': 'A', 'CUA': 'L', 'UCU': 'S', 'ACC': 'T', 'AUA': 'I', 'CCG': 'P', 'ACU': 'T', 'AAG': 'K', 'AUG': 'M', 'AGU': 'S', 'AGC': 'S', 'GGG': 'G', 'CCC': 'P', 'UAA': '*', 'UUC': 'F', 'GCU': 'A', 'CUC': 'L', 'CCA': 'P', 'GCG': 'A', 'AGG': 'R', 'CGC': 'R', 'UGA': '*', 'UAG': '*', 'AUU': 'I', 'GAU': 'D', 'UCC': 'S', 'CAC': 'H', 'UUU': 'F', 'CUG': 'L', 'UUA': 'L', 'GAC': 'D', 'UAU': 'Y', 'GGA': 'G', 'GAA': 'E', 'AAU': 'N', 'CCU': 'P', 'CAU': 'H', 'AAA': 'K', 'GAG': 'E', 'ACG': 'T', 'GUG': 'V', 'UGC': 'C', 'GUA': 'V', 'UGU': 'C', 'CGU': 'R', 'GGU': 'G', 'CGG': 'R'}

# Hash of codons, keys are one-letter abbreviations for AA or * for stop codons,
# values are lists of corresponding codons
Codons = {'V': ['GUU', 'GUA', 'GUC', 'GUG'], 'K': ['AAA', 'AAG'], 'D': ['GAU', 'GAC'], 'M': ['AUG'], 'F': ['UUC', 'UUU'], 'E': ['GAA', 'GAG'], 'P': ['CCG', 'CCC', 'CCU', 'CCA'], 'T': ['ACC', 'ACU', 'ACA', 'ACG'], 'Y': ['UAC', 'UAU'], 'H': ['CAU', 'CAC'], 'N': ['AAC', 'AAU'], 'W': ['UGG'], 'A': ['GCA', 'GCC', 'GCU', 'GCG'], 'S': ['UCU', 'UCA', 'UCC', 'UCG', 'AGU', 'AGC'], 'C': ['UGU', 'UGC'], 'I': ['AUU', 'AUC', 'AUA'], 'R': ['CGG', 'CGU', 'AGG', 'CGC', 'AGA', 'CGA'], 'L': ['CUA', 'CUU', 'UUG', 'CUG', 'CUC', 'UUA'], 'G': ['GGU', 'GGA', 'GGC', 'GGG'], '*': ['UGA', 'UAG', 'UAA'], 'Q': ['CAG', 'CAA']}

## Cut versions exclude one AA from the pair, which has identical masses
##(Mass(K) == Mass(Q) == 128)
# Hash of masses of corresponding AA. Key - AA, value - mass in Da (rounded)
Masses = {'K': 128, 'V': 99, 'S': 87, 'C': 103, 'L': 113, 'M': 131, 'G': 57, 'I': 113, 'N': 114, 'T': 101, 'D': 115, 'H': 137, 'Y': 163, 'F': 147, 'E': 129, 'P': 97, 'R': 156, 'Q': 128, 'W': 186, 'A': 71}
Masses_cut = {'V': 99, 'S': 87, 'C': 103, 'M': 131, 'G': 57, 'I': 113, 'N': 114, 'T': 101, 'D': 115, 'H': 137, 'Y': 163, 'F': 147, 'E': 129, 'P': 97, 'R': 156, 'Q': 128, 'W': 186, 'A': 71}

Mass_to_AA = {128: 'Q', 129: 'E', 163: 'Y', 131: 'M', 101: 'T', 71: 'A', 137: 'H', 103: 'C', 115: 'D', 113: 'I', 114: 'N', 147: 'F', 97: 'P', 99: 'V', 87: 'S', 57: 'G', 186: 'W', 156: 'R'}

# List of one-letter AA abreviations
AA = ['R', 'L', 'K', 'P', 'N', 'H', 'D', 'A', 'T', 'W', 'V', 'F', 'I', 'M', 'E', 'Q', 'C', 'Y', 'S', 'G']
AA_cut = ['R', 'P', 'N', 'H', 'D', 'A', 'T', 'W', 'V', 'F', 'I', 'M', 'E', 'Q', 'C', 'Y', 'S', 'G']

# List of AA masses
MassesList = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


import DNA
import string

from collections import Counter
from time import time as time
from random import shuffle as shuffle

shuffle(AA_cut)
# Helper functions, converts DNA to RNA and vice versa
def DtoR(dna):
    '''
    Converts DNA to RNA
    '''
    return dna.replace("T", "U")

def RtoD(rna):
    '''
    Converts RNA to DNA
    '''
    return rna.replace("U", "T")

# Converts peptide into a string, in which each AA is replaced with its mass
def AAtoInt(peptide):
    '''
    Converts peptide into a string, in which each AA is replaced with its mass
    NB: USES CUT VERSION OF AA LIST
    
    >>> AAtoInt("PQRS")
    '97-128-156-87'
    '''
    string = ""
    for i in range(len(peptide)):
        string += str(Masses_cut[peptide[i]]) + '-'
    return string.rstrip('-')

def IntToAA(int_peptide):
    peptide = ""
    intpep = list(map(int, list(int_peptide.split('-'))))
    for i in range(len(intpep)):
        peptide += Mass_to_AA[intpep[i]]
    return peptide

# Converts RNA sequence into protein
def ProteinTranslation(rna):
    '''
    (str) -> str

    >>> ProteinTranslation("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")
    'MAMAPRTEINSTRING'
    '''
    protein = ""
    for i in range(0, len(rna), 3):
        protein += Aminoacids[rna[i : i + 3]]

    # Deletes * at the end from stop codone
    return protein.rstrip('*')


# Searches for DNA sequences, which code for the protein sequence given. Take
# reverse comlement into account as well
def PeptideEncoding(dna, pattern):
    '''
    (str, str) -> list

    Searches for DNA sequences, which code for the protein sequence given. Take
    reverse comlement into account as well

    >>> PeptideEncoding("ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA", "MA")
    ['ATGGCC', 'GGCCAT', 'ATGGCC']
    '''

    # Array for storing found sequences
    sequences = []

    # Checks whether DNA regions code for peptide given, append the DNA sequence
    # to the list. Alternatively, check it's reverse complement
    for i in range(len(dna) - len(pattern) * 3 + 1):
        if ProteinTranslation(DtoR(dna[i : i + len(pattern) * 3])) == pattern:
            sequences.append(dna[i : i + len(pattern) * 3])
        elif ProteinTranslation(DtoR(DNA.ReverseComplement(dna[i : i + len(pattern) * 3]))) == pattern:
            sequences.append(dna[i : i + len(pattern) * 3])
            
    return sequences


# Counts peptides mass
def PeptideMass(peptide):
    '''
    Counts peptides mass (obvious, I guess))
    '''
    mass = 0

    for i in range(len(peptide)):
        mass += Masses[peptide[i]]

    return mass


# Counts mass for peptides entered in form of AA masses, i.e. 57-113-186
def PeptideMassInt(peptide):
    '''
    Counts mass for peptides entered in form of masses, i.e. 57-113-186
    '''
    
    return sum(list(map(int, list(peptide.split('-')))))

# Generates theoretical spectrum of linear peptide (masses are rounded)
def LinearSpectrum(peptide):
    '''
    (str) -> str
    
    Generates theoretical spectrum of linear peptide

    >>> LinearSpectrum("NQEL")
    '0 113 114 128 129 242 242 257 370 371 484'
    '''

    # Initializes array for keeping masses of parts of peptide peptide[:i]
    # I.e. in case of NQEL it will store masses for N, NQ, NQE, NQEL
    parts_mass = []

    # For each part of peptide (see previous comment), counts the mass with the
    # help of PeptideMass function and appends it to the array
    for i in range(1, len(peptide) + 1):
        parts_mass.append(PeptideMass(peptide[ : i]))

    # Initializes the array for spectrum output, 0 is the peak
    spectrum = [0]

    # Take each item in parts_mass, and subtract it from other masses
    for j in range(len(peptide)):
        for k in range(j + 1, len(peptide)):
            spectrum.append(parts_mass[k] - parts_mass[j])

    # Append masses from parts_mass to the output spectrum
    a = [spectrum.append(mass) for mass in parts_mass]

    # Sorting
    spectrum = sorted(spectrum)

    # Initializes string for storing spectrum for output
    string = ''

    # Converts each item from spectrum into string and adds it to the output
    # string
    for item in spectrum:
        string += str(item) + ' '

    # Strip off the trailing space
    return string.rstrip(' ')


# Generates theoretical spectrum of circular peptide (masses are rounded).
# Utilizes the same algorithm as in LinearSpectrum, with slight modifications.
# Commented only them.
def CyclicSpectrum(peptide):
    '''
    (str) -> str

    Generates theoretical spectrum of circular peptide

    >>> CyclicSpectrum("LEQN")
    '0 113 114 128 129 227 242 242 257 355 356 370 371 484'
    '''
    
    parts_mass = []
    for i in range(1, len(peptide) + 1):
        parts_mass.append(PeptideMass(peptide[ : i]))

    spectrum = [0]
    
    for j in range(len(peptide)):
        for k in range(j + 1, len(peptide)):
            spectrum.append(parts_mass[k] - parts_mass[j])

            # Modification for adding parts of peptide wrapping around the end
            # of the peptide. If peptide == NQEL, will add mass of ELN,
            # for instance
            if k < len(peptide) - 1:
                spectrum.append(PeptideMass(peptide) - (parts_mass[k]
                                                    - parts_mass[j]))
            
    a = [spectrum.append(mass) for mass in parts_mass]
    
    spectrum = sorted(spectrum)

    string = ''

    for item in spectrum:
        string += str(item) + ' '

    return string.rstrip(' ')

# Generates theoretical spectrum for peptide given in form of AA masses,
# algortihm is the same as in CyclicSpectrum
def CyclicSpectrumInt(peptide):
    '''
    (str) -> str

    Generates theoretical spectrum for peptide given in form of AA masses

    >>> CyclicSpectrumInt("113-114-128-129")
    '0 113 114 128 129 227 242 242 257 355 356 370 371 484'
    '''
    if type(peptide) == int:
        return peptide
    aa_mass = list(map(int, list(peptide.split('-'))))
    
    parts_mass = [aa_mass[0]]

    for i in range(1, len(aa_mass)):
        parts_mass.append(parts_mass[-1] + aa_mass[i])

    peptide_mass = max(parts_mass)
    
    spectrum = [0]
    
    for j in range(len(parts_mass)):
        for k in range(j + 1, len(parts_mass)):
            spectrum.append(parts_mass[k] - parts_mass[j])

            if k < len(parts_mass) - 1:
                spectrum.append(peptide_mass - (parts_mass[k] - parts_mass[j]))

    a = [spectrum.append(mass) for mass in parts_mass]

    spectrum = sorted(spectrum)

    output = ''
    for item in spectrum:
        output += str(item) + ' '

    return output.rstrip(' ')

# For given spectrum generates a list of possible peptides which results in the
# same spectrum
def CyclopeptideSequencing(spectrum):
    '''
    (str) -> str

    For given spectrum generates a list of possible peptides which results in
    the same spectrum. Outputed peptides is in form of integers (masses)
    instead of AA

    >>> CyclopeptideSequencing("0 113 128 186 241 299 314 427")
    '186-113-128 186-128-113 113-186-128 113-128-186 128-186-113 128-113-186'
    '''
    
    ## Part 1. First iteration, results in non-empty peptides list, which is
    ## crucial for while loop

    # Copies the cut list of AA
    peptides = list(AA_cut)

    # Converts string of given spectrum into a list. Map function converts each
    # item in the list from str into int
    spectrum_set = list(map(int, list(spectrum.split(' '))))

    # Mass of total peptide from spectrum
    max_mass = max(spectrum_set)

    # List of peptides consistent with the spectrum
    final_list = []

    # For while loop. At iteration, when generated peptides are 100% consistent
    # with the spectrum, changes its value to 1, which marks the end threshold
    # for while loop (see condition in while)
    switch = 0

    # For each item in peptides checks if item's LINEAR spectrum is consistent
    # with original spectrum. If True, appends it to the final_list.
    # Use import collections and help(collections.Counter) for
    # Counter function help. Breafly, if a = [1, 2, 2, 3, 4, 5, 6, 8],
    # b = [1, 3, 3, 5] and c = [1, 3, 5], not Counter(b) - Counter(a)
    # returns False meaning that not all items in b are in a (takes count
    # of each item into account), not Counter(c) - Counter(a) returns True,
    # meaning that all items in c are in a
    for item in peptides:
        if not Counter(list(map(int, list(LinearSpectrum(item).split(' '))))) - Counter(spectrum_set):
            final_list.append(item)

    peptides = list(final_list)

    final_list.clear()

    ## Part 2. During each while-iteration add each of 18 AA to each of peptides
    ## in peptides list, then check whether this new peptide has mass equal to
    ## max mass from spectrum, if so, changes switch value to 1, meaning that
    ## at this while-iteration desired peptides are generated, and appends peptide
    ## to the final list. If mass is not max mass, checks new_peptide's spectrum
    ## consisctency with the original spectrum. If True, appends it to final_list
    # See above about switch. During iteration
    while switch == 0:

        for item in peptides:
            for acid in AA_cut:
                new_peptide = item + acid                
                if max_mass == PeptideMass(new_peptide):
                    final_list.append(new_peptide)
                    switch = 1
                elif not Counter(list(map(int, list(LinearSpectrum(new_peptide).split(' '))))) - Counter(spectrum_set):
                    final_list.append(new_peptide)
                    
        peptides = list(final_list)
        
        final_list.clear()
        
    # Converts final_list into string
    string = ''
    
    for item in peptides:
        string += AAtoInt(item) + ' '
        
    return string.rstrip(' ')


# Counts number of masses, which are present both in theoretical spectrum,
# generated by CyclicSpectrum function based on peptide given, and
# experimental spectrum spectrum.
def CyclicScore(peptide, spectrum):
    '''
    (str, str) -> int

    Counts number of masses, which are present both in theoretical and
    experimental spectrum

    >>> CyclicScore("NQEL", "0 99 113 114 128 227 257 299 355 356 370 371 484")
    11
    '''
    
    # Generates theoretical spectrum in a form of list of masses
    
    theoretical_spectrum = list(map(int, list(CyclicSpectrum(peptide).split(' '))))
        
    # Experimental spectrum in same form as theoretical
    experimental_spectrum = list(map(int, list(spectrum.split(' '))))

    # Converts to the dicts
    a = dict(Counter(theoretical_spectrum))
    b = dict(Counter(experimental_spectrum))
    
    count = 0

    # For each key from one spectrum checks if there is same key in the second
    # spectrum. If True, checks their values: if equal, increase count by this
    # value, if not, increase by min of both
    for key in a.keys():
        if key in b:
            if a[key] == b[key]:
                count += a[key]
            else:
                count += min(a[key], b[key])
                
    return count


# Same as CyclicScore, but theoretical spectrum is considered to be of linear
# peptide.
def LinearScore(peptide, spectrum):
    '''
    (str, str) -> int

    Counts number of masses, which are present both in theoretical and
    experimental spectrum

    >>> LinearScore("NQEL", "0 99 113 114 128 227 257 299 355 356 370 371 484")
    8
    '''
    
    # Generates theoretical spectrum in a form of list of masses
    theoretical_spectrum = list(map(int, list(LinearSpectrum(peptide).split(' '))))

    # Experimental spectrum in same form as theoretical
    experimental_spectrum = list(map(int, list(spectrum.split(' '))))

    # Converts to the dicts
    a = dict(Counter(theoretical_spectrum))
    b = dict(Counter(experimental_spectrum))
    
    count = 0

    # For each key from one spectrum checks if there is same key in the second
    # spectrum. If True, checks their values: if equal, increase count by this
    # value, if not, increase by min of both
    for key in a.keys():
        if key in b:
            if a[key] == b[key]:
                count += a[key]
            else:
                count += min(a[key], b[key])
                
    return count

# Generates cyclopeptide which has theoretical spectrum most close to the
# given experimental spectrum. Based on leaderboard concept, i.e. after each
# iteration of expanding and trimming only top n candidates are chosen.
'''
In current implementation it is limited to 100 instead of n, cause with n it
exceeds 5-minutes limit
'''
# Output contains one of possible solutions
def LeaderboardCyclopeptideSequencing(spectrum, n):
    '''
    (str, int) -> str

    Generates cyclopeptide which has theoretical spectrum most close to the
    given experimental spectrum

    >>> LeaderboardCyclopeptideSequencing("0 71 113 129 147 200 218 260 313 331 347 389 460", 10)
    '113-147-71-129' #this may vary
    '''
    
    # For time stamps
    total_time = time()
    start_time = time()
    
    ## Part 1. First iteration, results in non-empty peptides list, which is
    ## crucial for while loop
    
    # Converts string of given spectrum into a list. Map function converts each
    # item in the list from str into int
    spectrum_list = list(map(int, list(spectrum.split(' '))))

    # Cuts the list of AA based on the spectrum. Pops out AA, which were not
    # presented in range between the most light and the most heavy AA
    AA_list = list(AA_cut)
    for item in AA_list:
        if not PeptideMass(item) in spectrum_list:
            a = AA_list.pop(AA_list.index(item))
    
    # Copies the cutted list of AA
    peptides = list(AA_list)

    # Mass of total peptide from spectrum
    max_mass = max(spectrum_list)
    
    final_list = []
    
    # Hash of peptides consistent with the spectrum
    final_hash = {}

    final_hash_inv = {}

    leader_peptide = ""

    # For while loop for selecting top n candidates
    added = 0
    
    # For each item in peptides checks if item's LINEAR spectrum is consistent
    # with original spectrum. If True, appends it to the final_hash.
    # Use import collections and help(collections.Counter) for
    # Counter function help. Breafly, if a = [1, 2, 2, 3, 4, 5, 6, 8],
    # b = [1, 3, 3, 5] and c = [1, 3, 5], not Counter(b) - Counter(a)
    # returns False meaning that not all items in b are in a (takes count
    # of each item into account), not Counter(c) - Counter(a) returns True,
    # meaning that all items in c are in a
    for item in peptides:
        final_hash[item] = CyclicScore(item, spectrum)
        if not leader_peptide:
            leader_peptide = item
        else:
            if CyclicScore(item, spectrum) > CyclicScore(leader_peptide, spectrum):
                leader_peptide = item

    # Inverse final_hash, keys are now counts, values are lists of peptides
    for key, value in final_hash.items():
        if not value in final_hash_inv:
            final_hash_inv[value] = [key]
        else:
            final_hash_inv[value].append(key)

    # Sorting, max counts go first        
    values = sorted(list(final_hash_inv.keys()), reverse = True)

    # Selecting top n peptides
    while added <= n and values:
        if values:
            final_list.extend(final_hash_inv[values[0]])
            added += len(final_hash_inv[values[0]])
            xyz = values.pop(0)
            
    peptides = list(final_list)

    final_list = []
    final_hash = {}
    final_hash_inv = {}

    ## Part 2. During each while-iteration add each of 18 AA to each of peptides
    ## in peptides list, then check whether this new peptide has mass equal to
    ## max mass from spectrum, if so, changes switch value to 1, meaning that
    ## at this while-iteration desired peptides are generated, and appends peptide
    ## to the final list. If mass is not max mass, checks new_peptide's spectrum
    ## consisctency with the original spectrum. If True, appends it to final_list
    # See above about switch. During iteration
    while peptides:
        
        new_peptide = []
        added = 0

        # Expanding each candidate per one AA
        for item in peptides:
            for acid in AA_list:
                new_peptide.append(item + acid)

        # Compares mass of each newly generated peptide with mass of whole peptide.
        # If less, adds to final_hash with key == peptide, value == score
        # If equal, calcualtes the score and, if it is better then for current
        # leader peptide, updates it. In other cases deletes this peptide from
        # the list
        for peptide in new_peptide:
            if PeptideMass(peptide) < max_mass:
                final_hash[peptide] = CyclicScore(peptide, spectrum)
            elif PeptideMass(peptide) == max_mass:
                if CyclicScore(peptide, spectrum) > CyclicScore(leader_peptide, spectrum):
                    leader_peptide = peptide
                    final_hash[peptide] = CyclicScore(peptide, spectrum)
            else:
                new_peptide.pop(new_peptide.index(peptide))
        # Inverse hash
        for key, value in final_hash.items():
            if not value in final_hash_inv:
                final_hash_inv[value] = [key]
            else:
                final_hash_inv[value].append(key)

        values = sorted(list(final_hash_inv.keys()), reverse = True)

        # Leaderboard trimming
        while added <= 100 and values:
            if values:
                final_list.extend(final_hash_inv[values[0]])
                added += len(final_hash_inv[values[0]])
                xyz = values.pop(0)
                
        peptides = list(final_list)
        print(str(len(peptides)) + " " + str(round(time() - start_time, 2)) + " " +
              str(round(time() - total_time, 2)))
        final_list = []
        final_hash = {}
        final_hash_inv = {}
        start_time = time()
        
    return leader_peptide


# Computes the convolution of given specrtum, i.e. difference between each pair
# of elements in the spectrum, and if this difference is > 0 adds it to the dict
# with keys = difference, values = how many times this difference appear
def Convolution(spectrum):
    '''
    (str) -> str

    Computes the convolution of given spectrum, i.e. difference between each
    pair of elements in the spectrum (only differences >0 are taken into account)

    >>> Convolution("0 137 186 323")
    '137 137 186 186 323 49'
    '''

    # Splits spectrum string into a list. Also converts each entry into int
    spectrum = list(map(int, list(spectrum.split(' '))))

    # Initializes dict for stroing diff : counts
    counts = {}

    # Subtract each item on the list from each item on the list, then add or
    # updates entry in the dict if difference is > 0
    for i in range(len(spectrum)):
        for j in range(len(spectrum)):
            convol = spectrum[j] - spectrum[i]
            if convol > 0:
                if not convol in counts:
                    counts[convol] = 1
                else:
                    counts[convol] += 1
                    
    # Initializes output string
    output = ""

    # Adds key to the string as many time as value
    for key, value in counts.items():
        for i in range(value):
            output += str(key) + ' '

    # azaza =)
    return output.rstrip(' ')


##def ConvolutionCyclopeptideSequencing(spectrum, m, n):
##
##    a = list(map(int, list(Convolution(spectrum).split(' '))))
##
##    b = Counter(a)
##
##    c = b.most_common(m)
##
##    AA_list = []
##
##    for i in range(len(c)):
##        AA_list.append(c[i][0])



##def LeaderboardCyclopeptideSequencing2(spectrum, n):
##
##    total_time = time()
##    start_time = time()
##    ## Part 1. First iteration, results in non-empty peptides list, which is
##    ## crucial for while loop
##
##    # Copies the cut list of AA
##    peptides = list(AA_cut)
##
##    # Converts string of given spectrum into a list. Map function converts each
##    # item in the list from str into int
##    spectrum_set = list(map(int, list(spectrum.split(' '))))
##
##    # Mass of total peptide from spectrum
##    max_mass = max(spectrum_set)
##    
##    final_list = []
##    
##    # Hash of peptides consistent with the spectrum
##    final_hash = {}
##
##    final_hash_inv = {}
##
##    leader_peptide = ""
##
##    # For while loop. At iteration, when generated peptides are 100% consistent
##    # with the spectrum, changes its value to 1, which marks the end threshold
##    # for while loop (see condition in while)
##    switch = 0
##    
##    added = 0
##    
##    # For each item in peptides checks if item's LINEAR spectrum is consistent
##    # with original spectrum. If True, appends it to the final_hash.
##    # Use import collections and help(collections.Counter) for
##    # Counter function help. Breafly, if a = [1, 2, 2, 3, 4, 5, 6, 8],
##    # b = [1, 3, 3, 5] and c = [1, 3, 5], not Counter(b) - Counter(a)
##    # returns False meaning that not all items in b are in a (takes count
##    # of each item into account), not Counter(c) - Counter(a) returns True,
##    # meaning that all items in c are in a
##    for item in peptides:
##        final_hash[item] = CyclicScore(item, spectrum)
##        if not leader_peptide:
##            leader_peptide = item
##        else:
##            if CyclicScore(item, spectrum) > CyclicScore(leader_peptide, spectrum):
##                leader_peptide = item
##    
##    for key, value in final_hash.items():
##        if not value in final_hash_inv:
##            final_hash_inv[value] = [key]
##        else:
##            final_hash_inv[value].append(key)
##            
##    values = sorted(list(final_hash_inv.keys()), reverse = True)
##    
##    while added <= n and values:
##        if values:
##            final_list.extend(final_hash_inv[values[0]])
##            added += len(final_hash_inv[values[0]])
##            xyz = values.pop(0)
##            print(str(added))
##            
##    peptides = list(final_list)
##
##    final_list.clear()
##    final_hash.clear()
##    final_hash_inv.clear()
##    i = 1
##    
####    print(CyclicScore(leader_peptide, spectrum))
##    ## Part 2. During each while-iteration add each of 18 AA to each of peptides
##    ## in peptides list, then check whether this new peptide has mass equal to
##    ## max mass from spectrum, if so, changes switch value to 1, meaning that
##    ## at this while-iteration desired peptides are generated, and appends peptide
##    ## to the final list. If mass is not max mass, checks new_peptide's spectrum
##    ## consisctency with the original spectrum. If True, appends it to final_list
##    # See above about switch. During iteration
##    while peptides:
##
##        N = int(n/i)
##        new_peptide = []
##        added = 0
##        
##        for item in peptides:
##            for acid in AA_cut:
##                new_peptide.append(item + acid)
##                
##        for peptide in new_peptide:
##            if PeptideMass(peptide) < max_mass:
##                final_hash[peptide] = CyclicScore(peptide, spectrum)
##                
##            elif PeptideMass(peptide) == max_mass:
##                if CyclicScore(peptide, spectrum) > CyclicScore(leader_peptide, spectrum):
##                    leader_peptide = peptide
##                    final_hash[peptide] = CyclicScore(peptide, spectrum)
##            else:
##                new_peptide.pop(new_peptide.index(peptide))
##                
##        values = sorted(list(final_hash.values()), reverse = True)
##        values = list(set(values[:N]))
##
##        
##        for key, value in final_hash.items():
##            if value in values and added < N:
##                final_list.append(key)
##                added += 1
##                
####                for i in range(len(values)):
####                    if value == values[i]:
####                        final_list.append(key)
####                        added += 1
##
####        for key, value in final_hash.items():
####            if not value in final_hash_inv:
####                final_hash_inv[value] = [key]
####            else:
####                final_hash_inv[value].append(key)
####
####        values = sorted(list(final_hash_inv.keys()), reverse = True)
####    
####        while added < n and values:
####            if values:
####                final_list.extend(final_hash_inv[values[0]])
####                added += len(final_hash_inv[values[0]])
####                xyz = values.pop(0)
####
##        #i *= 1.1        
##        peptides = list(final_list)
##        print(str(len(peptides)) + " " + str(round(time() - start_time, 2)) + " " +
##              str(round(time() - total_time, 2)))
##        final_list.clear()
##        final_hash.clear()
##        final_hash_inv.clear()
##        start_time = time()
##        
##    return leader_peptide




##def LeaderboardCyclopeptideSequencing(spectrum, n):
##
##    total_time = time()
##    start_time = time()
##    ## Part 1. First iteration, results in non-empty peptides list, which is
##    ## crucial for while loop
##
##    # Copies the cut list of AA
##    peptides = list(AA_cut)
##
##    # Converts string of given spectrum into a list. Map function converts each
##    # item in the list from str into int
##    spectrum_list = list(map(int, list(spectrum.split(' '))))
##
##    # Mass of total peptide from spectrum
##    max_mass = max(spectrum_list)
##
##    final_list = []
##    
##    # Hash of peptides consistent with the spectrum
##    final_hash = {}
##
##    final_hash_inv = {}
##
##    # For while loop. At iteration, when generated peptides are 100% consistent
##    # with the spectrum, changes its value to 1, which marks the end threshold
##    # for while loop (see condition in while)
##    switch = 0
##    
##    added = 0
##    
##    # For each item in peptides checks if item's LINEAR spectrum is consistent
##    # with original spectrum. If True, appends it to the final_hash.
##    # Use import collections and help(collections.Counter) for
##    # Counter function help. Breafly, if a = [1, 2, 2, 3, 4, 5, 6, 8],
##    # b = [1, 3, 3, 5] and c = [1, 3, 5], not Counter(b) - Counter(a)
##    # returns False meaning that not all items in b are in a (takes count
##    # of each item into account), not Counter(c) - Counter(a) returns True,
##    # meaning that all items in c are in a
##    for item in peptides:
##        final_hash[item] = CyclicScore(item, spectrum)
##    
##    for key, value in final_hash.items():
##        if not value in final_hash_inv:
##            final_hash_inv[value] = [key]
##        else:
##            final_hash_inv[value].append(key)
##            
##    values = sorted(list(final_hash_inv.keys()), reverse = True)
##    
##    while added <= n and values:
##        if values:
##            final_list.extend(final_hash_inv[values[0]])
##            added += len(final_hash_inv[values[0]])
##            xyz = values.pop(0)
##            print(str(added))
##            
##    peptides = list(final_list)
##
##    final_list.clear()
##    final_hash.clear()
##    final_hash_inv.clear()
##
##    ## Part 2. During each while-iteration add each of 18 AA to each of peptides
##    ## in peptides list, then check whether this new peptide has mass equal to
##    ## max mass from spectrum, if so, changes switch value to 1, meaning that
##    ## at this while-iteration desired peptides are generated, and appends peptide
##    ## to the final list. If mass is not max mass, checks new_peptide's spectrum
##    ## consisctency with the original spectrum. If True, appends it to final_list
##    # See above about switch. During iteration
##    """
##    СВИТЧ ХУЙНЯ
##    """
##    while switch == 0:
##        
##        added = 0
##        
##        for item in peptides:
##            for acid in AA_cut:
##                new_peptide = item + acid                
##                if max_mass == PeptideMass(new_peptide):
##                    switch = 1
##                final_hash[new_peptide] = CyclicScore(new_peptide, spectrum)
##
##        for key, value in final_hash.items():
##            if not value in final_hash_inv:
##                final_hash_inv[value] = [key]
##            else:
##                final_hash_inv[value].append(key)
##
##        values = sorted(list(final_hash_inv.keys()), reverse = True)
##    
##        while added <= n and values:
##            if values:
##                final_list.extend(final_hash_inv[values[0]])
##                added += len(final_hash_inv[values[0]])
##                xyz = values.pop(0)
##                
##        peptides = list(final_list)
##        print(str(len(peptides)) + " " + str(round(time() - start_time, 2)) + " " +
##              str(round(time() - total_time, 2)))
##        final_list.clear()
##        final_hash.clear()
##        final_hash_inv.clear()
##        start_time = time()
##        
##    # Converts final_list into string
##    string = ''
##    
##    for item in peptides:
##        string += AAtoInt(item) + ' '
##        
##    return string.rstrip(' ')
