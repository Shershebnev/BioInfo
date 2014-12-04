# dictionary of codones with corresponding aminoacids. Stop codones = ' '
codone_to_aminoacid_table = {'AUC': 'I', 'UAC': 'Y', 'GAG': 'E', 'GAA': 'E', 'CGU': 'R', 'UAA': ' ', 'GGG': 'G', 'AAA': 'K', 'UGG': 'W', 'GGC': 'G', 'GGA': 'G', 'CGC': 'R', 'CGA': 'R', 'UGC': 'C', 'CGG': 'R', 'AAU': 'N', 'UCU': 'S', 'GCU': 'A', 'UAG': ' ', 'GUU': 'V', 'GGU': 'G', 'GAU': 'D', 'UUA': 'L', 'AUA': 'I', 'GUA': 'V', 'CUG': 'L', 'CUA': 'L', 'GUG': 'V', 'CUC': 'L', 'AGU': 'S', 'ACA': 'T', 'AUG': 'M', 'CUU': 'L', 'GCA': 'A', 'GCG': 'A', 'CAU': 'H', 'AUU': 'I', 'UCG': 'S', 'UUC': 'F', 'CCG': 'P', 'UCA': 'S', 'UUU': 'F', 'GUC': 'V', 'GCC': 'A', 'UGA': ' ', 'UCC': 'S', 'UAU': 'Y', 'CAA': 'Q', 'CAC': 'H', 'CCC': 'P', 'CAG': 'Q', 'CCA': 'P', 'ACU': 'T', 'ACG': 'T', 'UGU': 'C', 'CCU': 'P', 'AGA': 'R', 'UUG': 'L', 'AGC': 'S', 'GAC': 'D', 'AAG': 'K', 'ACC': 'T', 'AAC': 'N', 'AGG': 'R'}

aminoacid_to_mass_table = {'H': 137, 'I': 113, 'K': 128, 'L': 113, 'M': 131, 'N': 114, 'A': 71, 'C': 103, 'D': 115, 'E': 129, 'F': 147, 'G': 57, 'Y': 163, 'P': 97, 'Q': 128, 'R': 156, 'S': 87, 'T': 101, 'V': 99, 'W': 186}

mass_to_aminoacid_table = {128: 'Q', 129: 'E', 131: 'M', 97: 'P', 163: 'Y', 71: 'A', 137: 'H', 103: 'C', 115: 'D', 113: 'L', 114: 'N', 147: 'F', 99: 'V', 87: 'S', 57: 'G', 186: 'W', 156: 'R', 101: 'T'}


def RNAtoProtein(rna_string):
    '''
    Converts given RNA_string into protein sequence

    >>> RNAtoProtein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    'MAMAPRTEINSTRING'
    '''
    if not rna_string.isupper():
        rna_string = rna_string.upper()
        
    protein = ''
    for i in range(0, len(rna_string), 3):
        protein += codone_to_aminoacid_table[rna_string[i : i + 3]]
    return protein.strip()

    


def DNAtoRNA(genome):
    '''
    Converts DNA genome to RNA by changing each T to U
    '''
    if not genome.isupper():
        genome = genome.upper()
        
    RNA_string = ''
    for i in range(len(genome)):
        if genome[i] in 'AGC':
            RNA_string += genome[i]
        else:
            RNA_string += 'U'

    return RNA_string.strip()


def PeptideEncoding(DNA_genome, protein):
    '''
    Finds substrings of a DNA_genome encoding a given aminoacid sequence (protein)
    Take into account pattern of len(protein) from DNA_genome and its' reverse complement
    Output contains either a pattern or a reverse complement of pattern if any of them
    codes for protein given

    >>> PeptideEncoding('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')
    ['ATGGCC', 'GGCCAT', 'ATGGCC']
    '''
    if not peptide.isupper():
        peptide = peptide.upper()

    if not DNA_genome.isupper():
        DNA_genome = DNA_genome.upper()
        
    substrings = []
    for i in range(len(DNA_genome) - len(protein)*3 + 1):
        pattern = RNAtoProtein(DNAtoRNA(DNA_genome[i : i + len(protein)*3]))
        reverse_complement = RNAtoProtein(DNAtoRNA(DNAtotal.reverse_complementary_strand(DNA_genome[i : i + len(protein)*3])))
        if pattern == protein:
            substrings.append(DNA_genome[i : i + len(protein)*3])
        if reverse_complement == protein:
            substrings.append(DNA_genome[i : i + len(protein)*3])
    return substrings

def PeptideMass(peptide):
    '''
    Calculates peptides' mass according to aminoacid_to_mass_table dictionary
    '''
    if not peptide.isupper():
        peptide = peptide.upper()
        
    mass = 0
    for i in range(len(peptide)):
        mass += aminoacid_to_mass_table[peptide[i]]
    return mass


# produces theoretical spectrum of peptide
def TheoreticalSpectrum(peptide):
    
    '''
    Generate masses of each subpeptides, in addition to the mass 0 and the mass of
    the entire peptide
    >>> TheoreticalSpectrum('LEQN')
    '0 113 114 128 129 227 242 242 257 355 356 370 371 484'
    >>> TheoreticalSpectrum('PTP')
    '0 97 97 101 194 198 198 295'
    '''
    
    if not peptide.isupper():
        peptide = peptide.upper()

    # create a linear representation of cyclopeptide by adding all but last char from
    # peptide, i.e. PTP -> PTPPT, TLKIL -> TLKILTLKI
    cyclopeptide = peptide + peptide[: - 1]
    
    spectrum = '0 '
    masses = []

    # generate all possible subpeptides except for the whole peptide. For PTP from
    # cyclo PTPPT: P, T, P, PT, TP, PP. i - length of subpeptide in range [1: len(peptide))
    # for subpeptide: j - index of the aminoacid from protein, j + i - index of (last AA + 1)
    # i = 1, j = 0 -> P, i = 1, j = 1 -> T... i = 2, j = 0 -> PT, i = 2, j = 1 -> TP...
    # appends masses of subpeptides into a list
    for i in range(1, len(peptide)):
        for j in range(len(peptide)):
            masses.append(PeptideMass(cyclopeptide[j : j + i]))

    # sort masses of subpeptides
    masses.sort()

    for item in masses:
        spectrum += str(item) + ' '

    # add mass of peptide
    spectrum += str(PeptideMass(peptide))

    return spectrum
        
        
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


def AcidToMass(peptide):

    '''
    Represents peptide given as masses of its' single aminoacids
    >>> AcidToMass('WQL')
    '186-128-113'
    '''
    
    if not peptide.isupper():
        peptide = peptide.upper()
        
    in_mass = ''
    for i in range(len(peptide)):
        in_mass += str(PeptideMass(peptide[i])) + '-'

    return in_mass.rstrip('-')


def CycloPeptideSequencing(spectrum):
    
    '''
    Generates a list of cyclopeptide sequence in form of aminoacid masses
    (e.g. 186-101-113) which have theoretical spectrum equal to given
    experimental spectrum

    >>> CycloPeptideSequencing('0 113 128 186 241 299 314 427')
    ['128-113-186', '113-186-128', '186-113-128', '113-128-186', '128-186-113', '186-128-113']
    '''

    # converts spectrum into a list with items = masses
    b = spectrum.split()

    # dict for iteration to proceed (if use same dict possible_variants in
    # while-loop, returns an error due to dict size changed
    d = {}

    # contains all possible sequences. Restriction - mass of newly generated
    # sequence must be presented in spectrum
    possible_variants = {}

    # list of single aminoacids of which peptide consists
    initial_list = []

    # list of variants, checked for equality to spectrum
    final_variant = []

    # generation = number of aminoacids in generated sequence for corresponding
    # step of while-loop
    generation = 1

    # produces list of single aminoacids
    for item in b:
        if int(item) in mass_to_aminoacid_table:
            initial_list.append(mass_to_aminoacid_table[int(item)])

    # fill in the possible_variants dictionary with initial pairs
    # single aminoacid : array
    for item in initial_list:
        array = list(initial_list)
        c = array.pop(initial_list.index(item))
        possible_variants[item] = array

    # generates a dict of all possible sequences. Restriction - mass of newly
    # generated sequence must be presented in spectrum
    while generation < len(initial_list):
        for item in possible_variants:
            for i in range(len(possible_variants[item])):
                string = item + possible_variants[item][i]

                # checks whether mass of newly generated sequence is in
                # spectrum. If True, adds it to dictionary
                if str(PeptideMass(string)) in b:
                    array = list(possible_variants[item])
                    f = array.pop(i)
                    d[string] = array
        possible_variants = d
        d = {}
        generation += 1

    # compare theoretical spectrum for each item in possible_variants dict
    # with spectrum given. If equal, appends it to the list of final variants.
    # Sequence is converted according to AcidToMass function
    for item in possible_variants:
        if TheoreticalSpectrum(item) == spectrum:
            final_variant.append(AcidToMass(item))

    return final_variant
