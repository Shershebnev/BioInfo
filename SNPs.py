_end = '_end_'

def make_trie(words):
    root = dict()
    for word in words:
        current_dict = root
        for letter in word:
            current_dict = current_dict.setdefault(letter, {})
        current_dict = current_dict.setdefault(_end, _end)
    return root

class Node:
    def __init__(self, node, start_i, end_i):
        self.node = node
        self.start_i = start_i
        self.end_i = end_i

    def get_node(self):
        return self.node

    def get_start(self):
        return self.start_i

    def get_end(self):
        return self.end_i

def make_trie2(words):
    root = dict()
    i = 1
    start_i = 0
    updated = 0
    for word in words:
        current_dict = root
        for j in range(len(word)):
            if j == 0:
                node = Node(word[j], j, i)
                if current_dict:
                    for key in current_dict.keys():
                        if word[j] == key.get_node():
                            current_dict = current_dict.setdefault(key, {})
                            updated = 1
                            break
                    if updated == 0:
                        i += 1
                        current_dict = current_dict.setdefault(Node(word[j], j, i), {})
                        #i += 1
                    updated = 0
                else:
                    current_dict = current_dict.setdefault(node, {})
            else:
                for key in current_dict.keys():
                    if word[j] == key.get_node():
                        current_dict = current_dict.setdefault(key, {})
                        updated = 1
                        start_i = key.get_end()
                        break
                if updated == 0:
                    if start_i:
                        i += 1
                        current_dict = current_dict.setdefault(Node(word[j], start_i, i), {})
                        start_i = 0
                    else:
                        current_dict = current_dict.setdefault(Node(word[j], i, i + 1), {})
                        i += 1
                updated = 0
        current_dict = current_dict.setdefault(_end, _end)
    return root

def adjacency_trie(trie, file):
    a = open(file, 'a+')
    if list(trie.keys())[0] == _end:
        return
    for key in trie.keys():
        string = str(key.get_start()) + '->' + str(key.get_end()) + ':' + key.get_node()
        a.write(string + '\n')
        adjacency_trie(trie[key], file)
    a.close()

    
def trie_matching(string, words):
    trie = make_trie(words)
    position = []
    start_len = len(string)
    while string:
        current_trie = trie
        string = string[1 : ]
        for j in range(len(string)):
            if string[j] in current_trie.keys():
                current_trie = current_trie[string[j]]
                if '_end_' in current_trie.keys():
                    position.append(start_len - len(string))
                    break
            else:
                break
    position = [str(pos) for pos in position]
    return ' '.join(position)


def with_file_reader(file, function, extra_dataset = True):
    if function == trie_matching:
        a = open(file, 'r')
        if extra_dataset:
            things = a.readlines()
            a.close()
            things = [thing.rstrip('\n') for thing in things]
            string = things.pop(1)
            words = things[1 : -2]
        else:
            things = a.readlines()
            a.close()
            things = [thing.rstrip('\n') for thing in things]
            string = things.pop(0)
            words = things
        return trie_matching(string, words)


# Creates suffix tree; final node is {'_end_' : '_end_'}
def suffix_tree(word):
    root = dict()
    for i in range(len(word)):
        new_word = word[i : ]
        current_dict = root
        for letter in new_word:
            current_dict = current_dict.setdefault(letter, {})
        current_dict = current_dict.setdefault(_end, _end)
    return root


##def concat_nodes(tree):
##    current_tree = tree
##    #for i in range(j):
##    for item in current_tree.keys():
##        if len(current_tree[item].keys()) == 1:
##            if list(current_tree[item].keys())[0] != _end:
##                popped = current_tree.pop(item)
##                current_tree[item + list(popped.keys())[0]] = popped[list(popped.keys())[0]]
##            else:
##                return
##        elif len(current_tree[item].keys()) > 1:
##            concat_nodes(current_tree[item])
##        print(current_tree)
##    return current_tree


##def concat_nodes(tree):
##    current_tree = tree
##    if list(current_tree.keys())[0] == _end:
##        return
##    for key in current_tree.keys():
##        if len(current_tree[key].keys()) == 1:
##            if list(current_tree[key].keys())[0] != _end:
##                popped = current_tree.pop(key)
##                current_tree[key + list(popped.keys())[0]] = popped[list(popped.keys())[0]]
##                concat_nodes(current_tree[key + list(popped.keys())[0]])
##            else:
##                #print('azaza')
##                continue
##        else:
##            concat_nodes(current_tree[key])
##    print(current_tree)
##    #concat_nodes(current_tree)
##    return current_tree


# Concatinates unbranched nodes
# There is a bug, so needs couple runs to produce the final suffix tree with
# truly concatenated branches. Function final_tree is for that.
def concat_nodes(tree):
    if list(tree.keys())[0] == _end:
        return
    for key in tree.keys():
        if len(tree[key].keys()) == 1: # if node is unbranched
            # bug should be somewhere in this if
            if list(tree[key].keys())[0] != _end: # if next node is not '_end_'
                popped = tree.pop(key) # nodes subsequent to the node to be concatenated
                # new concatenated node created
                tree[key + list(popped.keys())[0]] = popped[list(popped.keys())[0]]
                # on the nodes subsequent to the newly created node
                concat_nodes(tree[key + list(popped.keys())[0]])
            else:
                pass
        else:
            concat_nodes(tree[key])
    return tree


def final_tree(word):
    tree = suffix_tree(word)
    for i in range(len(word)):
        tree = concat_nodes(tree)
    return tree


# Writes all nodes into the file
def nodes(tree, file):
    a = open(file, 'a+')
    if list(tree.keys())[0] == _end:
        return
    for key in tree.keys():
        a.write(key + '\n')
        nodes(tree[key], file)
    a.close()



def longest_repeat(tree, repeat = '', predecessor = '', key = ''):
    for key in tree.keys():
        if not '$' in key:
            if len(predecessor + key) > len(repeat):
                repeat = predecessor + key
            predecessor += key
            repeat = longest_repeat(tree[key], repeat = repeat, predecessor = predecessor, key = key)
            predecessor = predecessor[ : -len(key)]
    return repeat

        
def generate_kmers(string, k):
    kmers = []
    for i in range(len(string) - k + 1):
        kmers.append(string[i : i + k])
    return kmers


# The following two can be done with the suffix trees, needs rewriting
def longest_shared_substring(string1, string2):
    substrings = []
    for k in range(2, int(len(string1) / 2)):
        kmers1 = set(generate_kmers(string1, k))
        kmers2 = set(generate_kmers(string2, k))
        intersect = kmers1.intersection(kmers2)
        if intersect:
            substrings = list(intersect)
    return substrings[0]

def shortest_non_shared_substring(string1, string2):
    substrings = []
    for k in range(2, int(len(string1) / 2)):
        kmers1 = set(generate_kmers(string1, k))
        kmers2 = set(generate_kmers(string2, k))
        difference = kmers1.difference(kmers2)
        if difference:
            return list(difference)[0]



def suffix_array(string):
    suffixes = {}
    for i in range(len(string)):
        suffixes[string[i : ]] = i
    sorted_array = sorted(suffixes)
    indexes = []
    for item in sorted_array:
        indexes.append(str(suffixes[item]))
    return ', '.join(indexes)

# Returns Burrows - Wheeler transform, i.e. the last column of the matrix, formed
# from all cyclic rotations of string
def BWT(string):
    
    rotations = []

    # generates all rotations
    for i in range(len(string)):
        rotations.append(string[i : ] + string[ : i])

    rotations =  sorted(rotations)

    # take last symbol from each string in rotations list, generating BWT
    word = ''    
    for item in rotations:
        word += item[-1]
        
    return word

# Returns first column after Burrows - Wheeler transform
def BWT_first_col(string):

    rotations = []

    # generates all rotations
    for i in range(len(string)):
        rotations.append(string[i : ] + string[ : i])

    rotations =  sorted(rotations)

    # take last symbol from each string in rotations list, generating BWT
    word = ''    
    for item in rotations:
        word += item[0]
        
    return word


# Returns original string from Burrows - Wheeler transform
def inverse_BWT(string):

    # returns first columns of matrix
    first_col = sorted(string)

    # add "subscripts" to each symbol based on the order of its appearance
    for i in range(len(first_col)):
        # if previous symbol is the same, increment subscrpit by one,
        # else set subscript = 1
        if i != 0:
            if first_col[i - 1][0] == first_col[i]:
                first_col[i] = first_col[i - 1][0] + str(int(first_col[i - 1][1 : ]) + 1)
            else:
                first_col[i] = first_col[i] + str(1)
        else:
            first_col[i] = '$1'

    # for each letter in string, if it is in letters, increment its subscript by 1
    # and add it with this subscript to the last_col, else use subsctrip = 1
    letters = {}
    last_col = []
    for letter in string:
        try:
            letters[letter] += 1
            last_col.append(letter + str(letters[letter]))
        except KeyError:
            letters[letter] = 1
            last_col.append(letter + str(letters[letter]))

    # take element from first_col at index, at which pointer appears in last_col
    original_string = ''
    pointer = first_col[last_col.index('$1')]
    for i in range(len(last_col)):
        original_string += pointer[0]
        pointer = first_col[last_col.index(pointer)]
    
    return original_string


# Given the original string, returns the last-to-first array, which gives
# a position of symbol in first column based on its position in the last column
def last_to_first(string):
    rotations = []
    for i in range(len(string)):
        rotations.append(string[i : ] + string[ : i])
    rotations =  sorted(rotations)
    first_col = []
    last_col = []
    for word in rotations:
        first_col.append(word[0])
        last_col.append(word[-1])

    for i in range(len(first_col)):
        # if previous symbol is the same, increment subscrpit by one,
        # else set subscript = 1
        if i != 0:
            if first_col[i - 1][0] == first_col[i]:
                first_col[i] = first_col[i - 1][0] + str(int(first_col[i - 1][1 : ]) + 1)
            else:
                first_col[i] = first_col[i] + str(1)
        else:
            first_col[i] = '$1'
            
    letters = {}
    for i in range(len(last_col)):
        # if previous symbol is the same, increment subscrpit by one,
        # else set subscript = 1
        try:
            letters[last_col[i]] += 1
            last_col[i] = last_col[i] + str(letters[last_col[i]])
        except KeyError:
            letters[last_col[i]] = 1
            last_col[i] = last_col[i] + str(letters[last_col[i]])
            
##    LastToFirst = {}
##    for item in last_col:
##        LastToFirst[item] = first_col.index(item)
    LastToFirst = []
    for item in last_col:
        LastToFirst.append(first_col.index(item))
    return LastToFirst


def BWMatching(last_col_string, patterns):
    lastToFirst = last_to_first(inverse_BWT(last_col_string))
    number_of_occurences = []
    for pattern in patterns:
        top = 0
        top_ind = 0
        bottom = len(last_col_string) - 1
        while top <= bottom:
            substring = last_col_string[top : bottom + 1]
            if pattern:
                symbol = pattern[-1]
                pattern = pattern[ : -1]
                if symbol in substring:
                    top_index = substring.index(symbol) + top_ind
                    bottom_index = substring.rindex(symbol) + top_ind
                    new_range = lastToFirst[top_index : bottom_index + 1]
                    top = new_range[0]
                    bottom = new_range[-1]
                    top_ind = top
                else:
                    number_of_occurences.append(str(0))
                    break
            else:
                number_of_occurences.append(str(bottom - top + 1))
                break
    return ' '.join(number_of_occurences)

from collections import Counter
def Count(symbol, index, last_col):
    c = Counter(last_col[:index])
    return c[symbol]
    

def BetterBWMatching(last_col_string, patterns):
    first_col = BWT_first_col(inverse_BWT(last_col_string))
    number_of_occurences = []
    for pattern in patterns:
        top = 0
        bottom = len(last_col_string) - 1
        while top <= bottom:
            if pattern:
                symbol = pattern[-1]
                pattern = pattern[ : -1]
                if symbol in last_col_string[top : bottom + 1]:
                    top = first_col.find(symbol) + last_col_string.count(symbol, 0, top)
                    bottom = first_col.find(symbol) + last_col_string.count(symbol, 0, bottom + 1) - 1
                else:
                    number_of_occurences.append(str(0))
                    break
            else:
                number_of_occurences.append(str(bottom - top + 1))
                break
    return ' '.join(number_of_occurences)

# something wrong, some positions are missing and some positions are wrong
def multiple_pattern_matching(string, patterns):
    last_col = BWT(string)
    first_col = BWT_first_col(string)
    array = suffix_array(string).split(', ')
    positions = []
    for pattern in patterns:
        top = 0
        bottom = len(last_col) - 1
        while top <= bottom:
            if pattern:
                symbol = pattern[-1]
                pattern = pattern[ : -1]
                if symbol in last_col[top : bottom + 1]:
                    top = first_col.find(symbol) + last_col.count(symbol, 0, top)
                    bottom = first_col.find(symbol) + last_col.count(symbol, 0, bottom + 1) - 1
                else:
                    break
            else:
                
                positions.extend(array[top : bottom + 1])
                break
    #positions = [str(item) for item in sorted([int(position) for position in positions])]
    return ' '.join(sorted(positions))

import DNA

def multiple_approx_matching(file):
    a = open(file, 'r')
    strings = a.readlines()
    a.close()
    string = strings[0].rstrip('\n')
    patterns = strings[1].rstrip('\n')
    d = int(strings[2].rstrip('\n'))
    patterns = patterns.split(' ')
    positions = []
    for pattern in patterns:
        b = DNA.ApproximatePatternMatching(string, pattern, d)
        b = b.split(' ')
        positions.extend(b)
    #positions = [str(position) for position in sorted([int(position) for position in positions])]
    return ' '.join(positions)
    
