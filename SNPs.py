_end = '_end_'

def make_trie(words):
    root = dict()
    for word in words:
        current_dict = root
        for letter in word:
            current_dict = current_dict.setdefault(letter, {})
        current_dict = current_dict.setdefault(_end, _end)
    return root

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
