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




class Node:
        def __init__(self, start, substr):
                self.start = start
                self.substr = substr
                self.branches = {}
              
def insert_into_tree(subroot, suffix, start):
        prefix_len = len(subroot.substr)
        new_suffix = str(suffix[prefix_len:])
        if(len(subroot.branches) == 0):
                left_child = Node(subroot.start, "")
                right_child = Node(start, new_suffix)
                subroot.branches[""] = left_child
                subroot.branches[new_suffix] = right_child
        else:
                for (substr, node) in subroot.branches.items():
                        if len(substr) > 0 and new_suffix.startswith(substr):
                                insert_into_tree(node, new_suffix, start)
                                break
                else:
                        new_child = Node(start, new_suffix)
                        subroot.branches[new_suffix] = new_child
              
def build_suffix_tree(t_str):
        len_str = len(t_str)
        i = len_str - 1
        root = Node(len_str, "")
        while i >= 0:
                insert_into_tree(root, str(t_str[i:]), i)
                i -= 1
        return root
              
def display_all_suffix(subroot, suffix_s_prefix, level = 0):
        if len(subroot.branches) == 0:
                print(suffix_s_prefix, level)
                return
        for (substr, node) in subroot.branches.items():
                display_all_suffix(node, suffix_s_prefix + substr, level + 1)
              
root = build_suffix_tree("ATAAATG$")
display_all_suffix(root, "")
