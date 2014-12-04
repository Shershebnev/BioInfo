def FastaToString(file):

    a = open(file, 'r')

    genome = ""

    string = a.readline().rstrip('\n')

    while string:

        genome += string
        string = a.readline().rstrip('\n')

    a.close()
    return genome
