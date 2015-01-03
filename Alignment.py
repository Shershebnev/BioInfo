from Scoring_matrixes import BLOSUM62 as BLOSUM62
from Scoring_matrixes import PAM250 as PAM250

def RecursiveChange(money, coins):
    if money == 0:
        return 0
    MinNumCoins = float('inf')
    for i in range(len(coins)):
        if money >= coins[i]:
            NumCoins = RecursiveChange(money - coins[i], coins)
            if NumCoins + 1 < MinNumCoins:
                MinNumCoins = NumCoins + 1
        
    return MinNumCoins


def DynamicChange(money, coins):
    MinNumCoins = [0]*(money+1)
    CoinsUsed = [[]]*(money+1)
    for i in range(1, money + 1):

        MinNumCoins[i] = float('inf')

        for j in range(len(coins)):
            if i >= coins[j]:
                if MinNumCoins[i - coins[j]] + 1 < MinNumCoins[i]:
                    MinNumCoins[i] = MinNumCoins[i - coins[j]] + 1

                    # These two strings are used to collect which coins exactly
                    # should be used for each amount of money
                    CoinsUsed[i] = list(CoinsUsed[i - coins[j]])
                    CoinsUsed[i].append(coins[j])
    return MinNumCoins[-1], CoinsUsed[-1]
        
# n == number of rows, m == number of columns
def ManhattanTourist(file):
    

    
    a = open(file, 'r')
    n_and_m = a.readline().rstrip('\n').split(' ')
    n = int(n_and_m[0])
    m = int(n_and_m[1])

    s = [[0]*(m+1) for k in range(n + 1)] # creates two-dimensional array
    down = [[0]*(m + 1) for k in range(n)]
    right = [[0] * m for k in range(n + 1)]
    
    weights = a.readline().rstrip('\n')
    row = 0
    while weights != '-':
        weights = weights.split(' ')
        weights = [int(x) for x in weights]
        down[row]= weights
        row += 1
        weights = a.readline().rstrip('\n')

    weights = a.readline().rstrip('\n')
    row = 0
    while weights:
        weights = weights.split(' ')
        weights = [int(x) for x in weights]
        right[row] = weights
        row += 1
        weights = a.readline().rstrip('\n')
        
    a.close()

    for i in range(1, n + 1):
        s[i][0] = s[i-1][0] + down[i-1][0] # for down index [i-1] for correction of i range

    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j-1] # for right same as for down

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])

    return s[n][m]

# len(v) == n umber of rows, len(w) == number of columns, 
# i - row, j - column
def LCSBacktrack(v, w):
    
    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack = [[0]*(len(w)) for k in range(len(v))]

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            if v[i - 1] == w[j - 1]:
                s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + 1)
            else:
                s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1])

            if s[i][j] == s[i-1][j]:
                backtrack[i-1][j-1] = 'down'
            elif s[i][j] == s[i][j-1]:
                backtrack[i-1][j-1] = 'right'
            elif s[i][j] == s[i-1][j-1] + 1:
                backtrack[i-1][j-1] = 'diagonal'
    return backtrack

### i and j should be len(sequence) - 1 if interested in alignment of whole sequence.
### This is due to indexes starting from 0. Causes IndexError in if/elif statements otherwise
##def OutputLCS(backtrack, v, i, j):
##    file = 'C:/Users/Alex/Desktop/Output.txt'
##    a = open(file, 'a')
##    if i == 0 or j == 0:
##        a.write(v[i])
##        return
##    if backtrack[i][j] == 'down':
##        OutputLCS(backtrack, v, i - 1, j)
##    elif backtrack[i][j] == 'right':
##        OutputLCS(backtrack, v, i, j - 1)
##    else:
##        OutputLCS(backtrack, v, i - 1, j - 1)
##        a.write(v[i])
##    a.close()

# i and j should be len(sequence) - 1 if interested in alignment of whole sequence.
# This is due to indexes starting from 0. Causes IndexError in if/elif statements otherwise
def OutputLCS(v, w):
    backtrack = LCSBacktrack(v, w)
    i = len(v) - 1
    j = len(w) - 1
    string = ''
    while i >= 0 and j >= 0:
        if backtrack[i][j] == 'down':
            i -= 1
        elif backtrack[i][j] == 'right':
            j -= 1
        else:
            string += v[i]
            i -= 1
            j -= 1
    return string[::-1] # [::-1] - reverses string (extended slices)



##def LCSBacktrackGlobal(v, w, Scoring_Matrix = BLOSUM62, sigma = -5):
##
##    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
##    backtrack = [[0]*(len(w)) for k in range(len(v))]
##
##    for i in range(1, len(s)):
##        s[i][0] = s[i - 1][0] + sigma
##    for i in range(1, len(s[0])):
##        s[0][i] = s[0][i-1] + sigma
##    
##    for i in range(1, len(v) + 1):
##        for j in range(1, len(w) + 1):
##            s[i][j] = max(s[i-1][j] + sigma, s[i][j-1] + sigma, s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]])
##
##            if s[i][j] == s[i-1][j] + sigma:
##                backtrack[i-1][j-1] = 'down'
##            elif s[i][j] == s[i][j-1] + sigma:
##                backtrack[i-1][j-1] = 'right'
##            elif s[i][j] == s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]:
##                backtrack[i-1][j-1] = 'diagonal'
##                
##    return backtrack#, s



def LCSBacktrackGlobal(v, w, Scoring_Matrix = BLOSUM62, sigma = -5):

    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack = [[0]*(len(w) + 1) for k in range(len(v) + 1)]

    for i in range(1, len(s)):
        s[i][0] = s[i - 1][0] + sigma
    for i in range(1, len(s[0])):
        s[0][i] = s[0][i-1] + sigma

    # adds additional row and column
    for i in range(1, len(backtrack)):
        backtrack[i][0] = 'down'
    for i in range(1, len(backtrack[0])):
        backtrack[0][i] = 'right'
        
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[i][j] = max(s[i-1][j] + sigma, s[i][j-1] + sigma, s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]])

            if s[i][j] == s[i-1][j] + sigma:
                backtrack[i][j] = 'down'
            elif s[i][j] == s[i][j-1] + sigma:
                backtrack[i][j] = 'right'
            elif s[i][j] == s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]:
                backtrack[i][j] = 'diagonal'
    return backtrack#, s


def GlobalAlignment(v, w, Scoring_Matrix = BLOSUM62, sigma = -5):
##    backtrack_s = LCSBacktrackGlobal(v, w)
##    backtrack = backtrack_s[0]
##    s = backtrack_s[1]
    backtrack = LCSBacktrackGlobal(v, w)
    i = len(v)# - 1
    j = len(w)# - 1
    score = 0
    v_string = ''
    w_string = ''
    while i > 0 or j > 0:
        if i > 0 and backtrack[i][j] == 'down':
            v_string += v[i - 1]
            w_string += '-'
            #s_score = s[i][j]
            i -= 1
            score += sigma
        elif j > 0 and backtrack[i][j] == 'right':
            v_string += '-'
            w_string += w[j - 1]
            #s_score = s[i][j]
            j -= 1
            score += sigma
        else:
            v_string += v[i - 1]
            w_string += w[j - 1]
            #s_score = s[i][j]
            score += Scoring_Matrix[v[i - 1]][w[j - 1]]
            i -= 1
            j -= 1

##    # For the case, when sequences start from different letter. This results in situation
##    # when during backtracking program reaches only one edge (either top or left) in the
##    # score table s (i.e. it reaches -5 on the top edge (see picture below)). 
##    #----------------
##    #| 0 | -5 | -10 |
##    #----------------
##    #| -5 | № |  №  |
##    #----------------
##    # In this case this code adds the rest of the corresponding string to one
##    # string and the same amount of '-' to another and updates the final score
##    # Maybe there is better option, but I couldn't find it right now.
##    # I tried adding extra row and column, which would correspond to first row
##    # and column in score table with all 'right' and 'down' respectively
##    # (see LCSBacktrackBlosum2 below), but the result was not the one I was
##    # expecting:
##    #>>> GlobalAlignment('PLEASANTLY', 'MEANLY')
##    #('PLEASANTLY', 'M-EAN--L-Y', -6)
##    # Maybe there is something wrong with indexing, but I don't have time to
##    # work this out now.
##    if s_score != 0:
##        if j != 0:
##            w_string += '-' * (s_score//sigma)
##            v_string += v[s_score//sigma - 1 ::-1]
##            score += s_score
##        elif i != 0:
##            v_string += '-' * (s_score//sigma)
##            w_string += w[s_score//sigma - 1 ::-1]
##            score += s_score

    return v_string[::-1], w_string[::-1], score



def LCSBacktrackLocal(v, w, Scoring_Matrix = PAM250, sigma = -5):

    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack = [[0]*(len(w)) for k in range(len(v))]

##    for i in range(1, len(s)):
##        s[i][0] = s[i - 1][0] + sigma
##    for i in range(1, len(s[0])):
##        s[0][i] = s[0][i-1] + sigma
    
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[i][j] = max(0, s[i-1][j] + sigma, s[i][j-1] + sigma, s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]])

            if s[i][j] == s[i-1][j] + sigma:
                backtrack[i-1][j-1] = 'down'
            elif s[i][j] == s[i][j-1] + sigma:
                backtrack[i-1][j-1] = 'right'
            elif s[i][j] == s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]:
                backtrack[i-1][j-1] = 'diagonal'
                
    return backtrack, s

##
##def LocalAlignment(v, w, sigma = -5, Scoring_matrix = PAM250):
##
##    
