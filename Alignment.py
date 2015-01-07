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



def LCSBacktrackGlobal(v, w, sigma = -5, Scoring_Matrix = BLOSUM62):

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


def GlobalAlignment(v, w, sigma = -5, Scoring_Matrix = BLOSUM62):
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



def LCSBacktrackLocal(v, w, sigma = -5, Scoring_Matrix = PAM250):

    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    
##    for i in range(1, len(s)):
##        s[i][0] = s[i - 1][0] + sigma
##    for i in range(1, len(s[0])):
##        s[0][i] = s[0][i-1] + sigma
    
    # adds additional row and column
    for i in range(1, len(backtrack)):
        backtrack[i][0] = 'down'
    for i in range(1, len(backtrack[0])):
        backtrack[0][i] = 'right'
        
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[i][j] = max(0, s[i-1][j] + sigma, s[i][j-1] + sigma, s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]])

            if s[i][j] == s[i-1][j] + sigma:
                backtrack[i][j] = 'down'
            elif s[i][j] == s[i][j-1] + sigma:
                backtrack[i][j] = 'right'
            elif s[i][j] == s[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]:
                backtrack[i][j] = 'diagonal'
                
    return backtrack, s


def LocalAlignment(v, w, sigma = -5, Scoring_Matrix = PAM250):

    backtrack_s = LCSBacktrackLocal(v, w)
    backtrack = backtrack_s[0]
    s = backtrack_s[1]

    max_i = 0
    max_j = 0
    maximum_score = -float('inf')
    for i in range(len(s)):
        for j in range(len(s[0])):
            if s[i][j] >= maximum_score:
                maximum_score = s[i][j]
                max_i = i
                max_j = j
                
    score = maximum_score
    i = max_i
    j = max_j
    v_string = ''
    w_string = ''
    
    while backtrack[i][j]:
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
            
    return v_string[::-1], w_string[::-1], maximum_score

# Edit distance is the minimum number of single character operations (deletion,
# insertion, substitution) required to change first string into second
def EditDistance(v, w):
    
    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    
    for i in range(1, len(s)):
        s[i][0] = i
    for i in range(1, len(s[0])):
        s[0][i] = i

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            if v[i - 1] == w[j - 1]:
                s[i][j] = s[i-1][j-1]
            else:
                s[i][j] = min(s[i-1][j] + 1, s[i][j-1] + 1, s[i-1][j-1] + 1)
                

    return s[i][j]


# Find substring v' in v which maximizes global alignment score with w
def FittingAlignment(v, w):
    
    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack = [[0]*(len(w) + 1) for k in range(len(v) + 1)]

    for i in range(1, len(s[0])):
        s[0][i] = -i

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            if v[i-1] == w[j-1]:
                s[i][j] = max(s[i-1][j] - 1, s[i][j-1] - 1, s[i-1][j-1] + 1)
            else:
                s[i][j] = max(s[i-1][j] - 1, s[i][j-1] - 1, s[i-1][j-1] - 1)
                
            if s[i][j] == s[i-1][j] - 1:
                backtrack[i][j] = 'down'
            elif s[i][j] == s[i][j-1] - 1:
                backtrack[i][j] = 'right'
            elif s[i][j] == s[i-1][j-1] + 1:
                backtrack[i][j] = 'diagonal'
          
    max_score = -float('inf')
    max_i = 0
    for i in range(1, len(s)):
        if s[i][-1] > max_score:
            max_score = s[i][-1]
            max_i = i

    v_string = ''
    w_string = ''
    i = max_i
    j = len(w)
    while j > 0:
        if backtrack[i][j] == 'down':
            v_string += v[i - 1]
            w_string += '-'
            i -= 1
        elif backtrack[i][j] == 'right':
            v_string += '-'
            w_string += w[j - 1]
            j -= 1
        else:
            v_string += v[i - 1]
            w_string += w[j - 1]
            i -= 1
            j -= 1

    return v_string[::-1], w_string[::-1], max_score


# Return score of optimal overlap alignment, suffix v and prefix w, which maximize
# this score
def OverlapAlignment(v, w):

    s = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            if v[i-1] == w[j-1]:
                s[i][j] = max(s[i-1][j] - 2, s[i][j-1] - 2, s[i-1][j-1] + 1)
            else:
                s[i][j] = max(s[i-1][j] - 2, s[i][j-1] - 2, s[i-1][j-1] - 2)
                
            if s[i][j] == s[i-1][j] - 2:
                backtrack[i][j] = 'down'
            elif s[i][j] == s[i][j-1] - 2:
                backtrack[i][j] = 'right'
            elif s[i][j] == s[i-1][j-1] + 1:
                backtrack[i][j] = 'diagonal'
  
    max_score = -float('inf')
    max_j = 0
    for j in range(1, len(s[0])):
        if s[-1][j] >= max_score:
            max_score = s[-1][j]
            max_j = j


    v_string = ''
    w_string = ''
    j = max_j
    i = len(v)
    while j > 0:
        if backtrack[i][j] == 'down':
            v_string += v[i - 1]
            w_string += '-'
            i -= 1
        elif backtrack[i][j] == 'right':
            v_string += '-'
            w_string += w[j - 1]
            j -= 1
        else:
            v_string += v[i - 1]
            w_string += w[j - 1]
            i -= 1
            j -= 1

    return v_string[::-1], w_string[::-1], max_score


def AffineGapAlignment(v, w, sigma = -11, epsilon = -1, Scoring_Matrix = BLOSUM62):

    lower = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    upper = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    middle = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    
    backtrack_l = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack_u = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack_m = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    
    for i in range(1, len(lower)):
        lower[i][0] = (i) * sigma
        upper[i][0] = (i) * sigma
        middle[i][0] = i  * sigma        
    for i in range(1, len(lower[0])):
        lower[0][i] = (i) * sigma
        upper[0][i] = (i) * sigma
        middle[0][i] = (i) * sigma
        
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            lower[i][j] = max(lower[i-1][j] + epsilon, middle[i-1][j] + sigma)
##            if lower[i][j] == lower[i-1][j] + epsilon:
##                backtrack_l[i][j] = 'down'
##            if lower[i][j] == middle[i-1][j] + sigma:
##                backtrack_l[i][j] = 'diagonal'
            
##            if lower[i-1][j] + epsilon > middle[i-1][j] + sigma:
##                backtrack[i][j] = 'lower'
##            else:
##                backtrack[i][j] = 'middle'
            
            upper[i][j] = max(upper[i][j-1] + epsilon, middle[i][j-1] + sigma)
##            if upper[i][j] == upper[i][j-1] + epsilon:
##                backtrack_u[i][j] = 'right'
##            if upper[i][j] == middle[i][j-1] + sigma:
##                backtrack_u[i][j] = 'diagonal'
            
##            if upper[i][j-1] + epsilon > middle[i][j-1] + sigma:
##                backtrack[i][j] = 'upper'
##            else:
##                backtrack[i][j] = 'middle'
            middle[i][j] = max(lower[i][j], upper[i][j], middle[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]])
            
##            if middle[i][j] == middle[i-1][j-1] + Scoring_Matrix[v[i-1]][w[j-1]]:
##                backtrack_m[i][j] = 'diagonal'
##            if middle[i][j] == lower[i][j]:
##                backtrack_m[i][j] = 'down'
##            if middle[i][j] == upper[i][j]:
##                backtrack_m[i][j] = 'right'
                
##            if lower[i][j] > middle[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]:
##                if lower[i][j] > upper[i][j]:
##                    backtrack[i][j] = 'lower'
##                else:
##                    backtrack[i][j] = 'upper'
##            else:
##                if middle[i][j] > upper[i][j]:
##                    backtrack[i][j] = 'middle'
##                else:
##                    backtrack[i][j] = 'upper'
            
##            if lower[i][j] > upper[i][j]:
##                if middle[i][j] > lower[i][j]:
##                    backtrack[i][j] = 'middle'
##                else:
##                    backtrack[i][j] = 'lower'
##            elif middle[i][j] > upper[i][j]:
##                backtrack[i][j] = 'middle'
##            else:
##                backtrack[i][j] = 'upper'
    #print(lower, upper, middle)
    i = len(v)
    j = len(w)
    v_string = ''
    w_string = ''
    max_score = max(lower[len(v)][len(w)], upper[len(v)][len(w)], middle[len(v)][len(w)])
##    if lower[len(v)][len(w)] == max_score:
##        t = 'down'
##    elif upper[len(v)][len(w)] == max_score:
##        t = 'right'
##    else:
##        t = 'diagonal'
    while i > 0 and j > 0:
##        print(t)
##        if t == 'down':
##            v_string += v[i-1]
##            w_string += '-'
##            i -= 1
##            t = backtrack_l[i][j]
##        elif t == 'right':
##            v_string += '-'
##            w_string += w[j-1]
##            j -= 1
##            t = backtrack_u[i][j]
##        else:
##            v_string += v[i-1]
##            w_string += w[j-1]
##            i -= 1
##            j -= 1
##            t = backtrack_m[i-1][j-1]
        
        
        #print(lower[i][j], upper[i][j], middle[i-1][j-1])
        a = max(lower[i][j], upper[i][j], middle[i-1][j-1])
        if lower[i][j] == a:# and lower[i][j] != middle[i][j]:
            v_string += v[i - 1]
            w_string += '-'
            i -= 1
        if upper[i][j] == a:#and upper[i][j] != middle[i][j]:
            v_string += '-'
            w_string += w[j - 1]
            j -= 1
        if middle[i-1][j-1] == a:
            v_string += v[i - 1]
            w_string += w[j - 1]
            i -= 1
            j -= 1
           
        #print(v_string, w_string)
    return v_string[::-1], w_string[::-1], max_score

# l <-> x, u <-> y, m <-> d
def AffineGapAlignment2(v, w, sigma = -11, epsilon = -1, Scoring_Matrix = BLOSUM62):

    lower = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    upper = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    middle = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    
    backtrack_l = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack_u = [[0]*(len(w) + 1) for k in range(len(v) + 1)]
    backtrack_m = [[0]*(len(w) + 1) for k in range(len(v) + 1)]

    lower[0][0] = sigma + (len(v) + len(w))*epsilon
    upper[0][0] = sigma + (len(v) + len(w))*epsilon
    
    for i in range(1, len(lower)):
        lower[i][0] = sigma + (len(v) + len(w)) * epsilon
        upper[i][0] = sigma + (i - 1) * epsilon
        middle[i][0] = sigma + (len(v) + len(w)) * epsilon
        backtrack_u[i][0] = 'upper'
        
    for i in range(1, len(lower[0])):
        lower[0][i] = sigma + (i - 1) * epsilon
        upper[0][i] = sigma + (len(v) + len(w)) * epsilon
        middle[0][i] = sigma + (len(v) + len(w)) * epsilon
        backtrack_l[0][i] = 'lower'
        
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            lower[i][j] = max(lower[i][j-1] + epsilon, middle[i][j-1] + sigma)
            if lower[i][j] == middle[i][j-1] + sigma:
                backtrack_l[i][j] = 'middle'
            else:
                backtrack_l[i][j] = 'lower'
            upper[i][j] = max(upper[i-1][j] + epsilon, middle[i-1][j] + sigma)
            if upper[i][j] == middle[i-1][j] + sigma:
                backtrack_u[i][j] = 'middle'
            else:
                backtrack_u[i][j] = 'upper'
            middle[i][j] = max(lower[i][j], upper[i][j], middle[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]) ###
            if middle[i][j] == middle[i-1][j-1] + Scoring_Matrix[v[i - 1]][w[j - 1]]:
                backtrack_m[i][j] = 'middle'
            if middle[i][j] == upper[i][j]:# + Scoring_Matrix[v[i - 1]][w[j - 1]]:
                backtrack_m[i][j] = 'upper'
            if middle[i][j] == lower[i][j]:# + Scoring_Matrix[v[i - 1]][w[j - 1]]:
                backtrack_m[i][j] = 'lower'
                
    max_score = max(lower[len(v)][len(w)], upper[len(v)][len(w)], middle[len(v)][len(w)])

    if middle[len(v)][len(w)] == max_score:
        t = 'middle'
    if lower[len(v)][len(w)] == max_score:
        t = 'lower'
    if upper[len(v)][len(w)] == max_score:
        t = 'upper'
        
    i = len(v)
    j = len(w)
    v_string = ''
    w_string = ''
    while i > 0 and j > 0:
        print(t)
        if t == 'middle':
            v_string += v[i-1]
            w_string += w[j-1]
            i -= 1
            j -= 1
            t = backtrack_m[i][j]
        elif t == 'upper':
            v_string += v[i-1]
            w_string += '-'
            i -= 1
            t = backtrack_u[i][j]
        elif t == 'lower':
            v_string += '-'
            w_string += w[j-1]
            j -= 1
            t = backtrack_l[i][j]
    return v_string[::-1], w_string[::-1]
##    i = len(v)
##    j = len(w)
##    v_string = ''
##    w_string = ''
##
##    
##    while i > 0 and j > 0:
##
##        print(lower[i][j], upper[i][j], middle[i][j])
##        a = max(lower[i][j], upper[i][j], middle[i][j])
##        if lower[i][j] == a:
##            v_string += v[i - 1]
##            w_string += '-'
##            i -= 1
##        elif upper[i][j] == a:
##            v_string += '-'
##            w_string += w[j - 1]
##            j -= 1
##        elif middle[i][j] == a:
##            v_string += v[i - 1]
##            w_string += w[j - 1]
##            i -= 1
##            j -= 1
##
##    return v_string[::1], w_string[::1], max_score


def MultipleLCS(v, w, u):

    s = [[[0]*(len(u) + 1) for k in range(len(w) + 1)] for i in range(len(v) + 1)]
    backtrack = [[[0]*(len(u) + 1) for k in range(len(w) + 1)] for i in range(len(v) + 1)]
    
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            for k in range(1, len(u) + 1):
                if v[i-1] == w[j-1] and v[i-1] == u[k-1]:
                    s[i][j][k] = max(s[i-1][j][k], s[i][j-1][k], s[i][j][k-1], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k-1], s[i-1][j-1][k-1] +1)
                else:
                    s[i][j][k] = max(s[i-1][j][k], s[i][j-1][k], s[i][j][k-1], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k-1], s[i-1][j-1][k-1])
##                if s[i][j][k] == s[i-1][j-1][k-1] and ((v[i-1] != w[j-1] and v[i-1] != u[k-1] and w[j-1] != u[k-1]) or
##                                                         (v[i-1] == w[j-1] and v[i-1] != u[k-1]) or
##                                                         (v[i-1] == u[k-1] and v[i-1] != w[j-1]) or
##                                                         (w[j-1] == u[k-1] and w[j-1] != v[i-1])):
##                    backtrack[i][j][k] = 'mismatch'
                if s[i][j][k] == s[i-1][j][k] and v[i-1] != w[j-1] and v[i-1] != u[k-1]:
                    backtrack[i][j][k] = 'second and third'
                elif s[i][j][k] == s[i][j-1][k] and w[j-1] != u[k-1] and v[i-1] != w[j-1]:
                    backtrack[i][j][k] = 'first and third'
                elif s[i][j][k] == s[i][j][k-1] and v[i-1] != u[k-1] and w[j-1] != u[k-1]:
                    backtrack[i][j][k] = 'first and second'
                elif s[i][j][k] == s[i-1][j-1][k-1] + 1:#and v[i-1] == w[j-1] and v[i-1] == u[k-1]:
                    backtrack[i][j][k] = 'all'
                elif s[i][j][k] == s[i-1][j-1][k] and v[i-1] == w[j-1]:
                    backtrack[i][j][k] = 'third'
                elif s[i][j][k] == s[i-1][j][k-1] and v[i-1] == u[k-1]:
                    backtrack[i][j][k] = 'second'
                elif s[i][j][k] == s[i][j-1][k-1] and w[j-1] == u[k-1]:
                    backtrack[i][j][k] = 'first'
                print(i, j, k)
                print(s[i][j][k], s[i-1][j][k], s[i][j-1][k], s[i][j][k-1], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k-1], s[i-1][j-1][k-1], s[i-1][j-1][k-1] + 1)
                print(backtrack[i][j][k])
    #print(backtrack)
    #return len(backtrack), len(backtrack[0]), len(backtrack[0][0])
    i = len(v)
    j = len(w)
    k = len(u)
    v_string = ''
    w_string = ''
    u_string = ''
    while i > 0 and j >0 and k > 0:
        print(i, j, k)
        print(backtrack[i][j][k])
        if backtrack[i][j][k] == 'second and third':
            v_string += v[i-1]
            w_string += '-'
            u_string += '-'
            i -= 1
        elif backtrack[i][j][k] == 'first and third':
            w_string += w[j-1]
            v_string += '-'
            u_string += '-'
            j -= 1
        elif backtrack[i][j][k] == 'first and second':
            u_string += u[k-1]
            v_string += '-'
            w_string += '-'
            k -= 1
        elif backtrack[i][j][k] == 'third':
            v_string += v[i-1]
            w_string += w[j-1]
            u_string += '-'
            i -= 1
            j -= 1
        elif backtrack[i][j][k] == 'second':
            v_string += v[i-1]
            u_string += u[k-1]
            w_string += '-'
            i -= 1
            k -= 1
        elif backtrack[i][j][k] == 'first':
            w_string += w[j-1]
            u_string += u[k-1]
            v_string += '-'
            j -= 1
            k -= 1
        else:
            v_string += v[i-1]
            w_string += w[j-1]
            u_string += u[k-1]
            i -= 1
            j -= 1
            k -= 1
    return v_string[::-1], w_string[::-1], u_string[::-1]
