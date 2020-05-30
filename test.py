from pprint import pprint

global S
global E
global match
global mismatch
global MIN
S   = -2.
E   = -1
match = 1.
mismatch = -1.
MIN = -float("inf")

#return match or mismatch score
def _match(s, t):
    if t == s:
        return match
    else:
        return mismatch

def generate(x, y):
        return [[0]*(y) for i in range(x)]

def distance_matrix(s, t):
    dim_i = len(t) + 1 # x
    dim_j = len(s) + 1 # y

    # DP行列の生成
    M = generate(1+len(t),1+len(s))
    X = generate(1+len(t),1+len(s))
    Y = generate(1+len(t),1+len(s))

    # 生成された行列の初期化  (affineベースケース)
    for i in range(1, dim_i):
        M[i][0] = -float('inf')
        X[i][0] = -float('inf')
        Y[i][0] = S + i*E

    for j in range(1, dim_j):
        M[0][j] = -float('inf')
        X[0][j] = S + j*E
        Y[0][j] = -float('inf')

    for i in range(1, dim_i):
        for j in range(1, dim_j):
            X[i][j] = max(S + E + M[i][j-1], E + X[i][j-1])
            Y[i][j] = max(S + E + M[i-1][j], E + Y[i-1][j])
            M[i][j] = max(_match(s[j-1], t[i-1]) + M[i-1][j-1], X[i][j], Y[i][j])


    af_score = max(M[dim_i-2][dim_j-2], X[dim_i-2][dim_j-2], Y[dim_i-2][dim_j-2])
    print(af_score)
    return [X, Y, M]

def backtrace(s, t, X, Y, M):
    sequ1 = ''
    sequ2 = ''
    i = len(t) # 2nd arg
    j = len(s) # 1st arg
    while (i>0 or j>0):
        if (i>0 and j>0 and M[i][j] == M[i-1][j-1] + _match(s[j-1], t[i-1])):
            sequ1 += s[j-1]
            sequ2 += t[i-1]
            i -= 1; j -= 1
        elif (i>0 and M[i][j] == Y[i][j]):
            sequ1 += '_'
            sequ2 += t[i-1]
            i -= 1
        elif (j>0 and M[i][j] == X[i][j]):
            sequ1 += s[j-1]
            sequ2 += '_'
            j -= 1

    sequ1r = ' '.join([sequ1[j] for j in range(-1, -(len(sequ1)+1), -1)])
    sequ2r = ' '.join([sequ2[j] for j in range(-1, -(len(sequ2)+1), -1)])

    return [sequ1r, sequ2r]

seq1 = "CTCCCCAGACGAAGTGGATCACGTCTAGTTAACAGAAGTTCACGATACACTAGGGCGGATTATCAGGACATGAAT"
seq2 = "TCTACGGGCTCCCCACTGTCGTTCGGTGTTATAACCTATAGTCTGACGCAGCG"

[X, Y, M] = distance_matrix(seq1, seq2)
#print(X, Y, M, "\n")
[str1, str2] = backtrace(seq1, seq2, X, Y, M)
print("-=Alignment=-")
print(str1)
print(str2)
print("\n")