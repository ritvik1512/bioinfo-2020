import sys
import re
import string
from random import seed
from random import randint
from random import choice

# スコア情報設定
class score:
    def __init__(self, match, mismatch, gap, gap_ex):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap # gap開始ペナルティー
        self.gap_ex = gap_ex # gap伸長ペナルティー

    def matching(self, x, y):
        # match確認
        if x == y:
            return self.match
        else:
            return self.mismatch

class matrix:
    def __init__(self, x, y):
        self.x = len(x)
        self.y = len(y)
    
    def generate(self, x, y):
        return [[0]*(y) for i in range(x)]

    # def initialize(self, x, y):
    #     for i in range(1, x+1):
    #         z1 = self.generate(x, y)
    #         z[i][0]  = -float('inf')

class pairing:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    # aligment 行列の作成とaffine-scoreの計算
    def align(self,x,y): 
        # X, Yはアラインメントの対象になるシークエンス
        mat = matrix(x,y)
        params = score(1, -1, -2, -1) #指定されたパラメーター

        # dimensions
        i_dim = len(y)+1 # t
        j_dim = len(x)+1 # s

        # DP行列の生成
        M = mat.generate(1+len(y),1+len(x))
        Ix = mat.generate(1+len(y),1+len(x))
        Iy = mat.generate(1+len(y),1+len(x))

        # 生成された行列の初期化  (affineベースケース)
        for i in range(1, i_dim):
            M[i][0] = -float('inf')
            Ix[i][0] = -float('inf')
            Iy[i][0] = params.gap + i*params.gap_ex

        for j in range(1, j_dim):
            M[0][j] = -float('inf')
            Ix[0][j] = params.gap + j*params.gap_ex
            Iy[0][j] = -float('inf')

        for i in range(1, i_dim):
            for j in range(1, j_dim):
                # score計算
                Ix[i][j] = max(params.gap + params.gap_ex + M[i][j-1], params.gap_ex + Ix[i][j-1])
                Iy[i][j] = max(params.gap + params.gap_ex + M[i-1][j], params.gap_ex + Iy[i-1][j])
                #M[i][j]  = params.matching(x[j-1], y[i-1]) + max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1])
                M[i][j]  = max(params.matching(x[j-1], y[i-1]) + M[i-1][j-1], Ix[i][j], Iy[i][j])
        
        af_score  = max(M[i_dim-1][j_dim-1], Ix[i_dim-1][j_dim-1], Iy[i_dim-1][j_dim-1])
        print(af_score)

        return [M, Ix, Iy], af_score

    def backtrace(self, x, y, M, Ix, Iy, params = score(1, -1, -2, -1)):
        align1, align2 = '', ''
        i = len(y) 
        j = len(x) # start from the bottom right cell

        while i>0 or j>0: # end toching the top or the left edge
            if (i>0 and j>0 and M[i][j] == M[i-1][j-1] + params.matching(x[j-1], y[i-1])):
                align1 += x[j-1]
                align2 += y[i-1]
                i -= 1
                j -= 1
            elif (i>0 and M[i][j] == Iy[i][j]):
                align1 += '_'
                align2 += y[i-1]
                i -= 1
            elif (j>0 and M[i][j] == Ix[i][j]):
                align1 += x[j-1]
                align2 += '_'
                j -= 1

        genseq1r = ' '.join([align1[j] for j in range(-1, -(len(align1)+1), -1)])
        genseq2r = ' '.join([align2[j] for j in range(-1, -(len(align2)+1), -1)])

        return [genseq1r, genseq2r]

class parse():
    # FASTA入力の解析
    def parse_fa(self, input):
        seq_l =  []

        for line in input:
            if line == "\n":
                continue
            else:
                line = line.strip()
                seq_l.append(line)
        return seq_l # 最後のsequenceのデータを出力


    # X_seq = str("First seq: "+ X +"\n")
    # Y_seq = str("Second seq: "+ Y +"\n\n")

    # with open("man_output.txt", "w") as out:
    #     out.write(X_seq)
    #     out.write(Y_seq)

def DNA(length):
        return ''.join(choice('CGTA') for _ in range(length))

if __name__ == "__main__":
    # 初期化
    header = []
    seq = []
    seed(1)

    # FASTAの標準入力
    #inFile = sys.argv[1]
    parser = parse()

    # DNAのFAファイル作成
    # with open(inFile, "w") as f:
    #     for s in range(10):
    #         f.write(DNA(randint(0, 100)))
    #         f.write("\n\n")

    # # FAファイル解析
    # with open(inFile, "r") as f:
    #     for sequence in parser.parse_fa(f):
    #         seq.append(sequence) # seqのリスト

    # x, y = seq[1], seq[2]
    # print(x, "\n\n", y)

    # 呼び出し
    x = "CTCCCCAGACGAAGTGGATCACGTCTAGTTAACAGAAGTTCACGATACACTAGGGCGGATTATCAGGACATGAAT"
    y = "TCTACGGGCTCCCCACTGTCGTTCGGTGTTATAACCTATAGTCTGACGCAGCG"

    sequence = pairing(x,y)
    [M, Ix, Iy], af_score = sequence.align(x,y)
    [res1, res2] = sequence.backtrace(x, y, M, Ix, Iy, params=score(1, -1, -2, -1))
    print(res1)
    print(res2)