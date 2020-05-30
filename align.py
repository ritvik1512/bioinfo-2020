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
        i_dim = len(x)+1
        j_dim = len(y)+1

        # DP行列の作成
        M = mat.generate(1+len(x),1+len(y))
        Ix = mat.generate(1+len(x),1+len(y))
        Iy = mat.generate(1+len(x),1+len(y))

        # traceback行列形成
        M_trace = mat.generate(1+len(x),1+len(y))
        Ix_trace = mat.generate(1+len(x),1+len(y))
        Iy_trace = mat.generate(1+len(x),1+len(y))

        # 生成された行列の初期化  (affineベースケース)
        for i in range(1, i_dim):
            M[i][0] = params.gap + params.gap_ex*(i-1)
            Ix[i][0] = -float('inf')
            Iy[i][0] = -float('inf')
        for j in range(1, j_dim):
            M[0][j] = params.gap + params.gap_ex*(j-1)
            Ix[0][j] = -float('inf')
            Iy[0][j] = -float('inf')


        for i in range(1, i_dim):
            for j in range(1, j_dim):
                # 生成された行列に値を入力 (DPとtraceback)
                
                Ix[i][j] = max([params.gap + M[i-1][j], params.gap_ex + Ix[i-1][j]])
                Ix_trace[i][j] = [params.gap + M[i-1][j], params.gap_ex + Ix[i-1][j]].index(Ix[i][j])

                Iy[i][j] = max([params.gap + M[i][j-1], params.gap_ex + Iy[i][j-1]])
                Iy_trace[i][j] = [params.gap + M[i][j-1], params.gap_ex + Iy[i][j-1]].index(Iy[i][j])

                M[i][j]  = max([params.matching(x[i-1], y[j-1]) + M[i-1][j-1], Ix[i][j], Iy[i][j]])
                M_trace[i][j] = [params.matching(x[i-1], y[j-1]) + M[i-1][j-1], Ix[i][j], Iy[i][j]].index(M[i][j])
                
        matrix_score = [Ix[i][j], Iy[i][j], M[i][j]]
        af_score  = max(matrix_score)
        print(af_score)

        self.backtrace = matrix_score.index(af_score)
        
        return [M_trace, Ix_trace, Iy_trace], af_score

    def backtracing(self, x, y, M_trace, Ix_trace, Iy_trace, params = score(1, -1, -2, -1)):
        # 与えられたsequenceから
        gen_x, gen_y = x, y

        i = len(x) 
        j = len(y) # start from the bottom right cell

        # backtracingから求めたアラインメントの作成
        while i>0 and j>0:
            if self.backtrace == 0:
                if Ix_trace[i][j] == 0:
                    self.backtrace = 2
                gen_y = gen_y[:j] + '_' +  gen_y[j:]
                i -= 1
            elif self.backtrace == 1:
                if Iy_trace[i][j] == 0:
                    self.backtrace = 2
                gen_x = gen_x[:i] + '_' + gen_x[i:]
                j -= 1
            elif self.backtrace == 2:
                if M_trace[i][j] == 1:
                    self.backtrace = 0
                elif M_trace[i][j] == 2:
                    self.backtrace = 1
                else:
                    i -= 1
                    j -= 1

        # 何かが残されたら分を追加
        for left in range(i):
            gen_y = gen_y[:0] + '_' + gen_y[0:]
        for left in range(j):
            gen_x = gen_x[:0] + '_' + gen_x[0:]

        return [gen_x, gen_y]

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

    # x = "TTTAGCAA"
    # y = "TAGC"

    sequence = pairing(x,y)
    [M_trace, Ix_trace, Iy_trace], af_score = sequence.align(x,y)
    [res1, res2] = sequence.backtracing(x, y, M_trace, Ix_trace, Iy_trace, params=score(1, -1, -2, -1))
    print(res1)
    print(res2)