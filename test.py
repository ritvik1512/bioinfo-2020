import sys
import re
import string
from random import seed
from random import randint

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

        # DP行列の生成
        M = mat.generate(1+len(x),1+len(y))
        Ix = mat.generate(1+len(x),1+len(y))
        Iy = mat.generate(1+len(x),1+len(y))

        # 生成された行列の初期化
        for i in range(1, i_dim):
            M[i][0] = -float('inf')
            Ix[i][0] = -float('inf')
            Iy[i][0] = params.gap + i*params.gap_ex # affine-gap計算

        for j in range(1, j_dim):
            M[0][j] = -float('inf')
            Ix[0][j] = params.gap + j*params.gap_ex # affine-gap計算
            Iy[0][j] = -float('inf')

        for i in range(1, i_dim):
            for j in range(1, j_dim):
                M[i][j]  = params.matching(x[i-1], y[j-1]) + max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1])
                Ix[i][j] = max(params.gap + params.gap_ex + M[i][j-1], params.gap_ex + Ix[i][j-1], params.gap + params.gap_ex + Iy[i][j-1])
                Iy[i][j] = max(params.gap + params.gap_ex + M[i-1][j], params.gap + params.gap_ex + Ix[i-1][j], params.gap_ex + Iy[i-1][j])
        
        af_score  = max(M[i_dim-1][j_dim-1], Ix[i_dim-1][j_dim-1], Iy[i_dim-1][j_dim-1])
        print(af_score)

        return [M, Ix, Iy], af_score

    def backtrace(self, x, y, M, Ix, Iy, params = score(1, -1, -2, -1)):
        genseq1, genseq2 = '', ''
        i = len(x)
        j = len(y)

        print(i, j)

        while (j>0 or i>0):
            if (j>0 and i>0 and M[i][j] == M[i-1][j-1] + params.matching(x[i-1], y[j-1])):
                genseq1 += x[i-1]
                genseq2 += y[j-1]
                i -= 1
                j -= 1
                print("first l", i, " ", j)
            elif (j>0 and M[i][j] == Iy[i][j]):
                genseq1 += '_'
                genseq2 += x[i-1]
                j -= 1
                print("second l", i, " ", j)
            elif (i>0 and M[i][j] == Ix[i][j]):
                genseq1 += y[j-1]
                genseq2 += '_'
                i -= 1
                print("third l", i, " ", j)

        genseq1r = ' '.join([genseq1[j] for j in range(-1, -(len(genseq1)+1), -1)])
        genseq2r = ' '.join([genseq2[j] for j in range(-1, -(len(genseq2)+1), -1)])

        return [genseq1r, genseq2r]

class parse():
    # FASTA入力の解析
    def parse_fa(self, input):
        head = None
        seq_l =  []

        for line in input:
            line = line.rstrip() # \n などをとる
            if regex.search(line) or line == "\n":
                if head: yield (head, "".join(seq_l)) # 溜まったsequenceリストをくっつく
                head, seq_l = line, [] # seq毎でリセット
            else:
                line = line.translate({ord(c): None for c in string.whitespace})
                seq_l.append(line)
        if head: yield(head, "".join(seq_l)) # 最後のsequenceのデータを出力



    # X_seq = str("First seq: "+ X +"\n")
    # Y_seq = str("Second seq: "+ Y +"\n\n")

    # with open("man_output.txt", "w") as out:
    #     out.write(X_seq)
    #     out.write(Y_seq)

if __name__ == "__main__":
    # 初期化
    header = []
    seq = []
    seed(1)
    regex = re.compile(r'>')

    # FASTAの標準入力
    inFile = sys.argv[1]
    parser = parse()

    with open(inFile, "r") as f:
        for head,sequence in parser.parse_fa(f):
            header.append(head) # headerのリスト
            seq.append(sequence) # seqのリスト
    x, y = seq[1], seq[2]
    print(x, "\n\n", y)

    # 呼び出し
    sequence = pairing(x,y)
    [M, Ix, Iy], af_score = sequence.align(x,y)
    [res1, res2] = sequence.backtrace(x, y, M, Ix, Iy, params=score(1, -1, -2, -1))
    print(res1, "\n\n")
    print(res2)