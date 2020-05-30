import sys
import re
import string
from random import seed
from random import randint
from random import choice

# スコア情報設定
class scoring:
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

class pairing:
    def __init__(self,x,y,params = scoring(1, -1, -2, -1)):
        self.x = x
        self.y = y
        self.gap = params.gap
        self.gap_ex = params.gap_ex
        self.matching = params.matching

    # aligment 行列の作成とaffine-scoreの計算
    def align(self,x,y): 
        # X, Yはアラインメントの対象になるシークエンス
        mat = matrix(x,y)
        #params = scoring(1, -1, -2, -1) #指定されたパラメーター

        # dimensions
        i_dim = len(x)+1
        j_dim = len(y)+1

        # DP行列の作成
        self.M = mat.generate(1+len(x),1+len(y))
        self.Ix = mat.generate(1+len(x),1+len(y))
        self.Iy = mat.generate(1+len(x),1+len(y))

        # traceback行列形成
        self.M_trace = mat.generate(1+len(x),1+len(y))
        self.Ix_trace = mat.generate(1+len(x),1+len(y))
        self.Iy_trace = mat.generate(1+len(x),1+len(y))

        # 生成された行列の初期化  (affineベースケース)
        for i in range(1, i_dim):
            self.M[i][0] = self.gap + self.gap_ex*(i-1)
            self.Ix[i][0] = -float('inf')
            self.Iy[i][0] = -float('inf')
        for j in range(1, j_dim):
            self.M[0][j] = self.gap + self.gap_ex*(j-1)
            self.Ix[0][j] = -float('inf')
            self.Iy[0][j] = -float('inf')


        for i in range(1, i_dim):
            for j in range(1, j_dim):
                # 生成された行列に値を入力 (DPとtraceback)
                
                self.Ix[i][j] = max([self.gap + self.M[i-1][j], self.gap_ex + self.Ix[i-1][j]])
                self.Ix_trace[i][j] = [self.gap + self.M[i-1][j], self.gap_ex + self.Ix[i-1][j]].index(self.Ix[i][j])

                self.Iy[i][j] = max([self.gap + self.M[i][j-1], self.gap_ex + self.Iy[i][j-1]])
                self.Iy_trace[i][j] = [self.gap + self.M[i][j-1], self.gap_ex + self.Iy[i][j-1]].index(self.Iy[i][j])

                self.M[i][j]  = max([self.matching(x[i-1], y[j-1]) + self.M[i-1][j-1], self.Ix[i][j], self.Iy[i][j]])
                self.M_trace[i][j] = [self.matching(x[i-1], y[j-1]) + self.M[i-1][j-1], self.Ix[i][j], self.Iy[i][j]].index(self.M[i][j])
                
        matrix_score = [self.Ix[i][j], self.Iy[i][j], self.M[i][j]]
        af_score  = max(matrix_score)
        print(af_score)

        # backtrace用スコア
        self.backtrace = matrix_score.index(af_score)
        
        return af_score

    def backtracing(self, x, y):
        # 元のsequenceからアラインメント開始
        gen_x, gen_y = x, y

        i, j = len(x), len(y)

        # backtracingから求めたアラインメントの作成
        while i>0 and j>0:
            if self.backtrace == 0:
                if self.Ix_trace[i][j] == 0:
                    self.backtrace = 2
                gen_y = gen_y[:j] + '_' +  gen_y[j:]
                i -= 1
            elif self.backtrace == 1:
                if self.Iy_trace[i][j] == 0:
                    self.backtrace = 2
                gen_x = gen_x[:i] + '_' + gen_x[i:]
                j -= 1
            elif self.backtrace == 2:
                if self.M_trace[i][j] == 1:
                    self.backtrace = 0
                elif self.M_trace[i][j] == 2:
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

class data():
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
    
    def DNA(self, length):
        # 引数に指定された長さのDNA sequenceを作成
        return ''.join(choice('CGTA') for _ in range(length))


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
    length = randint(0, 100)

    # FASTAの標準入力
    inFile = sys.argv[1]

    # dataメソッドでDNA配列の作成と受け取り (前処理)
    sequence = data()

    # randomでDNA sequenceのFAファイル作成
    with open(inFile, "w") as f:
        for s in range(10):
            f.write(sequence.DNA(length))
            f.write("\n\n")

    # FAファイル解析
    with open(inFile, "r") as f:
        for each_seq in sequence.parse_fa(f):
            seq.append(each_seq) # seqのリスト

    # 入力 sequence
    x, y = seq[randint(0,len(seq))], seq[randint(0,len(seq))]
    print(x, "\n\n", y)

    # pairingメソッドを呼び出して、アラインメント始める
    # paramsはスコアのパラメーター: scoring(match, mismatch, gap_start_penalty, gap_extend_penalty)
    sequence = pairing(x, y, params=scoring(1, -1, -2, -1)) 
    # affine アラインメントスコア
    af_score = sequence.align(x,y)
    # backtracingから作成されたアラインメント結果
    [res1, res2] = sequence.backtracing(x, y)

    print(res1)
    print(res2)