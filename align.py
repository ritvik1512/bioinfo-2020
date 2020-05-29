import sys
import re
from random import seed
from random import randint

# スコア情報設定
class score:
    def __init__(self, match, mismatch, gap, gap_ex):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap # gap開始ペナルティー
        self.gap_ex = gap_ex # gap伸長ペナルティー

    def matching(self, X, Y):
        # match確認
        if X == Y:
            return self.match
        else:
            return self.mismatch

# FASTA入力の解析
def parse_fa(input):
    head = None
    sequence =  []

    for line in input:
        line = line.rstrip() # \n などをとる
        if regex.search(line) or line == "\n":
            if head: yield (head, "".join(sequence)) # 溜まったsequenceリストをくっつく
            head, sequence = line, [] # seq毎でリセット
        else:
            sequence.append(line)
    if head: yield(head, "".join(sequence)) # 最後のsequenceのデータを出力

def align(X,Y):
    X_seq = str("First seq: "+ X +"\n")
    Y_seq = str("Second seq: "+ Y +"\n\n")

    with open("man_output.txt", "w") as out:
        out.write(X_seq)
        out.write(Y_seq)

if __name__ == "__main__":
    # 初期化
    header = []
    seq = []
    seed(1)
    regex = re.compile(r'>')

    # FASTAの標準入力
    inFile = sys.argv[1]  
    with open(inFile, "r") as f:
        for head,sequence in parse_fa(f):
            header.append(head) # headerのリスト
            seq.append(sequence) # seqのリスト
    
    align(seq[randint(0, len(seq))], seq[randint(0, len(seq))])