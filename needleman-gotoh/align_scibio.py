from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import re
from random import seed
from random import randint

regex = re.compile(r'>')

def parse_fa(input):
    seq_l =  []
    for line in input:
        if line == "\n":
            continue
        else:
            line = line.strip()
            seq_l.append(line)
    return seq_l # 最後のsequenceのデータを出力

def align(X,Y):
    X_seq = str("First seq: "+ X +"\n")
    Y_seq = str("Second seq: "+ Y +"\n\n")

    alignment = pairwise2.align.globalms(X, Y, 1, -1, -2, -1)

    with open("output.txt", "w") as out:
        out.write(X_seq)
        out.write(Y_seq)
        for al in alignment:
            #print(format_alignment(*al))
            out.write(format_alignment(*al))


if __name__ == "__main__":
    # 初期化
    header = []
    seq = []
    seed(1)

    # FASTAの標準入力
    inFile = sys.argv[1]  
    
    with open(inFile, "r") as f:
        for sequence in parse_fa(f):
            seq.append(sequence) # seqのリスト
    
    #align(seq[randint(0, len(seq))], seq[randint(0, len(seq))])
    align(seq[1], seq[2])