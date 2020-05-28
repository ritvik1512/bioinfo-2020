from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import re

regex = re.compile(r'>')

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


if __name__ == "__main__":
    # 初期化
    header = []
    seq = []

    # FASTAの標準入力
    inFile = sys.argv[1]  
    with open(inFile, "r") as f:
        for head,sequence in parse_fa(f):
            header.append(head) # headerのリスト
            seq.append(sequence) # seqのリスト
