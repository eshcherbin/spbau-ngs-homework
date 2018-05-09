#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from Bio import SeqIO
from collections import Counter, defaultdict


if __name__ == '__main__':
    assert len(sys.argv) == 4
    ref_seq = next(SeqIO.parse(sys.argv[2], 'fasta')).seq
    f = pysam.AlignmentFile(sys.argv[1], 'r')
    cnt = defaultdict(Counter)

    hp = dict()
    i = 0
    print(len(ref_seq))
    while i < len(ref_seq):
        j = i + 1
        while j < len(ref_seq) and ref_seq[i] == ref_seq[j]:
            j += 1
        if j - i >= 3:
            hp[i] = j - 1
        i = j

    print('Found all homopolymers')

    cr = 0
    for r in f:
        cr += 1
        if r.is_unmapped:
            continue
        ap = r.get_aligned_pairs()
        i = 0
        while i < len(ap):
            if ap[i][1] not in hp:
                i += 1
                continue
            j = i
            while j < len(ap) and ap[j][1] != hp[ap[i][1]]:
                j += 1
            if j == len(ap):
                break
            l = ap[j][1] - ap[i][1] + 1
            cur = 0
            while i <= j:
                if ap[i][0] is not None:
                    cur += 1
                i += 1
            cnt[l][cur] += 1

    for l in cnt:
        with open(sys.argv[3] + '.' + str(l) + '.tsv', 'w') as fout:
            for k in range(max(cnt[l]) + 1):
                fout.write(str(k) + '\t' + str(cnt[l][k]) + '\n')
