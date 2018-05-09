#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from Bio import SeqIO
from collections import Counter


BASES = 'ACGT_'


if __name__ == '__main__':
    assert len(sys.argv) == 5
    f = pysam.AlignmentFile(sys.argv[1], 'r')
    ref_seq = next(SeqIO.parse(sys.argv[2], 'fasta')).seq
    cnt = Counter()
    qual_cnt = (Counter(), Counter())
    indel_cnt = (Counter(), Counter())
    for r in f:
        if r.is_unmapped:
            continue
        # if not r.is_proper_pair:
            # continue
        i, j = 0, r.reference_start
        for op, l in r.cigartuples:
            if op == 0:
                for _ in range(l):
                    cnt[r.seq[i], ref_seq[j]] += 1
                    if r.seq[i] != ref_seq[j]:
                        qual_cnt[0][r.query_qualities[i]] += 1
                    i += 1
                    j += 1
            elif op == 1:
                indel_cnt[0][l] += 1
                for _ in range(l):
                    cnt[r.seq[i], '_'] += 1
                    qual_cnt[1][r.query_qualities[i]] += 1
                    i += 1
            elif op == 2:
                indel_cnt[1][l] += 1
                for _ in range(l):
                    cnt['_', ref_seq[j]] += 1
                    j += 1
            elif op == 4:
                i += l

    cnt = np.array([[cnt[(b_to, b_from)] for b_from in BASES] for b_to in BASES], dtype=np.int32)
    cnt = pd.DataFrame(data=cnt, index=list(BASES), columns=list(BASES))
    print(cnt)

    for i in (0, 1):
        for x in indel_cnt[i]:
            if indel_cnt[i][x] > 0:
                indel_cnt[i][x] = np.log10(indel_cnt[i][x])
    plt.bar(indel_cnt[0].keys(), indel_cnt[0].values(), color='b', label='Insertions')
    plt.bar([-x for x in indel_cnt[1].keys()], indel_cnt[1].values(), color='r', label='Deletions')
    plt.ylabel('log10 of cnt')
    plt.legend()
    plt.savefig(sys.argv[3])
    plt.clf()

    plt.scatter(qual_cnt[0].keys(), qual_cnt[0].values(), c='b', label='Mismatches')
    plt.scatter(qual_cnt[1].keys(), qual_cnt[1].values(), c='r', label='Insertions')
    plt.legend()
    plt.savefig(sys.argv[4])
    plt.clf()
