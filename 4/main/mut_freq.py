#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
import sys
from Bio import SeqIO
from collections import Counter


BASES = 'ACGT'


if __name__ == '__main__':
    assert len(sys.argv) == 3
    f = pysam.AlignmentFile(sys.argv[1], 'r')
    ref_seq = next(SeqIO.parse(sys.argv[2], 'fasta')).seq
    cnt = Counter()
    for r in f:
        if r.is_unmapped:
            continue
        # if not r.is_proper_pair:
            # continue
        i, j = 0, r.reference_start
        for op, l in r.cigartuples:
            if op == 0:
                for _ in range(l):
                    cnt[(r.seq[i], ref_seq[j])] += 1
                    i += 1
                    j += 1
            elif op == 1:
                i += l
            elif op == 2:
                j += l

    cnt = np.array([[cnt[(b_from, b_to)] for b_from in BASES] for b_to in BASES], dtype=np.int32)
    cnt = pd.DataFrame(data=cnt, index=list(BASES), columns=list(BASES))
    print(cnt)
