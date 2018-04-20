#!/usr/bin/env python3
import pysam
import sys


def is_mismatch(c):
    return c is None or str.islower(c)


def calc(ori_fn, cor_fn):
    ori = pysam.AlignmentFile(ori_fn, 'r')
    cor = pysam.AlignmentFile(cor_fn, 'r')
    fp, fn, tn = 0, 0, 0
    for r_cor in cor:
        r_ori = next(ori)
        while r_ori.qname != r_cor.qname:
            r_ori = next(ori)
        if r_ori.is_unmapped or r_cor.is_unmapped:
            continue
        m_ori = list(zip(*r_ori.get_aligned_pairs(with_seq=True)))[2]
        m_cor = list(zip(*r_cor.get_aligned_pairs(with_seq=True)))[2]
        if r_cor.is_reverse:
            m_ori = list(reversed(m_ori))
            m_cor = list(reversed(m_cor))
        for c_ori, c_cor in zip(m_ori, m_cor):
            if is_mismatch(c_ori) and is_mismatch(c_cor):
                fp += 1
            elif not is_mismatch(c_ori) and is_mismatch(c_cor):
                fn += 1
            elif is_mismatch(c_ori) and not is_mismatch(c_cor):
                tn += 1
    return fp, fn, tn


if __name__ == '__main__':
    ori1_fn, ori2_fn, cor1_fn, cor2_fn = sys.argv[1:]
    fp1, fn1, tn1 = calc(ori1_fn, cor1_fn)
    fp2, fn2, tn2 = calc(ori2_fn, cor2_fn)
    fp, fn, tn = fp1 + fp2, fn1 + fn2, tn1 + tn2
    print('False positives:', fp)
    print('False negatives:', fn)
    print('True negatives:', tn)
