#!/usr/bin/env python3
import pysam
import sys


def is_mismatch(c):
    return c is None or str.islower(c)


def calc(ori_fn, cor_fn):
    ori = pysam.AlignmentFile(ori_fn, 'r')
    cor = pysam.AlignmentFile(cor_fn, 'r')
    tlen, cnt, fp, fn, tn = 0, 0, 0, 0, 0
    for r_cor in cor:
        r_ori = next(ori)
        while r_ori.query_name != r_cor.query_name:
            r_ori = next(ori)
        if r_ori.is_unmapped or r_cor.is_unmapped:
            continue
        tlen += r_cor.query_length
        m_ori = r_ori.get_aligned_pairs(with_seq=True)
        m_cor = r_cor.get_aligned_pairs(with_seq=True)
        if r_cor.is_reverse:
            m_ori = [(r_ori.query_length - 1 - read_pos if read_pos is not None else None, ref_pos, c) 
                     for read_pos, ref_pos, c
                     in reversed(m_ori)]
            m_cor = [(r_cor.query_length - 1 - read_pos if read_pos is not None else None, ref_pos, c)
                     for read_pos, ref_pos, c
                     in reversed(m_cor)]
        i_ori = 0
        for read_pos_cor, ref_pos_cor, c_cor in m_cor:
            if read_pos_cor is None:
                continue
            while m_ori[i_ori][0] != read_pos_cor:
                i_ori += 1
                if i_ori == len(m_ori):  # debug output
                    print(r_ori.query_name)
                    print(m_ori)
                    print(m_cor)
                    print(read_pos_cor)
                    sys.exit(1)
            read_pos_ori, ref_pos_ori, c_ori = m_ori[i_ori]
            if read_pos_ori != read_pos_cor or ref_pos_ori != ref_pos_cor:
                continue
            cnt += 1
            if is_mismatch(c_ori) and is_mismatch(c_cor):
                fp += 1
            elif not is_mismatch(c_ori) and is_mismatch(c_cor):
                fn += 1
            elif is_mismatch(c_ori) and not is_mismatch(c_cor):
                tn += 1
    return tlen, cnt, fp, fn, tn


if __name__ == '__main__':
    ori1_fn, ori2_fn, cor1_fn, cor2_fn = sys.argv[1:]
    tlen1, cnt1, fp1, fn1, tn1 = calc(ori1_fn, cor1_fn)
    tlen2, cnt2, fp2, fn2, tn2 = calc(ori2_fn, cor2_fn)
    tlen, cnt, fp, fn, tn = tlen1 + tlen2, cnt1 + cnt2, fp1 + fp2, fn1 + fn2, tn1 + tn2
    print('Total mapped corrected reads\' length:', tlen)
    print(cnt, 'matching positions found, of which:')
    print('False positives:', fp)
    print('False negatives:', fn)
    print('True negatives:', tn)
