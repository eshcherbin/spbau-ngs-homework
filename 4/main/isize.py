#!/usr/bin/env python3
import pysam
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys


if __name__ == '__main__':
    assert len(sys.argv) == 3
    f = pysam.AlignmentFile(sys.argv[1], 'r')
    cnt, cnt_good = 0, 0
    isizes = []
    for r in f:
        cnt += 1
        if r.is_unmapped or not r.is_proper_pair:
            continue
        if r.reference_start >= r.next_reference_start and r.reference_end - r.next_reference_start <= 3000:
            cnt_good += 1
            isizes.append(r.reference_end - r.next_reference_start)
    print('Finished calculating insertion sizes')
    print(cnt_good, 'good read pairs out of', cnt)

    isizes = np.array(isizes)
    bincnt = np.bincount(isizes)
    with open(sys.argv[2][:-4], 'w') as fout:
        fout.write('\n'.join(map(str, bincnt)) + '\n')
    print('Mean insert size: {:.3f}'.format(isizes.mean()))
    print('Standard deviation: {:.3f}'.format(isizes.std()))
    target = 0.95 * sum(bincnt)
    l, r = len(bincnt), len(bincnt)
    cur = 0
    while cur < target:
        l -= 1
        cur += bincnt[l]
    conf_l, conf_r = l, r
    while l > 0:
        r -= 1
        cur -= bincnt[r]
        while cur < target and l > 0:
            l -= 1
            cur += bincnt[l]
        if cur >= target and r - l < conf_r - conf_l:
            conf_l, conf_r = l, r
    print('95% confidence interval: [{}, {}]'.format(conf_l, conf_r - 1))

    # threshold = 1000
    threshold = 100
    min_isize = 0
    max_isize = len(bincnt)
    while bincnt[min_isize] < threshold:
        min_isize += 1
    while bincnt[max_isize - 1] < threshold:
        max_isize -= 1
    plt.scatter(range(min_isize, conf_l), bincnt[min_isize:conf_l], c='r')
    plt.scatter(range(conf_l, conf_r), bincnt[conf_l:conf_r], c='b')
    plt.scatter(range(conf_r, max_isize), bincnt[conf_r:max_isize], c='r')
    # plt.axvline(conf_l, c='k', linestyle='dashed')
    # plt.axvline(conf_r - 1, c='k', linestyle='dashed')
    plt.xlabel('Insert size')
    plt.ylabel('Number of read pairs')
    plt.savefig(sys.argv[2])
