#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys


if __name__ == '__main__':
    assert len(sys.argv) == 4

    cov_df = pd.read_table(sys.argv[1], index_col=False, header=None)
    cov = cov_df.iloc[:, 2]

    bin_len = int(sys.argv[2])
    x = list(range(0, len(cov), bin_len))
    plt.plot(x, [cov[i:i + bin_len].mean() for i in x]) #, bin_len)
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.gca().set_xlim(left=0)
    plt.gca().set_ylim(bottom=0)
    plt.savefig(sys.argv[3])

    print('Mean coverage: {:.3f}'.format(cov.mean()))
    print('{:.2f}% of genome is covered'.format((cov > 0).sum() / len(cov) * 100))
