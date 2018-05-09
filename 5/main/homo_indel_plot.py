#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys
import os


if __name__ == '__main__':
    assert len(sys.argv) == 2
    p = sys.argv[1]
    for f in os.listdir('.'):
        if not (f.startswith(p) and f[-4:] == '.tsv'):
            continue
        df = pd.read_table(f, header=None)
        plt.bar(df.iloc[:, 0], np.where(df.iloc[:, 1] > 0, np.log10(df.iloc[:, 1] + 0.5), 0), 0.95)
        plt.savefig(f[:-4] + '.png')
        plt.clf()
