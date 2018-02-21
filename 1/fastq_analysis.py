#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
from Bio import SeqIO


DEBUG_STEPS = 100000
GC_PLOT_STEP = 0.01


class FastqAnalyzer:
    def __init__(self, max_read_length):
        self.gc_content = \
            np.zeros(round(1. / GC_PLOT_STEP) + 1, dtype=np.int64)
        self.error_cnt = np.zeros(max_read_length, dtype=np.int64)
        self.error_sum = np.zeros(max_read_length, dtype=np.float64)
        # self.error_max = np.zeros(max_read_length, dtype=np.float64)
        # self.error_min = np.zeros(max_read_length, dtype=np.float64)

    def process(self, filename):
        parser = SeqIO.parse(filename, 'fastq-illumina')
        start_time = timer()
        for i, record in enumerate(parser):
            if i and i % DEBUG_STEPS == 0:
                print('{} steps passed, {:.3f} seconds elapsed'
                      .format(i, timer() - start_time))
            gc = len([c for c in record.seq if c in 'GCgc']) / len(record.seq)
            self.gc_content[int(gc / GC_PLOT_STEP)] += 1
            # print(record.letter_annotations['phred_quality'])
            # print(list(map(lambda qual: 1. / (10.**(qual / 10.) + 1), record.letter_annotations['phred_quality'])))
            for pos, qual in enumerate(
                    record.letter_annotations['phred_quality']):
                p = 1. / (10.**(qual / 10.) + 1)
                print(qual, p)
                self.error_cnt[pos] += 1
                self.error_sum[pos] += qual

    def plot_gc_content(self, image_path):
        xs = np.linspace(0, 1, len(self.gc_content))
        plt.scatter(xs, self.gc_content)
        plt.savefig(image_path)

    def plot_error_distribution(self, image_path):
        xs = np.arange(0, len(self.error_cnt))
        ys = np.where(self.error_cnt == 0,
                      self.error_sum,
                      self.error_sum / self.error_cnt)
        print(self.error_cnt)
        print(self.error_sum)
        plt.plot(xs, ys, linestyle='solid')
        plt.savefig(image_path)


if __name__ == '__main__':
    assert len(sys.argv) == 3
    fastq_filename = sys.argv[1]
    max_read_length = int(sys.argv[2])
    analyzer = FastqAnalyzer(max_read_length)
    analyzer.process(fastq_filename)
    analyzer.plot_gc_content(fastq_filename + '.gc_content.png')
    analyzer.plot_error_distribution(fastq_filename + '.err_distribution.png')
