#!/usr/bin/env python3
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
from Bio import SeqIO
plt.switch_backend('agg')


DEBUG_STEPS = 100000
MAX_STEPS = 1000000000  # used to terminate early for debug purposes
GC_PLOT_STEP = 0.01

BP_QUALITY_THRESHOLD = 10
READ_QUALITY_THRESHOLD = 0.5


class FastqAnalyzer:
    def __init__(self, max_read_length):
        # counter of reads with GC content within
        # each interval of size GC_PLOT_STEP
        self.gc_content = \
            np.zeros(round(1. / GC_PLOT_STEP) + 1, dtype=np.int64)
        # sum of all GC content of all reads
        self.gc_sum = 0
        # number of reads whose GC content is accounted (good quality reads)
        self.gc_n_reads = 0
        # counter of all reads that include the given position
        self.error_cnt = np.zeros(max_read_length, dtype=np.int64)
        # sum of all error probabilities for the given position
        self.error_sum = np.zeros(max_read_length, dtype=np.float64)
        # sum of all quality values for the given position
        self.qual_sum = np.zeros(max_read_length, dtype=np.float64)

    def process(self, filename):
        parser = SeqIO.parse(filename, 'fastq-illumina')
        start_time = timer()
        for i, record in enumerate(parser):
            # debug information and termination
            if i and i % DEBUG_STEPS == 0:
                print('{} steps passed, {:.3f} seconds elapsed'
                      .format(i, timer() - start_time))
            if i == MAX_STEPS:
                break

            # number of positions with bad quality
            n_bad_qual_bps = \
                len([q for q in record.letter_annotations['phred_quality']
                     if q <= BP_QUALITY_THRESHOLD])
            # check overall quality of the read
            if n_bad_qual_bps < READ_QUALITY_THRESHOLD * len(record.seq):
                # filter out positions with bad quality
                good_qual_seq = \
                    [c for c, q
                     in zip(record.seq,
                            record.letter_annotations['phred_quality'])
                     if q > BP_QUALITY_THRESHOLD]
                gc = len([c for c in good_qual_seq
                          if c in 'GCgc']) / len(record.seq)
                self.gc_content[int(gc / GC_PLOT_STEP)] += 1
                self.gc_sum += gc
                self.gc_n_reads += 1

            for j, q in enumerate(record.letter_annotations['phred_quality']):
                p = 1. / (10.**(q / 10.) + 1)  # calculate error probability
                self.error_cnt[j] += 1
                self.error_sum[j] += p
                self.qual_sum[j] += q

    def plot_gc_content(self, image_path):
        print('Mean GC content: {:.2f}%'
              .format(100. * self.gc_sum / self.gc_n_reads))
        xs = np.linspace(0, 100, len(self.gc_content))
        plt.xlabel('GC content, %')
        plt.ylabel('Number of reads')
        plt.scatter(xs, self.gc_content, marker='+', color='red')
        plt.savefig(image_path)
        plt.clf()
        plt.cla()
        plt.close()

    def plot_error_distribution(self, image_path):
        xs = np.arange(0, len(self.error_cnt))

        plt.subplot(211)
        ys = np.where(self.error_cnt == 0,
                      self.error_sum,
                      self.error_sum / self.error_cnt)
        plt.ylabel('Average error probability')
        plt.plot(xs, ys)

        plt.subplot(212)
        ys = np.where(self.error_cnt == 0,
                      self.qual_sum,
                      self.qual_sum / self.error_cnt)
        plt.xlabel('Position')
        plt.ylabel('Average quality')
        plt.plot(xs, ys)

        plt.savefig(image_path)
        plt.clf()
        plt.cla()
        plt.close()


if __name__ == '__main__':
    assert len(sys.argv) == 3
    fastq_filename = sys.argv[1]
    max_read_length = int(sys.argv[2])
    analyzer = FastqAnalyzer(max_read_length)
    analyzer.process(fastq_filename)
    analyzer.plot_gc_content(os.getcwd() + '/gc_content.png')
    analyzer.plot_error_distribution(os.getcwd() + '/err_distribution.png')
