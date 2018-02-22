#!/usr/bin/env python3
import sys
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from timeit import default_timer as timer
from Bio import SeqIO


DEBUG_STEPS = 100000
MAX_STEPS = 1000000000
GC_PLOT_STEP = 0.01

BP_QUALITY_THRESHOLD = 20
READ_QUALITY_THRESHOLD = 0.8


class FastqAnalyzer:
    def __init__(self, max_read_length):
        self.gc_content = \
            np.zeros(round(1. / GC_PLOT_STEP) + 1, dtype=np.int64)
        self.gc_sum = 0
        self.gc_n_reads = 0
        self.error_cnt = np.zeros(max_read_length, dtype=np.int64)
        self.error_sum = np.zeros(max_read_length, dtype=np.float64)
        self.qual_sum = np.zeros(max_read_length, dtype=np.float64)
        # self.error_max = np.zeros(max_read_length, dtype=np.float64)
        # self.error_min = np.zeros(max_read_length, dtype=np.float64)

    def process(self, filename):
        parser = SeqIO.parse(filename, 'fastq-illumina')
        start_time = timer()
        for i, record in enumerate(parser):
            if i and i % DEBUG_STEPS == 0:
                print('{} steps passed, {:.3f} seconds elapsed'
                      .format(i, timer() - start_time))
            if i == MAX_STEPS:
                break
            
            n_bad_qual_bps = \
                len([q for q in record.letter_annotations['phred_quality'] if q <= BP_QUALITY_THRESHOLD]) 
            if n_bad_qual_bps < READ_QUALITY_THRESHOLD * len(record.seq):  # else too many bad quality bps
                good_qual_seq = [c for c, q in zip(record.seq, record.letter_annotations['phred_quality']) if q > BP_QUALITY_THRESHOLD]
                gc = len([c for c in good_qual_seq if c in 'GCgc']) / len(record.seq)
                self.gc_content[int(gc / GC_PLOT_STEP)] += 1
                self.gc_sum += gc
                self.gc_n_reads += 1

            # print(record.letter_annotations['phred_quality'])
            # print(list(map(lambda qual: 1. / (10.**(qual / 10.) + 1), record.letter_annotations['phred_quality'])))
            for pos, qual in enumerate(record.letter_annotations['phred_quality']):
                p = 1. / (10.**(qual / 10.) + 1)
                # print(qual, p)
                self.error_cnt[pos] += 1
                self.error_sum[pos] += p
                self.qual_sum[pos] += qual

    def plot_gc_content(self, image_path):
        print('Mean GC content: {:.2f}%'.format(100. * self.gc_sum / self.gc_n_reads))
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
        # print(self.error_cnt)
        # print(self.error_sum)
        # print(xs)
        # print(ys)
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
    analyzer.plot_gc_content(os.getcwd() + '/gc_content_v7.png')
    analyzer.plot_error_distribution(os.getcwd() + '/err_distribution_v7.png')
