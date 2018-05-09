#!/bin/bash
samtools view B22-730.sam -bh >B22-730.bam
samtools sort B22-730.bam -o B22-730.sorted.bam
samtools depth -a B22-730.sorted.bam >B22-730.cov
./cov_stats.py B22-730.cov 1000 B22-730.cov.png

samtools view B22-730.bt2.sam -bh >B22-730.bt2.bam
samtools sort B22-730.bt2.bam -o B22-730.bt2.sorted.bam
samtools depth -a B22-730.bt2.sorted.bam >B22-730.bt2.cov
./cov_stats.py B22-730.bt2.cov 1000 B22-730.bt2.cov.png
