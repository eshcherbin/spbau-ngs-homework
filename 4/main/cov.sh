#!/bin/bash
samtools view HTC80499.sam -bh >HTC80499.bam
samtools sort HTC80499.bam -o HTC80499.sorted.bam
samtools depth -a HTC80499.sorted.bam >HTC80499.cov
./cov_stats.py HTC10499.cov 1000 HTC10499.cov.png

samtools view HTC80508.sam -bh >HTC80508.bam
samtools sort HTC80508.bam -o HTC80508.sorted.bam
samtools depth -a HTC80508.sorted.bam >HTC80508.cov
./cov_stats.py HTC10508.cov 1000 HTC10508.cov.png
