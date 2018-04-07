#!/bin/bash
make clean && make
mkdir -p ../results
#for f in 'ECOLI_IS220_QUAKE_1K_paired_reads.fasta' 'ECOLI_IS220_QUAKE_1K_single_reads.fasta' 's_6.first100000.fastq' 's_6.first10000.fastq'  's_6.first1000.fastq'
for f in $(ls ../data)
do
    echo "$f"
    fn=${f%.*}
    echo "$fn"
    time ./debrujin "../data/$f" 55 && make graph.svg
    mv edges.fasta "../results/$fn.edges.fasta"
    mv graph.svg "../results/$fn.svg"
done
