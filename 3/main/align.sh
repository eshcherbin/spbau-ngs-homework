#!/bin/bash
minimap2 -ax sr -a --cs MG1655-K12.first400K.fasta ecoli_400K_err_1.fastq | samtools view -bh > ori1.bam
minimap2 -ax sr -a --cs MG1655-K12.first400K.fasta ecoli_400K_err_2.fastq | samtools view -bh > ori2.bam
minimap2 -ax sr -a --cs MG1655-K12.first400K.fasta ecoli_400K_err_1.cor.fastq | samtools view -bh > cor1.bam
minimap2 -ax sr -a --cs MG1655-K12.first400K.fasta ecoli_400K_err_2.cor.fastq | samtools view -bh > cor2.bam
