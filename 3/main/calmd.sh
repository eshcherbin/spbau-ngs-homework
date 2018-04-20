#!/bin/bash
samtools calmd ori1.bam MG1655-K12.first400K.fasta | samtools view -h > ori_md1.sam
samtools calmd ori2.bam MG1655-K12.first400K.fasta | samtools view -h > ori_md2.sam
samtools calmd cor1.bam MG1655-K12.first400K.fasta | samtools view -h > cor_md1.sam
samtools calmd cor2.bam MG1655-K12.first400K.fasta | samtools view -h > cor_md2.sam
