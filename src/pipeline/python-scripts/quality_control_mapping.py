#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 12:03:17 2016

@author: eva

tests the quality of an alignment:
    average coverage > 10 000
    percentage mapped > 99%
    average mapping quality > 55
prints True if this is all ok

usage: python alignment_quality.py infile_sorted.bam

"""
import pysam
import os
import numpy as np
import sys
import matplotlib.pyplot as plt

def good_alignment(samfile):
    coverage = np.array(samfile.count_coverage('NL43_ann_wk0virusPassRef_plasmid',start=1,end=10000))
    coverage_full = np.sum(coverage,0)
    plt.plot(coverage_full)


    average = np.mean(coverage_full)
    average_5 = np.mean(coverage_full[450:790])
    average_a = np.mean(coverage_full[450:2200])
    average_b = np.mean(coverage_full[3000:4000])
    average_c = np.mean(coverage_full[4800:5900])
    average_d = np.mean(coverage_full[6600:7400])
    average_e = np.mean(coverage_full[8400:9150])
    average_3 = np.mean(coverage_full[9500:9700])
    cover = average > 5000
    cover_all = all([i>= 1000 for i in [average_5,average_a,average_b,average_c,average_d,average_e,average_3]])
    
    mapped = 0 # for bam you can do this: samfile.mapped/float((samfile.mapped+samfile.unmapped))
    mapping = mapped > 0.99

    mapquality = 0
    for read in samfile.fetch():
        mapquality += read.mapping_quality
    mapquality = 0 # for bam you can do this mapquality/float(samfile.mapped)
    quality = mapquality > 55

    plt.title('{:.0f},{:.3f},{:.2f},{}'.format(average,mapped, mapquality, cover and mapping and quality))
    plt.savefig(snakemake.output[0])
    return average,cover,average_5,average_a,average_b,average_c,average_d,average_e,average_3,cover_all,mapped, mapquality and mapping and quality

try:
    samfile = pysam.AlignmentFile(snakemake.input[0])
except ValueError:
    samfile = None

if samfile is not None:
    result = good_alignment(samfile)
    samfile.close()
    print(result)
    with open(snakemake.output[1],'w') as f:
        f.write(snakemake.input[0]+'\t')
        f.write('\t'.join([str(i) for i in result]))
        f.write('\n')
else:
    with open(snakemake.output[1],'w') as f:
        f.write(snakemake.input[0]+'\t0\tFalse\t0\t0\t0\t0\t0\t0\t0\tFalse\t0\t0\n')
        plt.plot([],[])
        plt.savefig(snakemake.output[0])


#coverage,coverage_full = np.array(coverage_plot(samfile))
#coverage = samfile.count_coverage(start=1,end=10000)
