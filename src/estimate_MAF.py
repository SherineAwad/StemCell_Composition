#! /usr/bin/env python
import sys
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg') 
from pylab import *
import operator 

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',dest='infile')
    parser.add_argument('-v', '--vcf', dest='vcf')
    parser.add_argument('-n', type= int, dest='N')
    parser.add_argument('-o','--outfile',dest='outfile')
    args = parser.parse_args()
    N = args.N
    afile = open(args.outfile +"_alleles.txt", 'w+')
    gfile = open(args.outfile +"_genotypes.txt", "w+") 
    genotypes = {}
    print("UM14-1_P26_150115,UM4-6_P15_150116,WA09_P25_140528,WA14_P21_131206,WA17_P12_140530,WA18_P9_140530,WA19_P11_140529,WA26_P15_140606,WA27_P15_140613,WIBR2_P33_140617,WIBR6_P11_140617,rn", file=gfile)
    for line in open(args.vcf, 'r'):
        if line.startswith('#'): 
           continue
        if "./." in line: 
           continue 
        row = line.split('\t')
        index = row[0] +"_"+row[1] + "_"+row[4] 
        record=""
        for i in range(9,N+9) :
                if len(row[i]) >1:
                   geno =row[i].split(":")
                   nflag = 0 
                   if geno[0] == "0/0":
                             record+="0,"
                             nflag +=0
                   elif geno[0] =="1/1":
                             record+= "1,"
                             nflag +=1
                   elif geno[0] == "0/1":
                             record +="0.5,"
                             nflag +=0.5
        if nflag == 0:
                        continue
        if nflag == (1*N):
                       continue
        if nflag == (0.5*N):
                        continue
        else:
             genotypes[index] = record
    print("total_reads,A,C,T,G,b_estimate,ALT", file=afile)
    for line in open(args.infile, 'r'):
           row = line.split('\t')
           chrom = row[0] 
           pos = row[1] 
           ALT = row[2] 
           N_reads = int(row[3]) 
           index = chrom +"_"+pos +"_"+ALT
           if index in genotypes:
              if N_reads == 0: 
                  continue
              a_info = row[4].split(":") 
              c_info = row[5].split(":") 
              g_info = row[6].split(":")
              t_info = row[7].split(":") 
              a_reads = int(a_info[1])
              c_reads = int(c_info[1])
              g_reads = int(g_info[1])
              t_reads = int(t_info[1])
              if ALT =="A": 
                 b_estimate = float(a_reads/N_reads)
              elif ALT =="C": 
                 b_estimate = float(c_reads/N_reads)
              elif ALT =="G": 
                 b_estimate = float(g_reads/N_reads)
              elif ALT =="T":
                  b_estimate = float(t_reads/N_reads)
              alleles_record = str(N_reads)+','+str(a_reads)+','+str(c_reads)+','+str(t_reads)+','+str(g_reads)+','+str(b_estimate)+','+ALT
              print(alleles_record, file=afile)
              print(genotypes[index]+index, file =gfile)
           else: 
              continue 
if __name__ == '__main__':
    main()
