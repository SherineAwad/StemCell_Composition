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

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',dest='infile')
    parser.add_argument('-o','--outfile',dest='outfile')
    parser.add_argument('-n', dest='N')
    args = parser.parse_args()
    gfile = open(args.outfile +"_filtered.vcf", 'w+')
    genos =""
    N = int(args.N) 
    outname = args.outfile 
    #print("UM14-1_P26_150115,UM4-6_P15_150116,WA09_P25_140528,WA14_P21_131206,WA17_P12_140530,WA18_P9_140530,WA19_P11_140529,WA26_P15_140606,WA27_P15_140613,WIBR2_P33_140617,WIBR6_P11_140617,rn", file=gfile)
    region_list_file =open(outname+"_regions.txt", 'w+')
    for line in open(args.infile):
            genos =""
            if line.startswith('#'):
                print(line.strip(), file =gfile)
                continue 
            rows = (line.strip()).split('\t')
            alt = rows[4]
            if len(alt) > 1: 
                 continue
            flag = 0 
            count = 0
            for i in range(9,N+9) :
                if len(rows[i]) >1:
                   geno =rows[i].split(":")
                   if (geno[0] !='./.'):
                          flag +=1
                          genos+=geno[0]
                          genos+=","
            if flag > (N-1):
                    snp_type = rows[4] 
                    if "," in snp_type:
                        continue 
                    start = rows[1] 
                    end = (int(start) + len(rows[4])) -1  
                    genotypes = genos.split(",")
                    record = ""
                    nflag = 0
                    for i in range(0, N): 
                        if genotypes[i] == "0/0": 
                             record +="0,"
                             nflag +=0 
                        elif genotypes[i] =="1/1":
                             record+= "1,"
                             nflag +=1 
                        elif genotypes[i] == "0/1":
                             record +="0.5,"
                             nflag +=0.5 
                    if nflag == 0: 
                        continue 
                    if nflag == (1*N):
                        continue 
                    if nflag == (0.5*N):
                        continue
                    #record+= rows[0]
                    #record+= ","+str(start)+"_"+str(end)+"_"+alt
                    print(line.strip(), file=gfile)
                    print(rows[0], start, end, file= region_list_file)
            else: 
                 genos =""
if __name__ == '__main__':
    main()
