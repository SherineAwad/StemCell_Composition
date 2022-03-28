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
    parser.add_argument('infile')
    parser.add_argument('-n',dest='n', default =1)
    parser.add_argument('-o',dest='outfile')
    args = parser.parse_args()
    N = int(args.n) 
    gtsnps = {}
    MAFs = {}
    vid = [] 
    count =0 
    outname = args.outfile 
    outfile = open(outname+"_genotype.txt", 'w+')
    outfile2 = open(outname+"_maf.txt", "w+")
    print("D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,rn", file =outfile)
    for line in open(args.infile):
            if line.startswith('#'):
                continue
            if "./." in line: 
                continue 
            rows = (line.strip()).split('\t')
            qualities = rows[7].split(';')
            af = qualities[1]
            AF = af.split('=')[1]
            if (',') in AF:
                 continue
            MAF = 0
            if (float(AF) < 0.5):
                 MAF = float(AF)
            else:
                 MAF = 1 - float(AF)
            count += 1 
            sid = rows[0] +'_'+rows[1]+'_'+rows[3]+'_'+rows[4]
            GT = ""
            code =0 
            for i in range(9,9+N):
              if len(rows[i]) > 0:
                gt = rows[i].split(":")
                if gt[0] =="./.":
                    GT = "X"
                    break 
                elif gt[0] == "0/0":
                        GT += "0,"  
                        code += 0 
                elif gt[0] =="1/1":
                        GT += "1,"
                        code +=2 
                elif gt[0] == "0/1":
                        GT +="0.5,"
                else:
                    GT ="X"
                    break

            if GT != "X" and code !=0 and code != 22:
                gtsnps[sid]= GT
                MAFs[sid] = MAF      
    for i in gtsnps: 
        print(gtsnps[i],i, file =outfile) 
        print(i,",",MAFs[i], file=outfile2)    
if __name__ == '__main__':
    main()
