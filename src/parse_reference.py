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
    vid = [] 
    count =0 
    outname = args.outfile 
    outfile = open(outname+"_genotype.txt", 'w+')
    for line in open(args.infile):
            if line.startswith('#'):
                continue
            if "./." in line: 
                continue 
            rows = (line.strip()).split('\t')
            count += 1 
            sid = rows[0] +'.'+rows[1]
            GT = ""
            code =0 
            for i in range(9,9+N):
              if len(rows[i]) > 0:
                gt = rows[i].split(":")
                if gt[0] =="./.":
                    GT = "X"
                    break 
                elif gt[0] == "0/0":
                        GT += "0\t"  
                        code += 0 
                elif gt[0] =="1/1":
                        GT += "2\t"
                        code +=2 
                elif gt[0] == "0/1":
                        GT +="1\t"
                else:
                    GT ="X"
                    break

            if GT != "X" and code !=0 and code != 22:
                gtsnps[sid]= GT 
            
    for i in gtsnps: 
        print(i,'_', gtsnps[i], file =outfile) 

if __name__ == '__main__':
    main()
