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
    parser.add_argument('--mindp',type=int, default=0)
    parser.add_argument('--minqd', type=int, default=1)
    args = parser.parse_args()
    
    N = 11 
    mindp = int(args.mindp)
    minqd = int(args.minqd) 
    outname = args.outfile 
    vcf_file =open(outname, 'w+')
    x =0
    for line in open(args.infile):
            if line.startswith('#'):
                print(line.strip(), file=vcf_file)
                continue 
            rows = (line.strip()).split('\t')
            print(rows)
            x +=1
            qualities = rows[7].split(';')
            for i in qualities:
                if "QD" in i:
                    qd = str(i) #str(i).sptit('=')[1] 
                if "MQ=" in i: 
                    mq = str(i) 
                if "FS" in i: 
                    fs = str(i) 
                if "MQRankSum" in i: 
                    mqranksum = str(i) 
                if "ReadPosRankSum" in i: 
                    readposranksum = str(i) 
                    
            QD = qd.split('=')[1]
            MQ = mq.split('=')[1] 
            FS = fs.split('=')[1]
            MQRankSum = mqranksum.split('=')[1]
            ReadPosRankSum = readposranksum.split('=')[1] 
            flag = 0 
            count = 0
            for i in range(9,20) :
                if len(rows[i]) >1:
                   geno =rows[i].split(":")
                   if (geno[0] !='./.'):
                        dp = rows[i].split(":")[2]
                        if int(dp) > int(mindp) :
                   #    if   float(QD) > int(minqd) :
                   #        if float(MQ) > 30 and float(FS) < 60 and float(MQRankSum)  > -12.5 and float(ReadPosRankSum) > -8.0:
                                  flag +=1
            if flag > 10:
                    print(line.strip(), file= vcf_file)


if __name__ == '__main__':
    main()
