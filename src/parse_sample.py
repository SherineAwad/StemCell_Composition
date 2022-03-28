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
    parser.add_argument('infile1')
    parser.add_argument('infile2')
    parser.add_argument('-n',dest='n', default =1)
    parser.add_argument('-o',dest='outfile')
    parser.add_argument('-d', type =int, dest='mindp')  
    parser.add_argument('-q', type =int, dest ='minqd')
    args = parser.parse_args()
    N = args.n
    mindp = int(args.mindp)
    minqd = float(args.minqd)
    gtsnps = {}
    alsnps = {} 
    vid = [] 
    count =0 
    outname = args.outfile 
    gtfile =open(outname+"_genotype.txt", 'w+')
    print('D1\tD2\tD3\tD4\tD5\tD6\tD7\tD8\tD9\tD10\tD11',file=gtfile)
    acfile = open(outname+"_allelespecific.txt", "w+")
    print("Ref",'\t',"Alt", file =acfile)
    for line in open(args.infile1):
            record = line.split('_')
            gtsnps[record[0].strip()] = record[1].strip() 
    qd = "=0" 
    mq = "=0" 
    fs = "=0" 
    mqranksum = "=0" 
    readposranksum = "=0" 
    for line in open(args.infile2): 
            if line.startswith('#'):
                continue
            if "./." in line: 
                continue
            rows = (line.strip()).split('\t')
            qualities = rows[7].split(';')
            for i in qualities: 
                if "QD" in i: 
                    qd = str(i) 
                if "MQ=" in i:
                    mq = str(i)
                if "FS" in i:
                    fs = str(i)
                if "MQRankSum" in i:
                    mqranksum = str(i)
                if "ReadPosRankSum" in i:
                    readposranksum = str(i)
            if qd != "":        
                QD = float(qd.split('=')[1])
            if mq != "": 
                MQ = float(mq.split('=')[1])
            if fs !="": 
                FS = float(fs.split('=')[1])
            if mqranksum !="": 
                MQRankSum = float(mqranksum.split('=')[1])
            if readposranksum !="": 
                ReadPosRankSum = float(readposranksum.split('=')[1])
            sid = rows[0] +'.'+rows[1] 
            ref = 0
            alt = 0

            if sid in gtsnps.keys():
                gt = rows[9].split(":")
                dp = int(gt[2])
                if dp > mindp and QD > minqd :
                    #if float(MQ) > 30 and float(FS) < 60 and float(MQRankSum)  > -12.5 and float(ReadPosRankSum) > -8.0:
                        ad = gt[1].split(',')
                        if gt[0] == "0/0":
                           ref = int(ad[0])
                           alt = int(ad[1])
                        elif gt[0] =="1/1":
                           ref = int(ad[0])
                           alt = int(ad[1])
                        elif gt[0] == "0/1":
                           ref = int(ad[0])
                           alt = int(ad[1])
                        print (ref,alt,sep='\t', file=acfile)
                        print(str(gtsnps[sid]), file=gtfile)
    acfile.close()
    gtfile.close()

if __name__ == '__main__':
    main()
