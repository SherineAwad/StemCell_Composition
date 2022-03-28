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
    N = int(args.n)
    mindp = int(args.mindp)
    minqd = float(args.minqd)
    gtsnps = {}
    alsnps = {} 
    vid = [] 
    count =0 
    outname = args.outfile 
    acfile = open(outname+"_allelespecific.txt", "w+")
    print("rn,1",file=acfile)
    #gtfile =open(outname+"_genotype.txt", 'w+')
    #print('D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,rn',file=gtfile)
    for line in open(args.infile1):
            record = line.split(',')
            tmp_genotype =""
            for i in range(0,N):
                    tmp_genotype+=(record[i].strip())
                    tmp_genotype+=','
            gtsnps[record[N].strip()] = tmp_genotype 
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
            af = qualities[1] 
            AF = af.split('=')[1]
            if (',') in AF: 
                 continue 
            MAF = -1
            if (float(AF )< 0.5):
                 MAF = float(AF)
            else: 
                 MAF = 1 - float(AF)
            '''
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
            ''' 
            sid = rows[0] +'_'+rows[1]+'_'+rows[3]+'_'+rows[4]
            ref = 0
            alt = 0
            if sid in gtsnps.keys(): 
                print(str(sid).strip()+','+str(MAF), file=acfile)
                gt = rows[9].split(":")
                dp = int(gt[2])
                """
                #if dp > mindp and QD > minqd :
                    #if float(MQ) > 30 and float(FS) < 60 and float(MQRankSum)  > -12.5 and float(ReadPosRankSum) > -8.0:
                        
                        ad = gt[1].split(',')
                        if gt[0] == "0/0":
                           alt = min(int(ad[0]), int(ad[1]))
                           alt = alt / (int(ad[0]) + int(ad[1]))
                        elif gt[0] =="1/1":
                           alt = min(int(ad[0]), int(ad[1])) 
                           alt = alt / (int(ad[0]) + int(ad[1]))
                        elif gt[0] == "0/1":
                           alt = min(int(ad[0]), int(ad[1]))
                           alt = alt / (int(ad[0]) + int(ad[1]))
                        #print(str(sid).strip()+','+str(alt).strip(), file=acfile)
                        print(str(gtsnps[sid])+str(sid), file=gtfile)
                """
    acfile.close()

if __name__ == '__main__':
    main()
