#!/usr/bin/env python
import os,sys
import numpy as np
from collections import defaultdict

def get_control_dp(idfile):
    """ 
    Generating statistics read count from control samples
    C_i=median(c_ai1,c_ai2,…,c_aiN).  
    σ_i=std(c_ai1,c_ai2,…,c_aiN).   
    
    N: the number of  control samples.
    c: the read count of  a target region in a control sample
    """
    dic = defaultdict(lambda:[[],[]])
    with open(idfile,'r') as f:
        for i in f:
            sample = i.strip()
            read_norm_file = f'results/{sample}.count.tsv'
            regins = []
            sys.stderr.write(f'Teat:{read_norm_file}\n')
            with open(read_norm_file,'r') as f2:
                heads = next(f2).strip().split('\t')
                i_chr = heads.index('Chrom')
                i_s = heads.index('Start')
                i_e = heads.index('End')
                i_norm_count = heads.index('Read_Normalized')
                i_gc = heads.index('Read_GC_Normalized')
                for j in f2:
                    t = j.strip().split('\t')
                    chrom = t[i_chr]
                    start = t[i_s]
                    end = t[i_e]
                    ncount = t[i_norm_count]
                    gc_count = t[i_gc]
                    key = '\t'.join([chrom,start,end])
                    dic[key][0].append(float(ncount))
                    dic[key][1].append(float(gc_count))
                    regins.append(key)
    
    print("Chrom\tStart\tEnd\tRead_Counts_Median\tSd\tRead_Counts_Median_GC\tSd_GC")
    for key in regins:
        median = np.median(dic[key][0])  # 中位值
        std_dev = np.std(dic[key][0])    # 标准差
        gc_median = np.median(dic[key][1])  # 中位值
        gc_std_dev = np.std(dic[key][1])    # 标准差
        print(f'{key}\t{median}\t{std_dev}\t{gc_median}\t{gc_std_dev}')
                
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(sys.argv[0],'id.control')
        sys.exit(1)
        
    get_control_dp(sys.argv[1])