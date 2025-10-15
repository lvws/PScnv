#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np

def getCNR(infile,gtag,cnrtag,gene):
    sdf = pd.read_csv(infile,
                    sep="\t",        # 制表符分隔
                    header=0,        # 第一行作为列名；若文件无表头可设为 None
                    index_col=None)  # 不将任何列作为索引；可传列名或列号)
    mask = sdf[gtag].str.split(':').str[0] == gene
    df = sdf[mask] 
    cnr = np.median(df[cnrtag])
    return cnr

if len(sys.argv) < 5:
    print(sys.argv[0],'extron/genemetrics.txt\tgatg(Gene/gene)\tcnrtag(SFix_Ratio_log2\tlog2)\tgene')
    sys.exit(1)
    
cnr = getCNR(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
del_thd = -0.2
dup_thd = 0.246
sample = os.path.basename(sys.argv[1]).split('.')[0]
if cnr > dup_thd:
    cnv = 'dup'
elif cnr < del_thd:
    cnv = 'del'
else:
    cnv = 'normal'
print(f'{sample}\t{cnr}\t{cnv}')