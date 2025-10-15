#!/usr/bin/env python
import os,sys
from collections import defaultdict
import numpy as np

""" 
生成统计文件： gene + 外显子
"""

if len(sys.argv) < 5:
    print(sys.argv[0],"sample.log2.tsv\tsample.Dis_Fix.seg.tsv\tpanle.bed\toutfold")
    sys.exit(1)
    
count_file = sys.argv[1]
seg_file = sys.argv[2]
bed_file = sys.argv[3]
outfold = sys.argv[4]
sample = os.path.basename(count_file).split('.')[0]

# 获取位置->位置外显子信息
def getGenePos(bed_file):
    extron_dic = {}
    with open(bed_file,'r') as f:
        for i in f:
            t = i.strip().split('\t')
            chrom = t[0]
            if chrom in {'X','Y'}: # 不做性染色体
                continue
            start = t[1]
            end = t[2]
            extron = t[3].split('|')[0]
            key = '\t'.join([chrom,start,end])
            extron_dic[key] = extron
    return extron_dic

# 获取位置 -> seq_log2 信息
def getSegInfo(seg_file):
    seg_dic = {}
    with open(seg_file,'r') as f:
        heads = next(f).strip().split('\t')
        i_chr = heads.index('chrom')
        i_s = heads.index('loc.start')
        i_e = heads.index('loc.end')
        i_seg = heads.index('seg.mean')
        for i in f:
            t = i.strip().split('\t')
            chrom = t[i_chr]
            start = t[i_s]
            end = t[i_e]
            seg = t[i_seg]
            key = '\t'.join([chrom,start,end])
            seg_dic[key] = seg
    return seg_dic

# 获取位置 -> read_count 信息
def getRdInfo(count_file):
    rd_dic = {}
    with open(count_file,'r') as f:
        heads = next(f).strip().split('\t')
        i_chr = heads.index('Chrom')
        i_s = heads.index('Start')
        i_e = heads.index('End')
        i_ori = heads.index('Read_Ratio_Ori')
        i_fix = heads.index('Read_Ratio_Fix')
        i_sfix = heads.index('Read_Ratio_Fix_Smooth')
        for i in f:
            t = i.strip().split('\t')
            chrom = t[i_chr]
            start = t[i_s]
            end = t[i_e]
            ori_log2 = min(max(-4,float(t[i_ori])),4)
            ori = 2**float(ori_log2)
            fix_log2 =  min(max(-4,float(t[i_fix])),4)
            fix = 2**float(fix_log2)
            sfix_log2 = min(max(-4,float(t[i_sfix])),4)
            sfix = 2**float(sfix_log2)
            key = '\t'.join([chrom,start,end])
            rd_dic[key] = [ori,ori_log2,fix,fix_log2,sfix,sfix_log2]
    return rd_dic

def checkInRegin(key1,key2):
    c1,s1,e1 = key1.split('\t')
    c2,s2,e2 = key2.split('\t')
    if c1 != c2:
        return False
    if int(s1) >= int(s2) and int(e1) <= int(e2):
        return True
    return False

extron_dic = getGenePos(bed_file)
seg_dic = getSegInfo(seg_file) # CBS 后的 log2 ratio 值， 区域进行了合并
rd_dic = getRdInfo(count_file) # 根据bed 文件，计算的每个区域的 rd 标准化的值

# 以bed 区域为基础，注释出 基因，外显子， rd 、 seg 
extron_file = f'{outfold}/{sample}.extron.tsv'
w1 = open(extron_file,'w')
w1.write('ID\tChrom\tStart\tEnd\tGene\tGene.ID\tRead_Ratio\tRead_Ratio_log2\tFix_Ratio\tFix_Ratio_log2\tSFix_Ratio\tSFix_Ratio_log2\tSeg_Ratio\tSeg_Ratio_log2\n')
# gene_dic = defaultdict(lambda:[[],[],[],[]]) # rds \ segs \ regins
# 添加一个小的偏移量,不然会因为无限值报错
epsilon = 1e-12
for key in extron_dic:
    extron = extron_dic[key]
    ori,ori_log2,fix,fix_log2,sfix,sfix_log2 = rd_dic[key]
    for key2 in seg_dic:
        if checkInRegin(key,key2):
            seg_log2 = seg_dic[key2]
            break
    gene = extron.split(':')[0]
    # gene_dic[gene][0].append(float(rd))
    # gene_dic[gene][1].append(float(seg_log2))
    # gene_dic[gene][2].append(key)
    # 控制大小，防止溢出
    seg_log2 = min(float(seg_log2),4)
    seg_log2 = max(float(seg_log2),-4)
    seg_rate = 2**float(seg_log2)
    w1.write(f'{sample}\t{key}\t{gene}\t{extron}\t{ori}\t{ori_log2}\
        \t{fix}\t{fix_log2}\t{sfix}\t{sfix_log2}\t{seg_rate}\t{seg_log2}\n')
w1.close()

# gene_file = f'{outfold}/{sample}.gene.tsv'
# w2 = open(gene_file,'w')
# w2.write('ID\tChrom\tStart\tEnd\tGene\tRead_Ratio\tRead_Ratio_log2\tSeg_Ratio\tSeg_Ratio_log2\n')
# for gene in gene_dic:
#     rds,segs,rgs = gene_dic[gene]
#     rd = np.median(rds)
#     rd_log2 = np.log2(float(rd)+epsilon)
#     seg_log2 = np.median(segs)
#     chrom,start,_ = rgs[0].split('\t')
#     end = rgs[-1].split('\t')[-1]
#     # 控制大小，防止溢出
#     seg_log2 = min(float(seg_log2),4)
#     seg_log2 = max(float(seg_log2),-4)
#     seg_rate = 2**float(seg_log2)
#     w2.write(f'{sample}\t{chrom}\t{start}\t{end}\t{gene}\t{rd}\t{rd_log2}\t{seg_rate}\t{seg_log2}\n')
# w2.close()
