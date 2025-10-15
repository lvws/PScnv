#!/usr/bin/env python
import os,sys,pysam
from collections import defaultdict
import pyfaidx
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import pandas as pd

def getGC(bed_file):
    gc_dic = {}
    with open(bed_file) as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3]
            if chrom in {'X','Y'}: # 不做性染色体
                continue
            key = '\t'.join([chrom,start,end])
            sequence = genome[chrom][int(start):int(end)].seq.upper()
            # 计算GC含量
            gc_count = sequence.count('G') + sequence.count('C')
            total_bases = len(sequence)
            gc_content = gc_count / total_bases if total_bases > 0 else 0
            gc_dic[key] = gc_content
    return gc_dic

def count_reads_in_bed(bam_file, bed_file):
    """
    计算BED区间内的reads数量
    :param bam_file: 输入的BAM文件路径
    :param bed_file: BED格式的区间文件
    :return: BED区间-reads 数量 的字典
    """
    read_dic = defaultdict(lambda:[])
    bam = pysam.AlignmentFile(bam_file, "rb")
    regins = []
    
    with open(bed_file) as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3]
            key = '\t'.join([chrom,start,end])
            if chrom in {'X','Y'}: # 不做性染色体
                continue
            regins.append(key)
            for read in bam.fetch(chrom, int(start), int(end)):
                read_dic[key].append(read.query_name)
    dic = {}
    for key in regins:
        reads_num = len(set(read_dic[key]))
        dic[key] = reads_num
    return dic

# 构建dataframe 方便后续处理
def makeDataFrame(read_count_dic,gc_dic):
    chrom_lst = []
    start_lst = []
    end_lst = []
    reads_ori_lst = []   # 原始的reads 数量
    rd_dis_fix_lst = []  # reads 数量做区间长度的校正
    gc_lst = []
    
    for key in read_count_dic:
        chrom,start,end = key.split('\t')
        chrom_lst.append(chrom)
        start_lst.append(int(start))
        end_lst.append(int(end))
        reads_ori_lst.append(read_count_dic[key])
        rd_dis_fix_lst.append(read_count_dic[key]/(int(end)-int(start)))
        gc_lst.append(gc_dic[key])
        
    df = pd.DataFrame({'Chrom':chrom_lst,'Start':start_lst,'End':end_lst,\
        'Reads_Num_Ori':reads_ori_lst,'Read_Num_Length_Fix':rd_dis_fix_lst,'GC':gc_lst})
    return df

# GC校正 + 深度校正
# def normalize(df):
#     # GC 校正
#     # 使用 LOWESS 拟合 GC 含量与测序深度之间的关系
#     fitted = lowess(df['Read_Num_Length_Fix'], df['GC'], frac=0.3, it=3)  # frac 是平滑参数，可根据需要调整
#     # 直接添加拟合结果到原DataFrame，避免merge导致的重复
#     df = df.copy()  # 避免修改原DataFrame
#     df['Read_Num_GC'] = fitted[:, 1]  # 直接使用拟合结果
    
#     df['Read_Normalized'] = df['Read_Num_Length_Fix'] / np.median(df['Read_Num_Length_Fix'])
#     # 使用GC校正后中位值，再进行测序深度校正
#     read_median_dp = np.median(df['Read_Num_GC'])
#     df['Read_GC_Normalized'] = df['Read_Num_GC'] / read_median_dp
    
#     return df

def normalize(df):
    df = df.copy()
    
    # 1. 保持原始索引以确保顺序一致
    original_index = df.index
    
    # 2. 对GC排序以进行LOWESS拟合
    sorted_by_gc = df.sort_values('GC')
    gc_sorted = sorted_by_gc['GC'].values
    read_num_sorted = sorted_by_gc['Read_Num_Length_Fix'].values
    
    # 3. 在排序后的数据上执行LOWESS拟合
    fitted = lowess(read_num_sorted, gc_sorted, frac=0.3, it=3)
    
    # 4. 将拟合值重新映射回原始顺序
    # 创建Series，索引为排序后的索引，值为拟合值
    fitted_series = pd.Series(fitted[:, 1], index=sorted_by_gc.index)
    # 按原始索引重新排序
    gc_fit = fitted_series.reindex(original_index).values
    
    # 5. 更稳健的零值处理
    median_fitted = np.median(gc_fit)
    # 使用平滑过渡而非硬替换
    # 设置最小阈值，但保持相对关系
    min_threshold = median_fitted * 0.01
    # 使用对数空间平滑而不是直接替换
    gc_fit = np.where(gc_fit < min_threshold, 
                     min_threshold * (1 + np.log1p(gc_fit/min_threshold)), 
                     gc_fit)
    
    # 6. GC校正
    df['Read_Num_GC'] = (df['Read_Num_Length_Fix'] / gc_fit) * median_fitted
    
    # 后续步骤保持不变...
    read_median_dp = np.median(df['Read_Num_GC'])
    df['Read_GC_Normalized'] = df['Read_Num_GC'] / read_median_dp
    df['Read_Normalized'] = df['Read_Num_Length_Fix'] / np.median(df['Read_Num_Length_Fix'])
    
    return df


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(sys.argv[0],'bam\tbed\ths37d5.fa\tresult.tsv')
        sys.exit(1)
    
    bamfile = sys.argv[1]
    bedfile = sys.argv[2]
    fastafile = sys.argv[3]   
    outfile = sys.argv[4]
     
    genome = pyfaidx.Fasta(fastafile)
    read_count_dic = count_reads_in_bed(bamfile,bedfile)
    gc_dic = getGC(bedfile)
    df = makeDataFrame(read_count_dic,gc_dic)
    dfnorm = normalize(df)
    
    dfnorm.to_csv(outfile,sep='\t',index=False)
