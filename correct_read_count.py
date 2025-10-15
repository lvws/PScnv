#!/usr/bin/env python
import os,sys
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from matplotlib.font_manager import FontProperties
from sklearn.preprocessing import PolynomialFeatures
from scipy.stats import pearsonr, spearmanr
from sklearn.preprocessing import MinMaxScaler

def bed_check(bedfile):
    """ 
    # 计算bed 文件的区间大小，和数量，方便后续的染色体选择
    """
    # print('Chrom\tregins\tlength')
    dic = defaultdict(lambda:[])
    with open(bedfile,'r') as f:
        for i in f:
            chrom,start,end,*_ = i.strip().split('\t')
            if chrom in {'X','Y'}: # 不做性染色体
                continue
            dic[chrom].append(int(end)-int(start))

    ndic = {}
    for chrom in dic:
        vals = dic[chrom]
        # print(f"{chrom}\t{len(vals)}\t{sum(vals)}")
        ndic[chrom] = [len(vals),sum(vals)]
    return ndic

# 读取样本的Reads_Counts_Norm， 选择最稳定的染色体
def read_counts_norm(infile,chrom,normRD):
    dp_norms = []
    with open(infile,'r') as f:
        heads = next(f).strip().split('\t')
        i_chr = heads.index('Chrom')
        i_rcn = heads.index(normRD)
        for i in f:
            t = i.strip().split('\t')
            if t[i_chr] == chrom:
                dp_norms.append(float(t[i_rcn]))
    return dp_norms

# pandas 进行处理，方便后续
def choose_chrom(yfile,cfile,bed_dic,normRD,CnormRD,CnormSD):
    """ 
    利用公式，计算和对照样本测序标准化深度距离最接近的样本的染色体
    yfile: 样本的标准化后的bed 区域测序深度文件
    cfile: 对照样本标准化后的bed 区域测序深度文件，中位和标准差
    bed_dic: target 区域的 数量和 区域大小
    """
    sdf = pd.read_csv(yfile,
                    sep="\t",        # 制表符分隔
                    header=0,        # 第一行作为列名；若文件无表头可设为 None
                    index_col=None)  # 不将任何列作为索引；可传列名或列号)
    cdf = pd.read_csv(cfile,sep='\t',header=0,index_col=None)
    df = sdf.merge(cdf,on=["Chrom","Start","End"],how='left')
    # 添加一个小的偏移量,不然会因为无限值报错
    epsilon = 1e-12
    df[normRD] += epsilon
    df[CnormRD] += epsilon
    df[CnormSD] += epsilon
    df['Diff'] = abs(df[normRD]-df[CnormRD] )/df[CnormSD]
    
    # 1. 先算 log2
    df['logRD'] = np.log2(df[normRD] + epsilon)
    # 2. 再裁剪到 [-2, 2]
    df['logRD'] = df['logRD'].clip(-2, 2)   
    df['logCRD'] = np.log2(df[CnormRD] + epsilon)
    df['logCRD'] = df['logCRD'].clip(-2, 2) 
    df['logDiff'] = df['logRD'] - df['logCRD']
    
    select = ''
    diffp = 999999
    logdiff = 0
    pcor = 0
    # rag_sum = 0
    # for chrom in bed_dic:
    #     rag_sum += bed_dic[chrom][0]
    # rag_avg = rag_sum/22 
    results = []
    print("chrom\tregins\tregin_lenght\tdiff\tdiff2\tlog_diff\tpearson_corr")
    for chrom in bed_dic:
        if chrom in {'X','Y'}:
            continue
        rg,leg = bed_dic[chrom] # 染色体区域数量，区域大小
        diff = sum(df[df["Chrom"]==int(chrom)]['Diff'])/rg
        diff2 = sum(df[df["Chrom"]==int(chrom)]['Diff'])/rg/np.sqrt(rg)
        mean_log2 = sum(df[df["Chrom"]==int(chrom)]['logDiff'])/rg
        
        # 相关性
        sub = df.loc[df["Chrom"] == int(chrom), [normRD, CnormRD]]
        pearson_corr = sub[normRD].corr(sub[CnormRD])

        
        print(f'{chrom}\t{rg}\t{leg}\t{diff}\t{diff2}\t{mean_log2}\t{pearson_corr}')
        if diff2 < diffp:
            select = chrom
            diffp = diff2
            logdiff = mean_log2
        # if pearson_corr > pcor:
        #     select = chrom
        #     logdiff = mean_log2
        #     pcor = pearson_corr
        results.append({
        'chrom':          chrom,
        'regions':        rg,
        'region_length':  leg,
        'diff':           diff,
        'diff2':          diff2,
        'log_diff':       mean_log2,
        'pearson_corr':   pearson_corr
        })
    result_df = pd.DataFrame(results)
    # 1. 计算 abs_diff2
    result_df['abs_log_diff'] = result_df['log_diff'].abs()
    # 2. 用 Min–Max 归一化两列到 [0,1]
    scaler = MinMaxScaler()
    result_df[['n_log_diff','n_pearson']] = scaler.fit_transform(
        result_df[['abs_log_diff','pearson_corr']]
    )
    # 3. 计算到理想点 (0,1) 的欧氏距离
    #    注意要用 (n_diff2 - 0) 和 (n_pearson - 1)
    result_df['distance'] = np.sqrt(
        (result_df['n_log_diff'] - 0)**2 +
        (result_df['n_pearson'] - 1)**2
    )

    # 4. 排序并取最小距离
    best_df = result_df.sort_values('distance')
    return df,best_df['chrom'].iloc[0],best_df['log_diff'].iloc[0]
    
# 绘图展示拟合效果
def plot(real,pred,score):
    # 绘制拟合效果图
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.scatter(range(len(real)), real, color='blue', label='True', s=100)
    plt.plot(range(len(real)), pred, 'r--', marker='o', label='Predict')
    plt.title('回归模型拟合效果 (R2=%.2f)' % score,fontproperties=custom_font)
    plt.xlabel('样本编号', fontproperties=custom_font)
    plt.ylabel('目标值', fontproperties=custom_font)
    plt.legend()
    plt.grid(True)

    # 残差
    plt.subplot(2, 1, 2)
    residuals = real - pred
    plt.scatter(pred, residuals, color='purple', alpha=0.6)
    plt.axhline(y=0, color='black', linestyle='--')
    plt.xlabel('Predict')
    plt.ylabel('Diff')
    plt.title('残差分析', fontproperties=custom_font)
    # 设置y轴范围
    plt.ylim(-1.5, 1.5)  # 设置y轴下限为-1.5，上限为1.5
    # 设置y轴刻度
    plt.yticks(np.arange(-1.5, 1.5, 0.5))  # 从-1.5到2.0，步长为0.5

    plt.tight_layout()
    plt.show()
    
# 读取样本所有位置的Read_Normalized
def read_counts_norm_all(infile,normRD):
    dp_norms = []
    pos_lst = []
    with open(infile,'r') as f:
        heads = next(f).strip().split('\t')
        i_chr = heads.index('Chrom')
        i_s = heads.index('Start')
        i_e = heads.index('End')
        i_rcn = heads.index(normRD)
        for i in f:
            t = i.strip().split('\t')
            chrom = t[i_chr]
            s = t[i_s]
            e = t[i_e]
            if float(t[i_rcn]) < 0 :
                print(infile,chrom,s,e)
            dp_norms.append(float(t[i_rcn]))
            pos = '\t'.join([chrom,s,e])
            pos_lst.append(pos)
    return dp_norms,pos_lst
    
# 对染色体分段进行平滑降噪    
def assign_split_numbers(
    df: pd.DataFrame,
    chrom_col: str = 'Chrom',
    value_col: str = 'Read_Ratio_Fix',
    z_thresh: float = 3.0,
    min_sd: float = 0.1,
    split_col: str = 'Split_Num'
) -> pd.DataFrame:
    """
    在 DataFrame 中，按染色体分组，对每组 value_col 列做滑动 z‐score 分段，
    并在原表新增 split_col 列标注每条记录所属的段编号。

    参数
    ----
    df : pd.DataFrame
        原始数据表，必须包含 chrom_col 和 value_col 两列。
    chrom_col : str
        用于分组的“染色体”列名。
    value_col : str
        做分段计算的数值列名（如 Read_Ratio_Fix）。
    z_thresh : float
        z‐score 阈值，默认 3.0。新值与当前段均值的 z-score 超过此阈
        值时，启动新段。
    min_sd : float
        当当前段计算出的标准差小于此值时，替代为 min_sd，
        避免除以 0 或极小值导致的过度拆段。
    split_col : str
        在返回表中新增的分段编号列名。

    返回
    ----
    pd.DataFrame
        原始表的深拷贝，新增 split_col 列，表示每行的段 ID（从 1 开始计数）。
    """
    # 深拷贝输入，避免原表被修改
    out = df.copy()
    # 先创建一个全 0 的列，用于存放分段编号
    out[split_col] = 0

    # 按染色体分组，逐组分段
    for chrom, grp in out.groupby(chrom_col):
        idxs = grp.index.tolist()         # 该组所有行的原始索引
        values = grp[value_col].tolist()  # 该组所有点的数值
        split_ids = []                    # 存放对应每个值的段编号

        current_vals = []  # 缓存当前段所有值，用于实时计算均值/标准差
        seg_id = 1         # 从 1 开始的段编号

        for x in values:
            if not current_vals:
                # 第一条自动加入当前段
                current_vals.append(x)
                split_ids.append(seg_id)
                continue

            arr = np.array(current_vals, dtype=float)
            mu  = arr.mean()
            sd  = arr.std(ddof=1) if len(arr) > 1 else 0.0

            # 避免 sd 过小
            sd = max(sd, min_sd)

            z = abs(x - mu) / sd

            if z <= z_thresh:
                # 属于当前段
                current_vals.append(x)
                split_ids.append(seg_id)
            else:
                # 启动新段
                seg_id += 1
                current_vals = [x]
                min_sd = max(abs(x)*0.2,0.2)
                split_ids.append(seg_id)

        # 将本组的 split_ids 写回到输出表
        out.loc[idxs, split_col] = split_ids

    return out

# 高斯平滑
def smooth_exp_kernel_mixed(
    df: pd.DataFrame,
    alpha: float = 0.8,             # 平滑强度因子
    chrom_col: str = 'Chrom',
    seg_col:   str = 'Split_Num',
    coord_col: str = 'Start',
    value_col: str = 'Read_Ratio_Fix',
    smooth_col:str = 'Read_Ratio_Fix_Smooth'
) -> pd.DataFrame:
    """
    指数核平滑 + 原始值线性混合：
    1) 先用 exp(-|d|/length_scale) 做完全平滑
    2) 再用 alpha*(平滑值) + (1-alpha)*(原始值) 
       来控制平滑幅度
    """
    out = df.copy()
    out[smooth_col] = np.nan

    for (chrom, seg), grp in out.groupby([chrom_col, seg_col]):
        grp = grp.sort_values(coord_col)
        idxs  = grp.index.to_numpy()
        coords= grp[coord_col].to_numpy(float)
        vals  = grp[value_col].to_numpy(float)

        # 1) 计算距离矩阵
        dist = np.abs(coords[:,None] - coords[None,:])
        span   = coords.max() - coords.min()
        length_scale = ( (span/100)+1)*1
        # 2) 指数核权重
        W    = np.exp(-dist / length_scale)
        W   /= W.sum(axis=1, keepdims=True)
        # 3) 完全平滑结果
        sm0  = W.dot(vals)
        # 4) 与原始值按比例混合
        sm1  = alpha * sm0 + (1 - alpha) * vals

        out.loc[idxs, smooth_col] = sm1

    return out

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(sys.argv[0],'sample\tbed')
        sys.exit(1)
    
    sample = sys.argv[1]
    bed = sys.argv[2]
    path = os.path.abspath('.')
    if not os.path.exists(bed):
        bed = f'{path}/{bed}'
    # 指定字体文件路径
    font_path = '/usr/share/fonts/cjkuni-uming/uming.ttc'  # Windows系统黑体路径
    custom_font = FontProperties(fname=font_path)
    yfile = f'{path}/results/{sample}.count.tsv'
    cfile = f'{path}/control-reads-counts.tsv'
    normRD = 'Read_GC_Normalized'  # 降噪后的reads ratio
    CnormRD = 'Read_Counts_Median_GC' # 降噪后的基线样本 reads ratio 中位
    CnormSD = 'Sd_GC'  # 降噪后的基线样本 reads ratio 标准差
    bed_dic = bed_check(bed)
    df,chrom,logdiff = choose_chrom(yfile,cfile,bed_dic,normRD,CnormRD,CnormSD) # 选择与基线样本最接近的染色体，作为自校正的基板
    print("选择稳定染色体号：",chrom)
    # 根据上面的计算结果，排除 性染色体， 选择Differ_Per_Regin最小的染色体
    cdp_norms = []
    with open(f'{path}/id.control','r') as f:
        for i in f:
            s = i.strip()
            sfile = f'{path}/results/{s}.count.tsv'
            cdp_norms.append(read_counts_norm(sfile,chrom,normRD))
    y_norms = read_counts_norm(yfile,chrom,normRD)
    # 添加一个小的偏移量,不然会因为无限值报错
    epsilon = 1e-12
    X = np.log2(np.array(cdp_norms).T + epsilon)
    y = np.log2(np.array(y_norms).T + epsilon)  - logdiff
    
    # 构建标准化+Ridge回归模型
    model = make_pipeline(
        # PolynomialFeatures(degree=2),  # 二次多项式扩展
        StandardScaler(),
        Ridge(alpha=1.0)
    )
    model.fit(X, y)

    # # 预测结果
    # y_pred = model.predict(X)
    # score = model.score(X,y)
    # plot(y,y_pred,score)

    # 用模型预测正常的染色体的log2 RD 分布
    cdp_norms = []
    with open(f'{path}/id.control','r') as f:
        for i in f:
            s = i.strip()
            sfile = f'{path}/results/{s}.count.tsv'
            dp_norms,_ = read_counts_norm_all(sfile,normRD)
            cdp_norms.append(dp_norms)
    y1_norms,pos_lst = read_counts_norm_all(yfile,normRD)

    # 添加一个小的偏移量,不然会因为无限值报错
    epsilon = 1e-12
    X1 = np.log2(np.array(cdp_norms).T + epsilon)
    y1 = np.log2(np.array(y1_norms).T + epsilon)

    y1_pred = model.predict(X1)


    y_fix = y1 - y1_pred
    pos_dic = defaultdict(lambda:[])
    chrom_dic_ori = defaultdict(lambda:[])
    chrom_dic_fix = defaultdict(lambda:[])
    chrom_lst = []
    start_lst = []
    end_lst = []
    for i in range(len(pos_lst)):
        chrom,s,e = pos_lst[i].split('\t')
        pos_dic[chrom].append([int(s),int(e),y_fix[i]])
        chrom_dic_ori[chrom].append(y1[i])
        chrom_dic_fix[chrom].append(y_fix[i])
        chrom_lst.append(chrom)
        start_lst.append(s)
        end_lst.append(e)
    ndf = pd.DataFrame({'Chrom':chrom_lst,'Start':start_lst,'End':end_lst,'Read_Ratio_Ori':y1,'Read_Ratio_Fix':y_fix})
    
    
    # 测试,分隔
    sdf = assign_split_numbers(ndf, min_sd=0.05)
    sdf.to_csv(f'{path}/output/{sample}.split.tsv',sep='\t')
    
    # 基于距离的高斯平滑
    disdf = smooth_exp_kernel_mixed(sdf)
    disdf.to_csv(f'{path}/output/{sample}.smooth.tsv',sep='\t')
    cnvfile = f'{path}/output/{sample}.log2.tsv'
    disdf.to_csv(cnvfile,index=None,sep='\t')
    