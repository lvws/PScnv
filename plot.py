import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import sys,os

if len(sys.argv) < 3:
    print(sys.argv[0],'sample\tgene(MTAP,MET)')
    sys.exit(1)
    
sample = sys.argv[1]
highlight_genes = sys.argv[2].split(',')
path   = os.path.abspath('.')
# ----------------------------
# 一、全局绘图参数：提高 DPI、字体和线宽
# ----------------------------
mpl.rcParams.update({
    'figure.dpi':       300,    
    'savefig.dpi':      300,    
    'font.size':        10,    
    'axes.linewidth':   1.0,    
    'lines.linewidth':  1.0,    
    'xtick.major.size': 4,      
    'ytick.major.size': 4,      
    'legend.fontsize':  9,      
    'legend.frameon':   True,   
})

# ----------------------------
# 二、读取数据并计算 x 坐标（略，同前）
# ----------------------------
# path   = '/dssg07/InternalResearch07/cln02/20250729_CNV/PRIME_v7.1.0'
# sample = 'TB257N0627-S1APP92KXF1-L000KBUX'
# sample = "FB24750324-L1AAP9XXXFX-Y04G"
df     = pd.read_csv(f'{path}/output/{sample}.extron.tsv', sep='\t')
df     = df[df['Chrom'].isin(range(1,23))].copy()
df['Chrom_num'] = df['Chrom'].astype(int)
df.sort_values(['Chrom_num','Start'], inplace=True)
df.reset_index(drop=True, inplace=True)

# 等距铺点
df['order'] = df.groupby('Chrom_num').cumcount()
df['count'] = df.groupby('Chrom_num')['Start'].transform('count')
df['x_rel'] = df['order'] / (df['count'] - 1)
df.loc[df['count']==1, 'x_rel'] = 0.5
df['x'] = df['Chrom_num'] + df['x_rel']

# 目前没有高亮基因
# highlight_genes = ["MTAP"]
df['Highlight'] = df['Gene'].isin(highlight_genes)

# 配色字典（若 highlight_genes 非空可用）
colors    = cm.rainbow(np.linspace(0,1,len(highlight_genes)))
color_map = dict(zip(highlight_genes, colors))

# ----------------------------
# 三、四行一列子图
# ----------------------------
cols = [
    'Read_Ratio_log2',
    'Fix_Ratio_log2',
    'SFix_Ratio_log2',
    'Seg_Ratio_log2'
]

fig, axes = plt.subplots(
    nrows=4, ncols=1,
    sharex=True,
    figsize=(16, 10),  # 四行共 10 单位高，每行约 2.5
    dpi=300
)

labeled = set()
for ax, col in zip(axes, cols):
    # 普通点
    normal = df[~df['Highlight']]
    ax.scatter(
        normal['x'], normal[col],
        s=5, c='gray', alpha=1
    )

    # 高亮点（如果有的话）
    for _, row in df[df['Highlight']].iterrows():
        ax.scatter(
            row['x'], row[col],
            s=20,
            color=color_map[row['Gene']],
            alpha=1.0
        )
        # （可选）添加基因标签
        #     用相同方式为第一条出现的 gene 加 annotate
        if row['Gene'] not in labeled:
            labeled.add(row['Gene'])
            ax.annotate(
                row['Gene'],
                xy=(row['x'], row['Seg_Ratio_log2']),
                xytext=(5, 5), textcoords='offset points',
                fontsize=11, color=color_map[row['Gene']],
                bbox=dict(boxstyle="round,pad=0.2",
                        fc="white", ec=color_map[row['Gene']], alpha=0.8)
            )
        
    # 横向和纵向参考线
    ax.axhline(0, color='gray', linestyle='--', alpha=0.6)
    ax.set_ylim(-2,2)
    # 美化
    ax.set_ylabel(col)
    ax.grid(axis='y', alpha=0.2)

# 最下面子图设置 x 轴
axes[-1].set_xlim(1, 23)
axes[-1].set_xlabel('Chromosome (1–22)')
centers = np.arange(1.5, 23.5)
axes[-1].set_xticks(centers)
axes[-1].set_xticklabels([str(i) for i in range(1,23)])

# 染色体分隔线（在所有子图上画）
for ax in axes:
    for c in range(2,23):
        ax.axvline(c, color='gray', linestyle='--', alpha=0.3)

plt.tight_layout()

# ----------------------------
# 四、保存图像
# ----------------------------
out_png = f'{path}/output/{sample}.cnv_4panel.png'
out_pdf = f'{path}/output/{sample}.cnv_4panel.pdf'
plt.savefig(out_png, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')

