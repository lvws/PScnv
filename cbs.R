#!/usr/bin/env Rscript

#' 命令行 CNV 分析工具
#' 
#' 用法: Rscript a.R infile sid selfix output
#' 示例: Rscript a.R a.log2.tsv a Self_Fix out

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) < 4) {
  stop("用法: Rscript a.R infile sid selfix output\n",
       "参数说明:\n",
       "  infile: 输入文件路径 (TSV格式)\n",
       "  sid: 样本ID\n", 
       "  selfix: Self_Fix列名\n",
       "  output: 输出目录\n")
}

# 解析参数
infile <- args[1]
sid <- args[2]
selfix <- args[3]
output_dir <- args[4]

# 加载所需包
if (!requireNamespace("DNAcopy", quietly = TRUE)) {
  stop("需要安装DNAcopy包: 请先运行 BiocManager::install('DNAcopy')")
}
library(DNAcopy)

#' 执行DNA拷贝数变异(CNV)分析
#'
#' @param filename 输入文件的路径和名称
#' @param self_fix_col 包含校正后log2比值的列名
#' @param sample_id 样本标识符
#' @param outdir 输出目录路径
#' @param alpha 显著性阈值，值越小越严格（默认为0.01）
#' @param nperm 排列检验次数（默认为10000）
#' @param undo_splits 撤销不重要分割的方法（默认为"sdundo"）
#' @param undo_sd 撤销分割的标准差阈值（默认为3）
#' @param min_width 最小探针数量（默认为3）
#'
#' @return 包含分段结果和统计信息的列表
analyze_cnv <- function(filename, self_fix_col, sample_id, outdir,
                        alpha = 0.01, nperm = 10000, undo_splits = "sdundo",
                        undo_sd = 3, min_width = 3) {
  
  # 确保输出目录存在
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    cat("创建输出目录:", outdir, "\n")
  }
  
  # 读取数据
  data <- read.csv(filename, sep = "\t")
  
  # 检查所需的列是否存在
  if (!all(c("Chrom", "Start", self_fix_col) %in% colnames(data))) {
    stop(paste("数据必须包含'Chromosome', 'Position'和", self_fix_col, "列"))
  }
  
  # 准备染色体和位置信息
  chromosome <- ordered(as.character(data$Chrom), 
                        levels = c(1:22, "X", "Y"))
  position <- data$Start
  
  # 创建CNA对象
  cna_object <- CNA(genomdat = data[[self_fix_col]],
                    chrom = chromosome,
                    maploc = position,
                    data.type = "logratio",
                    sampleid = sample_id)
  
  cat("\nCNA 对象创建成功:\n")
  print(cna_object)
  
  # 数据平滑处理 (处理离群值)
  cna_smoothed <- smooth.CNA(cna_object,smooth.SD.scale=2, trim=0.025,
                             smooth.region = 10,
                             outlier.SD.scale = 5)
  
  # 执行分段分析
  segments_undo <- segment(cna_smoothed,
                           alpha = alpha,
                           undo.splits = undo_splits,
                           undo.SD = undo_sd,
                           min.width = min_width,
                           verbose = 1)
  
  # 获取分段结果数据框
  segment_results <- segments_undo$output
  
  cat("\n分段结果摘要:\n")
  print(segment_results)
  
  # 计算增益(gain)和缺失(loss)的数量
  gains <- sum(segment_results$seg.mean > 0.2)
  losses <- sum(segment_results$seg.mean < -0.2)
  
  cat("拷贝数增益(seg.mean > 0.2):", gains, "个片段\n")
  cat("拷贝数缺失(seg.mean < -0.2):", losses, "个片段\n")
  
  # 生成输出文件名
  output_filename <- paste0(sample_id, ".", self_fix_col, ".seg.tsv")
  output_path <- file.path(outdir, output_filename)
  
  # 保存分段结果到文件
  write.table(segment_results, 
              file = output_path, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  cat("分段结果已保存到:", output_path, "\n")
  
  # 返回结果
  result <- list(
    cna_object = cna_object,
    cna_smoothed = cna_smoothed,
    segments = segments_undo,
    segment_results = segment_results,
    gains = gains,
    losses = losses,
    summary = data.frame(
      Sample = sample_id,
      Gains = gains,
      Losses = losses,
      Total_Segments = nrow(segment_results)
    ),
    output_file = output_path
  )
  
  return(result)
}

# 执行分析
cat("开始CNV分析...\n")
cat("输入文件:", infile, "\n")
cat("样本ID:", sid, "\n")
cat("Self_Fix列:", selfix, "\n")
cat("输出目录:", output_dir, "\n")

tryCatch({
  result <- analyze_cnv(infile, selfix, sid, output_dir)
  cat("分析完成!\n")
}, error = function(e) {
  cat("分析过程中出错:", e$message, "\n")
  quit(status = 1)
})
