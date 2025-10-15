# PScnv

**PScnv** is a bioinformatics tool for detecting **Copy Number Variations (CNVs)** from **WGS, WES, and targeted panel sequencing data**.  
It leverages a baseline reference built from healthy control samples to improve CNV detection accuracy without requiring matched negative controls.

---

## ✨ Features

- Uses **50+ healthy control samples** as a baseline reference.
- Identifies the **most stable chromosome** in each test sample for self-normalization.
- Reduces experimental noise and batch effects by leveraging baseline stability.
- Works with **WGS, WES, and targeted panel sequencing** data.
- No need for matched negative control samples.

---

## 📦 Installation

```bash
git clone git@github.com:lvws/PScnv.git
cd PScnv
```
---
## Dependencies

- Python 3.x
- R (with DNAcopy package for CBS algorithm)
- SAMtools
---
## Usage
---
### Preparation Phase
---
Prepare your sample files:
- `id.control`: List of healthy control sample names
- `id.tumor`: List of test sample names  
- `BAM/`: Directory containing all BAM files
- `test.bed`: Target regions with gene information in 4th column (format: `chr\tstart\tend\tgene:transcript:exon`)
- `refer.fa`: Reference genome FASTA file

### Step 1: Calculate Read Counts and Normalization
---
```bash
for i in `cat id.control id.tumor`;do 
    echo "python read_count.py BAM/${i}.bam test.bed refer.fa results/${i}.count.tsv"
done | sh
```

### Step 2: Calculate Baseline Median and Standard Deviation
---
```bash
python control_median_sd.py id.control > control-reads-counts.tsv
```

### Step 3: Self-chromosome Correction Algorithm
---
```bash
for i in `cat id.tumor`;do 
    echo "python correct_read_count.py ${i} test.bed"
done > cnv-pre.sh
```

### Step 4: CBS Algorithm for Segmentation
---
```bash
for i in `cat id.tumor`;do 
    echo "Rscript cbs.R output/${i}.log2.tsv ${i} Read_Ratio_Fix_Smooth output"
done > cbs.sh
```

### Step 5: Gene-level Amplification Status
---
```bash
for i in `cat id.tumor`;do 
    echo "python extron_info.py output/${i}.log2.tsv output/${i}.Read_Ratio_Fix_Smooth.seg.tsv test.bed output/"
done | sh

for i in `cat id.tumor`;do 
    echo "python getGeneCNV.py output/${i}.extron.tsv Gene SFix_Ratio_log2 <Gene>"
done | sh
```

### Step 6: Visualization
---
```bash
for i in `cat id.tumor`;do 
    echo "python3.9 plot.py ${i}"
done | sh
```

## Output Files
---
- `results/*.count.tsv`: Normalized read counts
- `control-reads-counts.tsv`: Baseline statistics
- `output/*.log2.tsv`: Log2 ratio files
- `output/*.seg.tsv`: CBS segmentation results
- `output/*.extron.tsv`: Exon-level CNV information
- `output/*.png/pdf`: CNV visualization plots


## Citation
---
If you use PScnv in your research, please cite:

[Citation information to be added]

## License
---
This project is licensed under the MIT License.


