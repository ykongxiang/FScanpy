# FScanpy

## English

### Introduction
FScanpy is a Python package for predicting Programmed Ribosomal Frameshifting (PRF) sites in DNA sequences. It integrates machine learning models (Gradient Boosting and BiLSTM-CNN) to provide accurate PRF predictions.

### Installation

#### 1. Using pip
```bash
pip install FScanpy
```

#### 2. Clone from github
```bash
git clone https://github.com/
cd your_project_directory
pip install -e .
```

### Testing usage(if you want to test the package, you can use the following code)
```python
from FScanpy.data import get_test_data_path, list_test_data
from FScanpy import PRFPredictor
predictor = PRFPredictor()
predictor.predict_full(sequence='ATGCGTACGTATGCGTACGTATGCGTACGT',
                              window_size=3,
                              gb_threshold=0.1)
seq_30bp = 'ATGCGTACGT' * 3  
seq_300bp = 'ATGCGTACGT' * 30
result = predictor.predict_region(seq_30bp=seq_30bp, seq_300bp=seq_300bp)
print(result)
region_example = pd.read_excel(get_test_data_path('region_example.xlsx'))
results = predictor.predict_region(seq_30bp=region_example['30bp'], 
                                  seq_300bp=region_example['300bp'])
from FScanpy.utils import fscanr
blastx_output = pd.read_excel(get_test_data_path('blastx_example.xlsx'))
fscanr_result = fscanr(blastx_output, 
                      mismatch_cutoff=10,
                      evalue_cutoff=1e-5,
                      frameDist_cutoff=10)
from FScanpy.utils import extract_prf_regions
prf_regions = extract_prf_regions(mrna_file=get_test_data_path('mrna_example.fasta'),
                                prf_data=fscanr_result)
```

### Test data
test data can be found in FScanpy/data/test_data
```python
from FScanpy.data import get_test_data_path, list_test_data
blastx_file = get_test_data_path('blastx_example.xlsx')
mrna_file = get_test_data_path('mrna_example.fasta')
region_example = get_test_data_path('region_example.xlsx')
```

### Usage
#### Predict PRF sites in a full sequence
```python
from FScanpy import PRFPredictor
predictor = PRFPredictor()
predictor.predict_full(sequence='ATGCGTACGTATGCGTACGTATGCGTACGT',
                              window_size=3,
                              gb_threshold=0.1)
```

#### Predict PRF in specific regions
```python
from FScanpy import PRFPredictor
predictor = PRFPredictor()
```
#### Method 1: Single sequence(need to ensure + strand)
```python
seq_30bp = 'ATGCGTACGT' * 3  # Ensure length is 30
seq_300bp = 'ATGCGTACGT' * 30  # Ensure length is 300
result = predictor.predict_region(seq_30bp=seq_30bp, seq_300bp=seq_300bp)
print(result)
```
#### Method 2: Batch prediction
```python
region_example = pd.read_excel(get_test_data_path('region_example.xlsx'))
results = predictor.predict_region(seq_30bp=region_example['30bp'], 
                                  seq_300bp=region_example['300bp'])


#### Identify PRF sites from BLASTX output
```python
from FScanpy.utils import fscanr
blastx_output = pd.read_excel(get_test_data_path('blastx_example.xlsx'))
fscanr_result = fscanr(blastx_output, 
                      mismatch_cutoff=10,
                      evalue_cutoff=1e-5,
                      frameDist_cutoff=10)


#### Extract analysis window sequences
```python
from FScanpy.utils import extract_prf_regions
prf_regions = extract_prf_regions(mrna_file=get_test_data_path('mrna_example.fasta'),
                                prf_data=fscanr_result)


#Chinese

### 介绍
FScanpy 是一个用于预测 mRNA 序列中程序化核糖体框架移位 (PRF) 位点的 Python 包。它集成了机器学习模型（Gradient Boosting 和 CNN），以提供准确的 PRF 预测。
### 安装
#### 1. 使用 pip
```bash
pip install FScanpy

#### 2.从 github克隆
```bash
git clone https://github.com/
pip install -e .
```
### 测试用法(如果你想测试这个包，可以使用以下代码)
```python
from FScanpy.data import get_test_data_path, list_test_data
from FScanpy import PRFPredictor
predictor = PRFPredictor()
predictor.predict_full(sequence='ATGCGTACGTATGCGTACGTATGCGTACGT',
                              window_size=3,
                              gb_threshold=0.1)
seq_30bp = 'ATGCGTACGT' * 3  
seq_300bp = 'ATGCGTACGT' * 30
result = predictor.predict_region(seq_30bp=seq_30bp, seq_300bp=seq_300bp)
print(result)
from FScanpy.utils import fscanr
blastx_output = pd.read_excel(get_test_data_path('blastx_example.xlsx'))
fscanr_result = fscanr(blastx_output, 
                      mismatch_cutoff=10,
                      evalue_cutoff=1e-5,
                      frameDist_cutoff=10)
from FScanpy.utils import extract_prf_regions
prf_regions = extract_prf_regions(mrna_file=get_test_data_path('mrna_example.fasta'),
                                prf_data=fscanr_result)
```

### 测试数据
测试数据可以在 FScanpy/data/test_data 中找到
```python
from FScanpy.data import get_test_data_path, list_test_data
blastx_file = get_test_data_path('blastx_example.xlsx')
mrna_file = get_test_data_path('mrna_example.fasta')
region_example = get_test_data_path('region_example.xlsx')
```

### 使用
#### 预测全序列中的PRF位点
```python
from FScanpy import PRFPredictor
predictor = PRFPredictor()
predictor.predict_full_sequence(sequence='ATGCGTACGTATGCGTACGTATGCGTACGT',
                              window_size=3,
                              gb_threshold=0.1)
```

#### 预测特定区域中的PRF位点
```python
from FScanpy import PRFPredictor
predictor = PRFPredictor()
# 方法1：单个序列(需要确保是+链)
predictor.predict_region(seq_30bp='ATGCGTACGT'*3, 
                        seq_300bp='ATGCGTACGT'*30)
```

#### 方法2：批量预测
```python
region_example = pd.read_excel(get_test_data_path('region_example.xlsx'))
results = predictor.predict_region(seq_30bp=region_example['30bp'], 
                                  seq_300bp=region_example['300bp'])
```

#### 从 BLASTX 输出中识别 PRF 位点
```python
from FScanpy.utils import fscanr
blastx_output = pd.read_excel(get_test_data_path('blastx_example.xlsx'))
fscanr_result = fscanr(blastx_output, 
                      mismatch_cutoff=10,
                      evalue_cutoff=1e-5,
                      frameDist_cutoff=10)    
```
#### 提取分析窗口序列
```python
from FScanpy.utils import extract_prf_regions
prf_regions = extract_prf_regions(mrna_file=get_test_data_path('mrna_example.fasta'),
                                prf_data=fscanr_result)
```