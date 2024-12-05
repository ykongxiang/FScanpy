# FScanpy

<div align="right">
  <a href="README.md">English</a> |
  <b>中文</b>
</div>

### 介绍
FScanpy 是一个用于预测 mRNA 序列中程序化核糖体框架移位 (PRF) 位点的 Python 包。它集成了机器学习模型（Gradient Boosting 和 CNN），以提供准确的 PRF 预测。

### 安装

#### 1. 使用pip
```bash
pip install FScanpy
```

#### 2. 从GitHub克隆
```bash
git clone https://github.com/
cd your_project_directory
pip install -e .
```

### 测试用法
如果你想测试这个包，可以使用以下代码：

```python
# 基础预测
from FScanpy.data import get_test_data_path, list_test_data
from FScanpy import PRFPredictor

predictor = PRFPredictor()
predictor.predict_full(sequence='ATGCGTACGTATGCGTACGTATGCGTACGT',
                      window_size=3,
                      gb_threshold=0.1)

# 区域预测
seq_30bp = 'ATGCGTACGT' * 3  
seq_300bp = 'ATGCGTACGT' * 30
result = predictor.predict_region(seq_30bp=seq_30bp, seq_300bp=seq_300bp)
print(result)

# 批量预测
region_example = pd.read_excel(get_test_data_path('region_example.xlsx'))
results = predictor.predict_region(seq_30bp=region_example['30bp'], 
                                 seq_300bp=region_example['300bp'])

# BLASTX分析
from FScanpy.utils import fscanr
blastx_output = pd.read_excel(get_test_data_path('blastx_example.xlsx'))
fscanr_result = fscanr(blastx_output, 
                      mismatch_cutoff=10,
                      evalue_cutoff=1e-5,
                      frameDist_cutoff=10)

# 提取PRF区域
from FScanpy.utils import extract_prf_regions
prf_regions = extract_prf_regions(mrna_file=get_test_data_path('mrna_example.fasta'),
                                prf_data=fscanr_result)
```

### 测试数据
测试数据可以在 `FScanpy/data/test_data` 中找到：
```python
from FScanpy.data import get_test_data_path, list_test_data

blastx_file = get_test_data_path('blastx_example.xlsx')
mrna_file = get_test_data_path('mrna_example.fasta')
region_example = get_test_data_path('region_example.xlsx')
``` 