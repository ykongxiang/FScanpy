import numpy as np
import pandas as pd
import itertools
from typing import List, Dict, Union

class SequenceFeatureExtractor:
    """DNA序列特征提取器"""
    
    def __init__(self):
        """初始化特征提取器"""
        self.bases = ['A', 'T', 'G', 'C']
        self.valid_bases = set('ATGCN')
        self.feature_names = self._get_feature_names()
    
    def _get_feature_names(self) -> List[str]:
        """
        返回特征名称列表
        
        Returns:
            List[str]: 特征名称列表
        """
        features = []
        
        # 基础特征
        features.extend(['A', 'T', 'G', 'C', 'N'])
        
        # 2-mer特征
        kmers_2 = [''.join(p) for p in itertools.product(self.bases, repeat=2)]
        kmers_2.extend(['N' + base for base in self.bases])
        kmers_2.extend([base + 'N' for base in self.bases])
        features.extend(kmers_2)
        
        # 3-mer特征
        kmers_3 = [''.join(p) for p in itertools.product(self.bases, repeat=3)]
        kmers_3.extend(['N' + ''.join(p) for p in itertools.product(self.bases, repeat=2)])
        kmers_3.extend([''.join(p) + 'N' for p in itertools.product(self.bases, repeat=2)])
        kmers_3.extend(['NN' + base for base in self.bases])
        kmers_3.extend([base + 'NN' for base in self.bases])
        kmers_3.extend(['NNN'])
        features.extend(kmers_3)
        
        # 位置特征
        for i in range(3):
            for base in self.valid_bases:
                features.append(f'window_{i}_{base}')
        
        return features
    
    def _preprocess_sequence(self, sequence):
        """序列预处理"""
        # 如果是Series，获取单个值
        if isinstance(sequence, pd.Series):
            sequence = sequence.iloc[0]
        
        # 处理空值或非字符串
        if pd.isna(sequence) or not isinstance(sequence, str):
            sequence = str(sequence)
        
        # 转换为大写并验证序列
        sequence = sequence.upper()
        if not all(base in self.valid_bases for base in sequence):
            sequence = ''.join(base if base in self.valid_bases else 'N' 
                             for base in sequence)
        
        return sequence
    
    def extract_features(self, sequence: Union[str, float]) -> List[float]:
        """
        从单个序列中提取特征
        
        Args:
            sequence: DNA序列
            
        Returns:
            List[float]: 特征值列表
        """
        try:
            sequence = self._preprocess_sequence(sequence)
            features: Dict[str, float] = {}
            
            # 1. 碱基组成特征
            seq_length = max(len(sequence), 1)
            for base in self.valid_bases:
                features[base] = sequence.count(base) / seq_length
            
            # 2. k-mer特征
            for k in [2, 3]:
                for i in range(len(sequence) - k + 1):
                    kmer = sequence[i:i+k]
                    features[kmer] = features.get(kmer, 0) + 1
            
            # 3. 位置特征
            window_size = max(len(sequence) // 3, 1)
            for i in range(3):
                start = i * window_size
                end = min((i + 1) * window_size, len(sequence))
                window = sequence[start:end]
                window_len = max(len(window), 1)
                
                for base in self.valid_bases:
                    features[f'window_{i}_{base}'] = window.count(base) / window_len
            
            return [features.get(name, 0) for name in self.feature_names]
            
        except Exception as e:
            raise ValueError(f"特征提取失败: {str(e)}")
    
    def extract_features_batch(self, sequences: List[Union[str, float]]) -> np.ndarray:
        """
        批量提取特征
        
        Args:
            sequences: DNA序列列表
            
        Returns:
            np.ndarray: 特征矩阵
        """
        try:
            return np.array([self.extract_features(seq) for seq in sequences])
        except Exception as e:
            raise ValueError(f"批量特征提取失败: {str(e)}")
    
    def predict_region_batch(self, data: pd.DataFrame, gb_threshold: float = 0.1) -> pd.DataFrame:
        """
        批量预测区域序列
        
        Args:
            data: DataFrame包含'30bp'和'300bp'列
            gb_threshold: GB模型概率阈值（默认为0.1）
            
        Returns:
            DataFrame: 包含预测结果的DataFrame
        """
        results = []
        for idx, row in data.iterrows():
            try:
                # 确保序列是字符串
                seq_30bp = str(row['30bp'])
                seq_300bp = str(row['300bp'])
                
                # 确保序列长度正确
                seq_30bp = seq_30bp.ljust(30, 'N')[:30]
                seq_300bp = seq_300bp.ljust(300, 'N')[:300]
                
                # 预测
                result = self.predict_region(seq_30bp, seq_300bp, gb_threshold)
                
                # 添加原始数据的其他列
                for col in data.columns:
                    if col not in ['30bp', '300bp']:
                        result[col] = row[col]
                        
                results.append(result)
                
            except Exception as e:
                print(f"Error processing sequence at index {idx}: {str(e)}")
                continue
                
        return pd.DataFrame(results)