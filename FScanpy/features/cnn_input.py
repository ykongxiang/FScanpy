import numpy as np
from tensorflow.keras.preprocessing.sequence import pad_sequences
from typing import List, Union

class CNNInputProcessor:
    """CNN模型输入数据处理器"""
    
    def __init__(self, max_length: int = 300):
        self.max_length = max_length
        self.base_to_num = {'A': 1, 'T': 2, 'G': 3, 'C': 4, 'N': 0}
    
    def prepare_sequence(self, sequence: str) -> np.ndarray:
        """
        处理单个序列
        
        Args:
            sequence: DNA序列
            
        Returns:
            np.ndarray: 处理后的序列数组
        """
        try:
            # 序列验证和预处理
            sequence = str(sequence).upper()
            valid_bases = set('ATGCN')
            sequence = ''.join(base if base in valid_bases else 'N' for base in sequence)
            
            # 序列转换为数字
            seq_numeric = [self.base_to_num.get(base, 0) for base in sequence]
            
            # 填充序列
            if len(seq_numeric) < self.max_length:
                seq_numeric.extend([0] * (self.max_length - len(seq_numeric)))
            elif len(seq_numeric) > self.max_length:
                seq_numeric = seq_numeric[:self.max_length]
            
            return np.array(seq_numeric).reshape(1, self.max_length, 1)
            
        except Exception as e:
            raise ValueError(f"序列处理失败: {str(e)}")
    
    def prepare_sequences(self, sequences: List[str]) -> np.ndarray:
        """
        批量处理序列
        
        Args:
            sequences: DNA序列列表
            
        Returns:
            np.ndarray: 处理后的序列数组
        """
        try:
            # 序列验证和预处理
            sequences = [str(seq).upper() for seq in sequences]
            valid_bases = set('ATGCN')
            sequences = [
                ''.join(base if base in valid_bases else 'N' for base in seq)
                for seq in sequences
            ]
            
            # 转换所有序列
            numeric_seqs = [
                [self.base_to_num.get(base, 0) for base in seq]
                for seq in sequences
            ]
            
            # 填充序列
            padded_seqs = pad_sequences(
                numeric_seqs, 
                maxlen=self.max_length,
                padding='post',
                truncating='post'
            )
            
            return np.array(padded_seqs).reshape(-1, self.max_length, 1)
            
        except Exception as e:
            raise ValueError(f"批量序列处理失败: {str(e)}")