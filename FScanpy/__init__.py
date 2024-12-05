from .predictor import PRFPredictor
import pandas as pd
import numpy as np
from typing import Union, List, Dict

__version__ = '0.2.0'
__author__ = ''
__email__ = ''

__all__ = ['PRFPredictor', 'predict_prf', '__version__', '__author__', '__email__']

def predict_prf(
    sequence: Union[str, List[str], None] = None,
    data: Union[pd.DataFrame, None] = None,
    window_size: int = 3,
    gb_threshold: float = 0.1,
    model_dir: str = None
) -> pd.DataFrame:
    """
    PRF site prediction function
    
    Args:
        sequence: single or multiple DNA sequences, used for sliding window prediction
        data: DataFrame data, must contain '30bp' and '300bp' columns, used for region prediction
        window_size: sliding window size (default is 3)
        gb_threshold: GB model probability threshold (default is 0.1)
        model_dir: model file directory path (optional)
    
    Returns:
        pandas.DataFrame: prediction results
        
    Examples:
        # 1. single sequence sliding window prediction
        >>> from FScanR-V2 import predict_prf
        >>> sequence = "ATGCGTACGT..."
        >>> results = predict_prf(sequence=sequence)
        
        # 2. multiple sequences sliding window prediction
        >>> sequences = ["ATGCGTACGT...", "GCTATAGCAT..."]
        >>> results = predict_prf(sequence=sequences)
        
        # 3. DataFrame region prediction
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     '30bp': ['ATGCGT...', 'GCTATAG...'],
        ...     '300bp': ['ATGCGT...', 'GCTATAG...']
        ... })
        >>> results = predict_prf(data=data)
    """
    predictor = PRFPredictor(model_dir=model_dir)
    
    # 验证输入参数
    if sequence is None and data is None:
        raise ValueError("One of sequence or data parameters must be provided")
    if sequence is not None and data is not None:
        raise ValueError("sequence and data parameters cannot be provided simultaneously")
    
    # 滑动窗口预测模式
    if sequence is not None:
        if isinstance(sequence, str):
            # 单条序列预测
            return predictor.predict_sliding_window(
                sequence, window_size, gb_threshold)
        elif isinstance(sequence, (list, tuple)):
            # 多条序列预测
            results = []
            for i, seq in enumerate(sequence, 1):
                try:
                    result = predictor.predict_sliding_window(
                        seq, window_size, gb_threshold)
                    result['Sequence_ID'] = f'seq_{i}'
                    results.append(result)
                except Exception as e:
                    print(f"Warning:sequence {i} prediction failed - {str(e)}")
            return pd.concat(results, ignore_index=True)
        else:
            raise ValueError("sequence parameter must be a string or a sequence list")
    
    # 区域化预测模式
    else:
        if not isinstance(data, pd.DataFrame):
            raise ValueError("data parameter must be a pandas DataFrame type")
        
        if '30bp' not in data.columns or '300bp' not in data.columns:
            raise ValueError("DataFrame must contain '30bp' and '300bp' columns")
        
        results = []
        for i, row in data.iterrows():
            try:
                result = predictor.predict_region(
                    row['30bp'], row['300bp'], gb_threshold)
                result['Index'] = i
                result['30bp'] = row['30bp']
                result['300bp'] = row['300bp']
                for col in data.columns:
                    if col not in ['30bp', '300bp']:
                        result[col] = row[col]
                results.append(result)
            except Exception as e:
                print(f"Warning: sequence prediction at index {i} failed - {str(e)}")
                result = {
                    'Index': i,
                    'GB_Probability': 0.0,
                    'CNN_Probability': 0.0,
                    'Voting_Probability': 0.0,
                    '30bp': row['30bp'],
                    '300bp': row['300bp']
                }
                for col in data.columns:
                    if col not in ['30bp', '300bp']:
                        result[col] = row[col]
                results.append(result)
        
        return pd.DataFrame(results)