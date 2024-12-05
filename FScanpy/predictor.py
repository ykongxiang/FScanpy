import os
import pickle
import numpy as np
import pandas as pd
from tensorflow.keras.models import load_model
from .features.sequence import SequenceFeatureExtractor
from .features.cnn_input import CNNInputProcessor
from .utils import extract_window_sequences


class PRFPredictor:

    def __init__(self, model_dir=None):

        if model_dir is None:
            model_dir = os.path.join(os.path.dirname(__file__), 'pretrained')
        
        try:
            # 加载模型
            self.gb_model = self._load_pickle(os.path.join(model_dir, 'GradientBoosting_all.pkl'))
            self.cnn_model = load_model(os.path.join(model_dir, 'BiLSTM-CNN_all.keras'))
            self.voting_model = self._load_pickle(os.path.join(model_dir, 'Voting_all.pkl'))
            
            # 初始化特征提取器
            self.feature_extractor = SequenceFeatureExtractor()
            self.cnn_processor = CNNInputProcessor()
            
        except FileNotFoundError as e:
            raise FileNotFoundError(f"can't find model file: {str(e)}")
        except Exception as e:
            raise Exception(f"load model error: {str(e)}")
    
    def _load_pickle(self, path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    
    def predict_single_position(self, fs_period, full_seq, gb_threshold=0.1):
        '''
        Args:
            fs_period: 30bp sequence
            full_seq: 300bp sequence
            gb_threshold: GB model probability threshold (default is 0.1)
        Returns:
            dict: dictionary containing prediction probabilities
        '''
        try:
            # GB模型预测
            gb_features = self.feature_extractor.extract_features(fs_period)
            gb_prob = self.gb_model.predict_proba([gb_features])[0][1]
            
            # 如果GB概率低于阈值，直接返回
            if gb_prob < gb_threshold:
                return {
                    'GB_Probability': gb_prob,
                    'CNN_Probability': 0.0,
                    'Voting_Probability': 0.0
                }
            
            # CNN模型预测
            cnn_input = self.cnn_processor.prepare_sequence(full_seq)
            cnn_prob = self.cnn_model.predict(cnn_input, verbose=0)[0][0]
            
            # 投票模型预测
            voting_input = np.array([[gb_prob, cnn_prob]])
            voting_prob = self.voting_model.predict_proba(voting_input)[0][1]
            
            return {
                'GB_Probability': gb_prob,
                'CNN_Probability': cnn_prob,
                'Voting_Probability': voting_prob
            }
            
        except Exception as e:
            raise Exception(f"predict process error: {str(e)}")
    
    def predict_full(self, sequence, window_size=3, gb_threshold=0.1):
        '''
        predict whole sequence by sliding window
        
        Args:
            sequence: input DNA sequence
            window_size: sliding window size (default is 3)
            gb_threshold: GB model probability threshold (default is 0.1)

        Returns:
            DataFrame: containing prediction results for each position
        '''
        if window_size < 1:
            raise ValueError("window size must be greater than or equal to 1")
        if gb_threshold < 0:
            raise ValueError("GB threshold must be greater than or equal to 0")
        
        results = []
        
        try:
            for pos in range(0, len(sequence) - 2, window_size):
                fs_period, full_seq = extract_window_sequences(sequence, pos)
                
                if fs_period is None or full_seq is None:
                    continue
                
                pred = self.predict_single_position(fs_period, full_seq, gb_threshold)
                pred.update({
                    'Position': pos,
                    'Codon': sequence[pos:pos+3],
                    '30bp': fs_period,
                    '300bp': full_seq
                })
                results.append(pred)
            
            return pd.DataFrame(results)
            
        except Exception as e:
            raise Exception(f"sliding window predict error: {str(e)}")
    
    def predict_region(self, seq_30bp, seq_300bp, gb_threshold=0.1):
        '''
        predict region sequence
        
        Args:
            seq_30bp: 30bp sequence
            seq_300bp: 300bp sequence
            gb_threshold: GB model probability threshold (default is 0.1)
            
        Returns:
            dict: dictionary containing prediction probabilities
        ''' 
        try:
            return self.predict_single_position(seq_30bp, seq_300bp, gb_threshold)
            
        except Exception as e:
            raise Exception(f"region predict error: {str(e)}")