import numpy as np
import pandas as pd
from typing import Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq

def fscanr(blastx_output: pd.DataFrame, 
           mismatch_cutoff: float = 10,
           evalue_cutoff: float = 1e-5,
           frameDist_cutoff: float = 10) -> pd.DataFrame:
    """
    identify PRF sites from BLASTX output
    
    Args:
        blastx_output: BLASTX output DataFrame
        mismatch_cutoff: mismatch threshold
        evalue_cutoff: E-value threshold 
        frameDist_cutoff: frame distance threshold
        
    Returns:
        pd.DataFrame: DataFrame containing PRF site information
    """
    blastx = blastx_output.copy()
    
    blastx.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", 
                     "gapopen", "qstart", "qend", "sstart", "send", 
                     "evalue", "bitscore", "qframe", "sframe"]
    
    blastx = blastx[
        (blastx['evalue'] <= evalue_cutoff) & 
        (blastx['mismatch'] <= mismatch_cutoff)
    ].dropna()

    freq = blastx['qseqid'].value_counts()
    multi_hits = freq[freq > 1].index
    blastx = blastx[blastx['qseqid'].isin(multi_hits)]
    
    blastx = blastx.sort_values(['qseqid', 'sseqid', 'qstart'])
    
    prf_list = []
    for i in range(1, len(blastx)):
        curr = blastx.iloc[i]
        prev = blastx.iloc[i-1]
        
        if (curr['qseqid'] == prev['qseqid'] and 
            curr['sseqid'] == prev['sseqid'] and
            curr['qframe'] != prev['qframe'] and 
            curr['qframe'] * prev['qframe'] > 0):
            
            if curr['qframe'] > 0 and prev['qframe'] > 0:
                frame_start = prev['qend'] 
                frame_end = curr['qstart']
                pep_start = prev['send']
                pep_end = curr['sstart']
                strand = "+"
            elif curr['qframe'] < 0 and prev['qframe'] < 0:
                frame_start = prev['qstart']
                frame_end = curr['qend'] 
                pep_start = curr['send']
                pep_end = prev['sstart']
                strand = "-"
            else:
                continue
                
            q_dist = frame_end - frame_start - 1
            s_dist = pep_end - pep_start
            fs_type = q_dist + (1 - s_dist) * 3
            
            if (abs(q_dist) <= frameDist_cutoff and 
                abs(s_dist) <= frameDist_cutoff // 3 and
                -3 < fs_type < 3):
                
                prf_list.append({
                    'DNA_seqid': curr['qseqid'],
                    'FS_start': frame_start,
                    'FS_end': frame_end,
                    'Pep_seqid': curr['sseqid'],
                    'Pep_FS_start': prev['send'] + 1,
                    'Pep_FS_end': curr['sstart'],
                    'FS_type': fs_type,
                    'Strand': strand
                })
    
    if not prf_list:
        print("No PRF events detected!")
        return pd.DataFrame()
        
    prf = pd.DataFrame(prf_list)
    
    for col in ['DNA_seqid', 'Pep_seqid']:
        for pos in ['FS_start', 'FS_end']:
            loci = prf[col] + '_' + prf[pos].astype(str)
            prf = prf[~loci.duplicated()]
            
    return prf

def extract_prf_regions(mrna_file: str, 
                       prf_data: pd.DataFrame) -> pd.DataFrame:
    """
    extract PRF site sequences from mRNA sequences
    
    Args:
        mrna_file: mRNA sequence file path (FASTA format)
        prf_data: FScanR output PRF site data
        
    Returns:
        pd.DataFrame: DataFrame containing 30bp and 300bp sequences
    """
    mrna_dict = {rec.id: str(rec.seq) 
                 for rec in SeqIO.parse(mrna_file, "fasta")}
    
    results = []
    for _, row in prf_data.iterrows():
        seq_id = row['DNA_seqid']
        if seq_id not in mrna_dict:
            print(f"Warning: {seq_id} not found in mRNA file")
            continue
            
        sequence = mrna_dict[seq_id]
        strand = row['Strand']
        fs_start = int(row['FS_start'])
        
        try:
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            
            period_seq, full_seq = extract_window_sequences(sequence, fs_start)
            
            results.append({
                'DNA_seqid': seq_id,
                'FS_start': fs_start,
                'FS_end': int(row['FS_end']),
                'Strand': strand,
                '30bp': period_seq,
                '300bp': full_seq,
                'FS_type': row['FS_type']
            })
            
        except Exception as e:
            print(f"Error processing {seq_id}: {str(e)}")
            continue
            
    return pd.DataFrame(results)

def extract_window_sequences(seq: str, position: int) -> Tuple[Optional[str], Optional[str]]:
    """
    extract analysis window sequences at specified position
    
    Args:
        seq: input DNA sequence
        position: current analysis position (FS_start)
    
    Returns:
        Tuple[str, str]: (30bp序列, 300bp序列)
    """
    frame_position = position - (position % 3)

    start_30 = frame_position - 10
    end_30 = frame_position + 20
    extend_each_side = 135  # (300 - 30) // 2
    start_300 = start_30 - extend_each_side
    end_300 = end_30 + extend_each_side

    seq_30 = _extract_and_pad(seq, start_30, end_30, 30)
    seq_300 = _extract_and_pad(seq, start_300, end_300, 300)
    
    return seq_30, seq_300

def _extract_and_pad(seq: str, start: int, end: int, target_length: int) -> str:
    """extract sequence and pad with N"""
    if start < 0:
        prefix = 'N' * abs(start)
        extracted = prefix + seq[:end]
    elif end > len(seq):
        suffix = 'N' * (end - len(seq))
        extracted = seq[start:] + suffix
    else:
        extracted = seq[start:end]
    
    if len(extracted) < target_length:
        extracted += 'N' * (target_length - len(extracted))
    
    return extracted

def prepare_cnn_input(sequence: str) -> np.ndarray:
    """prepare CNN model input"""
    base_to_num = {'A': 1, 'T': 2, 'G': 3, 'C': 4, 'N': 0}
    seq_numeric = [base_to_num.get(base, 0) for base in sequence.upper()]
    return np.array(seq_numeric).reshape(1, len(sequence), 1) 