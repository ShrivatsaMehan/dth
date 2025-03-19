#!/usr/bin/env python3
import os
import subprocess
import argparse
import logging
import tempfile
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import re
import time
import shutil
import uuid
import json

# Set up logging with timestamp
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("msa_pipeline.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define 4-fold degenerate codons (where 3rd position can be any nucleotide)
FOURFOLD_CODONS = {
    'GC': True,  # Alanine (GCN)
    'CG': True,  # Arginine (CGN)
    'GG': True,  # Glycine (GGN)
    'CT': True,  # Leucine (CTN)
    'CC': True,  # Proline (CCN)
    'TC': True,  # Serine (TCN)
    'AC': True,  # Threonine (ACN)
    'GT': True   # Valine (GTN)
}

# Define genetic code table mapping codons to amino acids (including ambiguous codons)
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# Define amino acid similarity groups for better conservation scoring
AA_SIMILARITY_GROUPS = [
    {'G', 'A', 'S', 'T'},           # Small/polar
    {'C'},                          # Cysteine
    {'V', 'I', 'L', 'M'},           # Hydrophobic
    {'P'},                          # Proline
    {'F', 'Y', 'W'},                # Aromatic
    {'H', 'K', 'R'},                # Positive charge
    {'D', 'E', 'N', 'Q'},           # Negative charge/polar
]

def find_cds_files(directory='.'):
    """Find all .cds files in the given directory."""
    try:
        cds_files = [f for f in os.listdir(directory) if f.endswith('.cds')]
        logger.info(f"找到 {len(cds_files)} 个 .cds 文件在 {directory} 中")
        return sorted(cds_files)  # Sort for consistent processing order
    except Exception as e:
        logger.error(f"在 {directory} 中查找 CDS 文件时出错: {str(e)}")
        return []

def extract_species_name(header):
    """
    Extract species name from FASTA header using multiple methods.
    
    Args:
        header: FASTA header string
    
    Returns:
        Extracted species name
    """
    logger.debug(f"从头部提取物种名: {header}")
    
    # Method 1: Try to split by space and take second part
    try:
        if ' ' in header:
            parts = header.split(' ', 1)
            species = parts[1]
            logger.debug(f"方法 1 提取: {species}")
            return species
    except Exception as e:
        logger.debug(f"方法 1 失败: {str(e)}")
    
    # Method 2: Try regex to find pattern after '>' and before '/' or whitespace
    try:
        match = re.search(r'>?(.*?)(?:$|/|\s)', header)
        if match:
            species = match.group(1)
            logger.debug(f"方法 2 提取: {species}")
            return species
    except Exception as e:
        logger.debug(f"方法 2 失败: {str(e)}")
    
    # Method 3: Just use the whole header as species name
    logger.debug(f"使用完整的头部作为物种名")
    return header.strip('>')

def check_coding_sequence(sequence):
    """Check if sequence looks like a valid coding sequence.
    
    Returns:
        tuple: (is_valid, message)
    """
    if not sequence:
        return False, "空序列"
        
    if len(sequence) % 3 != 0:
        return False, f"序列长度 {len(sequence)} 不是 3 的倍数"
    
    # Check for stop codons in the middle
    for i in range(0, len(sequence) - 3, 3):
        codon = sequence[i:i+3].upper()
        if codon in ['TAA', 'TAG', 'TGA']:
            return False, f"序列中间发现终止密码子 {codon} 在位置 {i}"
    
    # Check for invalid characters
    valid_chars = set('ACGTRYKMSWBDHVN-')
    invalid_chars = set(sequence.upper()) - valid_chars
    if invalid_chars:
        return False, f"序列包含无效字符: {', '.join(invalid_chars)}"
    
    return True, "有效的编码序列"

def parse_fasta(file_path, duplicate_strategy='longest'):
    """
    Parse FASTA file and extract sequences with species info.
    
    Args:
        file_path: Path to FASTA file
        duplicate_strategy: How to handle duplicate species:
            'longest' - Keep the longest sequence
            'first' - Keep the first sequence encountered
            'rename' - Keep all sequences but append a suffix to duplicates
            'alignment_quality' - Will be handled separately in process_cds_file
    
    Returns:
        Dictionary mapping species to sequences
    """
    logger.debug(f"解析 FASTA 文件: {file_path}, 策略: {duplicate_strategy}")
    
    # If strategy is 'alignment_quality', we'll collect all sequences for now
    # and handle duplicates later after alignment
    if duplicate_strategy == 'alignment_quality':
        return parse_fasta_with_duplicates(file_path)
    
    species_seqs = {}
    species_counts = defaultdict(int)
    duplicate_found = False
    original_ids = {}  # Store original IDs for reference
    
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        logger.debug(f"在 {file_path} 中找到 {len(records)} 个序列")
        
        for record in records:
            # Extract species from header
            species = extract_species_name(record.description)
            sequence = str(record.seq)
            
            # Check if the sequence is a valid coding sequence
            is_valid, message = check_coding_sequence(sequence)
            if not is_valid:
                logger.warning(f"序列 {record.id} ({species}) 可能不是有效的编码序列: {message}")
            
            logger.debug(f"处理序列: ID={record.id}, 物种={species}, 长度={len(sequence)}")
            
            # Check if this species already exists
            if species in species_seqs:
                duplicate_found = True
                species_counts[species] += 1
                
                if duplicate_strategy == 'longest':
                    # Keep the longest sequence
                    if len(sequence) > len(species_seqs[species]):
                        logger.debug(f"用更长的序列替换 {species} 的较短序列 (长度: {len(species_seqs[species])} → {len(sequence)})")
                        species_seqs[species] = sequence
                        original_ids[species] = record.id
                
                elif duplicate_strategy == 'first':
                    # Keep the first sequence, ignore this one
                    logger.debug(f"忽略 {species} 的重复序列")
                    continue
                    
                elif duplicate_strategy == 'rename':
                    # Rename the duplicate by adding a suffix
                    new_species = f"{species}_{species_counts[species]}"
                    logger.debug(f"将重复的物种从 {species} 重命名为 {new_species}")
                    species_seqs[new_species] = sequence
                    original_ids[new_species] = record.id
            else:
                # First time seeing this species
                species_seqs[species] = sequence
                species_counts[species] = 1
                original_ids[species] = record.id
        
        if duplicate_found:
            logger.warning(f"在 {file_path} 中发现重复物种. 使用 '{duplicate_strategy}' 策略处理.")
            logger.info(f"{file_path} 中的物种计数: {dict(species_counts)}")
        
        return species_seqs
    except Exception as e:
        logger.error(f"解析 FASTA 文件 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def parse_fasta_with_duplicates(file_path):
    """
    Parse FASTA file preserving all duplicates with original IDs.
    Used for alignment_quality strategy.
    
    Returns:
        Dictionary mapping species to a dict of {id: sequence}
    """
    species_to_seqs = defaultdict(dict)
    
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        logger.debug(f"在 {file_path} 中找到 {len(records)} 个序列")
        
        for record in records:
            species = extract_species_name(record.description)
            record_id = record.id
            sequence = str(record.seq)
            
            # Check if the sequence is a valid coding sequence
            is_valid, message = check_coding_sequence(sequence)
            if not is_valid:
                logger.warning(f"序列 {record_id} ({species}) 可能不是有效的编码序列: {message}")
            
            # Store with original ID to track them
            species_to_seqs[species][record_id] = sequence
        
        # Log duplicates found
        duplicates = {species: len(ids) for species, ids in species_to_seqs.items() if len(ids) > 1}
        if duplicates:
            logger.info(f"在 {file_path} 中发现重复物种: {duplicates}")
            
        return species_to_seqs
    except Exception as e:
        logger.error(f"解析 FASTA 文件 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def format_sequences_for_codon_alignment(sequences, output_file):
    """
    Prepare sequences for alignment by ensuring they are valid coding sequences.
    
    Args:
        sequences: Dict of sequences or dict of {id: sequence}
        output_file: Path to write formatted sequences
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(output_file, 'w') as f:
            # Handle different input structures
            if isinstance(next(iter(sequences.values())), dict):
                # Format: {species: {id: sequence}}
                for species, id_to_seq in sequences.items():
                    for seq_id, seq in id_to_seq.items():
                        # Ensure length is multiple of 3
                        seq_len = len(seq)
                        if seq_len % 3 != 0:
                            padding = 'N' * (3 - seq_len % 3)
                            seq = seq + padding
                            logger.debug(f"序列 {seq_id} 长度不是3的倍数，添加 {len(padding)} 个N")
                        
                        f.write(f">{seq_id}\n{seq}\n")
            else:
                # Format: {species/id: sequence}
                for seq_id, seq in sequences.items():
                    # Ensure length is multiple of 3
                    seq_len = len(seq)
                    if seq_len % 3 != 0:
                        padding = 'N' * (3 - seq_len % 3)
                        seq = seq + padding
                        logger.debug(f"序列 {seq_id} 长度不是3的倍数，添加 {len(padding)} 个N")
                    
                    f.write(f">{seq_id}\n{seq}\n")
                    
        logger.debug(f"成功写入格式化序列到 {output_file}")
        return True
    except Exception as e:
        logger.error(f"格式化序列时发生错误: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def run_prank(input_file, output_prefix, prank_path, codon_aware=True, f=0.2, 
             gaprate=None, gapext=None, use_logs=False, penalize_terminal_gaps=False):
    """Run PRANK for multiple sequence alignment with optimization parameters.
    
    Args:
        input_file: Path to input file
        output_prefix: Prefix for output files
        prank_path: Path to PRANK executable
        codon_aware: Whether to use codon-aware alignment
        f: Parameter controlling insertion opening probability (default: 0.2)
        gaprate: Gap opening rate (if None, use PRANK defaults based on data type)
        gapext: Gap extension probability (if None, use PRANK defaults based on data type)
        use_logs: Use logarithm calculations for large datasets
        penalize_terminal_gaps: Whether to penalize terminal gaps
    """
    start_time = time.time()
    
    # Ensure prank_path exists and is executable
    if not os.path.exists(prank_path):
        logger.error(f"PRANK 可执行文件未在此路径找到: {prank_path}")
        return False, "可执行文件未找到"
    
    if not os.access(prank_path, os.X_OK):
        logger.error(f"PRANK 可执行文件没有执行权限: {prank_path}")
        return False, "可执行文件没有执行权限"
    
    cmd = [prank_path, '-d=' + input_file, '-o=' + output_prefix]
    
    # Add codon-aware alignment if requested
    if codon_aware:
        cmd.append('-codon')
    
    # Add optimization parameters
    cmd.append(f'-f={f}')  # Controls insertion opening probability
    
    # Set gap parameters only if specified, otherwise use PRANK defaults
    if gaprate is not None:
        cmd.append(f'-gaprate={gaprate}')
    
    if gapext is not None:
        cmd.append(f'-gapext={gapext}')
    
    # Add additional optimization parameters
    cmd.append('+F')  # Force insertions to always be skipped
    cmd.append('-prunetree')  # Prune guide tree branches with no sequence data
    cmd.append('-shortnames')  # Truncate names at first space for better processing
    
    if use_logs:
        cmd.append('-uselogs')  # Use logarithm calculations for large datasets
        
    if penalize_terminal_gaps:
        cmd.append('-termgap')  # Penalize terminal gaps normally
    
    # Always run once for consistency and speed
    cmd.append('-once')
    
    logger.info(f"运行 PRANK: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"PRANK 错误 (返回码 {result.returncode}): {result.stderr}")
            return False, f"返回码 {result.returncode}: {result.stderr}"
            
        execution_time = time.time() - start_time
        logger.info(f"PRANK 成功完成, 耗时 {execution_time:.1f} 秒")
        
        # Try to find alignment output with various extensions
        possible_outputs = [
            f"{output_prefix}.best.fas",
            f"{output_prefix}.fas",
            f"{output_prefix}.1.fas",
            f"{output_prefix}.aln"
        ]
        
        found_output = None
        for output_file in possible_outputs:
            if os.path.exists(output_file):
                found_output = output_file
                logger.info(f"找到PRANK输出文件: {found_output}")
                break
                
        if not found_output:
            # Check if PRANK created any files at all
            dir_path = os.path.dirname(output_prefix)
            base_name = os.path.basename(output_prefix)
            all_files = [f for f in os.listdir(dir_path) if f.startswith(base_name)]
            
            if all_files:
                # Find the most likely output file (e.g., largest .fas file)
                fas_files = [f for f in all_files if f.endswith('.fas')]
                if fas_files:
                    largest_fas = max(fas_files, key=lambda f: os.path.getsize(os.path.join(dir_path, f)))
                    found_output = os.path.join(dir_path, largest_fas)
                    logger.info(f"通过查找相关文件找到可能的输出: {found_output}")
            
        if not found_output:
            logger.error(f"PRANK 未创建预期的输出文件: {output_prefix}.best.fas")
            logger.error(f"PRANK 标准输出: {result.stdout}")
            logger.error(f"PRANK 标准错误: {result.stderr}")
            return False, "未创建输出文件"
            
        # Check if the file contains valid alignment
        valid, frame_info = check_alignment_validity(found_output)
        if valid:
            # Copy or rename the found output to the expected output name if different
            expected_output = f"{output_prefix}.best.fas"
            if found_output != expected_output:
                shutil.copy2(found_output, expected_output)
                logger.info(f"将找到的输出文件 {found_output} 复制到预期位置 {expected_output}")
            return True, expected_output
        else:
            logger.error(f"PRANK生成的文件 {found_output} 不包含有效的比对: {frame_info}")
            return False, f"输出文件不包含有效比对: {frame_info}"
            
    except Exception as e:
        logger.error(f"运行 PRANK 时发生异常: {str(e)}")
        logger.error(traceback.format_exc())
        return False, str(e)

def run_muscle(input_file, output_file, muscle_path, codon_aware=False):
    """Run MUSCLE for multiple sequence alignment with improved codon awareness.
    
    Args:
        input_file: Path to input file
        output_file: Path to output file
        muscle_path: Path to MUSCLE executable
        codon_aware: Whether to preserve codon structure
    """
    start_time = time.time()
    
    # Ensure muscle_path exists and is executable
    if not os.path.exists(muscle_path):
        logger.error(f"MUSCLE 可执行文件未在此路径找到: {muscle_path}")
        return False, "可执行文件未找到"
    
    if not os.access(muscle_path, os.X_OK):
        logger.error(f"MUSCLE 可执行文件没有执行权限: {muscle_path}")
        return False, "可执行文件没有执行权限"
    
    # If codon-aware alignment is requested, use a special pre-processing approach
    if codon_aware:
        logger.info("使用密码子感知方式运行MUSCLE")
        temp_dir = os.path.dirname(output_file)
        
        # 1. Convert sequences to aminoacid
        protein_input = os.path.join(temp_dir, os.path.basename(input_file).replace('.fasta', '.prot.fasta'))
        dna_to_protein_for_alignment(input_file, protein_input)
        
        # 2. Align proteins
        protein_output = os.path.join(temp_dir, os.path.basename(output_file).replace('.fasta', '.prot.aln'))
        protein_cmd = [muscle_path, '-in', protein_input, '-out', protein_output]
        
        logger.debug(f"运行蛋白质对齐: {' '.join(protein_cmd)}")
        protein_result = subprocess.run(protein_cmd, capture_output=True, text=True)
        
        if protein_result.returncode != 0:
            logger.error(f"MUSCLE蛋白质对齐错误: {protein_result.stderr}")
            return False, "蛋白质对齐失败"
            
        # 3. Map protein alignment back to DNA
        success = map_protein_alignment_to_dna(input_file, protein_output, output_file)
        if success:
            logger.info(f"成功完成基于蛋白质引导的密码子感知对齐，输出至 {output_file}")
            return True, output_file
        else:
            return False, "无法将蛋白质对齐映射回DNA"
    else:
        # Regular MUSCLE run for non-coding sequences
        cmd = [muscle_path, '-in', input_file, '-out', output_file]
        logger.info(f"运行 MUSCLE: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"MUSCLE 错误 (返回码 {result.returncode}): {result.stderr}")
                return False, f"返回码 {result.returncode}: {result.stderr}"
                
            execution_time = time.time() - start_time
            logger.info(f"MUSCLE 成功完成, 耗时 {execution_time:.1f} 秒")
            
            # Verify output
            valid, frame_info = check_alignment_validity(output_file)
            if valid:
                return True, output_file
            else:
                logger.error(f"MUSCLE生成的文件不包含有效的比对: {frame_info}")
                return False, f"输出文件不包含有效比对: {frame_info}"
        except Exception as e:
            logger.error(f"运行 MUSCLE 时发生异常: {str(e)}")
            logger.error(traceback.format_exc())
            return False, str(e)

def run_mafft(input_file, output_file, mafft_path, codon_aware=False):
    """Run MAFFT for multiple sequence alignment with improved codon awareness.
    
    Args:
        input_file: Path to input file
        output_file: Path to output file
        mafft_path: Path to MAFFT executable
        codon_aware: Whether to preserve codon structure
    """
    start_time = time.time()
    
    # Ensure mafft_path exists and is executable
    if not os.path.exists(mafft_path):
        logger.error(f"MAFFT 可执行文件未在此路径找到: {mafft_path}")
        return False, "可执行文件未找到"
    
    if not os.access(mafft_path, os.X_OK):
        logger.error(f"MAFFT 可执行文件没有执行权限: {mafft_path}")
        return False, "可执行文件没有执行权限"
    
    # If codon-aware alignment is requested, use protein guidance
    if codon_aware:
        logger.info("使用密码子感知方式运行MAFFT")
        temp_dir = os.path.dirname(output_file)
        
        # 1. Convert sequences to aminoacid
        protein_input = os.path.join(temp_dir, os.path.basename(input_file).replace('.fasta', '.prot.fasta'))
        dna_to_protein_for_alignment(input_file, protein_input)
        
        # 2. Align proteins
        protein_output = os.path.join(temp_dir, os.path.basename(output_file).replace('.fasta', '.prot.aln'))
        protein_cmd = [mafft_path, '--quiet', '--localpair', '--maxiterate', '1000', protein_input]
        
        logger.debug(f"运行蛋白质对齐: {' '.join(protein_cmd)}")
        
        try:
            with open(protein_output, 'w') as outfile:
                protein_result = subprocess.run(protein_cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
                
            if protein_result.returncode != 0:
                logger.error(f"MAFFT蛋白质对齐错误: {protein_result.stderr}")
                return False, "蛋白质对齐失败"
                
            # 3. Map protein alignment back to DNA
            success = map_protein_alignment_to_dna(input_file, protein_output, output_file)
            if success:
                logger.info(f"成功完成基于蛋白质引导的密码子感知对齐，输出至 {output_file}")
                return True, output_file
            else:
                return False, "无法将蛋白质对齐映射回DNA"
        except Exception as e:
            logger.error(f"运行MAFFT蛋白质对齐时出错: {str(e)}")
            logger.error(traceback.format_exc())
            return False, str(e)
    else:
        # Regular MAFFT run with general parameters for DNA
        cmd = [mafft_path]
        cmd.append('--auto')  # Auto algorithm selection
        cmd.extend(['--quiet', input_file])  # Redirect output
        
        logger.info(f"运行 MAFFT: {' '.join(cmd)}")
        
        try:
            with open(output_file, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"MAFFT 错误 (返回码 {result.returncode}): {result.stderr}")
                return False, f"返回码 {result.returncode}: {result.stderr}"
                
            execution_time = time.time() - start_time
            logger.info(f"MAFFT 成功完成, 耗时 {execution_time:.1f} 秒")
            
            # Verify output
            valid, frame_info = check_alignment_validity(output_file)
            if valid:
                return True, output_file
            else:
                logger.error(f"MAFFT生成的文件不包含有效的比对: {frame_info}")
                return False, f"输出文件不包含有效比对: {frame_info}"
        except Exception as e:
            logger.error(f"运行 MAFFT 时发生异常: {str(e)}")
            logger.error(traceback.format_exc())
            return False, str(e)

def dna_to_protein_for_alignment(dna_file, protein_file):
    """
    Convert DNA sequences to protein for alignment guidance.
    Handles codons and frame issues carefully.
    
    Args:
        dna_file: Path to DNA FASTA file
        protein_file: Path to output protein FASTA file
    
    Returns:
        bool: Success status
    """
    try:
        records = list(SeqIO.parse(dna_file, "fasta"))
        with open(protein_file, 'w') as out:
            for record in records:
                seq = str(record.seq).upper()
                
                # Ensure sequence length is multiple of 3
                remainder = len(seq) % 3
                if remainder > 0:
                    seq = seq + "N" * (3 - remainder)
                    logger.debug(f"序列 {record.id} 添加了 {3 - remainder} 个 N 以使其长度为3的倍数")
                
                # Translate to protein
                protein = ""
                for i in range(0, len(seq), 3):
                    codon = seq[i:i+3]
                    # Skip codons with gaps or ambiguous bases
                    if '-' in codon or 'N' in codon:
                        protein += 'X'
                    else:
                        try:
                            aa = GENETIC_CODE.get(codon, 'X')
                            protein += aa
                        except Exception:
                            protein += 'X'
                
                out.write(f">{record.id}\n{protein}\n")
        
        logger.debug(f"成功将 DNA 序列转换为蛋白质序列: {dna_file} -> {protein_file}")
        return True
    except Exception as e:
        logger.error(f"DNA转蛋白质时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def map_protein_alignment_to_dna(original_dna_file, protein_alignment_file, output_dna_alignment_file):
    """
    Map a protein alignment back to nucleotide sequences, preserving codon structure.
    
    Args:
        original_dna_file: Path to original unaligned DNA sequences
        protein_alignment_file: Path to aligned protein sequences
        output_dna_alignment_file: Path for output aligned DNA sequences
    
    Returns:
        bool: Success status
    """
    try:
        # Load original DNA sequences
        dna_records = SeqIO.to_dict(SeqIO.parse(original_dna_file, "fasta"))
        
        # Load aligned protein sequences
        aligned_proteins = {}
        for record in SeqIO.parse(protein_alignment_file, "fasta"):
            aligned_proteins[record.id] = str(record.seq)
        
        # Map each protein back to DNA
        aligned_dna = {}
        for seq_id, aligned_protein in aligned_proteins.items():
            if seq_id not in dna_records:
                logger.error(f"序列ID {seq_id} 未在原始DNA序列中找到")
                continue
                
            original_dna = str(dna_records[seq_id].seq).upper()
            
            # Ensure DNA length is multiple of 3
            remainder = len(original_dna) % 3
            if remainder > 0:
                original_dna = original_dna + "N" * (3 - remainder)
            
            # Build aligned DNA sequence
            aligned_dna_seq = ""
            dna_pos = 0
            
            for aa in aligned_protein:
                if aa == '-':
                    # Gap in protein = 3 gaps in DNA
                    aligned_dna_seq += "---"
                else:
                    # Regular amino acid = next codon from DNA
                    if dna_pos + 3 <= len(original_dna):
                        codon = original_dna[dna_pos:dna_pos+3]
                        aligned_dna_seq += codon
                        dna_pos += 3
                    else:
                        # Handle case where DNA sequence is shorter than expected
                        logger.warning(f"序列 {seq_id} DNA长度不足，在位置 {dna_pos} 添加 NNN")
                        aligned_dna_seq += "NNN"
            
            aligned_dna[seq_id] = aligned_dna_seq
        
        # Write aligned DNA sequences
        with open(output_dna_alignment_file, 'w') as out:
            for seq_id, seq in aligned_dna.items():
                out.write(f">{seq_id}\n{seq}\n")
        
        logger.info(f"成功将蛋白质比对映射回DNA: {len(aligned_dna)} 条序列")
        return True
    except Exception as e:
        logger.error(f"映射蛋白质比对到DNA时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def check_alignment_validity(alignment_file):
    """
    Check if the alignment file contains valid aligned sequences.
    
    Args:
        alignment_file: Path to alignment file
        
    Returns:
        tuple: (is_valid, frame_info)
    """
    try:
        # Check if file exists and is not empty
        if not os.path.exists(alignment_file) or os.path.getsize(alignment_file) == 0:
            logger.warning(f"比对文件不存在或为空: {alignment_file}")
            return False, "文件不存在或为空"
        
        # Parse the alignment file
        records = list(SeqIO.parse(alignment_file, "fasta"))
        
        # Check if any sequences were found
        if not records:
            logger.warning(f"比对文件没有包含序列: {alignment_file}")
            return False, "文件不包含序列"
        
        # Check for consistent sequence length (should be aligned)
        alignment_length = len(records[0].seq)
        if alignment_length == 0:
            logger.warning(f"比对序列长度为零: {alignment_file}")
            return False, "序列长度为零"
        
        # Check if the alignment length is a multiple of 3 (important for codon structure)
        is_codon_aligned = (alignment_length % 3 == 0)
        if not is_codon_aligned:
            logger.warning(f"比对长度 {alignment_length} 不是 3 的倍数，可能不保持密码子结构")
        
        # Check if all sequences have same length
        unequal_lengths = False
        for record in records:
            if len(record.seq) != alignment_length:
                logger.warning(f"比对序列长度不一致: {record.id} 的长度为 {len(record.seq)}, 而不是预期的 {alignment_length}")
                unequal_lengths = True
        
        if unequal_lengths:
            return False, "序列长度不一致"
        
        # Check for basic content
        # Ensure the file isn't just full of gaps or invalid characters
        valid_chars = set('ACGTRYKMSWBDHVN-')  # DNA + ambiguity codes + gap
        invalid_seqs = []
        
        for record in records:
            seq_str = str(record.seq).upper()
            # Check if sequence contains at least some valid DNA characters
            if not any(c in 'ACGT' for c in seq_str):
                logger.warning(f"序列 {record.id} 不包含有效的DNA字符")
                invalid_seqs.append(record.id)
                continue
                
            # Check if sequence has excessive invalid characters
            invalid_chars = set(seq_str) - valid_chars
            if invalid_chars and len(invalid_chars) > len(seq_str) * 0.1:  # >10% invalid chars
                logger.warning(f"序列 {record.id} 包含过多无效字符: {invalid_chars}")
                invalid_seqs.append(record.id)
        
        if invalid_seqs:
            return False, f"{len(invalid_seqs)} 条序列包含无效字符"
        
        # Calculate basic statistics
        total_gaps = sum(str(record.seq).count('-') for record in records)
        total_chars = sum(len(record.seq) for record in records)
        gap_percentage = (total_gaps / total_chars) * 100 if total_chars > 0 else 0
        
        # High gap percentage might indicate problems, but this is just informative
        if gap_percentage > 50:
            logger.warning(f"比对包含高比例的缺口: {gap_percentage:.1f}%")
        
        # Check for codon structure preservation by analyzing gap patterns
        if is_codon_aligned:
            codon_structure_preserved = True
            for record in records:
                seq_str = str(record.seq)
                # Check if gaps appear in multiples of 3
                for i in range(0, len(seq_str), 3):
                    if i+2 < len(seq_str):
                        codon = seq_str[i:i+3]
                        gap_count = codon.count('-')
                        if gap_count > 0 and gap_count < 3:
                            # Found a codon with some but not all positions as gaps
                            # This likely breaks codon structure
                            codon_structure_preserved = False
                            break
                            
                if not codon_structure_preserved:
                    break
                    
            if not codon_structure_preserved:
                logger.warning(f"比对可能不保持密码子结构（发现单/双碱基缺口）")
        
        frame_info = f"长度 {alignment_length}, 密码子对齐: {is_codon_aligned}, 缺口比例 {gap_percentage:.1f}%"
        logger.info(f"比对验证通过: {len(records)} 条序列, {frame_info}")
        
        return True, frame_info
        
    except Exception as e:
        logger.error(f"验证比对文件时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False, str(e)

def fix_alignment_frame(alignment_file, output_file=None):
    """
    Fix the alignment to ensure it maintains codon structure.
    
    Args:
        alignment_file: Path to alignment file to fix
        output_file: Path to write fixed alignment (if None, will use input path)
        
    Returns:
        bool: Success status
    """
    if output_file is None:
        output_file = alignment_file.replace('.fasta', '.framed.fasta')
        if output_file == alignment_file:
            output_file = alignment_file + '.framed'
    
    try:
        records = list(SeqIO.parse(alignment_file, "fasta"))
        if not records:
            logger.error(f"无法从 {alignment_file} 解析序列")
            return False
            
        alignment_length = len(records[0].seq)
        
        # Check if fix is needed
        if alignment_length % 3 == 0:
            logger.info(f"比对长度已经是3的倍数 ({alignment_length})，无需修复")
            return True
        
        # Calculate padding needed
        padding = 3 - (alignment_length % 3)
        logger.info(f"添加 {padding} 个缺口使比对长度为3的倍数")
        
        # Add padding to all sequences
        with open(output_file, 'w') as out:
            for record in records:
                fixed_seq = str(record.seq) + '-' * padding
                out.write(f">{record.id}\n{fixed_seq}\n")
        
        return True
    except Exception as e:
        logger.error(f"修复比对框架时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def calculate_similarity_score(aa1, aa2):
    """
    Calculate similarity score between two amino acids based on similarity groups.
    
    Args:
        aa1: First amino acid
        aa2: Second amino acid
        
    Returns:
        1.0 if identical, 0.5 if in same group, 0.0 if different groups
    """
    if aa1 == aa2:
        return 1.0
    
    for group in AA_SIMILARITY_GROUPS:
        if aa1 in group and aa2 in group:
            return 0.5
    
    return 0.0

def evaluate_alignment_quality(aligned_sequences, seq_id, is_protein=False):
    """
    Improved evaluation of alignment quality of a sequence compared to others.
    Lower score means better quality.
    
    Args:
        aligned_sequences: Dictionary of all aligned sequences
        seq_id: ID of the sequence to evaluate
        is_protein: Whether sequences are protein sequences
        
    Returns:
        Quality score (lower is better)
    """
    if seq_id not in aligned_sequences:
        logger.warning(f"序列ID {seq_id} 未在比对序列中找到")
        return float('inf')  # Return worst possible score if sequence not found
    
    sequence = aligned_sequences[seq_id]
    # Filter out the sequence being evaluated
    other_seqs = [seq for id, seq in aligned_sequences.items() if id != seq_id]
    
    if not other_seqs:
        logger.debug(f"序列 {seq_id} 没有其他序列可比较")
        return 0  # No other sequences to compare with
    
    # Initialize quality metrics with weights
    weights = {
        'gap_ratio': 0.3,           # Penalize gaps
        'conservation': 0.4,        # Reward conservation with other sequences
        'terminal_gaps': 0.1,       # Penalize terminal gaps less than internal gaps
        'consecutive_gaps': 0.2     # Penalize consecutive gaps (indels)
    }
    
    # Calculate gap metrics
    total_gaps = sequence.count('-')
    seq_length = len(sequence)
    gap_ratio = total_gaps / seq_length if seq_length > 0 else 1.0
    
    # Identify terminal gaps vs internal gaps
    left_term_gaps = 0
    for char in sequence:
        if char == '-':
            left_term_gaps += 1
        else:
            break
            
    right_term_gaps = 0
    for char in reversed(sequence):
        if char == '-':
            right_term_gaps += 1
        else:
            break
    
    internal_gaps = total_gaps - (left_term_gaps + right_term_gaps)
    terminal_gap_ratio = (left_term_gaps + right_term_gaps) / seq_length if seq_length > 0 else 0
    internal_gap_ratio = internal_gaps / seq_length if seq_length > 0 else 0
    
    # Count consecutive gap blocks (indels)
    in_gap = False
    gap_blocks = 0
    for char in sequence:
        if char == '-':
            if not in_gap:
                gap_blocks += 1
                in_gap = True
        else:
            in_gap = False
    
    consecutive_gap_score = gap_blocks / (seq_length / 10) if seq_length > 0 else 1.0
    # Normalize to 0-1 range (more blocks is worse but better than fewer large blocks)
    consecutive_gap_score = min(consecutive_gap_score, 1.0)
    
    # Calculate conservation with other sequences
    conservation_scores = []
    
    for pos in range(seq_length):
        # Skip positions where this sequence has a gap
        if sequence[pos] == '-':
            continue
        
        # Count matches with other sequences at this position
        match_scores = []
        
        for other_seq in other_seqs:
            if pos >= len(other_seq):
                continue
                
            if other_seq[pos] == '-':
                continue
                
            if is_protein:
                # For proteins, use similarity scoring
                match_scores.append(calculate_similarity_score(sequence[pos], other_seq[pos]))
            else:
                # For nucleotides, exact match only
                match_scores.append(1.0 if sequence[pos] == other_seq[pos] else 0.0)
        
        if match_scores:
            # Average conservation at this position
            conservation_scores.append(sum(match_scores) / len(match_scores))
    
    # Overall conservation score (higher is better, so invert for our scoring where lower is better)
    avg_conservation = 1.0 - (sum(conservation_scores) / len(conservation_scores) if conservation_scores else 0)
    
    # Combine all metrics (weighted)
    quality_score = (
        weights['gap_ratio'] * gap_ratio + 
        weights['conservation'] * avg_conservation +
        weights['terminal_gaps'] * terminal_gap_ratio +
        weights['consecutive_gaps'] * consecutive_gap_score
    )
    
    logger.debug(f"序列 {seq_id} 的比对质量得分: {quality_score:.4f} " + 
                f"(缺口比: {gap_ratio:.2f}, 保守性: {avg_conservation:.2f}, " +
                f"终端缺口: {terminal_gap_ratio:.2f}, 连续缺口块: {consecutive_gap_score:.2f})")
    
    return quality_score

def select_best_duplicate_sequences(species_to_seqs, aligned_file, is_protein=False):
    """
    Select the best sequence for each species based on alignment quality.
    
    Args:
        species_to_seqs: Dictionary mapping species to dict of {id: sequence}
        aligned_file: Path to alignment output file
        is_protein: Whether sequences are protein sequences
    
    Returns:
        Dictionary mapping species to their best sequence
    """
    logger.info(f"根据比对质量选择最佳重复序列 (蛋白质: {is_protein})")
    
    best_sequences = {}
    
    try:
        # Parse the aligned file
        aligned_seqs = {}
        for record in SeqIO.parse(aligned_file, "fasta"):
            aligned_seqs[record.id] = str(record.seq)
        
        # For each species, find the best sequence
        for species, id_to_seq in species_to_seqs.items():
            if len(id_to_seq) == 1:
                # Only one sequence, no need to evaluate
                seq_id = next(iter(id_to_seq))
                best_sequences[species] = id_to_seq[seq_id]
                logger.debug(f"物种 {species} 只有一个序列: {seq_id}")
                continue
                
            # Evaluate each duplicate sequence
            best_id = None
            best_score = float('inf')
            all_scores = {}
            
            for seq_id in id_to_seq.keys():
                if seq_id not in aligned_seqs:
                    logger.warning(f"序列ID {seq_id} 未在比对序列中找到，跳过评估")
                    continue
                    
                score = evaluate_alignment_quality(aligned_seqs, seq_id, is_protein=is_protein)
                all_scores[seq_id] = score
                
                if score < best_score:
                    best_score = score
                    best_id = seq_id
            
            if best_id:
                best_sequences[species] = id_to_seq[best_id]
                logger.info(f"为物种 {species} 选择了最佳序列 {best_id} (得分: {best_score:.4f})")
                
                # Log scores for all duplicates for this species
                for seq_id, score in all_scores.items():
                    logger.debug(f"  物种 {species} 序列 {seq_id} 得分: {score:.4f}" + 
                                (f" (已选择)" if seq_id == best_id else ""))
            else:
                # Fallback to first sequence if no best found
                first_id = next(iter(id_to_seq))
                best_sequences[species] = id_to_seq[first_id]
                logger.warning(f"无法为物种 {species} 找到最佳序列，使用第一个序列 {first_id}")
        
        return best_sequences
    except Exception as e:
        logger.error(f"选择最佳重复序列时出错: {str(e)}")
        logger.error(traceback.format_exc())
        
        # Fallback to first sequence for each species
        fallback = {}
        for species, seqs in species_to_seqs.items():
            fallback[species] = seqs[next(iter(seqs))]
            logger.warning(f"为物种 {species} 回退使用第一个序列")
        
        return fallback

def run_aligner(input_file, output_file, aligner_params):
    """
    Run the specified aligner to align sequences.
    
    Args:
        input_file: Path to input file
        output_file: Path or prefix for output file
        aligner_params: Dictionary of aligner parameters
        
    Returns:
        (success, output_file) tuple
    """
    aligner_type = aligner_params['aligner']
    
    if aligner_type == 'prank':
        output_prefix = output_file.replace('.best.fas', '')
        return run_prank(
            input_file, 
            output_prefix, 
            aligner_params['prank_path'],
            codon_aware=aligner_params['codon_aware'],
            f=aligner_params['f'],
            gaprate=aligner_params['gaprate'],
            gapext=aligner_params['gapext'],
            use_logs=aligner_params['use_logs'],
            penalize_terminal_gaps=aligner_params['penalize_terminal_gaps']
        )
    elif aligner_type == 'muscle':
        return run_muscle(
            input_file,
            output_file,
            aligner_params['muscle_path'],
            codon_aware=aligner_params['codon_aware']
        )
    elif aligner_type == 'mafft':
        return run_mafft(
            input_file,
            output_file,
            aligner_params['mafft_path'],
            codon_aware=aligner_params['codon_aware']
        )
    else:
        logger.error(f"未知的比对工具: {aligner_type}")
        return False, f"未知的比对工具: {aligner_type}"

def process_cds_file_with_alignment_quality(file_path, output_dir, aligner_params):
    """
    Process CDS file using alignment quality to select best duplicate sequences.
    
    Args:
        file_path: Path to CDS file
        output_dir: Output directory
        aligner_params: Parameters for alignment
    
    Returns:
        Gene info object with 4D sites
    """
    try:
        file_name = os.path.basename(file_path)
        gene_name = file_name.replace('.cds', '')
        temp_dir = os.path.join(output_dir, "temp", gene_name + "_" + str(uuid.uuid4())[:8])
        os.makedirs(temp_dir, exist_ok=True)
        
        logger.info(f"使用比对质量策略处理 {file_name}")
        
        # Step 1: Parse FASTA with duplicates preserved
        species_to_seqs = parse_fasta_with_duplicates(file_path)
        
        if not species_to_seqs:
            logger.error(f"无法从 {file_path} 解析序列")
            return None
        
        # Check if any duplicates exist
        has_duplicates = any(len(seqs) > 1 for seqs in species_to_seqs.values())
        
        if not has_duplicates:
            logger.info(f"{file_name} 没有重复物种，直接处理")
            # Create a regular dictionary with one sequence per species
            simplified_seqs = {species: next(iter(seqs.values())) for species, seqs in species_to_seqs.items()}
            return process_cds_file_standard(file_path, output_dir, aligner_params, simplified_seqs)
        
        # Step 2: Create a temporary FASTA file with all sequences including duplicates
        initial_fasta = os.path.join(temp_dir, f"{gene_name}_all_seqs.fasta")
        
        # Use codon-aware formatting for input sequences
        with open(initial_fasta, 'w') as f:
            for species, id_to_seq in species_to_seqs.items():
                for seq_id, sequence in id_to_seq.items():
                    # Ensure sequence length is multiple of 3
                    seq_len = len(sequence)
                    if seq_len % 3 != 0:
                        padding = 'N' * (3 - seq_len % 3)
                        sequence = sequence + padding
                        logger.debug(f"序列 {seq_id} 长度不是3的倍数，添加 {len(padding)} 个N")
                    f.write(f">{seq_id}\n{sequence}\n")
        
        # Step 3: Run initial alignment
        initial_output = os.path.join(temp_dir, f"{gene_name}_initial.best.fas")
        success, output_path = run_aligner(initial_fasta, initial_output, aligner_params)
        
        if not success:
            logger.error(f"初始比对 {file_name} 失败")
            # Fall back to longest sequence strategy
            logger.info(f"退回到使用最长序列策略")
            simplified_seqs = {}
            for species, id_to_seq in species_to_seqs.items():
                best_id = max(id_to_seq.items(), key=lambda x: len(x[1]))[0]
                simplified_seqs[species] = id_to_seq[best_id]
                logger.info(f"为物种 {species} 选择了最长序列 {best_id}")
            return process_cds_file_standard(file_path, output_dir, aligner_params, simplified_seqs)
        
        # Fix alignment frame if needed
        valid, frame_info = check_alignment_validity(output_path)
        if not valid or 'codon' in frame_info and '不保持' in frame_info:
            logger.warning(f"需要修复比对框架: {frame_info}")
            fixed_output = os.path.join(temp_dir, f"{gene_name}_initial_fixed.best.fas")
            fix_alignment_frame(output_path, fixed_output)
            output_path = fixed_output
        
        # Step 4: Select the best sequence for each species based on alignment quality
        best_sequences = select_best_duplicate_sequences(species_to_seqs, output_path)
        
        if not best_sequences:
            logger.error(f"无法选择最佳重复序列")
            return None
        
        # Step 5: Process with the selected best sequences
        return process_cds_file_standard(file_path, output_dir, aligner_params, best_sequences)
        
    except Exception as e:
        logger.error(f"处理 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None
    finally:
        # Clean up temporary directory
        try:
            if os.path.exists(temp_dir) and aligner_params['clean_temp']:
                shutil.rmtree(temp_dir)
                logger.debug(f"清理临时目录 {temp_dir}")
        except Exception as e:
            logger.warning(f"清理临时目录 {temp_dir} 时出错: {str(e)}")

def process_cds_file_standard(file_path, output_dir, aligner_params, sequences=None):
    """
    Process a single CDS file through alignment and 4D site extraction.
    
    Args:
        file_path: Path to CDS file
        output_dir: Output directory
        aligner_params: Parameters for alignment
        sequences: Optional pre-filtered sequences dictionary
    """
    try:
        file_name = os.path.basename(file_path)
        gene_name = file_name.replace('.cds', '')
        output_prefix = os.path.join(output_dir, "alignments", gene_name)
        
        # Create directory for alignments if it doesn't exist
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        
        # Use provided sequences or parse from file using traditional strategy
        if sequences is None:
            logger.debug(f"使用 {aligner_params['duplicate_strategy']} 策略解析序列")
            input_seqs = parse_fasta(file_path, aligner_params['duplicate_strategy'])
        else:
            logger.debug(f"使用预处理的序列集 (长度: {len(sequences)})")
            input_seqs = sequences
            
        if not input_seqs:
            logger.error(f"无法从 {file_path} 解析序列")
            return None
            
        logger.info(f"为 {file_name} 处理 {len(input_seqs)} 个序列")
        
        # Create temporary FASTA file with processed sequences
        temp_dir = os.path.join(output_dir, "temp")
        os.makedirs(temp_dir, exist_ok=True)
        temp_input = os.path.join(temp_dir, f"{gene_name}_processed.fasta")
        
        # Ensure sequences are properly formatted for codon-aware alignment
        with open(temp_input, 'w') as f:
            for species, seq in input_seqs.items():
                # Ensure sequence length is multiple of 3
                seq_len = len(seq)
                if seq_len % 3 != 0:
                    padding = 'N' * (3 - seq_len % 3)
                    seq = seq + padding
                    logger.debug(f"序列 {species} 长度不是3的倍数，添加 {len(padding)} 个N")
                f.write(f">{species}\n{seq}\n")
        
        # Align using the appropriate aligner
        aligned_output = f"{output_prefix}.best.fas"
        success, output_path = run_aligner(temp_input, aligned_output, aligner_params)
        
        if not success:
            logger.error(f"比对 {file_name} 失败")
            return None
        
        # Ensure we're using the correct path to the output file
        aligned_file = output_path if output_path else aligned_output
        
        # Parse aligned sequences
        if not os.path.exists(aligned_file):
            logger.error(f"比对输出文件 {aligned_file} 未找到")
            return None
        
        # Parse aligned file
        aligned_seqs = {}
        for record in SeqIO.parse(aligned_file, "fasta"):
            species = extract_species_name(record.description)
            aligned_seqs[species] = str(record.seq)
        
        logger.info(f"从 {aligned_file} 解析了 {len(aligned_seqs)} 个比对序列")
        
        # Fix alignment frame if needed to ensure codon structure
        aligned_length = len(next(iter(aligned_seqs.values())))
        if aligned_length % 3 != 0:
            logger.warning(f"比对长度 {aligned_length} 不是3的倍数，尝试修复")
            fixed_file = aligned_file.replace('.best.fas', '.framed.fas')
            if fix_alignment_frame(aligned_file, fixed_file):
                # Re-read the aligned sequences from fixed file
                aligned_seqs = {}
                for record in SeqIO.parse(fixed_file, "fasta"):
                    species = extract_species_name(record.description)
                    aligned_seqs[species] = str(record.seq)
                logger.info(f"使用修复后的比对 ({len(aligned_seqs)} 个序列)")
        
        # Identify and extract 4D sites
        fourfold_sites = identify_4d_sites(aligned_seqs)
        fourfold_seqs = extract_4d_sites(aligned_seqs, fourfold_sites)
        
        # Create gene info for debugging and tracking
        gene_info = {
            'name': gene_name,
            'species_count': len(fourfold_seqs),
            'site_count': len(fourfold_sites),
            'sequences': fourfold_seqs,
            'aligned_seqs': aligned_seqs  # Store full aligned sequences for amino acid translation
        }
        
        return gene_info
    except Exception as e:
        logger.error(f"处理 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def identify_4d_sites(aligned_seqs):
    """Identify 4-fold degenerate sites from codon-aligned sequences."""
    # Get alignment length and check if it's divisible by 3
    seq_values = list(aligned_seqs.values())
    if not seq_values:
        logger.warning("未提供序列给 identify_4d_sites")
        return []
    
    align_length = len(seq_values[0])
    logger.info(f"准备识别4D位点: 比对长度 {align_length}, 是否为3的倍数: {align_length % 3 == 0}")
    
    # Stop if alignment length is not a multiple of 3
    if align_length % 3 != 0:
        logger.warning(f"比对长度 {align_length} 不是 3 的倍数。无法识别4D位点。")
        return []
    
    # Count 4-fold sites
    fourfold_sites = []
    fourfold_candidates = 0
    gap_blocked = 0
    non_matching_prefix = 0
    
    # Check each codon position
    for i in range(0, align_length, 3):
        # Skip last incomplete codon if any
        if i+2 >= align_length:
            continue
        
        # Check if first two positions of the codon could form 4D site
        first_seq = seq_values[0]
        prefix = first_seq[i:i+2].upper()
        
        # Skip codons with gaps in any species (more strict check)
        has_gap = False
        for seq in seq_values:
            if i+2 >= len(seq) or '-' in seq[i:i+3]:
                has_gap = True
                break
        
        if has_gap:
            gap_blocked += 1
            continue
        
        # Is this a potential 4D codon?
        if prefix in FOURFOLD_CODONS:
            fourfold_candidates += 1
            
            # Check if all sequences have the same prefix at this position
            all_match = True
            for seq in seq_values[1:]:
                if seq[i:i+2].upper() != prefix:
                    all_match = False
                    break
            
            if all_match:
                fourfold_sites.append(i+2)  # Add the third position
            else:
                non_matching_prefix += 1
    
    logger.info(f"4D位点识别结果: 找到 {len(fourfold_sites)} 个4D位点, " +
               f"{fourfold_candidates} 个潜在4D密码子, " +
               f"{gap_blocked} 个密码子因缺口而跳过, " +
               f"{non_matching_prefix} 个密码子前缀不匹配")
    
    return fourfold_sites

def extract_4d_sites(aligned_seqs, fourfold_sites):
    """Extract 4-fold degenerate sites from aligned sequences."""
    if not fourfold_sites:
        logger.warning("没有 4D 位点可提取")
        return {}
        
    logger.debug(f"从 {len(aligned_seqs)} 个序列中提取 {len(fourfold_sites)} 个 4D 位点")
    fourfold_seqs = {}
    
    for species, seq in aligned_seqs.items():
        # Check if any 4D sites are out of bounds
        valid_sites = [pos for pos in fourfold_sites if pos < len(seq)]
        if len(valid_sites) != len(fourfold_sites):
            logger.warning(f"{species} 的一些 4D 位点超出范围: 序列长度 {len(seq)}, 最大位点 {max(fourfold_sites)}")
            
        fourfold_seq = ''.join([seq[pos] for pos in valid_sites])
        fourfold_seqs[species] = fourfold_seq
    
    return fourfold_seqs

def safe_translate_codon(codon):
    """
    Safely translate a codon to amino acid, handling ambiguous bases.
    
    Args:
        codon: Three-letter nucleotide codon
        
    Returns:
        Single letter amino acid or 'X' for unknown
    """
    # Convert to uppercase and remove whitespace
    codon = codon.upper().strip()
    
    # Check if codon exists in our lookup table
    if codon in GENETIC_CODE:
        return GENETIC_CODE[codon]
    
    # If codon contains N or other ambiguous nucleotides
    if 'N' in codon or any(base not in 'ACGT' for base in codon):
        return 'X'
        
    # This shouldn't happen with valid DNA sequences
    logger.warning(f"未知的密码子: '{codon}'")
    return 'X'

def translate_with_gaps(nucleotide_seq):
    """
    Translate a nucleotide sequence to amino acids, respecting gaps.
    Gaps (-) in the nucleotide sequence will be preserved as gaps in the amino acid sequence.
    """
    if not nucleotide_seq:
        return ""
    
    # 初始化padding变量为0
    needed_padding = 0
    
    # Ensure sequence length is multiple of 3 for proper translation
    original_len = len(nucleotide_seq)
    if original_len % 3 != 0:
        needed_padding = 3 - (original_len % 3)
        nucleotide_seq += 'N' * needed_padding
        logger.debug(f"序列长度 {original_len} 不是 3 的倍数，添加 {needed_padding} 个N进行翻译")
    
    amino_acids = []
    
    for i in range(0, len(nucleotide_seq), 3):
        codon = nucleotide_seq[i:i+3]
        
        if '-' in codon:
            # Count gaps in codon
            gap_count = codon.count('-')
            
            # Apply gap translation rules
            if gap_count == 3:
                # Complete gap codon
                amino_acids.append('-')
            elif gap_count > 0:
                # Partial gap codon - preserve frame by using X
                amino_acids.append('X')
            else:
                # This shouldn't happen
                amino_acids.append('X')
        else:
            try:
                amino_acid = safe_translate_codon(codon)
                amino_acids.append(amino_acid)
            except Exception as e:
                logger.warning(f"翻译密码子 '{codon}' 时出错: {str(e)}")
                amino_acids.append('X')  # Use X for unknown/problematic amino acid
    
    # Return only the portion corresponding to original sequence
    result = ''.join(amino_acids)
    if needed_padding > 0:
        # Trim any padding we might have added
        result = result[:-1] if needed_padding == 1 else result[:-2] if needed_padding == 2 else result
    
    return result

def process_cds_file(file_path, output_dir, aligner_params):
    """
    Process a single CDS file using the appropriate strategy.
    
    Args:
        file_path: Path to CDS file
        output_dir: Output directory
        aligner_params: Parameters for alignment
    
    Returns:
        Gene info object with 4D sites
    """
    if aligner_params['duplicate_strategy'] == 'alignment_quality':
        return process_cds_file_with_alignment_quality(file_path, output_dir, aligner_params)
    else:
        return process_cds_file_standard(file_path, output_dir, aligner_params)

def merge_sequences_by_species(all_gene_results, output_file, missing_species_strategy='gaps', min_coverage_pct=0):
    """
    Merge 4D sequences from all genes, organized by species.
    
    Args:
        all_gene_results: List of gene results with sequence data
        output_file: Base output file path (strategy will be appended)
        missing_species_strategy: Strategy for handling missing species
            'gaps' - Fill missing sequences with gaps (default)
            'exclude_species' - Exclude species below minimum coverage
            'exclude_genes' - Only use genes that have all species
        min_coverage_pct: Minimum percentage of genes a species must be present in
                          (only used with 'exclude_species' strategy)
    """
    logger.info(f"使用策略合并序列: {missing_species_strategy}")
    
    # Dictionary to store concatenated sequences for each species
    species_supergenes = defaultdict(str)
    
    # Dictionary to track which genes contribute to each species
    species_gene_coverage = defaultdict(list)
    
    # Keep track of all species encountered and valid genes
    all_species = set()
    valid_gene_results = [g for g in all_gene_results if g is not None]
    total_genes = len(valid_gene_results)
    
    logger.info(f"处理 {total_genes} 个有效基因，共 {len(all_gene_results)} 个总基因")
    
    # First pass: identify all species
    for gene_result in valid_gene_results:
        all_species.update(gene_result['sequences'].keys())
    
    logger.info(f"在所有基因中发现 {len(all_species)} 个不同的物种")
    
    # Calculate gene coverage for each species
    for gene_result in valid_gene_results:
        gene_name = gene_result['name']
        sequences = gene_result['sequences']
        
        for species in all_species:
            if species in sequences:
                species_gene_coverage[species].append(gene_name)
    
    # Apply filtering based on strategy
    included_species = set()
    included_genes = []
    
    if missing_species_strategy == 'exclude_species':
        # Only include species that meet the minimum coverage threshold
        min_genes_required = int(total_genes * min_coverage_pct / 100)
        for species, genes in species_gene_coverage.items():
            if len(genes) >= min_genes_required:
                included_species.add(species)
            else:
                logger.info(f"排除物种 {species}，覆盖率 {len(genes)}/{total_genes} 基因 ({len(genes)/total_genes*100:.1f}%) 低于阈值 {min_coverage_pct}%")
        
        logger.info(f"经过覆盖率过滤后包括 {len(included_species)}/{len(all_species)} 个物种")
        included_genes = valid_gene_results
        
    elif missing_species_strategy == 'exclude_genes':
        # Only include genes that have all species
        included_species = all_species
        for gene_result in valid_gene_results:
            gene_name = gene_result['name']
            sequences = gene_result['sequences']
            if set(sequences.keys()) == all_species:
                included_genes.append(gene_result)
            else:
                missing = all_species - set(sequences.keys())
                logger.info(f"排除基因 {gene_name}，缺失 {len(missing)}/{len(all_species)} 个物种")
        
        logger.info(f"包括 {len(included_genes)}/{total_genes} 个拥有所有物种的基因")
        
    else:  # 'gaps' strategy
        included_species = all_species
        included_genes = valid_gene_results
    
    # Check if we have any results after filtering
    if not included_species or not included_genes:
        logger.warning(f"策略 '{missing_species_strategy}' 过滤后没有剩余结果")
        return {}, []
    
    # Second pass: build the supergene
    logger.info(f"构建超基因，包含 {len(included_species)} 个物种和 {len(included_genes)} 个基因")
    gene_lengths = {}
    
    for gene_result in included_genes:
        sequences = gene_result['sequences']
        gene_name = gene_result['name']
        
        # Find the length of sequences for this gene
        if not sequences:
            logger.warning(f"基因 {gene_name} 没有序列")
            continue
            
        seq_length = len(next(iter(sequences.values())))
        gene_lengths[gene_name] = seq_length
        
        # Add sequence or gaps for each species
        for species in included_species:
            if species in sequences:
                species_supergenes[species] += sequences[species]
            else:
                # Add gaps for missing species
                species_supergenes[species] += '-' * seq_length
    
    logger.debug(f"基因长度: {gene_lengths}")
    
    # Calculate and log coverage statistics
    coverage_stats = {}
    for species in included_species:
        genes = species_gene_coverage[species]
        genes_in_supergene = [g for g in genes if g in [gene['name'] for gene in included_genes]]
        coverage_pct = (len(genes_in_supergene) / len(included_genes)) * 100 if included_genes else 0
        coverage_stats[species] = f"{len(genes_in_supergene)}/{len(included_genes)} 基因 ({coverage_pct:.1f}%)"
    
    logger.info(f"使用 '{missing_species_strategy}' 策略创建包含 {len(species_supergenes)} 个物种的超基因")
    logger.info(f"物种覆盖率统计: {coverage_stats}")
    
    # Create output directory structure
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    
    # Generate a detailed coverage matrix file
    output_dir = os.path.dirname(output_file)
    stats_dir = os.path.join(os.path.dirname(os.path.dirname(output_file)), "stats")
    os.makedirs(stats_dir, exist_ok=True)
    
    coverage_matrix = os.path.join(stats_dir, f"species_coverage_matrix_{missing_species_strategy}.tsv")
    with open(coverage_matrix, 'w') as f:
        # Write header with gene names
        f.write("Species\t" + "\t".join([g['name'] for g in included_genes]) + "\tCoverage(%)\n")
        
        # Write coverage for each species
        for species in sorted(included_species):
            row = [species]
            gene_present = []
            for gene in included_genes:
                if species in gene['sequences']:
                    row.append("1")
                    gene_present.append(1)
                else:
                    row.append("0")
                    gene_present.append(0)
            
            coverage = sum(gene_present) / len(gene_present) * 100 if gene_present else 0
            row.append(f"{coverage:.1f}")
            f.write("\t".join(row) + "\n")
    
    logger.info(f"物种覆盖率矩阵已写入 {coverage_matrix}")
    
    # Define output file name with strategy
    strategy_output_file = output_file.replace('.fasta', f'_{missing_species_strategy}.fasta')
    
    # Write the output file
    write_output(species_supergenes, strategy_output_file)
    
    return species_supergenes, included_genes

def merge_full_sequences(all_gene_results, output_file, missing_species_strategy='gaps', min_coverage_pct=0):
    """
    Merge full aligned sequences and translate to amino acids.
    
    Args:
        all_gene_results: List of gene results with sequence data
        output_file: Base output file path (strategy will be appended)
        missing_species_strategy: Strategy for handling missing species
        min_coverage_pct: Minimum percentage of genes a species must be present in
    """
    logger.info(f"使用策略 '{missing_species_strategy}' 合并和翻译完整序列")
    
    # Dictionary to store concatenated sequences for each species
    species_supergenes = defaultdict(str)
    
    # Dictionary to track which genes contribute to each species
    species_gene_coverage = defaultdict(list)
    
    # Keep track of all species encountered and valid genes
    all_species = set()
    valid_gene_results = [g for g in all_gene_results if g is not None]
    total_genes = len(valid_gene_results)
    
    # First pass: identify all species
    for gene_result in valid_gene_results:
        all_species.update(gene_result['aligned_seqs'].keys())
    
    # Calculate gene coverage for each species
    for gene_result in valid_gene_results:
        gene_name = gene_result['name']
        sequences = gene_result['aligned_seqs']
        
        for species in all_species:
            if species in sequences:
                species_gene_coverage[species].append(gene_name)
    
    # Apply filtering based on strategy
    included_species = set()
    included_genes = []
    
    if missing_species_strategy == 'exclude_species':
        # Only include species that meet the minimum coverage threshold
        min_genes_required = int(total_genes * min_coverage_pct / 100)
        for species, genes in species_gene_coverage.items():
            if len(genes) >= min_genes_required:
                included_species.add(species)
            else:
                logger.info(f"排除物种 {species}，覆盖率 {len(genes)}/{total_genes} 基因 ({len(genes)/total_genes*100:.1f}%) 低于阈值 {min_coverage_pct}%")
        
        included_genes = valid_gene_results
        
    elif missing_species_strategy == 'exclude_genes':
        # Only include genes that have all species
        included_species = all_species
        for gene_result in valid_gene_results:
            gene_name = gene_result['name']
            sequences = gene_result['aligned_seqs']
            if set(sequences.keys()) == all_species:
                included_genes.append(gene_result)
            else:
                missing = all_species - set(sequences.keys())
                logger.info(f"排除基因 {gene_name}，缺失 {len(missing)}/{len(all_species)} 个物种")
        
    else:  # 'gaps' strategy
        included_species = all_species
        included_genes = valid_gene_results
    
    # Check if we have any results after filtering
    if not included_species or not included_genes:
        logger.warning(f"策略 '{missing_species_strategy}' 过滤后没有剩余结果")
        return {}
    
    # Second pass: build the nucleotide supergene
    logger.info(f"构建完整序列超基因，包含 {len(included_species)} 个物种和 {len(included_genes)} 个基因")
    gene_lengths = {}
    
    for gene_result in included_genes:
        sequences = gene_result['aligned_seqs']
        gene_name = gene_result['name']
        
        # Find the length of sequences for this gene
        if not sequences:
            logger.warning(f"基因 {gene_name} 没有序列")
            continue
            
        seq_length = len(next(iter(sequences.values())))
        gene_lengths[gene_name] = seq_length
        
        # Add sequence or gaps for each species
        for species in included_species:
            if species in sequences:
                species_supergenes[species] += sequences[species]
            else:
                # Add gaps for missing species
                species_supergenes[species] += '-' * seq_length
    
    total_length = sum(gene_lengths.values()) if gene_lengths else 0
    logger.info(f"完整序列超基因长度: {total_length} bp")
    logger.debug(f"基因长度: {gene_lengths}")
    
    # Define output file names with strategy for nucleotide and amino acid files
    full_dir = os.path.join(os.path.dirname(os.path.dirname(output_file)), "full_cds")
    protein_dir = os.path.join(os.path.dirname(os.path.dirname(output_file)), "proteins")
    
    os.makedirs(full_dir, exist_ok=True)
    os.makedirs(protein_dir, exist_ok=True)
    
    strategy_nuc_file = os.path.join(full_dir, os.path.basename(output_file).replace('.fasta', f'_full_{missing_species_strategy}.fasta'))
    strategy_aa_file = os.path.join(protein_dir, os.path.basename(output_file).replace('.fasta', f'_protein_{missing_species_strategy}.fasta'))
    
    # Write nucleotide output
    write_output(species_supergenes, strategy_nuc_file)
    
    # Translate and write amino acid output
    amino_acid_seqs = {}
    translation_errors = 0
    
    for species, seq in species_supergenes.items():
        try:
            # Translate the sequence in frame 0
            aa_seq = translate_with_gaps(seq)
            amino_acid_seqs[species] = aa_seq
        except Exception as e:
            logger.error(f"翻译 {species} 序列时出错: {str(e)}")
            logger.error(traceback.format_exc())
            translation_errors += 1
            # Use placeholder with X's to maintain species in output
            amino_acid_seqs[species] = 'X' * (len(seq) // 3 + (1 if len(seq) % 3 else 0))
    
    if translation_errors:
        logger.warning(f"翻译过程中遇到 {translation_errors} 个错误")
    
    write_output(amino_acid_seqs, strategy_aa_file)
    logger.info(f"已为 '{missing_species_strategy}' 策略写入核苷酸和氨基酸文件")
    
    return species_supergenes

def run_trimal(input_file, output_file, trimal_path, automated=True, 
              gap_threshold=None, consistency_threshold=None, 
              conservation_threshold=None):
    """Run TrimAl to automatically trim MSA.
    
    Args:
        input_file: Path to input alignment file
        output_file: Path to output trimmed alignment
        trimal_path: Path to trimal executable
        automated: Whether to use automated trimming method
        gap_threshold: Minimum gap threshold (if automated=False)
        consistency_threshold: Minimum consistency threshold (if automated=False)
        conservation_threshold: Minimum conservation threshold (if automated=False)
    """
    start_time = time.time()
    
    # Ensure trimal_path exists and is executable
    if not os.path.exists(trimal_path):
        logger.error(f"TrimAl 可执行文件未在此路径找到: {trimal_path}")
        return False, "可执行文件未找到"
    
    if not os.access(trimal_path, os.X_OK):
        logger.error(f"TrimAl 可执行文件没有执行权限: {trimal_path}")
        return False, "可执行文件没有执行权限"
    
    # Build command
    cmd = [trimal_path, '-in', input_file, '-out', output_file]
    
    # Add method options
    if automated:
        cmd.append('-automated1')
    else:
        if gap_threshold is not None:
            cmd.extend(['-gt', str(gap_threshold)])
        if consistency_threshold is not None:
            cmd.extend(['-ct', str(consistency_threshold)])
        if conservation_threshold is not None:
            cmd.extend(['-st', str(conservation_threshold)])
    
    logger.info(f"运行 TrimAl: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"TrimAl 错误 (返回码 {result.returncode}): {result.stderr}")
            return False, f"返回码 {result.returncode}: {result.stderr}"
            
        execution_time = time.time() - start_time
        logger.info(f"TrimAl 成功完成, 耗时 {execution_time:.1f} 秒")
        
        # Verify output file was created
        if not os.path.exists(output_file):
            logger.error(f"TrimAl 未创建预期的输出文件: {output_file}")
            return False, "未创建输出文件"
            
        # Check if the output file has valid content
        if os.path.getsize(output_file) == 0:
            logger.error(f"TrimAl 创建了空的输出文件: {output_file}")
            return False, "输出文件为空"
            
        # Check if the trimming was successful by ensuring we have valid sequences
        valid, frame_info = check_alignment_validity(output_file)
        if not valid:
            logger.error(f"TrimAl 输出文件不包含有效比对: {output_file}")
            return False, f"输出文件不包含有效比对: {frame_info}"
            
        return True, output_file
    except Exception as e:
        logger.error(f"运行 TrimAl 时发生异常: {str(e)}")
        logger.error(traceback.format_exc())
        return False, str(e)

def align_protein_sequences(protein_file, output_file, aligner_params, trimal_params):
    """
    Align protein sequences and optionally run TrimAl.
    
    Args:
        protein_file: Path to protein sequences file
        output_file: Path for output aligned file
        aligner_params: Parameters for the aligner
        trimal_params: Parameters for TrimAl
    """
    logger.info(f"对蛋白质序列进行比对: {protein_file}")
    
    # Modify aligner parameters for protein alignment
    protein_aligner_params = aligner_params.copy()
    protein_aligner_params['codon_aware'] = False  # Protein sequences don't use codon awareness
    
    # Run the appropriate aligner
    success = False
    output_path = ""
    
    if aligner_params['aligner'] == 'prank':
        output_prefix = output_file.replace('.fasta', '')
        success, output_path = run_prank(
            protein_file,
            output_prefix,
            aligner_params['prank_path'],
            codon_aware=False,
            gaprate=0.005,      # Default protein gap opening rate for PRANK
            gapext=0.5          # Default protein gap extension for PRANK
        )
        if success:
            output_path = f"{output_prefix}.best.fas"
    
    elif aligner_params['aligner'] == 'muscle':
        success, output_path = run_muscle(
            protein_file,
            output_file,
            aligner_params['muscle_path'],
            codon_aware=False
        )
    
    elif aligner_params['aligner'] == 'mafft':
        success, output_path = run_mafft(
            protein_file,
            output_file,
            aligner_params['mafft_path'],
            codon_aware=False
        )
    
    if not success:
        logger.error(f"蛋白质序列对齐失败: {protein_file}")
        return False
    
    # Get the aligned output
    aligned_file = output_path if os.path.exists(output_path) else output_file
    
    # Run TrimAl if requested
    if trimal_params['use_trimal'] and trimal_params['trimal_path']:
        trimmed_file = output_file.replace('.fasta', '.trimmed.fasta')
        success, trimmed_path = run_trimal(
            aligned_file,
            trimmed_file,
            trimal_params['trimal_path'],
            automated=trimal_params['automated'],
            gap_threshold=trimal_params['gap_threshold'],
            consistency_threshold=trimal_params['consistency_threshold'],
            conservation_threshold=trimal_params['conservation_threshold']
        )
        
        if success:
            # Copy trimmed file to output file
            shutil.copy(trimmed_path, output_file)
            logger.info(f"蛋白质比对和修剪成功完成: {output_file}")
            return True
        else:
            logger.error(f"TrimAl修剪失败: {aligned_file}")
            # Still copy the untrimmed alignment
            shutil.copy(aligned_file, output_file)
            return False
    else:
        # Copy aligned file to output file
        shutil.copy(aligned_file, output_file)
        logger.info(f"蛋白质比对成功完成: {output_file}")
        return True

def create_protein_msa_supergene(all_gene_results, output_file, missing_species_strategy='gaps',
                                min_coverage_pct=0, aligner_params=None, trimal_params=None):
    """
    Create a supergene of protein sequences that have been aligned individually.
    
    Args:
        all_gene_results: List of gene results with sequence data
        output_file: Base output file path
        missing_species_strategy: Strategy for handling missing species
        min_coverage_pct: Minimum percentage of genes a species must be in
        aligner_params: Parameters for the aligner
        trimal_params: Parameters for TrimAl
    """
    logger.info(f"创建蛋白质MSA超基因，使用策略: {missing_species_strategy}")
    
    # First create nucleotide supergene to get filtering strategy results
    nucleotide_supergenes, included_genes = merge_sequences_by_species(
        all_gene_results,
        output_file,
        missing_species_strategy=missing_species_strategy,
        min_coverage_pct=min_coverage_pct
    )
    
    if not nucleotide_supergenes or not included_genes:
        logger.warning(f"无法创建核苷酸超基因，无法继续创建蛋白质MSA超基因")
        return {}
    
    # Process gene by gene
    species_protein_msa = defaultdict(str)
    temp_dir = os.path.join(os.path.dirname(os.path.dirname(output_file)), "temp", "protein_msa")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Create output directory for individual protein MSAs
    protein_msa_dir = os.path.join(os.path.dirname(os.path.dirname(output_file)), "protein_msa")
    os.makedirs(protein_msa_dir, exist_ok=True)
    
    for gene_result in included_genes:
        gene_name = gene_result['name']
        aligned_seqs = gene_result['aligned_seqs']
        
        # Translate each sequence
        protein_seqs = {}
        for species, seq in aligned_seqs.items():
            try:
                protein_seqs[species] = translate_with_gaps(seq)
            except Exception as e:
                logger.error(f"翻译基因 {gene_name} 的物种 {species} 时出错: {str(e)}")
                # Skip this species for this gene
                continue
        
        # Write protein sequences to temp file
        protein_file = os.path.join(temp_dir, f"{gene_name}_protein.fasta")
        with open(protein_file, 'w') as f:
            for species, seq in protein_seqs.items():
                f.write(f">{species}\n{seq}\n")
        
        # Align protein sequences
        protein_msa_file = os.path.join(protein_msa_dir, f"{gene_name}_protein_msa.fasta")
        align_success = align_protein_sequences(protein_file, protein_msa_file, aligner_params, trimal_params)
        
        if align_success:
            # Read aligned protein sequences
            aligned_proteins = {}
            try:
                for record in SeqIO.parse(protein_msa_file, "fasta"):
                    species = extract_species_name(record.description)
                    aligned_proteins[species] = str(record.seq)
                
                # Add to protein MSA supergene for species that were included
                included_species = set(nucleotide_supergenes.keys())
                
                seq_length = 0
                if aligned_proteins:
                    seq_length = len(next(iter(aligned_proteins.values())))
                    
                for species in included_species:
                    if species in aligned_proteins:
                        species_protein_msa[species] += aligned_proteins[species]
                    else:
                        # Add gaps for missing species
                        species_protein_msa[species] += '-' * seq_length
                
                logger.info(f"基因 {gene_name} 的蛋白质MSA成功添加到超基因 (长度: {seq_length})")
            except Exception as e:
                logger.error(f"处理基因 {gene_name} 的蛋白质MSA时出错: {str(e)}")
                logger.error(traceback.format_exc())
        else:
            logger.warning(f"基因 {gene_name} 的蛋白质序列比对失败，跳过此基因")
    
    # Write final protein MSA supergene
    protein_msa_supergene_file = os.path.join(protein_msa_dir, os.path.basename(output_file).replace('.fasta', f'_protein_msa_{missing_species_strategy}.fasta'))
    write_output(species_protein_msa, protein_msa_supergene_file)
    
    # Optionally trim the entire protein MSA supergene
    if trimal_params and trimal_params['use_trimal'] and trimal_params['trim_supergene']:
        trimmed_file = protein_msa_supergene_file.replace('.fasta', '.trimmed.fasta')
        success, _ = run_trimal(
            protein_msa_supergene_file,
            trimmed_file,
            trimal_params['trimal_path'],
            automated=trimal_params['automated'],
            gap_threshold=trimal_params['gap_threshold'],
            consistency_threshold=trimal_params['consistency_threshold'],
            conservation_threshold=trimal_params['conservation_threshold']
        )
        
        if success:
            logger.info(f"蛋白质MSA超基因修剪成功: {trimmed_file}")
    
    return species_protein_msa

def write_output(supergene_seqs, output_file):
    """Write supergene sequences to FASTA file."""
    if not supergene_seqs:
        logger.warning(f"没有序列可写入 {output_file}")
        return False
        
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        with open(output_file, 'w') as f:
            for species, seq in sorted(supergene_seqs.items()):
                f.write(f">{species}\n")
                # Write sequence in blocks of 60 characters for readability
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')
        
        seq_length = len(next(iter(supergene_seqs.values()))) if supergene_seqs else 0
        logger.info(f"成功写入 {len(supergene_seqs)} 条序列到 {output_file} (长度: {seq_length} bp)")
        return True
    except Exception as e:
        logger.error(f"写入 {output_file} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def setup_output_dirs(output_dir):
    """Set up the output directory structure."""
    # Create main output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subdirectories for different outputs
    dirs = {
        "alignments": os.path.join(output_dir, "alignments"),
        "4d_sites": os.path.join(output_dir, "4d_sites"),
        "full_cds": os.path.join(output_dir, "full_cds"),
        "proteins": os.path.join(output_dir, "proteins"),
        "protein_msa": os.path.join(output_dir, "protein_msa"),
        "stats": os.path.join(output_dir, "stats"),
        "temp": os.path.join(output_dir, "temp")
    }
    
    for dir_name, dir_path in dirs.items():
        os.makedirs(dir_path, exist_ok=True)
        logger.debug(f"创建目录: {dir_path}")
    
    return dirs

def main():
    parser = argparse.ArgumentParser(description="Multiple sequence alignment and 4D site extraction pipeline")
    parser.add_argument("--input_dir", default=".", help="Directory containing CDS files")
    parser.add_argument("--output_dir", default="./output", help="Directory for output files")
    parser.add_argument("--supergene_output", default="supergene_4d.fasta", help="Base output file for supergene")
    
    # Aligner selection
    parser.add_argument("--aligner", choices=["prank", "muscle", "mafft"], default="prank",
                      help="Alignment tool to use (default: prank)")
    parser.add_argument("--prank_path", help="Absolute path to PRANK executable")
    parser.add_argument("--muscle_path", help="Absolute path to MUSCLE executable")
    parser.add_argument("--mafft_path", help="Absolute path to MAFFT executable")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    
    # Common alignment parameters
    # 修复：将store_true和default=True的冲突修改为直接使用store_false和default=True
    parser.add_argument("--no_codon_aware", action="store_true", help="Disable codon-aware alignment")
    
    # PRANK specific parameters
    parser.add_argument("--f", type=float, default=0.2, help="PRANK insertion opening probability")
    parser.add_argument("--gaprate", type=float, default=None, 
                      help="PRANK gap opening rate (default: PRANK's default based on data type)")
    parser.add_argument("--gapext", type=float, default=None, 
                      help="PRANK gap extension probability (default: PRANK's default based on data type)")
    parser.add_argument("--use_logs", action="store_true", 
                      help="Use logarithm calculations in PRANK for large datasets")
    parser.add_argument("--penalize_terminal_gaps", action="store_true",
                      help="Penalize terminal gaps in PRANK normally")
    
    # TrimAl parameters
    parser.add_argument("--use_trimal", action="store_true", help="Use TrimAl for protein alignment trimming")
    parser.add_argument("--trimal_path", help="Absolute path to TrimAl executable")
    parser.add_argument("--trimal_automated", action="store_true", default=True, 
                      help="Use automated trimming method in TrimAl")
    parser.add_argument("--gap_threshold", type=float, default=None, help="TrimAl minimum gap threshold")
    parser.add_argument("--consistency_threshold", type=float, default=None, help="TrimAl consistency threshold")
    parser.add_argument("--conservation_threshold", type=float, default=None, help="TrimAl conservation threshold")
    parser.add_argument("--trim_supergene", action="store_true", help="Apply TrimAl to the final protein supergene")
    
    # Other parameters
    parser.add_argument("--duplicate_strategy", 
                      choices=['longest', 'first', 'rename', 'alignment_quality'], 
                      default='alignment_quality',
                      help="Strategy for handling duplicate species in a file")
    parser.add_argument("--skip_existing", action="store_true", help="Skip processing if alignment file exists")
    parser.add_argument("--min_coverage_pct", type=float, default=50.0,
                      help="Minimum percentage of genes a species must be present in (for exclude_species strategy)")
    parser.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO',
                      help="Set logging level")
    parser.add_argument("--clean_temp", action="store_true", default=True,
                      help="Clean temporary files after processing")
    parser.add_argument("--create_protein_msa", action="store_true", default=False,
                      help="Create multiple sequence alignments of protein sequences (abandoned)")
    
    args = parser.parse_args()
    
    # Validate aligner path based on selection
    if args.aligner == "prank" and not args.prank_path:
        parser.error("选择PRANK作为比对工具时，必须提供--prank_path参数")
    elif args.aligner == "muscle" and not args.muscle_path:
        parser.error("选择MUSCLE作为比对工具时，必须提供--muscle_path参数")
    elif args.aligner == "mafft" and not args.mafft_path:
        parser.error("选择MAFFT作为比对工具时，必须提供--mafft_path参数")
        
    # Validate TrimAl if enabled
    if args.use_trimal and not args.trimal_path:
        parser.error("启用TrimAl修剪时，必须提供--trimal_path参数")
    
    # Set logging level
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # Log start time and parameters
    logger.info(f"开始 MSA 流水线，参数: {vars(args)}")
    start_time = time.time()
    
    # Set up directory structure
    try:
        dirs = setup_output_dirs(args.output_dir)
    except Exception as e:
        logger.error(f"设置输出目录时出错: {str(e)}")
        return 1
    
    # Find CDS files
    cds_files = find_cds_files(args.input_dir)
    if not cds_files:
        logger.error("未找到 .cds 文件。退出。")
        return 1
        
    logger.info(f"将处理 {len(cds_files)} 个 CDS 文件")
    
    # Prepare parameters for aligners
    aligner_params = {
        'aligner': args.aligner,
        'prank_path': args.prank_path,
        'muscle_path': args.muscle_path,
        'mafft_path': args.mafft_path,
        'codon_aware': not args.no_codon_aware,  # 修复：默认启用密码子感知
        'duplicate_strategy': args.duplicate_strategy,
        'clean_temp': args.clean_temp,
        
        # PRANK specific
        'f': args.f,
        'gaprate': args.gaprate,
        'gapext': args.gapext,
        'use_logs': args.use_logs,
        'penalize_terminal_gaps': args.penalize_terminal_gaps
    }
    
    # TrimAl parameters
    trimal_params = {
        'use_trimal': args.use_trimal,
        'trimal_path': args.trimal_path,
        'automated': args.trimal_automated,
        'gap_threshold': args.gap_threshold,
        'consistency_threshold': args.consistency_threshold,
        'conservation_threshold': args.conservation_threshold,
        'trim_supergene': args.trim_supergene
    }
    
    # Process each CDS file in parallel
    all_gene_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {}
        
        for file_name in cds_files:
            file_path = os.path.join(args.input_dir, file_name)
            gene_name = file_name.replace('.cds', '')
            output_file = os.path.join(dirs['alignments'], f"{gene_name}.best.fas")
            
            # Skip existing files if requested
            if args.skip_existing and os.path.exists(output_file):
                logger.info(f"跳过 {file_name}，输出已存在")
                
                try:
                    # Still need to process the existing file for 4D sites
                    aligned_seqs = {}
                    for record in SeqIO.parse(output_file, "fasta"):
                        species = extract_species_name(record.description)
                        aligned_seqs[species] = str(record.seq)
                        
                    if not aligned_seqs:
                        logger.warning(f"无法解析现有的比对 {output_file}")
                        continue
                        
                    # Verify and fix codon alignment if needed
                    aligned_length = len(next(iter(aligned_seqs.values())))
                    if aligned_length % 3 != 0:
                        logger.warning(f"现有比对文件 {output_file} 长度 {aligned_length} 不是3的倍数，修复中")
                        fixed_file = output_file.replace('.best.fas', '.framed.fas')
                        if fix_alignment_frame(output_file, fixed_file):
                            # Re-read the aligned sequences from fixed file
                            aligned_seqs = {}
                            for record in SeqIO.parse(fixed_file, "fasta"):
                                species = extract_species_name(record.description)
                                aligned_seqs[species] = str(record.seq)
                            logger.info(f"使用修复后的比对")
                        
                    fourfold_sites = identify_4d_sites(aligned_seqs)
                    fourfold_seqs = extract_4d_sites(aligned_seqs, fourfold_sites)
                    
                    all_gene_results.append({
                        'name': gene_name,
                        'species_count': len(fourfold_seqs),
                        'site_count': len(fourfold_sites),
                        'sequences': fourfold_seqs,
                        'aligned_seqs': aligned_seqs
                    })
                    logger.info(f"处理现有比对 {file_name}: {len(fourfold_seqs)} 个物种, {len(fourfold_sites)} 个 4D 位点")
                except Exception as e:
                    logger.error(f"处理现有比对 {output_file} 时出错: {str(e)}")
                continue
            
            futures[executor.submit(
                process_cds_file, 
                file_path,
                args.output_dir,
                aligner_params
            )] = file_name
        
        # Collect results
        completed = 0
        total_futures = len(futures)
        for future in futures:
            file_name = futures[future]
            try:
                result = future.result()
                completed += 1
                progress = completed / total_futures * 100
                logger.info(f"进度: {completed}/{total_futures} 文件处理 ({progress:.1f}%)")
                
                if result:
                    all_gene_results.append(result)
                    logger.info(f"处理 {file_name}: {result['species_count']} 个物种, {result['site_count']} 个 4D 位点")
                else:
                    logger.error(f"处理 {file_name} 失败")
            except Exception as e:
                logger.exception(f"处理 {file_name} 时出错: {str(e)}")
    
    if not all_gene_results:
        logger.error("没有基因结果成功处理")
        return 1
        
    # Create output path for supergene
    supergene_4d_path = os.path.join(dirs['4d_sites'], args.supergene_output)
    
    # Process with all three strategies
    strategies = ['gaps', 'exclude_species', 'exclude_genes']
    
    # Create 4D site supergenes for each strategy
    for strategy in strategies:
        try:
            logger.info(f"使用 '{strategy}' 策略创建 4D 位点超基因")
            supergene_seqs, included_genes = merge_sequences_by_species(
                all_gene_results, 
                supergene_4d_path,
                missing_species_strategy=strategy,
                min_coverage_pct=args.min_coverage_pct
            )
            if not supergene_seqs:
                logger.warning(f"策略 '{strategy}' 未能生成有效的 4D 位点超基因")
        except Exception as e:
            logger.error(f"使用策略 '{strategy}' 创建 4D 位点超基因时出错: {str(e)}")
            logger.error(traceback.format_exc())
    
    # Create full and translated amino acid supergenes for each strategy
    for strategy in strategies:
        try:
            logger.info(f"使用 '{strategy}' 策略创建完整核苷酸和蛋白质超基因")
            merged_seqs = merge_full_sequences(
                all_gene_results,
                supergene_4d_path,
                missing_species_strategy=strategy,
                min_coverage_pct=args.min_coverage_pct
            )
            
            if not merged_seqs:
                logger.warning(f"策略 '{strategy}' 未能生成有效的完整/蛋白质超基因")
                
            # Create protein MSA supergene if requested
            if args.create_protein_msa and merged_seqs:
                logger.info(f"使用 '{strategy}' 策略创建蛋白质MSA超基因")
                create_protein_msa_supergene(
                    all_gene_results,
                    supergene_4d_path,
                    missing_species_strategy=strategy,
                    min_coverage_pct=args.min_coverage_pct,
                    aligner_params=aligner_params,
                    trimal_params=trimal_params
                )
                
        except Exception as e:
            logger.error(f"使用策略 '{strategy}' 创建完整/蛋白质超基因时出错: {str(e)}")
            logger.error(traceback.format_exc())
    
    # Clean up temporary files if requested
    if args.clean_temp:
        try:
            logger.info("清理临时文件")
            temp_dir = dirs['temp']
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        except Exception as e:
            logger.warning(f"清理临时文件时出错: {str(e)}")
    
    # Log execution time
    end_time = time.time()
    execution_time = end_time - start_time
    hours, remainder = divmod(execution_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    logger.info(f"流水线完成。耗时: {int(hours)}小时 {int(minutes)}分钟 {seconds:.1f}秒")
    logger.info(f"结果保存在: {args.output_dir}")
    
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logger.critical(f"未捕获的异常: {str(e)}")
        logger.critical(traceback.format_exc())
        sys.exit(1)
