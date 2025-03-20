#!/dellfsqd2/ST_OCEAN/USER/lishuo1/01_software/miniconda3/bin/python
import sys
import os
import getopt
sys.path.append('/dellfsqd2/ST_OCEAN/USER/lishuo1/01_software/miniconda3/lib/python3.7/site-packages/')
import subprocess
import argparse
import logging
import traceback
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from functools import lru_cache
import re
import time
import shutil
import uuid
import json
import itertools
import hashlib
import multiprocessing
# 提前导入可能需要的所有模块
from Bio.Blast import NCBIXML
import numpy as np
from functools import lru_cache

# 设置日志记录，带有时间戳
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("msa_pipeline.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 定义4倍简并密码子（第三位可以是任何核苷酸）
FOURFOLD_CODONS = {
    'GC': True,  # 丙氨酸 (GCN)
    'CG': True,  # 精氨酸 (CGN)
    'GG': True,  # 甘氨酸 (GGN)
    'CT': True,  # 亮氨酸 (CTN)
    'CC': True,  # 脯氨酸 (CCN)
    'TC': True,  # 丝氨酸 (TCN)
    'AC': True,  # 苏氨酸 (ACN)
    'GT': True   # 缬氨酸 (GTN)
}

# 定义遗传密码表，将密码子映射到氨基酸（包括模糊密码子）
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

# 定义氨基酸相似性组，用于更好的保守性评分
AA_SIMILARITY_GROUPS = [
    {'G', 'A', 'S', 'T'},           # 小/极性
    {'C'},                          # 半胱氨酸
    {'V', 'I', 'L', 'M'},           # 疏水性
    {'P'},                          # 脯氨酸
    {'F', 'Y', 'W'},                # 芳香族
    {'H', 'K', 'R'},                # 正电荷
    {'D', 'E', 'N', 'Q'},           # 负电荷/极性
]

# 最小报告重复次数 - 提高阈值减少误报
MIN_REPEAT_COUNT_TO_REPORT = 70  # 提高阈值，减少CDS中正常重复的警告

# 4D位点识别中允许的最大缺口比例
MAX_GAP_RATIO_FOR_4D = 0.2  # 允许20%的物种在4D位点位置有缺口

# 'exclude_genes'策略中允许缺失的最大物种比例
MAX_MISSING_SPECIES_RATIO = 0.1  # 允许10%的物种缺失

# 性能优化：BLAST数据库缓存
BLAST_DB_CACHE = {}
BLAST_DB_CACHE_SIZE = 5  # 最大缓存数据库数量

# 性能优化：缓存函数计算结果
@lru_cache(maxsize=1024)
def safe_translate_codon_cached(codon):
    """
    安全地将密码子翻译为氨基酸的缓存版本
    
    Args:
        codon: 三字母核苷酸密码子
        
    Returns:
        单字母氨基酸或'X'表示未知
    """
    return safe_translate_codon(codon)

def find_cds_files(directory='.'):
    """查找给定目录中的所有.cds文件"""
    try:
        cds_files = [f for f in os.listdir(directory) if f.endswith('.cds')]
        logger.info(f"找到 {len(cds_files)} 个 .cds 文件在 {directory} 中")
        return sorted(cds_files)  # 排序，保证处理顺序一致
    except Exception as e:
        logger.error(f"在 {directory} 中查找 CDS 文件时出错: {str(e)}")
        return []

def extract_species_name(header):
    """
    根据FASTA头部提取物种名称，支持">ID SPECIES_NAME"格式
    
    Args:
        header: FASTA头部字符串
    
    Returns:
        提取的物种名称
    """
    # 移除'>'前缀如果存在
    if header.startswith('>'):
        header = header[1:]
    
    # 按空格分割，第一部分是ID，后面的是物种名
    parts = header.split(' ', 1)
    
    if len(parts) > 1:
        # 返回空格后的部分作为物种名
        return parts[1].strip()
    else:
        # 如果没有空格，则使用整个标识符作为物种名
        return header.strip()

def check_coding_sequence(sequence):
    """
    检查序列是否看起来像有效的编码序列
    
    Returns:
        tuple: (is_valid, message)
    """
    if not sequence:
        return False, "空序列"
        
    if len(sequence) % 3 != 0:
        return False, f"序列长度 {len(sequence)} 不是 3 的倍数"
    
    # 检查中间的终止密码子
    for i in range(0, len(sequence) - 3, 3):  # 排除最后一个密码子
        codon = sequence[i:i+3].upper()
        if codon in ['TAA', 'TAG', 'TGA'] and i < len(sequence)-5:
            return False, f"序列中间发现终止密码子 {codon} 在位置 {i}"
    
    # 检查无效字符
    valid_chars = set('ACGTRYKMSWBDHVN-')  # DNA + 模糊代码 + 缺口
    invalid_chars = set(sequence.upper()) - valid_chars
    if invalid_chars:
        return False, f"序列包含无效字符: {', '.join(invalid_chars)}"
    
    return True, "有效的编码序列"

def analyze_codon_structure(sequence):
    """
    分析序列的密码子结构，识别可能的起始/终止位点和开放阅读框
    
    Args:
        sequence: 核苷酸序列
    
    Returns:
        包含分析结果的字典
    """
    # 性能优化：使用numpy数组加速处理
    # 去除缺口以进行更准确的ORF分析
    clean_seq = sequence.upper().replace('-', '').replace('N', '')
    
    # 分析3个可能的阅读框
    frames = []
    frame_orfs = []
    
    for frame in range(3):
        codons = []
        for i in range(frame, len(clean_seq)-2, 3):
            codons.append(clean_seq[i:i+3])
        frames.append(codons)
        
        # 查找这个阅读框中的所有ORFs
        orfs = []
        current_orf = []
        
        for i, codon in enumerate(codons):
            if codon in ['TAA', 'TAG', 'TGA']:
                if current_orf:
                    orfs.append({
                        'start': frame + 3*frames[frame].index(current_orf[0]),
                        'end': frame + 3*(frames[frame].index(current_orf[0])+len(current_orf)),
                        'codons': current_orf
                    })
                    current_orf = []
            else:
                current_orf.append(codon)
        
        # 添加最后一个ORF(如果没有终止密码子)
        if current_orf:
            orfs.append({
                'start': frame + 3*frames[frame].index(current_orf[0]),
                'end': frame + 3*(frames[frame].index(current_orf[0])+len(current_orf)),
                'codons': current_orf
            })
        
        frame_orfs.append(orfs)
    
    # 查找起始密码子
    start_positions = []
    for frame_idx, frame_codons in enumerate(frames):
        for i, codon in enumerate(frame_codons):
            if codon == 'ATG':
                # 记录可能的起始位置
                start_positions.append({
                    'frame': frame_idx,
                    'position': frame_idx + i*3,
                    'codon_idx': i
                })
    
    # 找最长的ORF
    longest_orf = None
    best_frame = 0
    
    for frame_idx, orfs in enumerate(frame_orfs):
        for orf in orfs:
            if not longest_orf or len(orf['codons']) > len(longest_orf['codons']):
                longest_orf = orf
                best_frame = frame_idx
    
    # 检查最长ORF是否以起始密码子开始
    has_start = False
    if longest_orf and longest_orf['codons']:
        has_start = (longest_orf['codons'][0] == 'ATG')
    
    result = {
        'best_frame': best_frame,
        'longest_orf': longest_orf,
        'has_start': has_start,
        'has_stop': False,  # 最长ORF可能没有终止子，会在后续逻辑中处理
        'start_positions': start_positions,
        'orf_length': len(longest_orf['codons'])*3 if longest_orf else 0
    }
    
    return result

def assess_sequence_quality(sequence, seq_id="Unknown"):
    """
    评估序列质量，检测潜在问题
    仅报告非常极端的重复情况，忽略CDS中正常的重复
    
    Args:
        sequence: 核苷酸序列
        seq_id: 序列ID (用于日志)
    
    Returns:
        问题列表
    """
    issues = []
    
    # 清理序列以便更准确的分析
    clean_seq = sequence.upper().replace('-', '')
    
    # 检查序列长度
    if len(clean_seq) < 30:  # 小于10个密码子
        issues.append("序列过短")
    
    # 检查GC含量
    gc_count = clean_seq.count('G') + clean_seq.count('C')
    gc_content = gc_count / len(clean_seq) if clean_seq else 0
    if gc_content < 0.25 or gc_content > 0.75:
        issues.append(f"异常GC含量: {gc_content:.2f}")
    
    # 检查未知碱基比例
    n_count = clean_seq.count('N')
    n_ratio = n_count / len(clean_seq) if clean_seq else 0
    if n_ratio > 0.1:
        issues.append(f"高比例未知碱基: {n_ratio:.2f}")
    
    # 检查终止密码子
    if len(clean_seq) % 3 == 0:  # 只有当序列是3的倍数时才检查密码子
        for i in range(0, len(clean_seq)-5, 3):  # 排除最后一个密码子
            codon = clean_seq[i:i+3]
            if codon in ['TAA', 'TAG', 'TGA'] and i < len(clean_seq)-5:
                issues.append(f"序列中部有终止密码子，位置 {i}")
    
    # 性能优化：只检查7和10碱基长度的重复，这不是3的倍数，可能暗示frame shift问题
    logger.debug(f"分析序列 {seq_id} 的重复模式 (长度: {len(clean_seq)})")

    repeat_checks = [7, 10]  # 非3倍数长度，可能表示框架问题
    significant_repeats = []
    
    # 性能优化：只在序列长度足够时才检查重复
    if len(clean_seq) > 100:
        for repeat_len in repeat_checks:
            # 优化：使用字典直接计数而不是构建中间列表
            repeat_counts = defaultdict(int)
            
            for i in range(len(clean_seq) - repeat_len + 1):
                fragment = clean_seq[i:i+repeat_len]
                if 'N' not in fragment:  # 忽略含N的片段
                    repeat_counts[fragment] += 1
                    
            # 计算总重复数，只统计出现超过1次的片段
            total_repeats = sum(1 for count in repeat_counts.values() if count > 1)
            
            # 使用更高的阈值，减少误报
            threshold = MIN_REPEAT_COUNT_TO_REPORT
            if repeat_len == 10:
                threshold = MIN_REPEAT_COUNT_TO_REPORT // 2  # 较长重复的阈值降低
                
            # 只有当重复数大幅超过阈值时才报告，这可能表示真正的问题
            if total_repeats >= threshold:
                significant_repeats.append(f"发现 {total_repeats} 个长度为 {repeat_len} 的重复 (非3倍数，可能影响密码子框架)")
                logger.debug(f"序列 {seq_id} 中{repeat_len}碱基重复统计: 总计 {total_repeats} 个")
    
    # 只有存在显著重复时才添加到问题列表
    if significant_repeats:
        issues.extend(significant_repeats)
    
    return issues

def handle_non_triplet_sequence(sequence, seq_id, position='end', context=None):
    """
    更智能地处理非3倍数的序列
    
    Args:
        sequence: 原始序列
        seq_id: 序列ID用于日志
        position: 填充位置 'start', 'end', 'smart'
        context: 提供额外上下文信息的词典
    
    Returns:
        处理后的序列
    """
    remainder = len(sequence) % 3
    if remainder == 0:
        return sequence
    
    padding_needed = 3 - remainder
    
    # 根据序列特征决定填充策略
    if position == 'smart' and context is not None:
        # 检查序列是否有起始/终止密码子
        has_start = context.get('has_start', False) or sequence[:3].upper() == 'ATG'
        has_stop = context.get('has_stop', False) or sequence[-3:].upper() in ['TAA', 'TAG', 'TGA']
        
        # 基于密码子情况决定填充位置
        if has_start and not has_stop:
            # 可能是3'端不完整
            padding = 'N' * padding_needed
            padded_seq = sequence + padding
            logger.debug(f"序列 {seq_id} 有起始密码子但无终止密码子，在3'端添加 {padding_needed} 个N")
        elif not has_start and has_stop:
            # 可能是5'端不完整
            padding = 'N' * padding_needed
            padded_seq = padding + sequence
            logger.debug(f"序列 {seq_id} 无起始密码子但有终止密码子，在5'端添加 {padding_needed} 个N")
        else:
            # 默认在末尾添加，但记录更详细的信息
            padding = 'N' * padding_needed
            padded_seq = sequence + padding
            logger.debug(f"序列 {seq_id} 不含完整密码子信息，在3'端添加 {padding_needed} 个N")
    elif position == 'start':
        padding = 'N' * padding_needed
        padded_seq = padding + sequence
        logger.debug(f"序列 {seq_id} 在5'端添加 {padding_needed} 个N")
    else:  # 'end' 或默认
        padding = 'N' * padding_needed
        padded_seq = sequence + padding
        logger.debug(f"序列 {seq_id} 在3'端添加 {padding_needed} 个N")
    
    return padded_seq

def cleanup_temp_files(temp_fasta, temp_output, db_prefix):
    """清理临时BLAST文件，但保留缓存的数据库"""
    try:
        # 删除临时FASTA和输出文件
        for file in [temp_fasta, temp_output]:
            if file and os.path.exists(file):
                os.remove(file)
                logger.debug(f"已删除临时文件: {file}")
        
        # 不删除可能被缓存的数据库文件
        if db_prefix and db_prefix not in BLAST_DB_CACHE:
            for ext in [".nhr", ".nin", ".nsq", ".fasta"]:
                db_file = db_prefix + ext
                if os.path.exists(db_file):
                    os.remove(db_file)
                    logger.debug(f"已删除BLAST数据库文件: {db_file}")
    except Exception as e:
        logger.debug(f"清理临时文件时出错: {str(e)}")

def hash_file_content(file_path):
    """
    计算文件内容的哈希值，用于数据库缓存
    """
    try:
        hasher = hashlib.md5()
        with open(file_path, 'rb') as f:
            buf = f.read(65536)  # 64kb chunks
            while buf:
                hasher.update(buf)
                buf = f.read(65536)
        return hasher.hexdigest()
    except Exception as e:
        logger.debug(f"计算文件哈希时出错: {str(e)}")
        return None

def create_or_get_blast_db(file_path, temp_dir, blast_path):
    """
    创建或从缓存获取BLAST数据库
    
    Args:
        file_path: CDS文件路径
        temp_dir: 临时目录
        blast_path: BLAST可执行文件路径
    
    Returns:
        临时数据库路径或None
    """
    # 计算文件内容哈希作为缓存键
    file_hash = hash_file_content(file_path)
    if not file_hash:
        return create_temp_blast_db(file_path, temp_dir, blast_path)
    
    # 检查缓存中是否已存在数据库
    if file_hash in BLAST_DB_CACHE:
        db_prefix = BLAST_DB_CACHE[file_hash]
        # 验证数据库文件是否仍然存在
        if os.path.exists(db_prefix + ".nsq"):
            logger.info(f"使用缓存的BLAST数据库: {db_prefix}")
            return db_prefix
        else:
            # 缓存项有效但文件已删除，从缓存中移除
            del BLAST_DB_CACHE[file_hash]
    
    # 创建新数据库
    db_prefix = create_temp_blast_db(file_path, temp_dir, blast_path)
    if db_prefix:
        # 管理缓存大小
        if len(BLAST_DB_CACHE) >= BLAST_DB_CACHE_SIZE:
            # 删除最旧的数据库
            oldest_key = next(iter(BLAST_DB_CACHE))
            oldest_db = BLAST_DB_CACHE[oldest_key]
            del BLAST_DB_CACHE[oldest_key]
            # 清理关联文件
            for ext in [".nhr", ".nin", ".nsq", ".fasta"]:
                old_file = oldest_db + ext
                if os.path.exists(old_file):
                    os.remove(old_file)
        
        # 添加新数据库到缓存
        BLAST_DB_CACHE[file_hash] = db_prefix
    
    return db_prefix

def create_temp_blast_db(file_path, temp_dir, blast_path):
    """
    从输入文件自动创建临时BLAST数据库
    
    Args:
        file_path: CDS文件路径
        temp_dir: 临时目录
        blast_path: BLAST可执行文件路径
    
    Returns:
        临时数据库路径或None
    """
    try:
        # 确保临时目录存在
        os.makedirs(temp_dir, exist_ok=True)
        
        logger.info(f"开始为 {os.path.basename(file_path)} 创建临时BLAST数据库")
        
        # 读取所有有效序列（3的倍数）- 性能优化：批量读取
        valid_seqs = []
        total_seqs = 0
        
        # 使用生成器而不是一次读取所有
        for record in SeqIO.parse(file_path, "fasta"):
            total_seqs += 1
            sequence = str(record.seq)
            if len(sequence) % 3 == 0:
                valid_seqs.append(record)
                # 当收集足够的序列后，停止读取
                if len(valid_seqs) >= 100:  # 限制使用的序列数量
                    break
        
        # 如果没有有效序列，返回None
        if not valid_seqs:
            logger.warning(f"文件 {file_path} 中没有长度为3倍数的有效序列, 无法创建BLAST数据库")
            return None
        
        logger.info(f"从 {file_path} 中找到 {len(valid_seqs)}/{total_seqs} 个有效序列用于创建BLAST数据库")
        
        # 创建临时FASTA文件
        db_prefix = os.path.join(temp_dir, f"temp_blastdb_{uuid.uuid4()}")
        temp_fasta = f"{db_prefix}.fasta"
        
        with open(temp_fasta, 'w') as f:
            for record in valid_seqs:
                f.write(f">{record.id}\n{record.seq}\n")
        
        # 构建makeblastdb路径
        blast_dir = os.path.dirname(blast_path)
        makeblastdb_path = os.path.join(blast_dir, "makeblastdb")
        
        # 尝试几种可能的路径
        makeblastdb_candidates = [
            makeblastdb_path,
            blast_path.replace("blastn", "makeblastdb"),
            os.path.join(blast_dir, "makeblastdb.exe"),  # Windows
            blast_path.replace("blastn.exe", "makeblastdb.exe")  # Windows
        ]
        
        makeblastdb_path = None
        for candidate in makeblastdb_candidates:
            if os.path.exists(candidate):
                makeblastdb_path = candidate
                break
                
        if not makeblastdb_path:
            logger.error(f"无法找到makeblastdb可执行文件，尝试了: {', '.join(makeblastdb_candidates)}")
            return None
            
        logger.debug(f"使用makeblastdb路径: {makeblastdb_path}")
        
        # 运行makeblastdb
        cmd = [
            makeblastdb_path,
            "-in", temp_fasta,
            "-dbtype", "nucl",
            "-out", db_prefix
        ]
        
        logger.debug(f"创建BLAST数据库命令: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"创建BLAST数据库失败: {result.stderr}")
            return None
        
        logger.info(f"成功创建临时BLAST数据库: {db_prefix}")
        return db_prefix
    
    except Exception as e:
        logger.error(f"创建临时BLAST数据库时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def handle_non_triplet_with_blast(sequence, seq_id, file_path, blast_path, temp_dir):
    """
    使用BLAST处理非3倍数序列，修改逻辑以明确仅在必要时调用
    性能优化：共享和缓存BLAST数据库
    
    Args:
        sequence: 原始序列
        seq_id: 序列ID
        file_path: 原始CDS文件路径
        blast_path: BLAST可执行文件路径
        temp_dir: 临时文件目录
    
    Returns:
        处理后的序列
    """
    # 只对非3倍数序列启动BLAST处理
    if len(sequence) % 3 == 0:
        return sequence
    
    logger.info(f"启动BLAST处理非3倍数序列: {seq_id}, 长度 {len(sequence)}")
        
    if not blast_path or not temp_dir or not file_path:
        logger.warning(f"未提供完整BLAST参数，使用默认填充方法处理序列 {seq_id}")
        return handle_non_triplet_sequence(sequence, seq_id, position='smart')
    
    try:
        # 获取或创建BLAST数据库
        db_prefix = create_or_get_blast_db(file_path, temp_dir, blast_path)
        
        if not db_prefix:
            logger.warning(f"无法创建临时BLAST数据库，使用默认填充方法处理序列 {seq_id}")
            return handle_non_triplet_sequence(sequence, seq_id, position='smart')
        
        # 创建临时序列文件
        temp_fasta = os.path.join(temp_dir, f"{seq_id}_{uuid.uuid4()}.fasta")
        
        with open(temp_fasta, 'w') as f:
            f.write(f">{seq_id}\n{sequence}\n")
        
        # 运行BLASTN
        temp_output = os.path.join(temp_dir, f"{seq_id}_{uuid.uuid4()}.xml")
        
        cmd = [
            blast_path,
            "-db", db_prefix,
            "-query", temp_fasta,
            "-outfmt", "5",  # XML格式输出
            "-out", temp_output,
            "-task", "blastn",  # 使用普通blastn而不是megablast以提高敏感度
            "-word_size", "7",  # 使用较小的word_size以便匹配更短的片段
            "-evalue", "0.01"   # 较严格的E值阈值
        ]
        
        logger.debug(f"BLASTN命令: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"BLASTN运行失败: {result.stderr}")
            # 清理临时文件
            cleanup_temp_files(temp_fasta, temp_output, None)  # 不删除数据库
            return handle_non_triplet_sequence(sequence, seq_id, position='smart')
        
        # 解析BLASTN结果
        best_hit = None
        with open(temp_output) as f:
            blast_records = NCBIXML.parse(f)
            for blast_record in blast_records:
                if blast_record.alignments:
                    best_hit = blast_record.alignments[0]
                    break
        
        if not best_hit:
            logger.warning(f"序列 {seq_id} 未找到BLAST匹配结果")
            # 清理临时文件
            cleanup_temp_files(temp_fasta, temp_output, None)
            return handle_non_triplet_sequence(sequence, seq_id, position='smart')
        
        # 分析比对结果
        hsp = best_hit.hsps[0]  # 最佳高分片段对
        
        # 检查框架和填充位置
        query_start = hsp.query_start  # 查询序列起始位置
        query_end = hsp.query_end  # 查询序列结束位置
        sbjct_start = hsp.sbjct_start  # 匹配序列起始位置
        sbjct_end = hsp.sbjct_end  # 匹配序列结束位置
        
        # 根据比对结果确定填充策略
        remainder = len(sequence) % 3
        padding_needed = 3 - remainder
        
        if query_start > 3:  # 序列前端有超过3个碱基的未比对部分
            padding = 'N' * padding_needed
            padded_seq = padding + sequence
            logger.info(f"序列 {seq_id} 根据BLAST结果在5'端添加 {padding_needed} 个N (起始位置: {query_start})")
        elif len(sequence) - query_end > 3:  # 序列后端有超过3个碱基的未比对部分
            padding = 'N' * padding_needed
            padded_seq = sequence + padding
            logger.info(f"序列 {seq_id} 根据BLAST结果在3'端添加 {padding_needed} 个N (结束位置: {query_end}/{len(sequence)})")
        else:
            # 根据主要比对对象的阅读框架推断
            frame_offset = (sbjct_start - 1) % 3
            current_offset = (query_start - 1) % 3
            
            if frame_offset != current_offset:
                # 需要调整框架使其与参考一致
                if current_offset < frame_offset:
                    padding = 'N' * (frame_offset - current_offset)
                    padded_seq = padding + sequence
                    logger.info(f"序列 {seq_id} 根据参考阅读框在5'端添加 {len(padding)} 个N (框架调整)")
                else:
                    padding = 'N' * (3 - (current_offset - frame_offset))
                    padded_seq = padding + sequence
                    logger.info(f"序列 {seq_id} 根据参考阅读框在5'端添加 {len(padding)} 个N (框架调整)")
            else:
                # 框架一致，根据序列内容分析添加填充
                # 检查是否有明确的起始或终止密码子
                has_start = sequence[:3].upper() == 'ATG'
                has_stop = sequence[-3:].upper() in ['TAA', 'TAG', 'TGA']
                
                if has_start and not has_stop:
                    # 有起始密码子但没有终止密码子，可能是3'端缺失
                    padding = 'N' * padding_needed
                    padded_seq = sequence + padding
                    logger.info(f"序列 {seq_id} 有起始密码子但无终止密码子，在3'端添加 {padding_needed} 个N")
                elif not has_start and has_stop:
                    # 没有起始密码子但有终止密码子，可能是5'端缺失
                    padding = 'N' * padding_needed
                    padded_seq = padding + sequence
                    logger.info(f"序列 {seq_id} 无起始密码子但有终止密码子，在5'端添加 {padding_needed} 个N")
                else:
                    # 默认在序列末尾添加填充
                    padding = 'N' * padding_needed
                    padded_seq = sequence + padding
                    logger.info(f"序列 {seq_id} 框架与参考一致，在3'端添加 {padding_needed} 个N (默认策略)")
        
        # 清理临时文件
        cleanup_temp_files(temp_fasta, temp_output, None)  # 不删除数据库
        
        return padded_seq
        
    except Exception as e:
        logger.error(f"使用BLAST处理序列 {seq_id} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return handle_non_triplet_sequence(sequence, seq_id, position='smart')

def preprocess_cds_sequence(sequence, seq_id, file_path=None, blast_params=None):
    """
    序列预处理流水线，包括质量检查、结构分析和必要的修复
    修改逻辑，只对非3倍数序列使用BLAST
    
    Args:
        sequence: 原始序列
        seq_id: 序列ID
        file_path: 原始CDS文件路径（用于创建BLAST数据库）
        blast_params: BLAST参数字典
    
    Returns:
        处理后的序列和结构信息
    """
    # 性能优化：仅在调试模式下执行详细的质量检查
    if logger.getEffectiveLevel() <= logging.DEBUG:
        issues = assess_sequence_quality(sequence, seq_id)
        if issues:
            logger.warning(f"序列 {seq_id} 存在质量问题: {', '.join(issues)}")
    
    # 分析密码子结构
    structure = analyze_codon_structure(sequence)
    
    # 如果需要，调整序列框架
    if structure['best_frame'] != 0:
        logger.info(f"序列 {seq_id} 最佳阅读框为 {structure['best_frame']}，调整中")
        sequence = sequence[structure['best_frame']:]
    
    # 处理非3倍数长度
    if len(sequence) % 3 != 0:
        # 只有非3倍数序列才需要处理
        logger.info(f"发现非3倍数序列: {seq_id}, 长度 {len(sequence)}")
        
        if blast_params and blast_params.get('use_blast', False) and file_path:
            # 使用BLAST处理
            logger.info(f"使用BLAST处理非3倍数序列 {seq_id}")
            sequence = handle_non_triplet_with_blast(
                sequence, 
                seq_id,
                file_path,
                blast_params.get('blast_path'),
                blast_params.get('temp_dir')
            )
        else:
            # 使用基于序列特征的处理
            context = {
                'has_start': structure.get('has_start', False),
                'has_stop': structure.get('has_stop', False),
                'best_frame': structure.get('best_frame', 0)
            }
            sequence = handle_non_triplet_sequence(sequence, seq_id, position='smart', context=context)
    
    return sequence, structure

def parse_fasta(file_path, duplicate_strategy='longest', blast_params=None):
    """
    从FASTA文件解析序列，使用指定的重复处理策略
    修复物种名称提取和BLAST处理逻辑
    性能优化：批量处理序列
    
    Args:
        file_path: FASTA文件路径
        duplicate_strategy: 重复处理策略
        blast_params: BLAST参数
    
    Returns:
        物种到序列的映射字典
    """
    logger.debug(f"解析 FASTA 文件: {file_path}, 策略: {duplicate_strategy}")
    
    # 如果策略是'alignment_quality'，使用专门的函数处理
    if duplicate_strategy == 'alignment_quality':
        return parse_fasta_with_duplicates(file_path, blast_params=blast_params)
    
    species_seqs = {}
    species_counts = defaultdict(int)
    duplicate_found = False
    original_ids = {}  # 存储原始ID作为参考
    
    try:
        # 性能优化：使用生成器而不是一次性加载所有记录
        record_count = 0
        batch_size = 100  # 每批处理的序列数
        batch = []
        
        for record in SeqIO.parse(file_path, "fasta"):
            batch.append(record)
            record_count += 1
            
            # 批量处理
            if len(batch) >= batch_size:
                process_fasta_batch(batch, species_seqs, species_counts, original_ids, duplicate_strategy, file_path, blast_params)
                batch = []
                logger.debug(f"已处理 {record_count} 条序列")
        
        # 处理最后一批
        if batch:
            process_fasta_batch(batch, species_seqs, species_counts, original_ids, duplicate_strategy, file_path, blast_params)
            
        logger.debug(f"在 {file_path} 中找到并处理了 {record_count} 个序列")
        
        # 检查是否有重复
        duplicate_found = any(count > 1 for count in species_counts.values())
        if duplicate_found:
            logger.info(f"在 {file_path} 中发现重复物种: {dict(species_counts)}")
        
        return species_seqs
    except Exception as e:
        logger.error(f"解析 FASTA 文件 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def process_fasta_batch(batch, species_seqs, species_counts, original_ids, duplicate_strategy, file_path, blast_params):
    """
    批量处理FASTA记录
    
    Args:
        batch: 要处理的记录批次
        species_seqs: 物种到序列的映射
        species_counts: 物种计数字典
        original_ids: 原始ID映射
        duplicate_strategy: 重复处理策略
        file_path: 原始文件路径
        blast_params: BLAST参数
    """
    for record in batch:
        # 从头部提取物种名
        species = extract_species_name(record.description)
        sequence = str(record.seq)
        
        # 预处理序列
        sequence, _ = preprocess_cds_sequence(sequence, record.id, file_path, blast_params)
        
        # 检查序列是否是有效的编码序列
        is_valid, message = check_coding_sequence(sequence)
        if not is_valid:
            logger.warning(f"序列 {record.id} ({species}) 可能不是有效的编码序列: {message}")
        
        # 检查这个物种是否已存在
        if species in species_seqs:
            species_counts[species] += 1
            
            if duplicate_strategy == 'longest':
                # 保留最长的序列
                if len(sequence) > len(species_seqs[species]):
                    logger.debug(f"用更长的序列替换 {species} 的较短序列 (长度: {len(species_seqs[species])} → {len(sequence)})")
                    species_seqs[species] = sequence
                    original_ids[species] = record.id
            
            elif duplicate_strategy == 'first':
                # 保留第一个序列，忽略这个
                logger.debug(f"忽略 {species} 的重复序列")
                continue
                
            elif duplicate_strategy == 'rename':
                # 通过添加后缀重命名重复物种
                new_species = f"{species}_{species_counts[species]}"
                logger.debug(f"将重复的物种从 {species} 重命名为 {new_species}")
                species_seqs[new_species] = sequence
                original_ids[new_species] = record.id
        else:
            # 第一次见到这个物种
            species_seqs[species] = sequence
            species_counts[species] = 1
            original_ids[species] = record.id

def parse_fasta_with_duplicates(file_path, blast_params=None):
    """
    解析FASTA文件，保留所有重复的序列和原始ID
    用于alignment_quality策略
    性能优化：批量处理序列
    
    Args:
        file_path: FASTA文件路径
        blast_params: BLAST参数
    
    Returns:
        物种到ID-序列字典的映射
    """
    species_to_seqs = defaultdict(dict)
    
    try:
        # 性能优化：使用生成器并批量处理
        record_count = 0
        batch_size = 100
        batch = []
        
        for record in SeqIO.parse(file_path, "fasta"):
            batch.append(record)
            record_count += 1
            
            if len(batch) >= batch_size:
                process_duplicate_batch(batch, species_to_seqs, file_path, blast_params)
                batch = []
        
        # 处理最后一批
        if batch:
            process_duplicate_batch(batch, species_to_seqs, file_path, blast_params)
            
        logger.debug(f"在 {file_path} 中找到并处理了 {record_count} 个序列")
        
        # 记录发现的重复
        duplicates = {species: len(ids) for species, ids in species_to_seqs.items() if len(ids) > 1}
        if duplicates:
            logger.info(f"在 {file_path} 中发现重复物种: {duplicates}")
            
        return species_to_seqs
    except Exception as e:
        logger.error(f"解析 FASTA 文件 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def process_duplicate_batch(batch, species_to_seqs, file_path, blast_params):
    """
    处理包含可能重复的FASTA记录批次
    
    Args:
        batch: 要处理的记录批次
        species_to_seqs: 物种到ID-序列字典的映射
        file_path: 原始文件路径
        blast_params: BLAST参数
    """
    for record in batch:
        # 从头部提取物种名
        species = extract_species_name(record.description)
        record_id = record.id
        sequence = str(record.seq)
        
        # 预处理序列
        sequence, _ = preprocess_cds_sequence(sequence, record_id, file_path, blast_params)
        
        # 使用原始ID存储以跟踪它们
        species_to_seqs[species][record_id] = sequence

def dna_to_protein_for_alignment(dna_file, protein_file):
    """
    转换DNA序列为蛋白质用于比对引导
    改进ID处理，避免重复ID问题
    性能优化：批量处理序列
    
    Args:
        dna_file: DNA FASTA文件路径
        protein_file: 输出蛋白质FASTA文件路径
    
    Returns:
        bool: 成功状态
    """
    try:
        # 记录原始序列ID和它们的索引，避免重复ID问题
        id_mapping = {}  # 存储新ID到原始ID的映射
        
        # 性能优化：使用生成器和批处理
        batch_size = 100
        records_batch = []
        batch_count = 0
        
        with open(protein_file, 'w') as out:
            for i, record in enumerate(SeqIO.parse(dna_file, "fasta")):
                records_batch.append((i, record))
                
                if len(records_batch) >= batch_size:
                    process_dna_to_protein_batch(records_batch, out, id_mapping)
                    batch_count += 1
                    logger.debug(f"DNA到蛋白质批处理 {batch_count}，处理了 {len(records_batch)} 条序列")
                    records_batch = []
            
            # 处理最后一批
            if records_batch:
                process_dna_to_protein_batch(records_batch, out, id_mapping)
        
        # 保存ID映射关系以便后续处理
        mapping_file = protein_file + ".id_map.json"
        with open(mapping_file, 'w') as f:
            json.dump(id_mapping, f)
            
        logger.debug(f"成功将DNA序列转换为蛋白质: {dna_file} -> {protein_file}")
        return True
    except Exception as e:
        logger.error(f"DNA转蛋白质时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def process_dna_to_protein_batch(records_batch, output_file, id_mapping):
    """
    批量处理DNA到蛋白质的转换
    
    Args:
        records_batch: (索引, 记录)元组的列表
        output_file: 输出文件句柄
        id_mapping: ID映射字典
    """
    for idx, record in records_batch:
        seq = str(record.seq).upper()
        original_id = record.id
        
        # 生成唯一ID
        unique_id = f"{original_id}_prot_{idx}"
        id_mapping[unique_id] = original_id
        
        # 确保序列长度是3的倍数
        remainder = len(seq) % 3
        if remainder > 0:
            seq = seq + "N" * (3 - remainder)
        
        # 翻译为蛋白质
        protein = improved_translate_with_gaps(seq)
        output_file.write(f">{unique_id}\n{protein}\n")

def map_protein_alignment_to_dna(original_dna_file, protein_alignment_file, output_dna_alignment_file):
    """
    将蛋白质比对映射回核苷酸序列，保持密码子结构
    改进ID处理，解决重复ID问题
    性能优化：批量处理序列
    
    Args:
        original_dna_file: 原始未比对DNA序列文件路径
        protein_alignment_file: 比对的蛋白质序列文件路径
        output_dna_alignment_file: 输出比对DNA序列文件路径
    
    Returns:
        bool: 成功状态
    """
    try:
        # 检查ID映射文件是否存在
        mapping_file = protein_alignment_file + ".id_map.json"
        id_mapping = {}
        
        if os.path.exists(mapping_file):
            try:
                with open(mapping_file, 'r') as f:
                    id_mapping = json.load(f)
                logger.debug(f"已加载序列ID映射，包含 {len(id_mapping)} 个映射关系")
            except Exception as e:
                logger.warning(f"无法加载ID映射文件 {mapping_file}: {str(e)}")
        
        # 性能优化：分批加载原始DNA序列
        dna_records = {}
        seq_count = 0
        batch_size = 100
        
        for record in SeqIO.parse(original_dna_file, "fasta"):
            dna_records[record.id] = str(record.seq).upper()
            seq_count += 1
            
            # 定期释放内存，避免大文件处理时内存占用过高
            if seq_count % 1000 == 0:
                gc.collect()  # 强制垃圾回收
        
        logger.debug(f"从 {original_dna_file} 加载了 {seq_count} 条DNA序列")
        
        # 加载比对的蛋白质序列
        aligned_proteins = {}
        for record in SeqIO.parse(protein_alignment_file, "fasta"):
            # 使用映射回原始ID，如果存在
            record_id = record.id
            target_id = record_id
            
            if record_id in id_mapping:
                target_id = id_mapping[record_id]
                logger.debug(f"将蛋白质ID {record_id} 映射回原始ID {target_id}")
            
            aligned_proteins[target_id] = str(record.seq)
        
        # 性能优化：分批映射蛋白质回DNA
        aligned_dna = {}
        missing_ids = []
        batch_count = 0
        protein_batch = []
        
        for seq_id, aligned_protein in aligned_proteins.items():
            protein_batch.append((seq_id, aligned_protein))
            
            if len(protein_batch) >= batch_size:
                process_protein_to_dna_batch(protein_batch, dna_records, aligned_dna, missing_ids)
                batch_count += 1
                logger.debug(f"蛋白质到DNA映射批处理 {batch_count}，处理了 {len(protein_batch)} 条序列")
                protein_batch = []
        
        # 处理最后一批
        if protein_batch:
            process_protein_to_dna_batch(protein_batch, dna_records, aligned_dna, missing_ids)
        
        # 如果有缺失的ID，尝试处理
        if missing_ids:
            process_missing_ids(missing_ids, aligned_proteins, dna_records, aligned_dna)
        
        # 写入比对后的DNA序列
        with open(output_dna_alignment_file, 'w') as out:
            for seq_id, seq in aligned_dna.items():
                out.write(f">{seq_id}\n{seq}\n")
        
        logger.info(f"成功将蛋白质比对映射回DNA: {len(aligned_dna)} 条序列")
        return True
    except Exception as e:
        logger.error(f"映射蛋白质比对到DNA时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def process_protein_to_dna_batch(protein_batch, dna_records, aligned_dna, missing_ids):
    """
    批量处理蛋白质到DNA的映射
    
    Args:
        protein_batch: (ID, 蛋白质序列)元组的列表
        dna_records: DNA序列字典
        aligned_dna: 输出的比对DNA字典
        missing_ids: 缺失ID的列表
    """
    for seq_id, aligned_protein in protein_batch:
        if seq_id not in dna_records:
            # 记录并稍后尝试找到备选ID
            missing_ids.append(seq_id)
            continue
        
        original_dna = dna_records[seq_id]
        
        # 确保DNA长度是3的倍数
        remainder = len(original_dna) % 3
        if remainder > 0:
            original_dna = original_dna + "N" * (3 - remainder)
        
        # 构建比对后的DNA序列
        aligned_dna_seq = ""
        dna_pos = 0
        
        for aa in aligned_protein:
            if aa == '-':
                # 蛋白质中的缺口对应DNA中的三个缺口
                aligned_dna_seq += "---"
            else:
                # 常规氨基酸对应DNA中的下一个密码子
                if dna_pos + 3 <= len(original_dna):
                    codon = original_dna[dna_pos:dna_pos+3]
                    aligned_dna_seq += codon
                    dna_pos += 3
                else:
                    # 处理DNA序列长度不足的情况
                    logger.warning(f"序列 {seq_id} DNA长度不足，在位置 {dna_pos} 添加 NNN")
                    aligned_dna_seq += "NNN"
        
        aligned_dna[seq_id] = aligned_dna_seq

def process_missing_ids(missing_ids, aligned_proteins, dna_records, aligned_dna):
    """
    处理在DNA记录中找不到的ID
    
    Args:
        missing_ids: 缺失ID列表
        aligned_proteins: 比对蛋白质序列字典
        dna_records: DNA序列字典
        aligned_dna: 输出的比对DNA字典
    """
    for missing_id in missing_ids:
        found = False
        
        # 检查是否是蛋白质ID而非原始ID
        if "_prot_" in missing_id:
            base_id = missing_id.split("_prot_")[0]
            if base_id in dna_records:
                logger.info(f"找到ID的匹配: {missing_id} -> {base_id}")
                
                # 使用上面的相同逻辑
                original_dna = dna_records[base_id]
                remainder = len(original_dna) % 3
                if remainder > 0:
                    original_dna = original_dna + "N" * (3 - remainder)
                
                aligned_dna_seq = ""
                dna_pos = 0
                
                for aa in aligned_proteins[missing_id]:
                    if aa == '-':
                        aligned_dna_seq += "---"
                    else:
                        if dna_pos + 3 <= len(original_dna):
                            codon = original_dna[dna_pos:dna_pos+3]
                            aligned_dna_seq += codon
                            dna_pos += 3
                        else:
                            aligned_dna_seq += "NNN"
                
                aligned_dna[base_id] = aligned_dna_seq
                found = True
        
        # 如果上面方法失败，尝试前缀匹配
        if not found:
            for orig_id in dna_records.keys():
                # 检查前缀匹配
                if orig_id.startswith(missing_id.split('_prot_')[0]):
                    logger.info(f"找到ID的前缀匹配: {missing_id} -> {orig_id}")
                    
                    # 重复上面的映射过程
                    original_dna = dna_records[orig_id]
                    remainder = len(original_dna) % 3
                    if remainder > 0:
                        original_dna = original_dna + "N" * (3 - remainder)
                    
                    aligned_dna_seq = ""
                    dna_pos = 0
                    
                    for aa in aligned_proteins[missing_id]:
                        if aa == '-':
                            aligned_dna_seq += "---"
                        else:
                            if dna_pos + 3 <= len(original_dna):
                                codon = original_dna[dna_pos:dna_pos+3]
                                aligned_dna_seq += codon
                                dna_pos += 3
                            else:
                                aligned_dna_seq += "NNN"
                    
                    aligned_dna[orig_id] = aligned_dna_seq
                    found = True
                    break
        
        if not found:
            logger.warning(f"无法找到序列ID {missing_id} 的匹配，将从输出中排除")

def improved_translate_with_gaps(nucleotide_seq):
    """
    优化的翻译函数，更好地处理含有缺口的密码子
    性能优化：使用缓存减少重复计算
    
    Args:
        nucleotide_seq: 核苷酸序列
    
    Returns:
        翻译后的氨基酸序列
    """
    if not nucleotide_seq:
        return ""
    
    # 确保序列长度是3的倍数
    original_len = len(nucleotide_seq)
    needed_padding = 0
    
    if original_len % 3 != 0:
        needed_padding = 3 - (original_len % 3)
        nucleotide_seq += 'N' * needed_padding
        logger.debug(f"序列长度 {original_len} 不是 3 的倍数，添加 {needed_padding} 个N进行翻译")
    
    # 性能优化：预分配列表大小
    amino_acids = [""] * (len(nucleotide_seq) // 3)
    
    for i in range(0, len(nucleotide_seq), 3):
        pos = i // 3
        codon = nucleotide_seq[i:i+3]
        
        if '-' in codon:
            # 计算缺口数量
            gap_count = codon.count('-')
            
            # 应用不同的翻译规则
            if gap_count == 3:
                # 完全缺口密码子
                amino_acids[pos] = '-'
            elif gap_count == 1 and codon[1:3] == '--':
                # 第一个碱基有但后两个是缺口，可能是框架位移
                amino_acids[pos] = 'X'
            elif gap_count == 1 and codon[0] == '-' and codon[2] == '-':
                # 仅中间位置有碱基的情况
                amino_acids[pos] = 'X'
            elif gap_count == 2 and codon[0] == '-':
                # 第一个位置是缺口
                amino_acids[pos] = 'X'
            elif gap_count == 1:
                # 仅一个缺口，尝试推断氨基酸
                non_gap_pos = [i for i, c in enumerate(codon) if c != '-']
                if len(non_gap_pos) == 2 and all(codon[p] in 'ACGT' for p in non_gap_pos):
                    # 检查是否可以唯一确定氨基酸
                    possible_aas = set()
                    for base in 'ACGT':
                        temp_codon = list(codon)
                        gap_pos = codon.index('-')
                        temp_codon[gap_pos] = base
                        temp_codon = ''.join(temp_codon)
                        possible_aas.add(safe_translate_codon_cached(temp_codon))
                    
                    if len(possible_aas) == 1:
                        amino_acids[pos] = next(iter(possible_aas))
                    else:
                        # 无法唯一确定，使用模糊氨基酸符号
                        amino_acids[pos] = 'X'
                else:
                    amino_acids[pos] = 'X'
            else:
                amino_acids[pos] = 'X'
        elif 'N' in codon:
            # 处理含N的密码子
            n_count = codon.count('N')
            if n_count == 3:
                amino_acids[pos] = 'X'
            else:
                # 尝试推断可能的氨基酸
                possible_aas = set()
                n_positions = [i for i, c in enumerate(codon) if c == 'N']
                
                # 生成所有可能的替换组合
                replacements = ['A', 'C', 'G', 'T']
                for combo in itertools.product(replacements, repeat=len(n_positions)):
                    temp_codon = list(codon)
                    for idx, pos in enumerate(n_positions):
                        temp_codon[pos] = combo[idx]
                    possible_aas.add(safe_translate_codon_cached(''.join(temp_codon)))
                
                if len(possible_aas) == 1:
                    amino_acids[pos] = next(iter(possible_aas))
                else:
                    amino_acids[pos] = 'X'
        else:
            try:
                # 使用缓存版本提高性能
                amino_acid = safe_translate_codon_cached(codon)
                amino_acids[pos] = amino_acid
            except Exception as e:
                logger.warning(f"翻译密码子 '{codon}' 时出错: {str(e)}")
                amino_acids[pos] = 'X'
    
    # 返回对应于原始序列的部分
    result = ''.join(amino_acids)
    if needed_padding > 0:
        # 去除由于填充添加的部分
        if needed_padding < 3:  # 只有当填充少于1个完整密码子时才需调整
            result = result[:-1]
    
    return result

def safe_translate_codon(codon):
    """
    安全地将密码子翻译为氨基酸，处理模糊碱基
    
    Args:
        codon: 三字母核苷酸密码子
        
    Returns:
        单字母氨基酸或'X'表示未知
    """
    # 转为大写并移除空白
    codon = codon.upper().strip()
    
    # 检查密码子是否存在于查找表中
    if codon in GENETIC_CODE:
        return GENETIC_CODE[codon]
    
    # 如果密码子包含N或其他模糊核苷酸
    if 'N' in codon or any(base not in 'ACGT' for base in codon):
        return 'X'
        
    # 这种情况不应该出现在有效的DNA序列中
    logger.warning(f"未知的密码子: '{codon}'")
    return 'X'

def run_mafft(input_file, output_file, mafft_path, codon_aware=False):
    """
    运行MAFFT进行序列比对，改进密码子感知处理
    性能优化：调整MAFFT参数以优化性能
    
    Args:
        input_file: 输入文件路径
        output_file: 输出文件路径
        mafft_path: MAFFT可执行文件路径
        codon_aware: 是否保持密码子结构
    """
    start_time = time.time()
    
    # 确保mafft_path存在且可执行
    if not os.path.exists(mafft_path):
        logger.error(f"MAFFT 可执行文件未在此路径找到: {mafft_path}")
        return False, "可执行文件未找到"
    
    if not os.access(mafft_path, os.X_OK):
        logger.error(f"MAFFT 可执行文件没有执行权限: {mafft_path}")
        return False, "可执行文件没有执行权限"
    
    # 如果请求密码子感知比对，使用蛋白质引导方法
    if codon_aware:
        logger.info("使用密码子感知方式运行MAFFT")
        temp_dir = os.path.dirname(output_file)
        
        # 1. 将序列转换为氨基酸
        protein_input = os.path.join(temp_dir, os.path.basename(input_file).replace('.fasta', '.prot.fasta'))
        dna_to_protein_for_alignment(input_file, protein_input)
        
        # 2. 比对蛋白质 - 性能优化：使用快速算法
        protein_output = os.path.join(temp_dir, os.path.basename(output_file).replace('.fasta', '.prot.aln'))
        
        # 检查序列数量，为大数据集使用快速算法
        seq_count = sum(1 for _ in SeqIO.parse(protein_input, "fasta"))
        
        if seq_count > 200:
            # 使用FFT-NS-2算法处理大型数据集
            logger.info(f"检测到大型数据集（{seq_count}个序列），使用FFT-NS-2算法")
            protein_cmd = [mafft_path, '--quiet', '--retree', '2', '--maxiterate', '0', protein_input]
        else:
            # 使用较精确的算法处理小型数据集
            protein_cmd = [mafft_path, '--quiet', '--localpair', '--maxiterate', '1000', protein_input]
        
        logger.debug(f"运行蛋白质对齐: {' '.join(protein_cmd)}")
        
        try:
            with open(protein_output, 'w') as outfile:
                protein_result = subprocess.run(protein_cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
                
            if protein_result.returncode != 0:
                logger.error(f"MAFFT蛋白质对齐错误: {protein_result.stderr}")
                return False, "蛋白质对齐失败"
                
            # 3. 将蛋白质比对映射回DNA
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
        # 普通MAFFT运行，带有适用于DNA的通用参数
        cmd = [mafft_path]
        
        # 性能优化：根据序列数量自动选择算法
        seq_count = sum(1 for _ in SeqIO.parse(input_file, "fasta"))
        
        if seq_count > 200:
            # 大型数据集：使用快速算法
            cmd.extend(['--retree', '2', '--maxiterate', '0'])
            logger.info(f"检测到大型数据集（{seq_count}个序列），使用FFT-NS-2算法")
        else:
            # 小型数据集：使用较精确的算法
            cmd.append('--auto')
            
        cmd.extend(['--quiet', input_file])  # 重定向输出
        
        logger.info(f"运行 MAFFT: {' '.join(cmd)}")
        
        try:
            with open(output_file, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"MAFFT 错误 (返回码 {result.returncode}): {result.stderr}")
                return False, f"返回码 {result.returncode}: {result.stderr}"
                
            execution_time = time.time() - start_time
            logger.info(f"MAFFT 成功完成, 耗时 {execution_time:.1f} 秒")
            
            # 验证输出
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

def run_muscle(input_file, output_file, muscle_path, codon_aware=False):
    """
    运行MUSCLE进行序列比对，改进密码子感知处理
    性能优化：调整MUSCLE参数
    
    Args:
        input_file: 输入文件路径
        output_file: 输出文件路径
        muscle_path: MUSCLE可执行文件路径
        codon_aware: 是否保持密码子结构
    """
    start_time = time.time()
    
    # 确保muscle_path存在且可执行
    if not os.path.exists(muscle_path):
        logger.error(f"MUSCLE 可执行文件未在此路径找到: {muscle_path}")
        return False, "可执行文件未找到"
    
    if not os.access(muscle_path, os.X_OK):
        logger.error(f"MUSCLE 可执行文件没有执行权限: {muscle_path}")
        return False, "可执行文件没有执行权限"
    
    # 如果请求密码子感知比对，使用蛋白质引导方法
    if codon_aware:
        logger.info("使用密码子感知方式运行MUSCLE")
        temp_dir = os.path.dirname(output_file)
        
        # 1. 将序列转换为氨基酸
        protein_input = os.path.join(temp_dir, os.path.basename(input_file).replace('.fasta', '.prot.fasta'))
        dna_to_protein_for_alignment(input_file, protein_input)
        
        # 2. 比对蛋白质
        protein_output = os.path.join(temp_dir, os.path.basename(output_file).replace('.fasta', '.prot.aln'))
        
        # 检查MUSCLE版本并设置参数
        muscle_version = get_muscle_version(muscle_path)
        
        if muscle_version >= 5:
            # MUSCLE v5+使用新的参数集
            protein_cmd = [muscle_path, '-align', protein_input, '-output', protein_output]
            logger.debug(f"检测到MUSCLE v5+，使用新参数格式")
        else:
            # MUSCLE v3/v4使用传统参数
            protein_cmd = [muscle_path, '-in', protein_input, '-out', protein_output]
        
        logger.debug(f"运行蛋白质对齐: {' '.join(protein_cmd)}")
        protein_result = subprocess.run(protein_cmd, capture_output=True, text=True)
        
        if protein_result.returncode != 0:
            logger.error(f"MUSCLE蛋白质对齐错误: {protein_result.stderr}")
            return False, "蛋白质对齐失败"
            
        # 3. 将蛋白质比对映射回DNA
        success = map_protein_alignment_to_dna(input_file, protein_output, output_file)
        if success:
            logger.info(f"成功完成基于蛋白质引导的密码子感知对齐，输出至 {output_file}")
            return True, output_file
        else:
            return False, "无法将蛋白质对齐映射回DNA"
    else:
        # 普通MUSCLE运行，自动检测版本和调整参数
        muscle_version = get_muscle_version(muscle_path)
        
        if muscle_version >= 5:
            # MUSCLE v5+
            cmd = [muscle_path, '-align', input_file, '-output', output_file]
        else:
            # MUSCLE v3/v4
            cmd = [muscle_path, '-in', input_file, '-out', output_file]
            
        logger.info(f"运行 MUSCLE: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"MUSCLE 错误 (返回码 {result.returncode}): {result.stderr}")
                return False, f"返回码 {result.returncode}: {result.stderr}"
                
            execution_time = time.time() - start_time
            logger.info(f"MUSCLE 成功完成, 耗时 {execution_time:.1f} 秒")
            
            # 验证输出
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

def get_muscle_version(muscle_path):
    """
    检测MUSCLE的版本
    
    Args:
        muscle_path: MUSCLE可执行文件路径
        
    Returns:
        版本号(主版本)，例如3, 4, 5
    """
    try:
        # 尝试运行带-version参数
        result = subprocess.run([muscle_path, '-version'], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout:
            version_text = result.stdout
        else:
            # 可能是MUSCLE v5，尝试不带参数运行
            result = subprocess.run([muscle_path], capture_output=True, text=True)
            version_text = result.stderr if result.stderr else result.stdout
        
        # 从输出中提取版本号
        if 'MUSCLE v5' in version_text or 'MUSCLE v 5' in version_text:
            return 5
        elif 'MUSCLE v3' in version_text or 'MUSCLE v 3' in version_text:
            return 3
        elif 'MUSCLE v4' in version_text or 'MUSCLE v 4' in version_text:
            return 4
        else:
            # 默认假设是v3
            return 3
    except Exception as e:
        logger.debug(f"检测MUSCLE版本时出错: {str(e)}")
        return 3  # 默认假设是v3

def run_prank(input_file, output_prefix, prank_path, codon_aware=True, f=0.2, 
             gaprate=None, gapext=None, use_logs=False, penalize_terminal_gaps=False):
    """
    运行PRANK进行多序列比对，带有优化参数
    性能优化：添加线程和内存优化参数
    
    Args:
        input_file: 输入文件路径
        output_prefix: 输出文件前缀
        prank_path: PRANK可执行文件路径
        codon_aware: 是否使用密码子感知比对
        f: 控制插入打开概率的参数 (默认: 0.2)
        gaprate: 缺口打开率 (如果为None，使用PRANK基于数据类型的默认值)
        gapext: 缺口扩展概率 (如果为None，使用PRANK基于数据类型的默认值)
        use_logs: 对大型数据集使用对数计算
        penalize_terminal_gaps: 是否正常惩罚终端缺口
    """
    start_time = time.time()
    
    # 确保prank_path存在且可执行
    if not os.path.exists(prank_path):
        logger.error(f"PRANK 可执行文件未在此路径找到: {prank_path}")
        return False, "可执行文件未找到"
    
    if not os.access(prank_path, os.X_OK):
        logger.error(f"PRANK 可执行文件没有执行权限: {prank_path}")
        return False, "可执行文件没有执行权限"
    
    cmd = [prank_path, '-d=' + input_file, '-o=' + output_prefix]
    
    # 如果请求，添加密码子感知比对
    if codon_aware:
        cmd.append('-codon')
    
    # 添加优化参数
    cmd.append(f'-f={f}')  # 控制插入打开概率
    
    # 仅在指定时设置缺口参数，否则使用PRANK默认值
    if gaprate is not None:
        cmd.append(f'-gaprate={gaprate}')
    
    if gapext is not None:
        cmd.append(f'-gapext={gapext}')
    
    # 添加额外的优化参数
    cmd.append('+F')  # 强制总是跳过插入
    cmd.append('-prunetree')  # 修剪没有序列数据的指导树分支
    cmd.append('-shortnames')  # 在第一个空格处截断名称以便更好地处理
    
    # 性能优化：添加线程参数
    # 获取系统可用的CPU核心数
    cpu_count = multiprocessing.cpu_count()
    thread_count = max(1, min(cpu_count - 1, 4))  # 使用最多4个线程或系统最大值-1
    cmd.append(f'-threads={thread_count}')
    
    if use_logs:
        cmd.append('-uselogs')  # 对大型数据集使用对数计算
        
    if penalize_terminal_gaps:
        cmd.append('-termgap')  # 正常惩罚终端缺口
    
    # 总是运行一次以保持一致性和速度
    cmd.append('-once')
    
    logger.info(f"运行 PRANK: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"PRANK 错误 (返回码 {result.returncode}): {result.stderr}")
            return False, f"返回码 {result.returncode}: {result.stderr}"
            
        execution_time = time.time() - start_time
        logger.info(f"PRANK 成功完成, 耗时 {execution_time:.1f} 秒")
        
        # 尝试找到不同扩展名的比对输出
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
            # 检查PRANK是否创建了任何文件
            dir_path = os.path.dirname(output_prefix)
            base_name = os.path.basename(output_prefix)
            all_files = [f for f in os.listdir(dir_path) if f.startswith(base_name)]
            
            if all_files:
                # 查找最可能的输出文件（例如，最大的.fas文件）
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
            
        # 检查文件是否包含有效比对
        valid, frame_info = check_alignment_validity(found_output)
        if valid:
            # 如果不同，将找到的输出复制或重命名为预期的输出名称
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

def check_alignment_validity(alignment_file):
    """
    检查比对文件是否包含有效的比对序列
    性能优化：快速验证大型文件
    
    Args:
        alignment_file: 比对文件路径
        
    Returns:
        tuple: (is_valid, frame_info)
    """
    try:
        # 检查文件是否存在且非空
        if not os.path.exists(alignment_file) or os.path.getsize(alignment_file) == 0:
            logger.warning(f"比对文件不存在或为空: {alignment_file}")
            return False, "文件不存在或为空"
        
        # 性能优化：对于大文件，只检查前100条序列
        seq_iterator = SeqIO.parse(alignment_file, "fasta")
        records = []
        
        # 限制读取的记录数量
        for i, record in enumerate(seq_iterator):
            records.append(record)
            if i >= 99:  # 读取前100条序列
                break
        
        # 检查是否找到任何序列
        if not records:
            logger.warning(f"比对文件没有包含序列: {alignment_file}")
            return False, "文件不包含序列"
        
        # 检查序列长度是否一致（应该是比对的）
        alignment_length = len(records[0].seq)
        if alignment_length == 0:
            logger.warning(f"比对序列长度为零: {alignment_file}")
            return False, "序列长度为零"
        
        # 检查比对长度是否是3的倍数（对密码子结构很重要）
        is_codon_aligned = (alignment_length % 3 == 0)
        if not is_codon_aligned:
            logger.warning(f"比对长度 {alignment_length} 不是 3 的倍数，可能不保持密码子结构")
        
        # 性能优化：只检查序列长度，不检查所有序列内容
        # 检查所有序列是否长度相同
        unequal_lengths = False
        for record in records:
            if len(record.seq) != alignment_length:
                logger.warning(f"比对序列长度不一致: {record.id} 的长度为 {len(record.seq)}, 而不是预期的 {alignment_length}")
                unequal_lengths = True
                break  # 一旦发现不一致就可以停止
        
        if unequal_lengths:
            return False, "序列长度不一致"
        
        # 性能优化：仅对少量记录进行详细检查
        # 选择最多10条记录进行详细验证
        sample_records = records[:min(10, len(records))]
        
        # 检查基本内容
        # 确保文件不是只有缺口或无效字符
        valid_chars = set('ACGTRYKMSWBDHVN-')  # DNA + 模糊碱基代码 + 缺口
        invalid_seqs = []
        
        for record in sample_records:
            seq_str = str(record.seq).upper()
            # 检查序列是否包含至少一些有效的DNA字符
            if not any(c in 'ACGT' for c in seq_str):
                logger.warning(f"序列 {record.id} 不包含有效的DNA字符")
                invalid_seqs.append(record.id)
                continue
                
            # 计算无效字符比例
            invalid_chars = set(seq_str) - valid_chars
            if invalid_chars and sum(seq_str.count(c) for c in invalid_chars) > len(seq_str) * 0.1:  # >10% 无效字符
                logger.warning(f"序列 {record.id} 包含过多无效字符: {invalid_chars}")
                invalid_seqs.append(record.id)
        
        if invalid_seqs:
            return False, f"{len(invalid_seqs)} 条序列包含无效字符"
        
        # 快速估算缺口百分比
        total_gaps = sum(str(record.seq).count('-') for record in sample_records)
        total_chars = sum(len(record.seq) for record in sample_records)
        gap_percentage = (total_gaps / total_chars) * 100 if total_chars > 0 else 0
        
        # 高缺口百分比可能表示问题，但这只是信息性的
        if gap_percentage > 50:
            logger.warning(f"比对包含高比例的缺口: {gap_percentage:.1f}%")
        
        # 通过分析缺口模式检查密码子结构保存（只在必要时进行）
        codon_structure_preserved = True
        if is_codon_aligned and logger.getEffectiveLevel() <= logging.DEBUG:
            for record in sample_records:
                seq_str = str(record.seq)
                # 只检查序列的子集以提高性能
                positions_to_check = min(len(seq_str), 300)  # 检查最多300个位置
                
                for i in range(0, positions_to_check, 3):
                    if i+2 < len(seq_str):
                        codon = seq_str[i:i+3]
                        gap_count = codon.count('-')
                        if gap_count > 0 and gap_count < 3:
                            # 找到一个密码子，其中一些但不是所有位置都是缺口
                            # 这可能会破坏密码子结构
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

def run_aligner(input_file, output_file, aligner_params):
    """
    运行指定的比对工具进行序列比对
    
    Args:
        input_file: 输入文件路径
        output_file: 输出文件路径或前缀
        aligner_params: 比对参数字典
        
    Returns:
        (success, output_file) 元组
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

def calculate_similarity_score(aa1, aa2):
    """
    基于相似性组计算两个氨基酸之间的相似性评分
    
    Args:
        aa1: 第一个氨基酸
        aa2: 第二个氨基酸
        
    Returns:
        如果相同则为1.0，如果在同一组中则为0.5，如果在不同组中则为0.0
    """
    if aa1 == aa2:
        return 1.0
    
    for group in AA_SIMILARITY_GROUPS:
        if aa1 in group and aa2 in group:
            return 0.5
    
    return 0.0

def evaluate_alignment_quality(aligned_sequences, seq_id, is_protein=False):
    """
    改进的比对质量评估，评估序列与其他序列相比的比对质量
    分数越低代表质量越好
    性能优化：使用numpy加速计算
    
    Args:
        aligned_sequences: 所有比对序列的字典
        seq_id: 要评估的序列ID
        is_protein: 序列是否是蛋白质序列
        
    Returns:
        质量分数（越低越好）
    """
    if seq_id not in aligned_sequences:
        logger.warning(f"序列ID {seq_id} 未在比对序列中找到")
        return float('inf')  # 如果找不到序列，返回最差的可能分数
    
    sequence = aligned_sequences[seq_id]
    # 过滤掉正在评估的序列
    other_seqs = [seq for id, seq in aligned_sequences.items() if id != seq_id]
    
    if not other_seqs:
        logger.debug(f"序列 {seq_id} 没有其他序列可比较")
        return 0  # 没有其他序列可比较
    
    # 使用权重初始化质量指标
    weights = {
        'gap_ratio': 0.3,           # 惩罚缺口
        'conservation': 0.4,        # 奖励与其他序列的保守性
        'terminal_gaps': 0.1,       # 惩罚终端缺口的程度低于内部缺口
        'consecutive_gaps': 0.2     # 惩罚连续缺口（插入/删除）
    }
    
    # 计算缺口指标
    total_gaps = sequence.count('-')
    seq_length = len(sequence)
    gap_ratio = total_gaps / seq_length if seq_length > 0 else 1.0
    
    # 识别终端缺口与内部缺口
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
    
    # 计算连续缺口块（插入/删除）
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
    # 归一化到0-1范围（更多块更糟，但比更少的大块好）
    consecutive_gap_score = min(consecutive_gap_score, 1.0)
    
    # 性能优化：使用numpy加速保守性计算
    sequence_array = np.array(list(sequence))
    other_arrays = [np.array(list(seq)) for seq in other_seqs]
    
    # 只分析非缺口位置
    non_gap_positions = np.where(sequence_array != '-')[0]
    
    # 如果序列只有缺口，返回最差分数
    if len(non_gap_positions) == 0:
        return 1.0
    
    conservation_scores = []
    
    # 批量处理保守性计算
    batch_size = 100  # 每批处理的位置数
    for batch_start in range(0, len(non_gap_positions), batch_size):
        batch_end = min(batch_start + batch_size, len(non_gap_positions))
        batch_positions = non_gap_positions[batch_start:batch_end]
        
        for pos in batch_positions:
            # 获取当前位置的字符
            current_char = sequence_array[pos]
            
            # 计算与其他序列在这个位置的匹配情况
            match_scores = []
            
            for other_array in other_arrays:
                if pos >= len(other_array) or other_array[pos] == '-':
                    continue
                    
                other_char = other_array[pos]
                
                if is_protein:
                    # 对于蛋白质，使用相似性评分
                    match_scores.append(calculate_similarity_score(current_char, other_char))
                else:
                    # 对于核苷酸，仅精确匹配
                    match_scores.append(1.0 if current_char == other_char else 0.0)
            
            if match_scores:
                # 这个位置的平均保守性
                conservation_scores.append(sum(match_scores) / len(match_scores))
    
    # 总体保守性分数（更高更好，所以对我们的评分系统进行反转，其中更低更好）
    avg_conservation = 1.0 - (sum(conservation_scores) / len(conservation_scores) if conservation_scores else 0)
    
    # 组合所有指标（加权）
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
    根据比对质量为每个物种选择最佳序列
    性能优化：批量处理评估
    
    Args:
        species_to_seqs: 物种到{id: 序列}字典的映射
        aligned_file: 比对输出文件路径
        is_protein: 序列是否是蛋白质序列
    
    Returns:
        物种到最佳序列的字典
    """
    logger.info(f"根据比对质量选择最佳重复序列 (蛋白质: {is_protein})")
    
    best_sequences = {}
    
    try:
        # 解析比对文件
        aligned_seqs = {}
        for record in SeqIO.parse(aligned_file, "fasta"):
            aligned_seqs[record.id] = str(record.seq)
        
        # 并行处理物种评估
        species_batch = []
        batch_size = 10  # 每批处理的物种数
        
        for species, id_to_seq in species_to_seqs.items():
            species_batch.append((species, id_to_seq))
            
            if len(species_batch) >= batch_size:
                process_species_batch(species_batch, aligned_seqs, best_sequences, is_protein)
                species_batch = []
        
        # 处理最后一批
        if species_batch:
            process_species_batch(species_batch, aligned_seqs, best_sequences, is_protein)
            
        return best_sequences
    except Exception as e:
        logger.error(f"选择最佳重复序列时出错: {str(e)}")
        logger.error(traceback.format_exc())
        
        # 退回到每个物种的第一个序列
        fallback = {}
        for species, seqs in species_to_seqs.items():
            fallback[species] = seqs[next(iter(seqs))]
            logger.warning(f"为物种 {species} 回退使用第一个序列")
        
        return fallback

def process_species_batch(species_batch, aligned_seqs, best_sequences, is_protein):
    """
    批量处理物种序列评估
    
    Args:
        species_batch: (物种名, id_to_seq字典)元组的列表
        aligned_seqs: 比对序列字典
        best_sequences: 输出的最佳序列字典
        is_protein: 是否是蛋白质序列
    """
    for species, id_to_seq in species_batch:
        if len(id_to_seq) == 1:
            # 只有一个序列，无需评估
            seq_id = next(iter(id_to_seq))
            best_sequences[species] = id_to_seq[seq_id]
            logger.debug(f"物种 {species} 只有一个序列: {seq_id}")
            continue
            
        # 评估每个重复序列
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
            
            # 记录这个物种所有重复序列的分数
            for seq_id, score in all_scores.items():
                logger.debug(f"  物种 {species} 序列 {seq_id} 得分: {score:.4f}" + 
                            (f" (已选择)" if seq_id == best_id else ""))
        else:
            # 如果找不到最佳序列，退回到第一个序列
            first_id = next(iter(id_to_seq))
            best_sequences[species] = id_to_seq[first_id]
            logger.warning(f"无法为物种 {species} 找到最佳序列，使用第一个序列 {first_id}")

def identify_4d_sites_with_tolerance(aligned_seqs):
    """
    从密码子比对序列中识别4倍简并位点，允许一定比例的缺口
    性能优化：使用numpy数组加速处理
    
    Args:
        aligned_seqs: 比对序列字典（物种->序列）
    
    Returns:
        4D位点位置列表
    """
    # 获取比对长度并检查是否是3的倍数
    seq_values = list(aligned_seqs.values())
    if not seq_values:
        logger.warning("未提供序列给 identify_4d_sites_with_tolerance")
        return []
    
    align_length = len(seq_values[0])
    logger.info(f"准备识别4D位点: 比对长度 {align_length}, 是否为3的倍数: {align_length % 3 == 0}")
    
    # 如果比对长度不是3的倍数，停止
    if align_length % 3 != 0:
        logger.warning(f"比对长度 {align_length} 不是 3 的倍数。无法识别4D位点。")
        return []
    
    # 性能优化：将序列转换为numpy数组进行快速处理
    seqs_array = np.array([list(seq) for seq in seq_values])
    
    # 物种总数，用于计算缺口比例
    total_species = len(aligned_seqs)
    max_allowed_gaps = int(total_species * MAX_GAP_RATIO_FOR_4D)
    
    logger.info(f"4D位点识别：共有 {total_species} 个物种，允许最多 {max_allowed_gaps} 个物种有缺口 ({MAX_GAP_RATIO_FOR_4D*100:.0f}%)")
    
    # 初始化计数器
    fourfold_sites = []
    fourfold_candidates = 0
    gap_blocked = 0
    non_matching_prefix = 0
    
    # 性能优化：批量处理密码子
    for i in range(0, align_length, 3):
        # 如果有不完整的密码子，跳过
        if i+2 >= align_length:
            continue
        
        # 提取当前位置的密码子为所有序列
        codons = seqs_array[:, i:i+3]
        
        # 检查每个密码子是否包含缺口
        has_gap = np.any(codons == '-', axis=1)
        gap_count = np.sum(has_gap)
        
        # 如果缺口数超过允许的阈值，跳过
        if gap_count > max_allowed_gaps:
            gap_blocked += 1
            continue
        
        # 提取没有缺口的密码子的前两个位置
        valid_codons = codons[~has_gap]
        prefixes = valid_codons[:, :2]
        
        # 检查前缀是否是4D密码子
        is_4d = np.zeros(len(prefixes), dtype=bool)
        
        for j, prefix in enumerate(prefixes):
            prefix_str = ''.join(prefix).upper()
            is_4d[j] = prefix_str in FOURFOLD_CODONS
        
        # 统计有效4D前缀的数量
        valid_4d_count = np.sum(is_4d)
        
        # 如果有足够的4D前缀
        if valid_4d_count > 0 and (len(prefixes) - valid_4d_count) <= max_allowed_gaps:
            fourfold_candidates += 1
            fourfold_sites.append(i+2)  # 添加第三个位置
        else:
            non_matching_prefix += 1
    
    logger.info(f"4D位点识别结果: 找到 {len(fourfold_sites)} 个4D位点, " +
               f"{fourfold_candidates} 个潜在4D密码子, " +
               f"{gap_blocked} 个密码子因缺口过多而跳过, " +
               f"{non_matching_prefix} 个密码子前缀不匹配")
    
    return fourfold_sites

def extract_4d_sites(aligned_seqs, fourfold_sites):
    """
    从比对序列中提取4倍简并位点
    性能优化：使用列表推导式
    """
    if not fourfold_sites:
        logger.warning("没有 4D 位点可提取")
        return {}
        
    logger.debug(f"从 {len(aligned_seqs)} 个序列中提取 {len(fourfold_sites)} 个 4D 位点")
    
    # 优化：预先检查所有位点是否在范围内
    max_pos = max(fourfold_sites)
    fourfold_seqs = {}
    
    for species, seq in aligned_seqs.items():
        seq_len = len(seq)
        
        # 快速检查是否所有位点都在范围内
        if seq_len > max_pos:
            # 所有位点在范围内，使用列表推导式直接提取
            fourfold_seqs[species] = ''.join(seq[pos] for pos in fourfold_sites)
        else:
            # 可能有超出范围的位点，需要过滤
            valid_sites = [pos for pos in fourfold_sites if pos < seq_len]
            if len(valid_sites) != len(fourfold_sites):
                logger.warning(f"{species} 的一些 4D 位点超出范围: 序列长度 {seq_len}, 最大位点 {max_pos}")
            
            fourfold_seqs[species] = ''.join(seq[pos] for pos in valid_sites)
    
    return fourfold_seqs

def process_cds_file_standard(file_path, output_dir, aligner_params, sequences=None):
    """
    通过比对和4D位点提取处理单个CDS文件
    性能优化：检测文件大小选择合适的处理策略
    
    Args:
        file_path: CDS文件路径
        output_dir: 输出目录
        aligner_params: 比对参数
        sequences: 可选的预过滤序列字典
    """
    try:
        file_name = os.path.basename(file_path)
        gene_name = file_name.replace('.cds', '')
        output_prefix = os.path.join(output_dir, "alignments", gene_name)
        
        # 如果不存在，创建比对目录
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        
        # 准备BLAST参数
        blast_params = None
        if aligner_params.get('use_blast', False) and aligner_params.get('blast_path'):
            logger.info(f"为 {file_name} 启用BLAST参数")
            blast_params = {
                'use_blast': True,
                'blast_path': aligner_params.get('blast_path'),
                'temp_dir': os.path.join(output_dir, "temp")
            }
        
        # 使用提供的序列或使用传统策略从文件解析
        if sequences is None:
            logger.debug(f"使用 {aligner_params['duplicate_strategy']} 策略解析序列")
            # 传递文件路径给blast处理
            input_seqs = parse_fasta(file_path, aligner_params['duplicate_strategy'], blast_params)
        else:
            logger.debug(f"使用预处理的序列集 (长度: {len(sequences)})")
            input_seqs = sequences
            
        if not input_seqs:
            logger.error(f"无法从 {file_path} 解析序列")
            return None
            
        logger.info(f"为 {file_name} 处理 {len(input_seqs)} 个序列")
        
        # 创建带处理序列的临时FASTA文件
        temp_dir = os.path.join(output_dir, "temp")
        os.makedirs(temp_dir, exist_ok=True)
        temp_input = os.path.join(temp_dir, f"{gene_name}_processed.fasta")
        
        # 序列已在解析阶段预处理，直接写入临时文件
        with open(temp_input, 'w') as f:
            for species, seq in input_seqs.items():
                f.write(f">{species}\n{seq}\n")
        
        # 使用适当的比对工具进行比对
        aligned_output = f"{output_prefix}.best.fas"
        success, output_path = run_aligner(temp_input, aligned_output, aligner_params)
        if not success:
            logger.error(f"比对 {file_name} 失败")
            return None
        
        # 确保我们使用正确的输出文件路径
        aligned_file = output_path if output_path else aligned_output
        
        # 解析比对序列
        if not os.path.exists(aligned_file):
            logger.error(f"比对输出文件 {aligned_file} 未找到")
            return None
        
        # 解析比对文件
        aligned_seqs = {}
        for record in SeqIO.parse(aligned_file, "fasta"):
            species = extract_species_name(record.description)
            aligned_seqs[species] = str(record.seq)
        
        logger.info(f"从 {aligned_file} 解析了 {len(aligned_seqs)} 个比对序列")
        
        # 如果需要，修复比对框架以确保密码子结构
        aligned_length = len(next(iter(aligned_seqs.values())))
        if aligned_length % 3 != 0:
            logger.warning(f"比对长度 {aligned_length} 不是3的倍数，尝试修复")
            fixed_file = aligned_file.replace('.best.fas', '.framed.fas')
            if advanced_fix_alignment_frame(aligned_file, fixed_file):
                # 从修复的文件重新读取比对序列
                aligned_seqs = {}
                for record in SeqIO.parse(fixed_file, "fasta"):
                    species = extract_species_name(record.description)
                    aligned_seqs[species] = str(record.seq)
                logger.info(f"使用修复后的比对 ({len(aligned_seqs)} 个序列)")
        
        # 识别并提取4D位点 - 性能优化版本
        fourfold_sites = identify_4d_sites_with_tolerance(aligned_seqs)
        fourfold_seqs = extract_4d_sites(aligned_seqs, fourfold_sites)
        
        # 创建用于调试和跟踪的基因信息
        gene_info = {
            'name': gene_name,
            'species_count': len(fourfold_seqs),
            'site_count': len(fourfold_sites),
            'sequences': fourfold_seqs,
            'aligned_seqs': aligned_seqs  # 存储完整比对序列用于氨基酸翻译
        }
        
        return gene_info
    except Exception as e:
        logger.error(f"处理 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def advanced_fix_alignment_frame(alignment_file, output_file=None):
    """
    高级比对框架修复，考虑缺口模式和序列特征
    性能优化：使用numpy分析缺口模式
    
    Args:
        alignment_file: 输入比对文件路径
        output_file: 输出修复后比对文件路径
        
    Returns:
        bool: 成功状态
    """
    if output_file is None:
        output_file = alignment_file.replace('.fasta', '.framed.fasta')
        if output_file == alignment_file:
            output_file = alignment_file + '.framed'
    
    try:
        # 性能优化：限制加载的序列数量
        records = []
        max_records = 100  # 限制分析的序列数
        
        for i, record in enumerate(SeqIO.parse(alignment_file, "fasta")):
            records.append(record)
            if i >= max_records - 1:
                break
                
        if not records:
            logger.error(f"无法从 {alignment_file} 解析序列")
            return False
            
        alignment_length = len(records[0].seq)
        
        # 检查是否需要修复
        if alignment_length % 3 == 0:
            logger.info(f"比对长度已经是3的倍数 ({alignment_length})，无需修复")
            # 如果输出文件与输入不同，复制文件
            if output_file != alignment_file:
                shutil.copy2(alignment_file, output_file)
            return True
        
        # 性能优化：使用numpy分析缺口模式
        # 将序列转换为二维numpy数组，其中'-'为1，其他为0
        seq_matrix = np.zeros((len(records), alignment_length), dtype=np.int8)
        
        for i, record in enumerate(records):
            seq_str = str(record.seq)
            for j, char in enumerate(seq_str):
                if char == '-':
                    seq_matrix[i, j] = 1
        
        # 计算每列的缺口比例
        gap_columns = np.mean(seq_matrix, axis=0)
        
        # 找出缺口最多的列作为候选插入位置
        padding_needed = 3 - (alignment_length % 3)
        candidate_positions = []
        
        # 使用argsort找出缺口比例最高的列
        sorted_indices = np.argsort(-gap_columns)  # 降序排列
        
        # 先尝试找缺口比例高于阈值的位置
        threshold = 0.5
        high_gap_positions = np.where(gap_columns >= threshold)[0]
        
        if len(high_gap_positions) >= padding_needed:
            # 有足够的高缺口比例位置
            candidate_positions = high_gap_positions[:padding_needed].tolist()
        else:
            # 退回到使用排序后的最高缺口比例位置
            candidate_positions = sorted_indices[:padding_needed].tolist()
        
        # 确保位置是排序的
        insert_positions = sorted(candidate_positions[:padding_needed])
        logger.info(f"在位置 {insert_positions} 插入缺口以修复框架")
        
        # 执行插入
        with open(output_file, 'w') as out:
            for record in records:
                seq = list(str(record.seq))
                # 从后向前插入，避免位置偏移
                for pos in reversed(insert_positions):
                    seq.insert(pos, '-')
                fixed_seq = ''.join(seq)
                out.write(f">{record.id}\n{fixed_seq}\n")
                
            # 处理剩余的序列（如果有限制的话）
            if max_records < sum(1 for _ in SeqIO.parse(alignment_file, "fasta")):
                remaining_count = 0
                for record in SeqIO.parse(alignment_file, "fasta"):
                    # 跳过已处理的序列
                    if remaining_count < max_records:
                        remaining_count += 1
                        continue
                        
                    seq = list(str(record.seq))
                    # 从后向前插入，避免位置偏移
                    for pos in reversed(insert_positions):
                        seq.insert(pos, '-')
                    fixed_seq = ''.join(seq)
                    out.write(f">{record.id}\n{fixed_seq}\n")
        
        return True
    except Exception as e:
        logger.error(f"高级修复比对框架时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def process_cds_file(file_path, output_dir, aligner_params):
    """
    使用适当的策略处理单个CDS文件
    
    Args:
        file_path: CDS文件路径
        output_dir: 输出目录
        aligner_params: 比对参数
    
    Returns:
        带有4D位点的基因信息对象
    """
    # 性能优化：检查文件大小，对大文件使用特殊处理
    file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
    
    # 大文件使用分块处理
    if file_size > 10 * 1024 * 1024:  # > 10MB
        logger.info(f"检测到大型CDS文件: {file_path} ({file_size/1024/1024:.1f}MB)，使用优化处理")
        # 继续使用正常处理路径，但是内部会有优化
    
    if aligner_params.get('duplicate_strategy') == 'alignment_quality':
        return process_cds_file_with_alignment_quality(file_path, output_dir, aligner_params)
    else:
        return process_cds_file_standard(file_path, output_dir, aligner_params)

def process_cds_file_with_alignment_quality(file_path, output_dir, aligner_params):
    """
    使用比对质量处理CDS文件以选择最佳重复序列
    性能优化：使用分批处理
    
    Args:
        file_path: CDS文件路径
        output_dir: 输出目录
        aligner_params: 比对参数
    
    Returns:
        带有4D位点的基因信息对象
    """
    try:
        file_name = os.path.basename(file_path)
        gene_name = file_name.replace('.cds', '')
        temp_dir = os.path.join(output_dir, "temp", gene_name + "_" + str(uuid.uuid4())[:8])
        os.makedirs(temp_dir, exist_ok=True)
        
        logger.info(f"使用比对质量策略处理 {file_name}")
        
        # 准备BLAST参数
        blast_params = None
        if aligner_params.get('use_blast', False) and aligner_params.get('blast_path'):
            blast_params = {
                'use_blast': True,
                'blast_path': aligner_params.get('blast_path'),
                'temp_dir': temp_dir
            }
            logger.info(f"为 {file_name} 启用BLAST参数: {blast_params}")
        
        # 性能优化：检查文件大小，大文件使用分块处理
        file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
        
        # 步骤1: 保留重复序列解析FASTA
        species_to_seqs = parse_fasta_with_duplicates(file_path, blast_params=blast_params)
        
        if not species_to_seqs:
            logger.error(f"无法从 {file_path} 解析序列")
            return None
        
        # 检查是否存在重复
        has_duplicates = any(len(seqs) > 1 for seqs in species_to_seqs.values())
        
        if not has_duplicates:
            logger.info(f"{file_name} 没有重复物种，直接处理")
            # 创建一个每个物种一个序列的常规字典
            simplified_seqs = {species: next(iter(seqs.values())) for species, seqs in species_to_seqs.items()}
            return process_cds_file_standard(file_path, output_dir, aligner_params, simplified_seqs)
        
        # 步骤2: 创建一个包含所有序列（包括重复）的临时FASTA文件
        initial_fasta = os.path.join(temp_dir, f"{gene_name}_all_seqs.fasta")
        
        # 性能优化：批量写入
        batch_size = 100  # 每批写入的序列数
        
        with open(initial_fasta, 'w') as f:
            seq_count = 0
            for species, id_to_seq in species_to_seqs.items():
                for seq_id, sequence in id_to_seq.items():
                    # 序列已在解析阶段预处理
                    f.write(f">{seq_id}\n{sequence}\n")
                    seq_count += 1
                    
                    # 定期刷新缓冲区，避免内存问题
                    if seq_count % batch_size == 0:
                        f.flush()
        
        # 步骤3: 运行初始比对 - 根据序列数量选择适当的参数
        initial_output = os.path.join(temp_dir, f"{gene_name}_initial.best.fas")
        
        # 计算序列总数决定使用哪种比对策略
        total_seqs = sum(len(seqs) for seqs in species_to_seqs.values())
        if total_seqs > 200:
            logger.info(f"检测到大数据集 ({total_seqs} 序列)，调整比对参数以优化性能")
            # 复制参数并修改
            optimized_params = aligner_params.copy()
            
            if aligner_params['aligner'] == 'mafft':
                # 添加MAFFT快速模式参数
                optimized_params['mafft_fast'] = True
            elif aligner_params['aligner'] == 'muscle':
                # 使用MUSCLE的快速模式
                optimized_params['muscle_fast'] = True
                
            success, output_path = run_aligner(initial_fasta, initial_output, optimized_params)
        else:
            success, output_path = run_aligner(initial_fasta, initial_output, aligner_params)
        
        if not success:
            logger.error(f"初始比对 {file_name} 失败")
            # 退回到最长序列策略
            logger.info(f"退回到使用最长序列策略")
            simplified_seqs = {}
            for species, id_to_seq in species_to_seqs.items():
                best_id = max(id_to_seq.items(), key=lambda x: len(x[1]))[0]
                simplified_seqs[species] = id_to_seq[best_id]
                logger.info(f"为物种 {species} 选择了最长序列 {best_id}")
            return process_cds_file_standard(file_path, output_dir, aligner_params, simplified_seqs)
        
        # 如果需要，修复比对框架
        valid, frame_info = check_alignment_validity(output_path)
        if not valid or 'codon' in frame_info and '不保持' in frame_info:
            logger.warning(f"需要修复比对框架: {frame_info}")
            fixed_output = os.path.join(temp_dir, f"{gene_name}_initial_fixed.best.fas")
            advanced_fix_alignment_frame(output_path, fixed_output)
            output_path = fixed_output
        
        # 步骤4: 基于比对质量为每个物种选择最佳序列
        best_sequences = select_best_duplicate_sequences(species_to_seqs, output_path)
        
        if not best_sequences:
            logger.error(f"无法选择最佳重复序列")
            return None
        
        # 步骤5: 使用选定的最佳序列进行处理
        return process_cds_file_standard(file_path, output_dir, aligner_params, best_sequences)
        
    except Exception as e:
        logger.error(f"处理 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None
    finally:
        # 清理临时目录
        try:
            if os.path.exists(temp_dir) and aligner_params.get('clean_temp', True):
                shutil.rmtree(temp_dir)
                logger.debug(f"清理临时目录 {temp_dir}")
        except Exception as e:
            logger.warning(f"清理临时目录 {temp_dir} 时出错: {str(e)}")

def merge_sequences_by_species(all_gene_results, output_file, missing_species_strategy='gaps', min_coverage_pct=0):
    """
    合并所有基因的4D序列，按物种组织
    性能优化：批量处理和并行计算
    
    Args:
        all_gene_results: 带有序列数据的基因结果列表
        output_file: 基本输出文件路径（将添加策略）
        missing_species_strategy: 处理缺失物种的策略
            'gaps' - 用缺口填充缺失序列（默认）
            'exclude_species' - 排除低于最小覆盖率的物种
            'exclude_genes' - 只使用包含符合最小物种覆盖率的基因
        min_coverage_pct: 物种必须存在的最小基因百分比
                          （仅在'exclude_species'策略中使用）
    """
    logger.info(f"使用策略合并序列: {missing_species_strategy}")
    
    # 计时性能
    start_time = time.time()
    
    # 存储每个物种级联序列的字典
    species_supergenes = defaultdict(str)
    
    # 跟踪每个物种贡献哪些基因的字典
    species_gene_coverage = defaultdict(list)
    
    # 跟踪遇到的所有物种和有效基因
    all_species = set()
    valid_gene_results = [g for g in all_gene_results if g is not None]
    total_genes = len(valid_gene_results)
    
    logger.info(f"处理 {total_genes} 个有效基因，共 {len(all_gene_results)} 个总基因")
    
    # 第一遍：识别所有物种 - 性能优化：并行收集
    # 使用并发处理加速
    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as executor:
        futures = []
        for gene_result in valid_gene_results:
            futures.append(executor.submit(lambda g: set(g['sequences'].keys()), gene_result))
        
        # 收集所有物种
        for future in futures:
            all_species.update(future.result())
    
    logger.info(f"在所有基因中发现 {len(all_species)} 个不同的物种")
    
    # 计算每个物种的基因覆盖率 - 性能优化：批量处理
    batch_size = 50  # 每批处理的基因数
    for i in range(0, len(valid_gene_results), batch_size):
        batch = valid_gene_results[i:i+batch_size]
        
        for gene_result in batch:
            gene_name = gene_result['name']
            sequences = gene_result['sequences']
            
            for species in all_species:
                if species in sequences:
                    species_gene_coverage[species].append(gene_name)
    
    # 根据策略应用过滤
    included_species = set()
    included_genes = []
    
    if missing_species_strategy == 'exclude_species':
        # 只包括满足最小覆盖率阈值的物种
        min_genes_required = int(total_genes * min_coverage_pct / 100)
        for species, genes in species_gene_coverage.items():
            if len(genes) >= min_genes_required:
                included_species.add(species)
            else:
                logger.info(f"排除物种 {species}，覆盖率 {len(genes)}/{total_genes} 基因 ({len(genes)/total_genes*100:.1f}%) 低于阈值 {min_coverage_pct}%")
        
        logger.info(f"经过覆盖率过滤后包括 {len(included_species)}/{len(all_species)} 个物种")
        included_genes = valid_gene_results
        
    elif missing_species_strategy == 'exclude_genes':
        # 修改策略：包括符合最低物种覆盖率的基因，而不是要求全部物种
        included_species = all_species
        
        # 计算每个基因包含的物种比例
        min_species_required = int(len(all_species) * (1.0 - MAX_MISSING_SPECIES_RATIO))
        
        for gene_result in valid_gene_results:
            gene_name = gene_result['name']
            sequences = gene_result['sequences']
            species_count = len(sequences)
            
            # 如果基因的物种覆盖率达到要求
            if species_count >= min_species_required:
                included_genes.append(gene_result)
                logger.debug(f"包含基因 {gene_name}，具有 {species_count}/{len(all_species)} 个物种 ({species_count/len(all_species)*100:.1f}%)")
            else:
                logger.info(f"排除基因 {gene_name}，仅有 {species_count}/{len(all_species)} 个物种 ({species_count/len(all_species)*100:.1f}%) 低于要求的 {min_species_required} 个物种")
        
        logger.info(f"包括 {len(included_genes)}/{total_genes} 个满足物种覆盖率要求的基因（要求至少包含 {min_species_required}/{len(all_species)} 个物种）")
        
    else:  # 'gaps' 策略
        included_species = all_species
        included_genes = valid_gene_results
    
    # 检查过滤后是否有任何结果
    if not included_species or not included_genes:
        logger.warning(f"策略 '{missing_species_strategy}' 过滤后没有剩余结果")
        return {}, []
    
    # 第二遍：构建超基因 - 性能优化：并行构建
    logger.info(f"构建超基因，包含 {len(included_species)} 个物种和 {len(included_genes)} 个基因")
    gene_lengths = {}
    
    # 使用并发处理加速构建
    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as executor:
        # 为每个物种预分配一个超基因字典
        species_supergenes = {species: [] for species in included_species}
        gene_order = []  # 保持基因顺序
        
        # 收集每个基因的长度
        future_to_gene = {}
        for gene_result in included_genes:
            future = executor.submit(
                lambda g: (g['name'], len(next(iter(g['sequences'].values()))) if g['sequences'] else 0), 
                gene_result
            )
            future_to_gene[future] = gene_result['name']
        
        # 收集基因长度
        for future in future_to_gene:
            gene_name, length = future.result()
            gene_lengths[gene_name] = length
            gene_order.append(gene_name)
            
        # 并行处理每个基因的序列分配
        future_to_gene = {}
        for gene_result in included_genes:
            future = executor.submit(
                process_gene_sequences, 
                gene_result, 
                included_species,
                gene_order.index(gene_result['name'])
            )
            future_to_gene[future] = gene_result['name']
        
        # 收集结果并填充超基因字典
        for future in future_to_gene:
            gene_name = future_to_gene[future]
            species_sequences = future.result()
            
            for species, (sequence, position) in species_sequences.items():
                # 使用位置信息保持正确的顺序
                while len(species_supergenes[species]) <= position:
                    species_supergenes[species].append("")
                species_supergenes[species][position] = sequence
    
    # 合并每个物种的序列段
    final_supergenes = {}
    for species, sequence_parts in species_supergenes.items():
        final_supergenes[species] = ''.join(sequence_parts)
    
    # 计算并记录覆盖率统计
    coverage_stats = {}
    for species in included_species:
        genes = species_gene_coverage[species]
        genes_in_supergene = [g for g in genes if g in [gene['name'] for gene in included_genes]]
        coverage_pct = (len(genes_in_supergene) / len(included_genes)) * 100 if included_genes else 0
        coverage_stats[species] = f"{len(genes_in_supergene)}/{len(included_genes)} 基因 ({coverage_pct:.1f}%)"
    
    elapsed_time = time.time() - start_time
    logger.info(f"使用 '{missing_species_strategy}' 策略创建包含 {len(final_supergenes)} 个物种的超基因 (耗时: {elapsed_time:.2f}秒)")
    logger.info(f"物种覆盖率统计: {coverage_stats}")
    
    # 创建输出目录结构
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    
    # 生成详细的覆盖率矩阵文件
    output_dir = os.path.dirname(output_file)
    stats_dir = os.path.join(os.path.dirname(os.path.dirname(output_file)), "stats")
    os.makedirs(stats_dir, exist_ok=True)
    
    coverage_matrix = os.path.join(stats_dir, f"species_coverage_matrix_{missing_species_strategy}.tsv")
    with open(coverage_matrix, 'w') as f:
        # 写入带有基因名的标题
        f.write("Species\t" + "\t".join([g['name'] for g in included_genes]) + "\tCoverage(%)\n")
        
        # 写入每个物种的覆盖率
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
    
    # 定义带策略的输出文件名
    strategy_output_file = output_file.replace('.fasta', f'_{missing_species_strategy}.fasta')
    
    # 写入输出文件
    write_output(final_supergenes, strategy_output_file)
    
    return final_supergenes, included_genes

def process_gene_sequences(gene_result, included_species, position):
    """
    处理单个基因的序列，用于并行计算
    
    Args:
        gene_result: 基因结果字典
        included_species: 包含的物种集合
        position: 基因在超基因中的位置
        
    Returns:
        物种到(序列,位置)的字典
    """
    sequences = gene_result['sequences']
    results = {}
    
    if not sequences:
        return results
        
    seq_length = len(next(iter(sequences.values())))
    
    # 为每个物种添加序列或缺口
    for species in included_species:
        if species in sequences:
            results[species] = (sequences[species], position)
        else:
            # 为缺失物种添加缺口
            results[species] = ('-' * seq_length, position)
    
    return results

def create_protein_msa_supergene(all_gene_results, output_file_base, 
                               missing_species_strategy='gaps', min_coverage_pct=0, 
                               aligner_params=None, trimal_params=None):
    """
    创建蛋白质多序列比对超基因
    性能优化：并行处理蛋白质转换
    
    Args:
        all_gene_results: 带有序列数据的基因结果列表
        output_file_base: 基本输出文件路径
        missing_species_strategy: 处理缺失物种的策略
        min_coverage_pct: 最小覆盖率百分比
        aligner_params: 比对参数
        trimal_params: TrimAl参数
    
    Returns:
        成功状态
    """
    try:
        logger.info(f"创建蛋白质MSA超基因，使用策略 '{missing_species_strategy}'")
        start_time = time.time()
        
        # 性能优化：使用并行处理转换蛋白质序列
        protein_gene_results = []
        
        # 多线程处理蛋白质转换
        with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as executor:
            futures = []
            
            for gene_result in all_gene_results:
                if not gene_result or 'aligned_seqs' not in gene_result:
                    continue
                
                futures.append(executor.submit(
                    convert_gene_to_protein,
                    gene_result
                ))
            
            # 收集结果
            for future in futures:
                try:
                    protein_gene_result = future.result()
                    if protein_gene_result:
                        protein_gene_results.append(protein_gene_result)
                except Exception as e:
                    logger.error(f"处理蛋白质转换时出错: {str(e)}")
        
        logger.info(f"成功转换 {len(protein_gene_results)}/{len(all_gene_results)} 个基因为蛋白质序列")
        
        # 合并蛋白质序列
        output_file = output_file_base.replace('_4d', '_protein')
        protein_supergenes, included_genes = merge_sequences_by_species(
            protein_gene_results,
            output_file,
            missing_species_strategy=missing_species_strategy,
            min_coverage_pct=min_coverage_pct
        )
        
        if not protein_supergenes:
            logger.warning("未能创建蛋白质超基因")
            return False
            
        # 如果需要，应用TrimAl
        if trimal_params and trimal_params.get('trim_supergene', False):
            logger.info("对蛋白质超基因应用TrimAl修剪")
            strategy_suffix = f"_{missing_species_strategy}"
            input_file = output_file.replace('.fasta', f'{strategy_suffix}.fasta')
            output_trimmed = input_file.replace('.fasta', '.trimmed.fasta')
            
            success = run_trimal(input_file, output_trimmed, trimal_params)
            if success:
                logger.info(f"成功修剪蛋白质超基因，保存至 {output_trimmed}")
            else:
                logger.warning("修剪蛋白质超基因失败")
        
        elapsed_time = time.time() - start_time
        logger.info(f"创建蛋白质MSA超基因完成，耗时: {elapsed_time:.2f}秒")
        return True
    except Exception as e:
        logger.error(f"创建蛋白质MSA超基因时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def convert_gene_to_protein(gene_result):
    """
    将单个基因的DNA序列转换为蛋白质序列
    用于并行处理
    
    Args:
        gene_result: 基因结果字典
        
    Returns:
        蛋白质序列的基因结果字典或None
    """
    try:
        if not gene_result or 'aligned_seqs' not in gene_result:
            return None
            
        aligned_dna = gene_result['aligned_seqs']
        protein_seqs = {}
        
        # 转换每个物种的比对DNA为蛋白质
        for species, dna_seq in aligned_dna.items():
            if len(dna_seq) % 3 != 0:
                logger.warning(f"基因 {gene_result['name']} 物种 {species} 的DNA长度 {len(dna_seq)} 不是3的倍数，将跳过")
                continue
                
            # 按密码子提取并翻译
            protein_seq = ""
            codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
            
            for codon in codons:
                # 如果整个密码子是缺口，输出缺口
                if codon == '---':
                    protein_seq += '-'
                # 对于含有部分缺口的密码子，输出X
                elif '-' in codon:
                    protein_seq += 'X'
                else:
                    # 翻译有效密码子 - 使用缓存版本
                    aa = safe_translate_codon_cached(codon)
                    protein_seq += aa
            
            protein_seqs[species] = protein_seq
        
        # 创建新的基因结果对象，使用蛋白质序列
        protein_gene_result = {
            'name': gene_result['name'],
            'species_count': len(protein_seqs),
            'sequences': protein_seqs
        }
        
        return protein_gene_result
    except Exception as e:
        logger.error(f"转换基因 {gene_result.get('name', 'unknown')} 为蛋白质时出错: {str(e)}")
        return None

def run_trimal(input_file, output_file, trimal_params):
    """
    运行TrimAl修剪比对
    性能优化：添加多线程参数
    
    Args:
        input_file: 输入比对文件
        output_file: 输出修剪后文件
        trimal_params: TrimAl参数字典
        
    Returns:
        成功状态
    """
    if not trimal_params:
        logger.error("未提供TrimAl参数")
        return False
        
    trimal_path = trimal_params.get('trimal_path')
    if not trimal_path:
        logger.error("未提供TrimAl路径")
        return False
        
    # 验证TrimAl路径
    if not os.path.exists(trimal_path):
        logger.error(f"TrimAl可执行文件未找到: {trimal_path}")
        return False
        
    # 构建基本命令
    cmd = [trimal_path, '-in', input_file, '-out', output_file]
    
    # 添加修剪方法
    if trimal_params.get('automated', True):
        cmd.append('-automated1')
    else:
        # 添加手动参数
        if trimal_params.get('gap_threshold') is not None:
            cmd.extend(['-gt', str(trimal_params['gap_threshold'])])
            
        if trimal_params.get('consistency_threshold') is not None:
            cmd.extend(['-ct', str(trimal_params['consistency_threshold'])])
            
        if trimal_params.get('conservation_threshold') is not None:
            cmd.extend(['-st', str(trimal_params['conservation_threshold'])])
    
    # 性能优化：添加多线程支持参数（如果TrimAl版本支持）
    cpu_count = multiprocessing.cpu_count()
    thread_count = max(1, min(cpu_count - 1, 4))  # 使用最多4个线程或系统最大值-1
    
    # TrimAl某些版本支持-threads参数
    cmd.extend(['-threads', str(thread_count)])
    
    logger.debug(f"TrimAl命令: {' '.join(cmd)}")
    
    try:
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 如果带threads参数失败，可能是旧版本，尝试不带此参数
        if result.returncode != 0 and "-threads" in result.stderr:
            logger.warning("TrimAl不支持多线程参数，尝试不带此参数再次运行")
            cmd = [c for c in cmd if c != '-threads' and c != str(thread_count)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
        if result.returncode != 0:
            logger.error(f"TrimAl运行失败: {result.stderr}")
            return False
            
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            logger.error("TrimAl运行完成但未生成有效输出文件")
            return False
            
        # 验证输出
        records = list(SeqIO.parse(output_file, "fasta"))
        if not records:
            logger.error("TrimAl输出文件不包含有效序列")
            return False
            
        elapsed_time = time.time() - start_time
        logger.info(f"TrimAl成功修剪比对: {len(records)} 条序列，长度 {len(records[0].seq)} bp，耗时: {elapsed_time:.2f}秒")
        return True
    except Exception as e:
        logger.error(f"运行TrimAl时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def write_output(supergene_seqs, output_file):
    """
    将超基因序列写入FASTA文件
    性能优化：批量写入
    """
    if not supergene_seqs:
        logger.warning(f"没有序列可写入 {output_file}")
        return False
        
    try:
        # 如果不存在则创建目录
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # 批量写入以提高性能
        batch_size = 10  # 每批处理的物种数
        species_list = sorted(supergene_seqs.keys())
        
        with open(output_file, 'w') as f:
            for i in range(0, len(species_list), batch_size):
                batch = species_list[i:i+batch_size]
                
                for species in batch:
                    seq = supergene_seqs[species]
                    f.write(f">{species}\n")
                    # 以60个字符一行写入序列，提高可读性
                    for j in range(0, len(seq), 60):
                        f.write(seq[j:j+60] + '\n')
                
                # 定期刷新缓冲区
                f.flush()
        
        seq_length = len(next(iter(supergene_seqs.values()))) if supergene_seqs else 0
        logger.info(f"成功写入 {len(supergene_seqs)} 条序列到 {output_file} (长度: {seq_length} bp)")
        return True
    except Exception as e:
        logger.error(f"写入 {output_file} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def create_full_cds_supergene(all_gene_results, output_file_base, 
                            missing_species_strategy='gaps', min_coverage_pct=0):
    """
    创建完整CDS多序列比对超基因
    
    Args:
        all_gene_results: 带有序列数据的基因结果列表
        output_file_base: 基本输出文件路径
        missing_species_strategy: 处理缺失物种的策略
        min_coverage_pct: 最小覆盖率百分比
    
    Returns:
        成功状态
    """
    try:
        logger.info(f"创建完整CDS超基因，使用策略 '{missing_species_strategy}'")
        start_time = time.time()
        
        # 修改输出文件路径 - 指向full_cds目录而不是4d_sites
        output_file = output_file_base.replace('_4d', '_full_cds')
        output_file = output_file.replace('/4d_sites/', '/full_cds/')
        
        # 收集每个基因的完整CDS比对序列而不是4D位点
        full_cds_gene_results = []
        
        for gene_result in all_gene_results:
            if not gene_result or 'aligned_seqs' not in gene_result:
                continue
                
            # 创建新的基因结果对象，使用完整比对序列
            full_cds_gene_result = {
                'name': gene_result['name'],
                'species_count': len(gene_result['aligned_seqs']),
                'sequences': gene_result['aligned_seqs']  # 使用比对序列而不是4D位点
            }
            
            full_cds_gene_results.append(full_cds_gene_result)
        
        logger.info(f"处理 {len(full_cds_gene_results)}/{len(all_gene_results)} 个基因的完整CDS序列")
        
        # 合并完整CDS序列
        full_cds_supergenes, included_genes = merge_sequences_by_species(
            full_cds_gene_results,
            output_file,
            missing_species_strategy=missing_species_strategy,
            min_coverage_pct=min_coverage_pct
        )
        
        if not full_cds_supergenes:
            logger.warning("未能创建完整CDS超基因")
            return False
            
        elapsed_time = time.time() - start_time
        logger.info(f"创建完整CDS超基因完成，耗时: {elapsed_time:.2f}秒")
        return True
    except Exception as e:
        logger.error(f"创建完整CDS超基因时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def main():
    """主函数：解析参数并执行流程"""
    parser = argparse.ArgumentParser(description="多序列比对和4D位点提取流水线（性能优化版）")
    parser.add_argument("--input_dir", default=".", help="包含CDS文件的目录")
    parser.add_argument("--output_dir", default="./output", help="输出文件目录")
    parser.add_argument("--supergene_output", default="supergene_4d.fasta", help="超基因输出文件的基础名称")
    
    # 比对工具选择
    parser.add_argument("--aligner", choices=["prank", "muscle", "mafft"], default="prank",
                      help="使用的比对工具 (默认: prank)")
    parser.add_argument("--prank_path", help="PRANK可执行文件的绝对路径")
    parser.add_argument("--muscle_path", help="MUSCLE可执行文件的绝对路径")
    parser.add_argument("--mafft_path", help="MAFFT可执行文件的绝对路径")
    parser.add_argument("--threads", type=int, default=4, help="用于并行处理的线程数")
    
    # BLAST参数
    parser.add_argument("--use_blast", action="store_true", help="使用BLAST帮助修复非3倍数序列")
    parser.add_argument("--blast_path", help="BLASTN可执行文件的绝对路径")
    parser.add_argument("--blast_db_cache_size", type=int, default=5, help="BLAST数据库缓存大小")
    
    # 通用比对参数
    parser.add_argument("--no_codon_aware", action="store_true", default=False, 
                       help="禁用密码子感知比对 (默认: 启用)")
    
    # PRANK特定参数
    parser.add_argument("--f", type=float, default=0.2, help="PRANK插入开放概率")
    parser.add_argument("--gaprate", type=float, default=None, 
                      help="PRANK缺口开放率 (默认: 根据数据类型使用PRANK的默认值)")
    parser.add_argument("--gapext", type=float, default=None, 
                      help="PRANK缺口扩展概率 (默认: 根据数据类型使用PRANK的默认值)")
    parser.add_argument("--use_logs", action="store_true", 
                      help="对大型数据集在PRANK中使用对数计算")
    parser.add_argument("--penalize_terminal_gaps", action="store_true",
                      help="在PRANK中正常惩罚终端缺口")
    
    # TrimAl参数
    parser.add_argument("--use_trimal", action="store_true", help="使用TrimAl进行蛋白质比对修剪")
    parser.add_argument("--trimal_path", help="TrimAl可执行文件的绝对路径")
    parser.add_argument("--trimal_automated", action="store_true", default=True, 
                      help="在TrimAl中使用自动修剪方法")
    parser.add_argument("--gap_threshold", type=float, default=None, help="TrimAl最小缺口阈值")
    parser.add_argument("--consistency_threshold", type=float, default=None, help="TrimAl一致性阈值")
    parser.add_argument("--conservation_threshold", type=float, default=None, help="TrimAl保守性阈值")
    parser.add_argument("--trim_supergene", action="store_true", help="对最终蛋白质超基因应用TrimAl")
    
    # 性能优化参数
    parser.add_argument("--batch_size", type=int, default=100, help="批处理大小，用于优化内存使用")
    parser.add_argument("--memory_limit", type=int, default=0, 
                      help="内存限制（MB），0表示无限制。当处理大文件时有用。")
    
    # 其他参数
    parser.add_argument("--duplicate_strategy", 
                      choices=['longest', 'first', 'rename', 'alignment_quality'], 
                      default='alignment_quality',
                      help="处理文件中重复物种的策略")
    parser.add_argument("--skip_existing", action="store_true", help="如果比对文件存在则跳过处理")
    parser.add_argument("--min_coverage_pct", type=float, default=30.0,
                      help="物种必须存在的最小基因百分比 (用于exclude_species策略)")
    parser.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO',
                      help="设置日志级别")
    parser.add_argument("--clean_temp", action="store_true", default=True,
                      help="处理后清理临时文件")
    parser.add_argument("--create_full_cds", action="store_true", default=True,
                      help="创建完整CDS序列的超基因 (默认: 启用)")
    parser.add_argument("--create_protein_msa", action="store_true", default=True,
                      help="Create multiple sequence alignments of protein sequences")
    parser.add_argument("--min_repeat_report", type=int, default=70,
                      help="报告为质量问题所需的最小重复数量 (默认: 70)")
    parser.add_argument("--max_gap_ratio", type=float, default=0.2,
                      help="4D位点允许的最大缺口物种比例 (默认: 0.2)")
    parser.add_argument("--max_missing_species", type=float, default=0.1,
                      help="基因包含允许的最大缺失物种比例 (默认: 0.1)")
    
    args = parser.parse_args()
    
    # 更新基于命令行参数的全局常量
    global MIN_REPEAT_COUNT_TO_REPORT, MAX_GAP_RATIO_FOR_4D, MAX_MISSING_SPECIES_RATIO, BLAST_DB_CACHE_SIZE
    MIN_REPEAT_COUNT_TO_REPORT = args.min_repeat_report
    MAX_GAP_RATIO_FOR_4D = args.max_gap_ratio
    MAX_MISSING_SPECIES_RATIO = args.max_missing_species
    BLAST_DB_CACHE_SIZE = args.blast_db_cache_size
    
    # 验证选择的比对工具路径
    if args.aligner == "prank" and not args.prank_path:
        parser.error("选择PRANK作为比对工具时，必须提供--prank_path参数")
    elif args.aligner == "muscle" and not args.muscle_path:
        parser.error("选择MUSCLE作为比对工具时，必须提供--muscle_path参数")
    elif args.aligner == "mafft" and not args.mafft_path:
        parser.error("选择MAFFT作为比对工具时，必须提供--mafft_path参数")
        
    # 如果启用，验证TrimAl
    if args.use_trimal and not args.trimal_path:
        parser.error("启用TrimAl修剪时，必须提供--trimal_path参数")
        
    # 如果启用，验证BLAST
    if args.use_blast and not args.blast_path:
        parser.error("启用BLAST处理时，必须提供--blast_path参数")
    
    # 设置日志级别
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # 配置内存限制（如果指定）
    if args.memory_limit > 0:
        import resource
        # 将软限制设置为用户请求的限制（以字节为单位）
        resource.setrlimit(
            resource.RLIMIT_AS, 
            (args.memory_limit * 1024 * 1024, resource.RLIM_INFINITY)
        )
        logger.info(f"设置内存限制为 {args.memory_limit} MB")
    
    # 记录启动时间和参数
    logger.info(f"启动MSA流水线 (性能优化版)，参数: {vars(args)}")
    start_time = time.time()
    
    # 设置目录结构
    try:
        # 确保输出目录存在
        os.makedirs(args.output_dir, exist_ok=True)
        
        # 创建必要的子目录
        subdirs = ["alignments", "4d_sites", "full_cds", "stats", "temp"]
        for subdir in subdirs:
            os.makedirs(os.path.join(args.output_dir, subdir), exist_ok=True)
    except Exception as e:
        logger.error(f"设置输出目录时出错: {str(e)}")
        return 1
    
    # 查找CDS文件
    cds_files = find_cds_files(args.input_dir)
    if not cds_files:
        logger.error("未找到.cds文件。退出。")
        return 1
        
    logger.info(f"将处理 {len(cds_files)} 个CDS文件")
    
    # 准备比对工具参数
    aligner_params = {
        'aligner': args.aligner,
        'prank_path': args.prank_path,
        'muscle_path': args.muscle_path,
        'mafft_path': args.mafft_path,
        'codon_aware': not args.no_codon_aware,
        'duplicate_strategy': args.duplicate_strategy,
        'clean_temp': args.clean_temp,
        'batch_size': args.batch_size,
        
        # PRANK特定参数
        'f': args.f,
        'gaprate': args.gaprate,
        'gapext': args.gapext,
        'use_logs': args.use_logs,
        'penalize_terminal_gaps': args.penalize_terminal_gaps
    }
    
    # 添加BLAST参数
    if args.use_blast:
        aligner_params['use_blast'] = args.use_blast
        aligner_params['blast_path'] = args.blast_path
    
    # TrimAl参数
    trimal_params = None
    if args.use_trimal:
        trimal_params = {
            'use_trimal': args.use_trimal,
            'trimal_path': args.trimal_path,
            'automated': args.trimal_automated,
            'gap_threshold': args.gap_threshold,
            'consistency_threshold': args.consistency_threshold,
            'conservation_threshold': args.conservation_threshold,
            'trim_supergene': args.trim_supergene
        }
    
    # 使用线程池并行处理每个CDS文件
    all_gene_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {}
        
        for file_name in cds_files:
            file_path = os.path.join(args.input_dir, file_name)
            gene_name = file_name.replace('.cds', '')
            output_file = os.path.join(args.output_dir, "alignments", f"{gene_name}.best.fas")
            
            # 如果请求，跳过已存在的文件
            if args.skip_existing and os.path.exists(output_file):
                logger.info(f"跳过 {file_name}，输出已存在")
                
                try:
                    # 仍需处理现有文件以提取4D位点
                    aligned_seqs = {}
                    for record in SeqIO.parse(output_file, "fasta"):
                        species = extract_species_name(record.description)
                        aligned_seqs[species] = str(record.seq)
                        
                    if not aligned_seqs:
                        logger.warning(f"无法解析现有的比对 {output_file}")
                        continue
                        
                    # 验证并修复密码子比对（如需要）
                    aligned_length = len(next(iter(aligned_seqs.values())))
                    if aligned_length % 3 != 0:
                        logger.warning(f"现有比对文件 {output_file} 长度 {aligned_length} 不是3的倍数，修复中")
                        fixed_file = output_file.replace('.best.fas', '.framed.fas')
                        if advanced_fix_alignment_frame(output_file, fixed_file):
                            # 从修复的文件重新读取比对序列
                            aligned_seqs = {}
                            for record in SeqIO.parse(fixed_file, "fasta"):
                                species = extract_species_name(record.description)
                                aligned_seqs[species] = str(record.seq)
                            logger.info(f"使用修复后的比对")
                        
                    fourfold_sites = identify_4d_sites_with_tolerance(aligned_seqs)
                    fourfold_seqs = extract_4d_sites(aligned_seqs, fourfold_sites)
                    
                    all_gene_results.append({
                        'name': gene_name,
                        'species_count': len(fourfold_seqs),
                        'site_count': len(fourfold_sites),
                        'sequences': fourfold_seqs,
                        'aligned_seqs': aligned_seqs
                    })
                    logger.info(f"处理现有比对 {file_name}: {len(fourfold_seqs)} 个物种, {len(fourfold_sites)} 个4D位点")
                except Exception as e:
                    logger.error(f"处理现有比对 {output_file} 时出错: {str(e)}")
                    logger.error(traceback.format_exc())
                continue
            
            futures[executor.submit(
                process_cds_file, 
                file_path,
                args.output_dir,
                aligner_params
            )] = file_name
        
        # 收集结果并显示进度
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
                    logger.info(f"处理 {file_name}: {result['species_count']} 个物种, {result['site_count']} 个4D位点")
                else:
                    logger.error(f"处理 {file_name} 失败")
            except Exception as e:
                logger.exception(f"处理 {file_name} 时出错: {str(e)}")
    
    if not all_gene_results:
        logger.error("没有基因结果成功处理")
        return 1
        
    # 创建超基因输出路径
    supergene_4d_path = os.path.join(args.output_dir, "4d_sites", args.supergene_output)
    
    # 使用三种策略处理
    strategies = ['gaps', 'exclude_species', 'exclude_genes']
    
    for strategy in strategies:
        try:
            logger.info(f"使用 '{strategy}' 策略创建4D位点超基因")
            supergene_seqs, merged_genes = merge_sequences_by_species(
                all_gene_results, 
                supergene_4d_path,
                missing_species_strategy=strategy,
                min_coverage_pct=args.min_coverage_pct
            )
            if not supergene_seqs:
                logger.warning(f"策略 '{strategy}' 未能生成有效的4D位点超基因")
                
            # 如果请求，创建蛋白质MSA超基因
            if args.create_protein_msa and merged_genes:
                logger.info(f"使用 '{strategy}' 策略创建蛋白质MSA超基因")
                create_protein_msa_supergene(
                    all_gene_results,
                    supergene_4d_path,
                    missing_species_strategy=strategy,
                    min_coverage_pct=args.min_coverage_pct,
                    aligner_params=aligner_params,
                    trimal_params=trimal_params
                )
        
            if args.create_full_cds:
                logger.info(f"使用 '{strategy}' 策略创建完整CDS超基因")
                create_full_cds_supergene(
                    all_gene_results,
                    supergene_4d_path,
                    missing_species_strategy=strategy,
                    min_coverage_pct=args.min_coverage_pct
                )
            
        except Exception as e:
            logger.error(f"使用策略 '{strategy}' 创建超基因时出错: {str(e)}")
            logger.error(traceback.format_exc())
    
    # 如果请求，清理临时文件
    if args.clean_temp:
        try:
            logger.info("清理临时文件")
            temp_dir = os.path.join(args.output_dir, "temp")
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        except Exception as e:
            logger.warning(f"清理临时文件时出错: {str(e)}")
    
    # 记录执行时间
    end_time = time.time()
    execution_time = end_time - start_time
    hours, remainder = divmod(execution_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    logger.info(f"流水线完成。耗时: {int(hours)}小时 {int(minutes)}分钟 {seconds:.1f}秒")
    logger.info(f"结果保存在: {args.output_dir}")
    
    return 0

if __name__ == "__main__":
    try:
        # 启用垃圾收集诊断，帮助调试内存问题
        import gc
        gc.set_debug(gc.DEBUG_STATS)
        
        sys.exit(main())
    except Exception as e:
        logger.critical(f"未捕获的异常: {str(e)}")
        logger.critical(traceback.format_exc())
        sys.exit(1)
