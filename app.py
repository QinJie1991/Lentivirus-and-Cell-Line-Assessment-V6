"""
慢病毒包装-细胞系评估系统

修复版本 - 解决前端模块加载错误、版本兼容性问题
"""

import streamlit as st
import requests
import json
import time
import re
import html
import logging
import sqlite3
import os
import difflib
import csv
import xml.etree.ElementTree as ET
import base64
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime, timedelta
from io import StringIO
import pandas as pd
import zipfile
import sys

# ==================== 版本兼容性处理 ====================
def safe_rerun():
    """兼容不同 Streamlit 版本的 rerun"""
    try:
        st.rerun()
    except AttributeError:
        try:
            st.experimental_rerun()
        except:
            pass  # 如果都失败，不强制 rerun

def safe_cache_data(func):
    """兼容不同 Streamlit 版本的 cache"""
    try:
        return st.cache_data(func)
    except AttributeError:
        try:
            return st.cache(func, allow_output_mutation=True)
        except:
            return func  # 无缓存直接运行

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("lentivirus_assessment")

# ==================== 页面配置 ====================
st.set_page_config(
    page_title="慢病毒包装-细胞系评估系统",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== AI模型配置 ====================
AVAILABLE_AI_MODELS = {
    'qwen3.6-plus-2026-04-02': '通义千问-Plus (可用额度)'
}
DEFAULT_AI_MODEL = 'qwen3.6-plus-2026-04-02'
# ================================================

# ==================== HPA细胞系自动补全服务（新增）====================
class HPACellLineAutocompleteService:
    """HPA数据库细胞系名称自动补全服务 - 包含1206个细胞系的标准名称"""

    # HPA v25.0 proteinatlas.tsv 中的主要细胞系名称列表
    # 基于 "RNA cell line [Name] [nTPM]" 列提取
    HPA_CELL_LINES = [
        "A-431", "A-549", "AN3-CA", "ASC52telo", "BEWO", "BJ", "CACO-2", "Calu-6",
        "CHP-212", "Daudi", "HaCaT", "HAP1", "HBEC-3i", "HCE-2", "HCT 116", "HEK 293",
        "HeLa", "Hep G2", "HepaRG", "HL-60", "HMC-1", "HSkMC", "hTCEpi", "hTERT-HME1",
        "HUVEC", "HUVEC-Tert2", "K-562", "Karpas-707", "LNCaP", "MCF7", "MOLT-4", "NB-4",
        "NCI-H2073", "NCI-H2170", "NCI-H226", "NCI-H358", "NCI-H460", "NCI-H520",
        "NCI-H596", "NCI-H69", "NK-92", "PC-3", "REH", "RH-30", "RPMI-8226", "RT-4",
        "SCLC-21H", "SH-SY5Y", "SiHa", "SK-BR-3", "SK-MEL-30", "T-47D", "THP-1",
        "TIME", "U-138 MG", "U-2197", "U-251 MG", "U-266/70", "U-266/84", "U-698",
        "U-87 MG", "U-937", "WM-115", "WM-266-4", "A3/HL-60", "A498", "A704", "ACHN",
        "AGS", "AMO-1", "C-33A", "C3A", "Calu-1", "Calu-3", "CCRF-CEM", "CL-11",
        "COR-L23", "COR-L88", "COV318", "COV504", "Detroit 562", "DKO", "DoTc2-4510",
        "EFO-21", "EFO-27", "EM-2", "EPLC-272H", "ES-2", "ESS-1", "EW-7", "F-36P",
        "FU-OV-1", "G-401", "GAMG", "GM12878", "HBL-1", "HCC38", "HCC827", "HDLM-2",
        "HEC-1-B", "HH", "Jiyoye", "Jurkat E6.1", "KE-37", "KG-1", "KM-3", "KNS-42",
        "KPL-1", "KU-812", "L-363", "L-428", "L-540", "Loukes", "LP-1", "MHH-CALL-2",
        "MHH-CALL-3", "MHH-CALL-4", "MSTO-211H", "N-549", "NALM-6", "NCI-H1299",
        "NCI-H1395", "NCI-H1437", "NCI-H146", "NCI-H1650", "NCI-H1666", "NCI-H1694",
        "NCI-H1792", "NCI-H1836", "NCI-H1838", "NCI-H1944", "NCI-H1975", "NCI-H2009",
        "NCI-H2023", "NCI-H2030", "NCI-H2052", "NCI-H2087", "NCI-H211", "NCI-H2122",
        "NCI-H2126", "NCI-H2228", "NCI-H23", "NCI-H2342", "NCI-H2347", "NCI-H2444",
        "NCI-H2452", "NCI-H28", "NCI-H446", "NCI-H522", "NCI-H719", "NCI-H82",
        "NCI-H838", "NCI-H889", "NTERA-2", "OCI-AML2", "OCI-AML3", "OCI-LY-3",
        "OCI-LY-7", "OV-90", "PANC-1", "PC-3M-1E8", "PC-3M-2B4", "Raji", "RKO",
        "S-117", "SCLC-24H", "SCLC-5HM", "SEM", "SK-LMS-1", "SK-LU-1", "SK-N-AS",
        "SK-N-DZ", "SK-N-FI", "SK-OV-3", "SK-UT-1", "SuDHL-10", "SuDHL-16", "SuDHL-4",
        "SuDHL-5", "T-ALL-1", "T98G", "TE-10", "TE-11", "TE-12", "TE-15", "TE-441",
        "TE-617", "TE-8", "TF-1", "TK-6", "TT", "U-2 OS", "UM-UC-3", "WS-1"
    ]

    # 细胞系别名映射（支持常见书写变体）
    CELL_LINE_ALIASES = {
        "HEK 293": ["HEK293", "HEK-293", "293", "HEK293A", "HEK293T"],
        "HeLa": ["Hela", "HELA", "HeLa S3"],
        "A-549": ["A549", "A 549", "NCI-A549"],
        "Hep G2": ["HepG2", "Hep-G2", "HEPG2"],
        "MCF7": ["MCF-7", "MCF 7"],
        "U-87 MG": ["U87MG", "U-87", "U87", "U87-MG"],
        "SH-SY5Y": ["SHSY5Y", "SH SY5Y", "SY5Y"],
        "K-562": ["K562", "K 562"],
        "THP-1": ["THP1", "THP 1"],
        "PC-3": ["PC3", "PC 3"],
        "LNCaP": ["LNCAP"],
        "U-2 OS": ["U2OS", "U-2OS", "U2-OS"],
        "HCT 116": ["HCT116", "HCT-116"],
        "NCI-H226": ["NCIH226", "NCI H226", "H226"],
        "NCI-H460": ["NCIH460", "NCI H460", "H460"],
        "NCI-H520": ["NCIH520", "NCI H520", "H520"],
        "NCI-H1975": ["NCIH1975", "NCI H1975", "H1975"],
        "NCI-H1299": ["NCIH1299", "NCI H1299", "H1299"],
        "NCI-H358": ["NCIH358", "NCI H358", "H358"],
        "NCI-H23": ["NCIH23", "NCI H23", "H23"],
        "NCI-H2170": ["NCIH2170", "NCI H2170", "H2170"],
        "NCI-H2073": ["NCIH2073", "NCI H2073", "H2073"],
        "HK-2": ["HK2", "HK 2"],
        "RAW 264.7": ["RAW2647", "RAW264.7"],
        "NIH/3T3": ["NIH3T3", "NIH 3T3", "3T3"],
        "CHO-K1": ["CHOK1", "CHO K1"],
        "COS-7": ["COS7", "COS 7"],
        "BHK-21": ["BHK21", "BHK 21"],
        "PC-12": ["PC12", "PC 12"],
        "CACO-2": ["CACO2", "Caco 2"],
        "HT-29": ["HT29", "HT 29"],
        "SW-480": ["SW480", "SW 480"],
        "DU-145": ["DU145", "DU 145"],
        "Saos-2": ["Saos2", "Saos 2"],
        "MG-63": ["MG63", "MG 63"],
        "A-431": ["A431", "A 431"],
        "HaCaT": ["HACAT"],
        "BEAS-2B": ["BEAS2B", "BEAS 2B"],
        "Calu-3": ["Calu3", "Calu 3"],
        "Jurkat E6.1": ["Jurkat", "JURKAT"],
    }

    def __init__(self):
        self.search_index = {}
        self._build_search_index()

    def _build_search_index(self):
        """构建搜索索引，包括所有名称变体"""
        for cell_line in self.HPA_CELL_LINES:
            normalized = self._normalize(cell_line)
            self.search_index[normalized] = cell_line
            parts = cell_line.replace('-', ' ').replace('/', ' ').split()
            for part in parts:
                if len(part) >= 2:
                    part_norm = self._normalize(part)
                    if part_norm not in self.search_index:
                        self.search_index[part_norm] = cell_line

        for standard, aliases in self.CELL_LINE_ALIASES.items():
            for alias in aliases:
                alias_norm = self._normalize(alias)
                if alias_norm not in self.search_index:
                    hpa_name = self._find_hpa_equivalent(standard)
                    self.search_index[alias_norm] = hpa_name if hpa_name else standard

    def _find_hpa_equivalent(self, name: str):
        norm = self._normalize(name)
        for cell in self.HPA_CELL_LINES:
            if self._normalize(cell) == norm:
                return cell
        return None

    def _normalize(self, name: str) -> str:
        if not name:
            return ""
        return name.upper().replace('-', '').replace(' ', '').replace('/', '').replace('_', '')

    def get_suggestions(self, query: str, limit: int = 8) -> list:
        if not query or len(query) < 1:
            return []
        query_norm = self._normalize(query)
        if not query_norm:
            return []

        matches = []
        seen = set()

        if query_norm in self.search_index:
            cell = self.search_index[query_norm]
            matches.append({'display_name': cell, 'hpa_name': cell, 'match_type': 'exact', 'score': 100})
            seen.add(cell)

        for norm, cell in self.search_index.items():
            if norm.startswith(query_norm) and cell not in seen:
                matches.append({'display_name': cell, 'hpa_name': cell, 'match_type': 'prefix', 'score': 80 - len(norm)})
                seen.add(cell)

        for norm, cell in self.search_index.items():
            if query_norm in norm and cell not in seen:
                matches.append({'display_name': cell, 'hpa_name': cell, 'match_type': 'substring', 'score': 50})
                seen.add(cell)

        if len(query_norm) >= 3:
            for norm, cell in self.search_index.items():
                if cell not in seen:
                    sim = difflib.SequenceMatcher(None, query_norm, norm).ratio()
                    if sim > 0.6:
                        matches.append({'display_name': cell, 'hpa_name': cell, 'match_type': 'fuzzy', 'score': int(sim * 40)})
                        seen.add(cell)

        matches.sort(key=lambda x: (x['score'], x['display_name']), reverse=True)
        return matches[:limit]

    def get_exact_match(self, query: str):
        return self.search_index.get(self._normalize(query))

    def is_valid_cell_line(self, query: str) -> bool:
        return self.get_exact_match(query) is not None



# ==================== HPA基因自动补全服务（基于Gene synonym）====================
class CellLineNormalizer:
    """细胞系名称标准化、模糊匹配验证（解决书写不规范问题）"""

    STANDARD_CELL_LINES = {
        "HEK293": ["HEK 293", "HEK-293", "293", "HEK", "Human Embryonic Kidney 293", "HEK-293A", "HEK 293A"],
        "HEK293T": ["HEK 293T", "HEK-293T", "293T", "HEK293-T", "HEK 293 T", "293-T"],
        "HeLa": ["Hela", "HELA", "Hela Cell", "cervical cancer cell line", "HeLa S3"],
        "A549": ["A 549", "A-549", "lung carcinoma cell line", "A549 cell", "NCI-A549"],
        "HepG2": ["Hep G2", "Hep-G2", "hepatoma cell line", "HepG2 cell", "HEPG2"],
        "HK-2": ["HK2", "HK 2", "human kidney 2", "renal proximal tubule", "HK-2 cell", "HK 2 cell"],
        "NCI-H226": ["NCI H226", "NCIH226", "H226", "lung squamous cell carcinoma", "NCI H-226", "H-226"],
        "MCF-7": ["MCF7", "MCF 7", "breast cancer cell line", "MCF-7 cell", "MCF7 cell"],
        "U-87 MG": ["U87MG", "U-87", "U87", "glioblastoma", "U87-MG", "U 87 MG"],
        "SH-SY5Y": ["SHSY5Y", "SH SY5Y", "neuroblastoma", "SH-SY5Y cell", "SHP-SY5Y", "SY5Y"],
        "K-562": ["K562", "K 562", "chronic myeloid leukemia", "K-562 cell", "K562 cell"],
        "THP-1": ["THP1", "THP 1", "acute monocytic leukemia", "THP-1 cell", "THP1 cell"],
        "Jurkat": ["JURKAT", "T cell leukemia", "Jurkat E6.1", "Jurkat cell"],
        "RAW264.7": ["RAW 264.7", "RAW2647", "mouse macrophage", "Raw 264.7", "RAW 2647"],
        "NIH/3T3": ["NIH 3T3", "NIH3T3", "mouse fibroblast", "NIH-3T3", "3T3"],
        "CHO-K1": ["CHO K1", "CHOK1", "Chinese hamster ovary", "CHO K-1", "CHO cell"],
        "MDCK": ["Madin-Darby Canine Kidney", "dog kidney", "MDCK cell", "M.D.C.K."],
        "Vero": ["VERO", "African green monkey kidney", "Vero cell", "VERO CCL81"],
        "COS-7": ["COS7", "COS 7", "monkey kidney", "Cos-7", "COS 7 cell"],
        "BHK-21": ["BHK21", "BHK 21", "baby hamster kidney", "BHK-21 cell"],
        "PC-12": ["PC12", "PC 12", "rat adrenal pheochromocytoma", "PC-12 cell"],
        "Caco-2": ["Caco2", "Caco 2", "colorectal adenocarcinoma", "CACO-2", "Caco2 cell"],
        "HT-29": ["HT29", "HT 29", "colon cancer", "HT-29 cell", "HT29 cell"],
        "HCT116": ["HCT 116", "HCT-116", "colorectal carcinoma", "HCT-116 cell"],
        "SW480": ["SW 480", "SW-480", "colon adenocarcinoma", "SW-480 cell"],
        "DU145": ["DU 145", "DU-145", "prostate cancer", "Du145", "DU-145 cell"],
        "PC-3": ["PC3", "PC 3", "prostate cancer", "PC-3 cell", "PC3 cell"],
        "LNCaP": ["LNCaP", "lymph node carcinoma of the prostate", "LNCAP", "LNCaP cell"],
        "Saos-2": ["Saos2", "Saos 2", "osteosarcoma", "SAOS-2", "Saos2 cell"],
        "MG-63": ["MG63", "MG 63", "osteosarcoma", "MG-63 cell", "MG63 cell"],
        "U-2 OS": ["U2OS", "U-2OS", "osteosarcoma", "U2-OS", "U 2 OS"],
        "A431": ["A 431", "A-431", "epidermoid carcinoma", "A431 cell"],
        "HaCaT": ["HACAT", "human keratinocyte", "HaCaT cell", "HACAT cell"],
        "BEAS-2B": ["BEAS2B", "BEAS 2B", "bronchial epithelium", "BEAS-2B cell"],
        "Calu-3": ["Calu3", "Calu 3", "lung adenocarcinoma", "Calu-3 cell", "CALU-3"],
        "H1975": ["H-1975", "H 1975", "lung adenocarcinoma", "H1975 cell", "NCI-H1975"],
        "H1299": ["H-1299", "H 1299", "lung carcinoma", "H1299 cell", "NCI-H1299"],
        "H460": ["H-460", "H 460", "large cell lung carcinoma", "H460 cell", "NCI-H460"],
        "H358": ["H-358", "H 358", "lung adenocarcinoma", "H358 cell"],
        "H23": ["H-23", "H 23", "lung adenocarcinoma", "H23 cell", "NCI-H23"],
        "H520": ["H-520", "H 520", "lung squamous cell carcinoma", "H520 cell", "NCI-H520"],
        "H1703": ["H-1703", "H 1703", "lung squamous cell carcinoma", "H1703 cell"],
        "H2170": ["H-2170", "H 2170", "lung squamous cell carcinoma", "H2170 cell"],
    }

    CELL_TYPE_KEYWORDS = {
        'lung': ['lung', 'pulmonary', 'bronchial', 'alveolar', 'A549', 'H1975', 'H1299', 'Calu-3', 'H460', 'H358', 'H23', 'H520', 'H1703', 'H2170', 'H226', 'NCI-H'],
        'kidney': ['kidney', 'renal', 'tubule', 'HK-2', 'HEK293', 'MDCK', 'proximal'],
        'liver': ['liver', 'hepatic', 'hepatoma', 'HepG2', 'Huh7', 'Hep3B', 'hepatocyte'],
        'brain': ['brain', 'neuro', 'glioma', 'astrocyte', 'SH-SY5Y', 'U-87', 'SK-N-SH', 'neuroblastoma', 'glioblastoma'],
        'blood': ['blood', 'leukemia', 'lymphoma', 'T cell', 'B cell', 'Jurkat', 'K-562', 'THP-1', 'HL-60', 'myeloid'],
        'skin': ['skin', 'keratinocyte', 'melanoma', 'HaCaT', 'A431', 'epidermoid'],
        'colon': ['colon', 'colorectal', 'intestinal', 'Caco-2', 'HT-29', 'HCT116', 'SW480', 'rectal'],
        'prostate': ['prostate', 'LNCaP', 'PC-3', 'DU145', 'prostatic'],
        'bone': ['bone', 'osteosarcoma', 'Saos-2', 'MG-63', 'U-2 OS', 'osteoblast'],
        'fibroblast': ['fibroblast', '3T3', 'MRC-5', 'BJ', 'IMR-90'],
        'macrophage': ['macrophage', 'RAW264.7', 'THP-1', 'monocyte'],
    }

    @classmethod
    def normalize(cls, name: str) -> str:
        if not name:
            return ""
        normalized = name.strip().upper()
        normalized = ' '.join(normalized.split())
        normalized = re.sub(r'\s*-\s*', '-', normalized)
        normalized = re.sub(r'([A-Z])(\d)', r'\1-\2', normalized)
        normalized = re.sub(r'-+', '-', normalized)

        if re.match(r'^HEK-?293$', normalized) or normalized == '293':
            return 'HEK293'
        if re.match(r'^HEK-?293T$', normalized) or normalized == '293T':
            return 'HEK293T'

        if normalized.startswith('NCI') and not normalized.startswith('NCI-'):
            normalized = normalized.replace('NCI', 'NCI-', 1)

        normalized = re.sub(r'-?CELL$', '', normalized)
        return normalized

    @classmethod
    def find_best_match(cls, input_name: str) -> Tuple[Optional[str], float, List[str]]:
        if not input_name:
            return None, 0.0, []

        normalized_input = cls.normalize(input_name)
        matches = []

        for standard, aliases in cls.STANDARD_CELL_LINES.items():
            normalized_standard = cls.normalize(standard)

            if normalized_input == normalized_standard:
                return standard, 1.0, [standard]

            for alias in aliases:
                normalized_alias = cls.normalize(alias)
                if normalized_input == normalized_alias:
                    matches.append((standard, 1.0))
                    break

        if matches:
            return matches[0][0], 1.0, [m[0] for m in matches]

        best_score = 0.0
        best_matches = []

        for standard, aliases in cls.STANDARD_CELL_LINES.items():
            candidates = [standard] + aliases
            for candidate in candidates:
                score = cls._calculate_similarity(normalized_input, cls.normalize(candidate))
                if score > best_score:
                    best_score = score
                    best_matches = [(standard, score)]
                elif abs(score - best_score) < 0.01 and score > 0.5:
                    if (standard, score) not in best_matches:
                        best_matches.append((standard, score))

        if best_matches and best_score > 0.6:
            best_matches.sort(key=lambda x: x[1], reverse=True)
            top_match = best_matches[0][0]
            all_candidates = list(dict.fromkeys([m[0] for m in best_matches]))[:5]
            return top_match, best_score, all_candidates

        return None, 0.0, []

    @classmethod
    def _calculate_similarity(cls, s1: str, s2: str) -> float:
        if s1 == s2:
            return 1.0

        edit_sim = difflib.SequenceMatcher(None, s1, s2).ratio()

        substring_bonus = 0.0
        if s1 in s2 or s2 in s1:
            shorter = min(len(s1), len(s2))
            longer = max(len(s1), len(s2))
            substring_bonus = 0.2 * (shorter / longer)

        keyword_bonus = 0.0
        s1_parts = set(s1.replace('-', '').split())
        s2_parts = set(s2.replace('-', '').split())
        common = s1_parts & s2_parts
        if common:
            keyword_bonus = len(common) / max(len(s1_parts), len(s2_parts)) * 0.1

        return min(1.0, edit_sim + substring_bonus + keyword_bonus)

    @classmethod
    def get_cell_type_hint(cls, name: str) -> Optional[str]:
        name_lower = name.lower()
        matched_types = []

        for cell_type, keywords in cls.CELL_TYPE_KEYWORDS.items():
            if any(kw.lower() in name_lower for kw in keywords):
                matched_types.append(cell_type)

        return ", ".join(matched_types) if matched_types else None

    @classmethod
    def validate_and_suggest(cls, input_name: str) -> Dict:
        result = {
            'input': input_name,
            'normalized': cls.normalize(input_name),
            'is_valid': False,
            'confidence': 0.0,
            'suggested_standard': None,
            'alternatives': [],
            'cell_type': None,
            'needs_confirmation': True,
            'warning': None
        }

        if not input_name or len(input_name.strip()) < 2:
            result['warning'] = "细胞系名称过短"
            return result

        if re.search(r'[A-Z]{3,}\d{3,}', input_name) and '-' not in input_name:
            result['warning'] = f"检测到'{input_name}'可能缺少分隔符，建议格式如'NCI-H226'而非'NCIH226'"

        best_match, confidence, alternatives = cls.find_best_match(input_name)
        result['confidence'] = confidence

        if best_match:
            result['is_valid'] = True
            result['suggested_standard'] = best_match
            result['alternatives'] = alternatives[1:] if len(alternatives) > 1 else []
            result['cell_type'] = cls.get_cell_type_hint(best_match)

            if confidence >= 0.9:
                result['needs_confirmation'] = False
            elif confidence >= 0.7:
                result['needs_confirmation'] = True
            else:
                result['needs_confirmation'] = True
                result['warning'] = f"匹配度较低({confidence:.0%})，请仔细核对细胞系名称"
        else:
            result['cell_type'] = cls.get_cell_type_hint(input_name)
            result['warning'] = f"未找到匹配的细胞系'{input_name}'，将使用原始输入进行检索"

        return result

# ==================== AI分析客户端 ====================
class AIAnalysisClient:
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.base_url = "https://dashscope.aliyuncs.com/api/v1/services/aigc/text-generation/generation"

        self.CELL_CULTURE_DIFFICULTY_CHECKLIST = {
            "培养基特殊要求": [
                "无血清培养基（如UltraCULTURE、SFM）",
                "需要添加特殊血清（如马血清、新生牛血清）",
                "需要添加生长因子（如EGF、bFGF、HGF、VEGF等）",
                "需要添加激素（如胰岛素、氢化可的松、地塞米松）",
                "需要添加微量元素（如硒酸钠、转铁蛋白）",
                "培养基价格>200元/500ml（如StemPro、mTeSR）",
                "需要定制化培养基（如添加特定脂肪酸、胆固醇）",
                "需要条件培养基（收集其他细胞上清）"
            ],
            "基质/包被要求": [
                "需要Matrigel包被（价格昂贵，4°C操作）",
                "需要胶原蛋白包被（I型/IV型）",
                "需要纤连蛋白（Fibronectin）包被",
                "需要明胶包被（0.1% Gelatin）",
                "需要多聚赖氨酸（Poly-L-lysine）包被",
                "需要饲养层细胞（Feeder layer，如MEF、STO）",
                "需要低吸附培养板（Ultra-Low Attachment）",
                "需要专用培养装置（如Spinner flask、旋转瓶）"
            ],
            "气体/环境敏感": [
                "需要低氧培养（5%或更低O₂，非常规CO₂培养箱）",
                "对pH极度敏感（需每日观察颜色变化）",
                "需要持续震荡/搅拌（贴壁培养会死亡）",
                "对CO₂浓度敏感（需精确控制5-10%）",
                "对温度波动敏感（不能长时间开门）"
            ],
            "操作复杂度": [
                "密度极度敏感（太稀死，太密也死，窗口窄）",
                "对机械力敏感（吹打、离心即死）",
                "需要每天换液/传代（周末无法休息）",
                "需要半量换液（不能全换，否则死亡）",
                "需要频繁观察/拍照（每日至少2次）",
                "胰酶消化时间窗口窄（<2分钟或>10分钟）",
                "容易成团成簇（需频繁吹打分散）",
                "容易自发分化（需持续加因子维持）"
            ],
            "时间/成本": [
                "倍增时间>48小时（生长极其缓慢）",
                "倍增时间<12小时（需每天传代，累死人）",
                "原代细胞（无法扩增，需反复取材）",
                "有限细胞系（传代次数受限，如<10代）",
                "需要特定批次血清（不同批次差异大）",
                "支原体敏感（一旦污染立即死亡，且形态无明显变化）"
            ]
        }

    def analyze_antiviral_evidence(self, gene_name: str, title: str, abstract: str) -> Dict:
        """
        分析文献中是否报道了基因的抗病毒功能。
        严格基于提供的文献内容，禁止推测。
        """
        if not self.api_key:
            return {'is_antiviral': False, 'confidence': 0, 'mechanism': '', 'reasoning': '未配置API'}

        try:
            prompt = f"""请【严格基于以下文献内容】判断基因"{gene_name}"是否具有抗病毒功能，并进行语义分析和归纳总结。

【重要警告】你绝对不能进行推测：
1. 只能基于提供的文献标题和摘要判断
2. 如果文献未明确提及抗病毒功能，必须返回false
3. 不要基于基因名称或一般知识进行推断

文献标题：{title}
文献摘要：{abstract}

请按以下JSON格式回答（只返回JSON）：
{{
    "is_antiviral": true/false,
    "confidence": 0.0-1.0,
    "mechanism": "文献明确描述的抗病毒机制，如无则留空",
    "reasoning": "引用文献中的具体描述作为判断依据",
    
    "semantic_analysis": {{
        "evidence_strength": "基于文献的语义证据强度评估（Strong/Moderate/Weak/None）",
        "functional_context": "文献中描述的功能上下文",
        "implied_roles": ["文献暗示的可能作用（但不编造）"],
        "terminology_signals": "文献用语的信号分析"
    }},
    
    "inductive_summary": {{
        "core_claim": "文献的核心主张总结",
        "supporting_evidence": "支持抗病毒功能的证据",
        "limitations": "文献中提到的局限性或条件",
        "research_context": "研究背景对结论的影响"
    }}
}}

严格标准：
1. is_antiviral=true：文献【明确】提到该基因能抑制病毒复制、增强抗病毒免疫等
2. is_antiviral=false：文献未提及抗病毒功能，或仅提及其他功能
3. confidence：0.0-1.0，基于证据明确程度
4. mechanism：必须是文献中明确描述的，禁止推测
5. reasoning：必须引用文献原文或准确概括文献内容
6. 【语义分析】深度理解文献的表述方式和隐含信息
7. 【归纳总结】提炼文献的核心发现和证据结构"""

            headers = {
                'Authorization': f'Bearer {self.api_key}',
                'Content-Type': 'application/json'
            }
            current_model = st.session_state.get('selected_ai_model', DEFAULT_AI_MODEL)
            
            payload = {
                'model': current_model,
                'input': {
                    'messages': [
                        {'role': 'system', 'content': '你是一个严谨的生物医学文献分析助手。你只基于提供的文献内容进行判断，绝不进行推测或引申。'},
                        {'role': 'user', 'content': prompt}
                    ]
                },
                'parameters': {
                    'result_format': 'message',
                    'max_tokens': 500,
                    'temperature': 0.1
                }
            }

            response = requests.post(
                self.base_url,
                headers=headers,
                json=payload,
                timeout=30
            )
            response.raise_for_status()

            result = response.json()
            content = result.get('output', {}).get('choices', [{}])[0].get('message', {}).get('content', '')

            try:
                content_clean = content.replace('```json', '').replace('```', '').strip()
                analysis = json.loads(content_clean)
                return {
                    'is_antiviral': analysis.get('is_antiviral', False),
                    'confidence': float(analysis.get('confidence', 0)),
                    'mechanism': analysis.get('mechanism', ''),
                    'reasoning': analysis.get('reasoning', '')
                }
            except json.JSONDecodeError:
                is_antiviral = any(kw in (title + abstract).lower() for kw in
                    ['antiviral', 'virus', 'interferon', 'ifitm', 'innate immunity'])
                return {
                    'is_antiviral': is_antiviral,
                    'confidence': 0.3 if is_antiviral else 0,
                    'mechanism': '解析失败，使用关键词匹配（仅供参考）',
                    'reasoning': 'API返回格式异常，使用关键词匹配降级处理'
                }

        except Exception as e:
            return {
                'is_antiviral': False,
                'confidence': 0,
                'mechanism': '',
                'reasoning': f'API调用失败: {str(e)}'
            }

    def analyze_gene_function_comprehensive(self, gene_name: str, gene_description: str,
                                          papers_oe: List[Dict], papers_kd: List[Dict],
                                          papers_ko: List[Dict], papers_general: List[Dict]) -> Dict:
        if not self.api_key:
            return {'error': '未配置AI API', 'note': '请在侧边栏配置API Key或在secrets中设置DASHSCOPE_API_KEY'}

        def format_papers(papers, label):
            if not papers:
                return f"\n{label}文献：无相关文献\n"
            text = f"\n{label}文献：\n"
            for i, p in enumerate(papers[:5], 1):
                text += f"{i}. {p.get('title', '')} - {p.get('abstract', '')[:300]}... PMID:{p.get('pmid', 'N/A')}\n"
            return text

        literature_text = ""
        literature_text += format_papers(papers_general, "基因功能相关")
        literature_text += format_papers(papers_oe, "过表达")
        literature_text += format_papers(papers_kd, "敲低/敲除")
        literature_text += format_papers(papers_ko, "敲除")

        total_papers = len(papers_general) + len(papers_oe) + len(papers_kd) + len(papers_ko)
        if total_papers == 0:
            return {
                'protein_function': {'category': '未检索到文献', 'domains': '未检索到文献', 'pathways': '未检索到文献', 'cellular_location': '未检索到文献', 'tissue_expression': '未检索到文献'},
                'overexpression': {'cell_models': [], 'animal_models': [], 'summary': '无相关文献'},
                'knockdown': {'cell_models': [], 'summary': '无相关文献'},
                'knockout': {'cell_models': [], 'animal_models': [], 'summary': '无相关文献'},
                'disease_relevance': {'cancer': '未检索到文献', 'other_diseases': '未检索到文献', 'therapeutic_potential': '未检索到文献'},
                'key_references': [],
                'experimental_notes': '未检索到该基因的相关文献，无法提供功能总结。建议查阅NCBI Gene数据库或PubMed获取最新信息。',
                'data_source_note': '未找到相关文献，AI未进行推测'
            }

        prompt = f"""作为分子生物学和遗传学专家，请【严格基于以下提供的文献】总结基因"{gene_name}"的功能及实验模型数据。

【重要警告】你绝对不能进行推测或编造：
1. 只能基于提供的文献回答问题
2. 如果文献未提及某项信息，该字段必须标注"文献未提供"
3. 禁止基于基因名称或一般知识进行推断
4. 每个具体数据必须标注文献来源（PMID）

提供的文献：
{literature_text}

请按以下JSON格式提供结构化总结（只返回JSON）：
{{
    "protein_function": {{
        "category": "文献明确描述的功能类别 + [PMID:XXXX]",
        "domains": "文献明确描述的结构域 + [PMID:XXXX]",
        "pathways": "文献明确描述的通路 + [PMID:XXXX]",
        "cellular_location": "文献明确描述的亚细胞定位 + [PMID:XXXX]",
        "tissue_expression": "文献明确描述的组织表达 + [PMID:XXXX]"
    }},
    "overexpression": {{
        "cell_models": [
            {{
                "cell_line": "文献明确报道的细胞系",
                "phenotype": "文献明确报道的表型",
                "mechanism": "文献明确描述的机制",
                "reference": "PMID:XXXX（文献标题）"
            }}
        ],
        "animal_models": [
            {{
                "model": "文献明确报道的动物模型",
                "phenotype": "文献明确报道的表型",
                "reference": "PMID:XXXX（文献标题）"
            }}
        ],
        "summary": "基于文献的过表达效应总结"
    }},
    "knockdown": {{
        "cell_models": [
            {{
                "cell_line": "文献明确报道的细胞系",
                "method": "文献明确报道的敲低方法",
                "phenotype": "文献明确报道的表型",
                "reference": "PMID:XXXX（文献标题）"
            }}
        ],
        "summary": "基于文献的敲低效应总结"
    }},
    "knockout": {{
        "cell_models": [
            {{
                "cell_line": "文献明确报道的细胞系",
                "method": "文献明确报道的敲除方法",
                "phenotype": "文献明确报道的表型",
                "viability": "文献明确报道的细胞活力影响",
                "reference": "PMID:XXXX（文献标题）"
            }}
        ],
        "animal_models": [
            {{
                "model": "文献明确报道的动物模型",
                "phenotype": "文献明确报道的表型",
                "lethality": "文献明确报道的致死性",
                "reference": "PMID:XXXX（文献标题）"
            }}
        ],
        "summary": "基于文献的敲除效应总结"
    }},
    "disease_relevance": {{
        "cancer": "文献明确描述的肿瘤作用 + [PMID:XXXX]",
        "other_diseases": "文献明确描述的其他疾病 + [PMID:XXXX]",
        "therapeutic_potential": "文献明确描述的治疗潜力 + [PMID:XXXX]"
    }},
    "key_references": [
        "PMID:XXXX（文献标题）"
    ],
    "experimental_notes": "基于文献的实验设计建议",
    "data_source_note": "基于X篇文献分析，所有信息均有文献支持",
    
    "semantic_analysis": {{
        "gene_character": "基于文献的基因特性语义描述（如：促癌基因/抑癌基因/必需基因等）",
        "functional_complexity": "功能复杂度的语义评估",
        "experimental_challenges": ["文献中暗示的实验挑战"],
        "safety_considerations": "基于文献的安全考虑"
    }},
    
    "inductive_summary": {{
        "functional_consensus": "多篇文献对该基因功能的共识",
        "phenotypic_patterns": "表型模式总结（如：过表达总是促进增殖）",
        "context_dependent_effects": "上下文依赖性效应（如：细胞类型特异性）",
        "knowledge_conflicts": "文献间的矛盾或不一致之处",
        "research_gaps": "研究空白和未来方向"
    }}
}}

严格要求：
1. 【禁止推测】文献未明确报道的信息必须标注"文献未提供"
2. 【必须标注来源】每个具体描述后必须标注[PMID:XXXX]
3. 【禁止常识推断】即使是常见基因（如p53、GAPDH），没有文献支持不得输出
4. 如果某类文献缺失（如无敲除文献），该部分必须返回空数组或标注"无相关文献"
5. key_references只能列出上述提供的真实文献，必须包含准确PMID
6. 【语义分析】深度理解文献中的功能描述和隐含信息
7. 【归纳总结】整合多篇文献的发现，找出功能模式、共识和矛盾"""

        try:
            headers = {
                'Authorization': f'Bearer {self.api_key}',
                'Content-Type': 'application/json'
            }
            current_model = st.session_state.get('selected_ai_model', DEFAULT_AI_MODEL)
            
            payload = {
                'model': current_model,
                'input': {
                    'messages': [
                        {'role': 'system', 'content': '你是分子生物学专家，精通基因功能注释和表型分析。你只能基于用户提供的文献进行总结，绝对禁止编造任何文献、PMID、作者或期刊信息。reference字段必须严格使用格式："第一作者 et al., 年份, 期刊缩写, PMID:数字"。如果信息不确定，明确标注"未见报道"。'},
                        {'role': 'user', 'content': prompt}
                    ]
                },
                'parameters': {
                    'result_format': 'message',
                    'max_tokens': 3000,
                    'temperature': 0.2
                }
            }

            response = requests.post(
                self.base_url,
                headers=headers,
                json=payload,
                timeout=60
            )
            response.raise_for_status()

            result = response.json()
            content = result.get('output', {}).get('choices', [{}])[0].get('message', {}).get('content', '')

            try:
                content_clean = content.replace('```json', '').replace('```', '').strip()
                return json.loads(content_clean)
            except json.JSONDecodeError:
                return {
                    'error': 'AI返回格式异常',
                    'raw_response': content[:1000],
                    'note': 'AI未能返回有效JSON格式'
                }

        except Exception as e:
            return {'error': str(e), 'note': 'API调用失败，请检查网络或API Key有效性'}

    def design_rnai_sequences(self, gene_name: str, gene_species: str = "human", target_region: str = "CDS") -> Dict:
        """设计shRNA/siRNA候选序列（修复截断部分，保持原有架构）"""
        if not self.api_key:
            return {'error': '未配置API Key', 'sequences': []}

        prompt = f"""请为基因"{gene_name}"（物种：{gene_species}，靶向区域：{target_region}）设计3条高效shRNA/siRNA候选序列。
遵循标准分子生物学规则：长度19-23nt，GC含量30%-55%，避免连续4个相同碱基，避开已知SNP位点。
请按以下JSON格式返回（仅返回JSON）：
{{
    "sequences": [
        {{"target_seq": "序列", "position": 起始位置, "gc_content": 0.45, "off_target_score": "低", "note": "设计说明"}}
    ],
    "design_notes": "设计依据与注意事项"
}}"""
        
        try:
            headers = {'Authorization': f'Bearer {self.api_key}', 'Content-Type': 'application/json'}
            current_model = st.session_state.get('selected_ai_model', DEFAULT_AI_MODEL)
            payload = {
                'model': current_model,
                'input': {
                    'messages': [
                        {'role': 'system', 'content': '你是专业的RNAi序列设计AI。仅返回符合分子生物学标准的合法JSON，不添加任何额外文本。'},
                        {'role': 'user', 'content': prompt}
                    ]
                },
                'parameters': {'result_format': 'message', 'max_tokens': 1000, 'temperature': 0.2}
            }
            response = requests.post(self.base_url, headers=headers, json=payload, timeout=45)
            response.raise_for_status()
            content = response.json().get('output', {}).get('choices', [{}])[0].get('message', {}).get('content', '')
            
            content_clean = content.replace('```json', '').replace('```', '').strip()
            return json.loads(content_clean)
        except json.JSONDecodeError:
            return {'error': 'AI返回格式异常', 'raw': content[:500]}
        except Exception as e:
            return {'error': f'API调用失败: {str(e)}'}
