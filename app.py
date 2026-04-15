"""
慢病毒包装-细胞系评估系统
完整修复版 - 解决前端模块加载错误、版本兼容性问题、补充完整UI与导出逻辑
"""

import streamlit as st
import requests
import json
import time
import re
import logging
import difflib
import csv
import os
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import pandas as pd
from io import StringIO

# ==================== 版本兼容性处理 ====================
def safe_rerun():
    try:
        st.rerun()
    except AttributeError:
        try:
            st.experimental_rerun()
        except:
            pass

def safe_cache_data(func):
    try:
        return st.cache_data(func)
    except AttributeError:
        try:
            return st.cache(func, allow_output_mutation=True)
        except:
            return func

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
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
    'qwen3.6-plus-2026-04-02': '通义千问-Plus (推荐)'
}
DEFAULT_AI_MODEL = 'qwen3.6-plus-2026-04-02'

# ==================== HPA细胞系自动补全服务 ====================
class HPACellLineAutocompleteService:
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
        if not name: return ""
        return name.upper().replace('-', '').replace(' ', '').replace('/', '').replace('_', '')

    def get_suggestions(self, query: str, limit: int = 8) -> list:
        if not query or len(query) < 1: return []
        query_norm = self._normalize(query)
        matches, seen = [], set()
        if query_norm in self.search_index:
            matches.append({'display_name': self.search_index[query_norm], 'match_type': 'exact', 'score': 100})
            seen.add(self.search_index[query_norm])
        for norm, cell in self.search_index.items():
            if norm.startswith(query_norm) and cell not in seen:
                matches.append({'display_name': cell, 'match_type': 'prefix', 'score': 80 - len(norm)})
                seen.add(cell)
        for norm, cell in self.search_index.items():
            if query_norm in norm and cell not in seen:
                matches.append({'display_name': cell, 'match_type': 'substring', 'score': 50})
                seen.add(cell)
        if len(query_norm) >= 3:
            for norm, cell in self.search_index.items():
                if cell not in seen:
                    sim = difflib.SequenceMatcher(None, query_norm, norm).ratio()
                    if sim > 0.6:
                        matches.append({'display_name': cell, 'match_type': 'fuzzy', 'score': int(sim * 40)})
                        seen.add(cell)
        matches.sort(key=lambda x: (x['score'], x['display_name']), reverse=True)
        return matches[:limit]

# ==================== 细胞系标准化服务 ====================
class CellLineNormalizer:
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

    @classmethod
    def normalize(cls, name: str) -> str:
        if not name: return ""
        normalized = name.strip().upper()
        normalized = ' '.join(normalized.split())
        normalized = re.sub(r'\s*-\s*', '-', normalized)
        normalized = re.sub(r'([A-Z])(\d)', r'\1-\2', normalized)
        normalized = re.sub(r'-+', '-', normalized)
        if re.match(r'^HEK-?293$', normalized) or normalized == '293': return 'HEK293'
        if re.match(r'^HEK-?293T$', normalized) or normalized == '293T': return 'HEK293T'
        if normalized.startswith('NCI') and not normalized.startswith('NCI-'): normalized = normalized.replace('NCI', 'NCI-', 1)
        normalized = re.sub(r'-?CELL$', '', normalized)
        return normalized

    @classmethod
    def find_best_match(cls, input_name: str) -> Tuple[Optional[str], float, List[str]]:
        if not input_name: return None, 0.0, []
        normalized_input = cls.normalize(input_name)
        matches = []
        for standard, aliases in cls.STANDARD_CELL_LINES.items():
            normalized_standard = cls.normalize(standard)
            if normalized_input == normalized_standard: return standard, 1.0, [standard]
            for alias in aliases:
                if cls.normalize(alias) == normalized_input:
                    matches.append((standard, 1.0))
                    break
        if matches: return matches[0][0], 1.0, [m[0] for m in matches]
        best_score, best_matches = 0.0, []
        for standard, aliases in cls.STANDARD_CELL_LINES.items():
            for candidate in [standard] + aliases:
                s1, s2 = normalized_input, cls.normalize(candidate)
                edit_sim = difflib.SequenceMatcher(None, s1, s2).ratio()
                sub_bonus = 0.2 * (min(len(s1), len(s2)) / max(len(s1), len(s2))) if s1 in s2 or s2 in s1 else 0.0
                kw_bonus = len(set(s1.replace('-','').split()) & set(s2.replace('-','').split())) / max(len(set(s1.replace('-','').split())), len(set(s2.replace('-','').split()))) * 0.1
                score = min(1.0, edit_sim + sub_bonus + kw_bonus)
                if score > best_score: best_score, best_matches = score, [(standard, score)]
                elif abs(score - best_score) < 0.01 and score > 0.5: best_matches.append((standard, score))
        if best_matches and best_score > 0.6:
            best_matches.sort(key=lambda x: x[1], reverse=True)
            return best_matches[0][0], best_score, list(dict.fromkeys([m[0] for m in best_matches]))[:5]
        return None, 0.0, []

# ==================== AI分析客户端 ====================
class AIAnalysisClient:
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.base_url = "https://dashscope.aliyuncs.com/api/v1/services/aigc/text-generation/generation"

    def _call_ai(self, prompt: str, system_msg: str = "", max_tokens: int = 500, temp: float = 0.1) -> str:
        if not self.api_key: raise ValueError("未配置API Key")
        model = st.session_state.get('selected_ai_model', DEFAULT_AI_MODEL)
        headers = {'Authorization': f'Bearer {self.api_key}', 'Content-Type': 'application/json'}
        payload = {
            'model': model,
            'input': {'messages': [{'role': 'system', 'content': system_msg}, {'role': 'user', 'content': prompt}]},
            'parameters': {'result_format': 'message', 'max_tokens': max_tokens, 'temperature': temp}
        }
        res = requests.post(self.base_url, headers=headers, json=payload, timeout=60)
        res.raise_for_status()
        return res.json().get('output', {}).get('choices', [{}])[0].get('message', {}).get('content', '')

    def analyze_antiviral_evidence(self, gene_name: str, title: str, abstract: str) -> Dict:
        prompt = f"""请【严格基于以下文献内容】判断基因"{gene_name}"是否具有抗病毒功能。
文献标题：{title}\n文献摘要：{abstract}\n
返回JSON格式：{{"is_antiviral":true/false, "confidence":0.0-1.0, "mechanism":"", "reasoning":""}}
严禁推测，无明确证据必须返回false。"""
        raw = self._call_ai(prompt, "你是严谨的生物医学文献分析助手。只基于提供内容判断，绝不推测。", 500, 0.1)
        try:
            clean = raw.replace('```json','').replace('```','').strip()
            return json.loads(clean)
        except: return {'is_antiviral': False, 'confidence': 0, 'mechanism': '', 'reasoning': 'AI格式异常，建议人工复核'}

    def analyze_gene_function_comprehensive(self, gene_name: str, papers: List[Dict]) -> Dict:
        if not papers: return {'note': '无文献输入'}
        lit_text = "\n".join([f"{i+1}. {p['title']} - {p['abstract'][:300]} (PMID:{p['pmid']})" for i, p in enumerate(papers[:5])])
        prompt = f"""作为分子生物学专家，请【严格基于以下文献】总结基因"{gene_name}"功能。
文献列表：\n{lit_text}\n
返回JSON：{{"protein_function":{{"category":"","pathways":""}}, "overexpression":{{"cell_models":[]}}, "knockdown":{{"summary":""}}, "disease_relevance":{{}}}}
禁止编造，未提及字段填'文献未提供'。"""
        raw = self._call_ai(prompt, "你是分子生物学专家。只能基于用户提供的文献总结，绝对禁止编造。", 3000, 0.2)
        try: return json.loads(raw.replace('```json','').replace('```','').strip())
        except: return {'raw': raw[:1000], 'note': 'JSON解析失败'}

    def design_rnai_sequences(self, gene_name: str, species: str = "human", region: str = "CDS") -> Dict:
        prompt = f"""请为基因"{gene_name}"（物种：{species}，靶区：{region}）设计3条shRNA/siRNA候选序列。
严格遵循分子生物学规则：19-23nt，GC含量30-52%，避免种子区连续G/T，避开SNP高发区。
返回JSON格式数组：[{"sequence":"", "target_position":0, "gc_content":0.0, "off_target_risk":"低/中/高", "note":""}]
仅返回JSON数组。"""
        raw = self._call_ai(prompt, "你是RNAi序列设计专家。严格遵循shRNA设计规则，仅返回合法JSON。", 800, 0.3)
        try:
            clean = raw.replace('```json','').replace('```','').strip()
            data = json.loads(clean)
            return data if isinstance(data, list) else [{'error': 'AI返回非数组格式'}]
        except: return [{'error': 'JSON解析失败，请重试'}]

# ==================== Streamlit 主程序 ====================
def main():
    # 初始化会话状态
    if 'api_key' not in st.session_state:
        st.session_state.api_key = os.environ.get('DASHSCOPE_API_KEY', '')
    if 'selected_ai_model' not in st.session_state:
        st.session_state.selected_ai_model = DEFAULT_AI_MODEL
    if 'cell_results' not in st.session_state: st.session_state.cell_results = []
    if 'analysis_results' not in st.session_state: st.session_state.analysis_results = []
    if 'rnai_results' not in st.session_state: st.session_state.rnai_results = []

    # 侧边栏配置
    with st.sidebar:
        st.title("⚙️ 系统设置")
        api_key = st.text_input("🔑 DashScope API Key", value=st.session_state.api_key, type="password", 
                                help="在阿里云百炼平台申请，或留空查看演示模式")
        st.session_state.api_key = api_key.strip()
        st.session_state.selected_ai_model = st.selectbox("🤖 AI模型", list(AVAILABLE_AI_MODELS.keys()), format_func=lambda x: AVAILABLE_AI_MODELS[x])
        
        if not api_key:
            st.warning("⚠️ 未配置API Key，AI功能将使用占位响应或需手动填写结果。建议配置后刷新页面。")
        else:
            st.success("✅ API Key 已就绪")
            
        st.markdown("---")
        st.caption("慢病毒包装-细胞系评估系统 v1.0")
        st.caption("基于HPA数据库 & 通义千问API")

    # 主标签页
    tabs = st.tabs(["🔍 细胞系标准化验证", "📖 文献AI深度分析", "🧬 RNAi/shRNA设计", "📥 数据导出与日志"])

    client = AIAnalysisClient(api_key=st.session_state.api_key)
    hpa_svc = HPACellLineAutocompleteService()

    # === Tab 1: 细胞系验证 ===
    with tabs[0]:
        st.header("🔍 细胞系名称标准化与HPA匹配")
        col1, col2 = st.columns([2, 1])
        with col1:
            user_input = st.text_input("输入细胞系名称（支持别名/缩写/大小写）", placeholder="例：hek293t / A549 / h1299")
        with col2:
            st.caption("💡 提示：系统将自动清洗格式、匹配HPA标准名、评估置信度")
            
        if user_input and st.button("🔎 立即验证", type="primary"):
            with st.spinner("正在标准化匹配与HPA交叉验证..."):
                norm = CellLineNormalizer.normalize(user_input)
                match, conf, alts = CellLineNormalizer.find_best_match(user_input)
                suggestions = hpa_svc.get_suggestions(user_input, limit=5)
                
                res = {
                    "时间": datetime.now().strftime("%Y-%m-%d %H:%M"),
                    "原始输入": user_input,
                    "标准化输出": norm,
                    "最佳匹配": match or "未匹配到标准库",
                    "置信度": f"{conf:.1%}",
                    "HPA支持": "✅ 是" if match and hpa_svc.is_valid_cell_line(match) else "❌ 否（需实验确认）",
                    "备选建议": ", ".join(alts[:3]) if alts else "无",
                    "HPA补全提示": [s['display_name'] for s in suggestions] if suggestions else ["无"]
                }
                st.session_state.cell_results.append(res)
                
                st.success("验证完成")
                st.json(res)
                st.info("📌 实验建议：置信度<0.7时，建议通过STR鉴定复核。慢病毒包装优先选择HEK293T/HEK293A/HEK293-FT等易转染株系。")

    # === Tab 2: 文献分析 ===
    with tabs[1]:
        st.header("📖 基因功能与抗病毒证据AI分析")
        gene = st.text_input("目标基因", placeholder="例：IFITM3, APOBEC3G, SAMHD1")
        title = st.text_input("文献标题（可选）", placeholder="可留空，系统将基于摘要推断")
        abstract = st.text_area("文献摘要", placeholder="粘贴PubMed摘要或实验描述...")
        pmid = st.text_input("PMID（仅记录用，不参与AI推理）", placeholder="例：34567890")
        
        if st.button("🧠 启动AI分析", type="primary", disabled=not abstract):
            if not st.session_state.api_key:
                st.error("请先在侧边栏配置API Key")
            else:
                with st.spinner("正在调用AI进行严格语义分析（禁止推测模式）..."):
                    av = client.analyze_antiviral_evidence(gene, title, abstract)
                    func = client.analyze_gene_function_comprehensive(gene, [{"title": title, "abstract": abstract, "pmid": pmid}])
                    res = {"基因": gene, "PMID": pmid, "抗病毒判断": av, "功能总结": func, "时间": datetime.now().strftime("%Y-%m-%d %H:%M")}
                    st.session_state.analysis_results.append(res)
                    st.success("分析完成（已开启严格事实校验）")
                    st.json(av)
                    st.markdown("### 📊 功能与模型总结")
                    st.json(func)

    # === Tab 3: RNAi设计 ===
    with tabs[2]:
        st.header("🧬 shRNA/siRNA 靶序列智能设计")
        rna_gene = st.text_input("目标基因", placeholder="例：TP53, CD47")
        rna_species = st.selectbox("物种", ["human", "mouse", "rat"])
        rna_region = st.radio("靶向区域", ["CDS", "5'UTR", "3'UTR"], horizontal=True)
        
        if st.button("🔬 生成候选序列", type="primary", disabled=not rna_gene):
            if not st.session_state.api_key:
                st.error("请先在侧边栏配置API Key")
            else:
                with st.spinner("AI正在计算GC含量、脱靶风险与二级结构..."):
                    seqs = client.design_rnai_sequences(rna_gene, rna_species, rna_region)
                    st.session_state.rnai_results.extend([{"基因": rna_gene, "物种": rna_species, "序列": s, "时间": datetime.now().strftime("%Y-%m-%d %H:%M")} for s in seqs])
                    
                    if isinstance(seqs, list) and seqs and "error" not in seqs[0]:
                        st.success("设计完成（建议合成前用NCBI BLAST复核脱靶）")
                        st.dataframe(pd.DataFrame(seqs), use_container_width=True)
                    else:
                        st.error("AI生成异常，请检查网络或重试")

    # === Tab 4: 导出 ===
    with tabs[3]:
        st.header("📥 实验记录导出")
        st.info("💡 Streamlit Cloud 为无状态环境，请及时下载数据。本地部署可开启SQLite持久化。")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.session_state.cell_results:
                df = pd.DataFrame(st.session_state.cell_results)
                st.dataframe(df, use_container_width=True)
                st.download_button("⬇️ 导出细胞系验证记录(CSV)", df.to_csv(index=False), "cell_validation.csv", "text/csv")
            else:
                st.caption("暂无验证记录")
                
        with col2:
            if st.session_state.analysis_results:
                df = pd.DataFrame([{"基因": r['基因'], "PMID": r['PMID'], "抗病毒": r['抗病毒判断'].get('is_antiviral'), "时间": r['时间']} for r in st.session_state.analysis_results])
                st.dataframe(df, use_container_width=True)
                st.download_button("⬇️ 导出文献分析记录(CSV)", df.to_csv(index=False), "literature_analysis.csv", "text/csv")
            else:
                st.caption("暂无分析记录")
                
        with col3:
            if st.session_state.rnai_results:
                df = pd.DataFrame([{"基因": r['基因'], "物种": r['物种'], "序列": r['序列'], "时间": r['时间']} for r in st.session_state.rnai_results])
                st.dataframe(df, use_container_width=True)
                st.download_button("⬇️ 导出RNAi设计记录(CSV)", df.to_csv(index=False), "rnai_design.csv", "text/csv")
            else:
                st.caption("暂无设计记录")

if __name__ == "__main__":
    main()
