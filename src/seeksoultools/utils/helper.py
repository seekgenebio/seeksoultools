import os
import sys
import re
import json
import subprocess
from loguru import logger
from xopen import xopen
from cutadapt.adapters import BackAdapter, Sequence
from .chemistry import CHEMISTRY
from ._version import __version__

# setup logger
logger.remove()
if os.environ.get("seeksoultools_debug", "False")=="True":
    logger.add(sys.stderr, level="DEBUG")
else:
    logger.add(sys.stderr, level="INFO")

def include_introns_callback(ctx, param:str, value:bool):
    # --include-introns
    if value:
        return "transcript"
    else:
        return "exon"

def abs_path_callback(ctx, param:str, value:bool):
    return os.path.abspath(value)

def get_file_path(filename:str):
    return os.path.dirname(os.path.abspath(filename))

def hamming_distance(s1, s2):
    return len([(i, j) for i, j in zip(s1, s2) if i != j])

def parse_structure(string:str) -> tuple:
    """解析接头结构

    使用字母B、L、U、X和T以及数字表示reads结构。
    B表示barcode部分碱基；
    L表示linker部分碱基；
    U表示umi部分碱基；
    X表示任意碱基，用于占位；
    T表示T碱基；
    字母后数字表示碱基长度。

    Args:
        string: 接头结构描述

    Returns:
        返回二维tuple,内容为按顺序的各部分结构和长度。
        例如：
            当string是B8L8B8L10B8U8,返回:
            (('B', 8), ('L', 8), ('B', 8), ('L', 10), ('B', 8), ('U', 8))
    """
    regex = re.compile(r'([BLUXT])(\d+)')
    groups = regex.findall(string)
    return tuple([(_[0], int(_[1])) for _ in groups])

def read_file(file_list: list) -> dict:
    """准备白名单set
        Args:
            file_list: 每段白名单文件的路径
    """
    wl_dict = dict()
    for i, wl_file in enumerate(file_list):
        white_list = set()
        with xopen(wl_file, "r") as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                la = l.strip()
                if not la:
                    continue
                white_list.add(la)
        wl_dict[i] = white_list
    return wl_dict


def get_new_bc(bc:str, white_list:set)->set:
    """返回原始barcode各位置错配后的set与白名单set的交集"""

    BASE_LIST = ["T", "C", "G", "A"]
    mm_dict = dict()
    for i, c in enumerate(bc):
        if c == "N":
            mm_dict = { bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST }
            break  
        else:
            mm_dict.update({ bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST if base!=c })
            
    bc_set = set(mm_dict.keys()).intersection(white_list)
    # return {k: mm_dict[k] for k in bc_set}
    return bc_set

class AdapterFilter:
    """过滤接头"""
    def __init__(self, adapters:list):
        self.adapters = [BackAdapter(sequence=_) for _ in adapters]
    
    def filter(self, sequence:Sequence)->tuple:
        flag = False
        for _ in self.adapters:
            m = _.match_to(sequence.sequence)
            if m:
                flag = True
                sequence =  m.trimmed(sequence)
        return flag, sequence

class QcStat:
    """汇总统计"""
    def __init__(self):
        self.data = {}

    def update(self, **d):
        if not self.data:
            self.data = d
        else:
            for k, v in d.items():
                self.data[k] += v

    @staticmethod
    def _sort_gc(d):
        idx_max = max([k[0] for k in d])
        return {
            b: [d.get((i, b), 0) for i in range(idx_max+1)] for b in 'ATCGN'
        }

    @staticmethod
    def _sort_q(d, phred=33):
        idx_max = max([k[0] for k in d])
        q_max = max([ord(k[1])-phred for k in d])
        return {
            i: [d.get((i, chr(q+phred)), 0) for q in range(q_max+1)] for i in range(idx_max+1)
        }

    def save(self, path='summary.json'):
        tmp = {'__version__': __version__}
        for k in self.data:
            if k.endswith('_gc'):
                tmp[k] = self._sort_gc(self.data[k])
            elif k.endswith('_q'):
                tmp[k] = self._sort_q(self.data[k])
            else:
                tmp[k] = dict(self.data[k])
        with open(path, 'w') as fh:
            json.dump(tmp, fh, indent=4)

def cut_fq(fq:str, outdir_tmp:str, reads_num:int=100000):
    
    fq_tmp = os.path.join(outdir_tmp, os.path.basename(fq))
    cmd_cut = f"zcat {fq}|head -n {reads_num*4} | gzip > {fq_tmp}"
    p = subprocess.run(cmd_cut, shell=True)
    return fq_tmp
    
def call_process(fq1, fq2, outdir_tmp, chem, barcode:list=[], shift:str=True, shift_pattern:str='A',
                 structure:str='B8L8B8L10B8U12T15', linker: list=[], adapters:list=["AAAAAAAAAAAA", ],
                 core:int=4, do_B_correction=True, do_L_correction=True,
                 use_multi=True, **kwargs):

    from .barcode import process_barcode

    #create test dir/file
    outdir_chem = os.path.join(outdir_tmp, chem)
    os.makedirs(outdir_chem, exist_ok=True)
    fq_out = os.path.join(outdir_chem, "test_r2.fq.gz")
    fqout_multi = os.path.join(outdir_chem, "test_multi.fq.gz")

    #params needed
    r1_structure = parse_structure(structure)
    barcode_wl_dict = read_file(barcode)
    linker_wl_dict = read_file(linker)

    summary_dict = process_barcode(
        fq1, fq2, fq_out, fqout_multi, r1_structure, shift, shift_pattern, barcode_wl_dict,
        linker_wl_dict, adapters, do_B_correction, do_L_correction, use_multi
    )
    return summary_dict["stat"]["valid"]/summary_dict["stat"]["total"] * 100, fq_out
    
def r1_structure_validation(fq1:list, fq2:list, samplename: str, outdir:str, **kwargs):
    """
    validate read1 strucutre, if valid_rate less than 50.0%, test for other pre defined chemistry.
    """

    logger.info("Starting r1 structure validation")
    os.makedirs(f"{outdir}/step1/_fq_test", exist_ok=True)
    outdir_tmp = os.path.join(outdir, f"step1/_fq_test")
    _valid_rate = {}

    for _fq1, _fq2 in zip(fq1, fq2):
        fq1_tmp = cut_fq(_fq1, outdir_tmp)
        fq2_tmp = cut_fq(_fq2, outdir_tmp)

        if not kwargs['chemistry']:
            kwargs['chemistry'] = "custom"
            valid_rate, fq_out = call_process(fq1_tmp, fq2_tmp, outdir_tmp, "custom", **kwargs)
            _valid_rate[f"{_fq1},{_fq2}"] = valid_rate
            logger.info(f"{_fq1},{_fq2}\nValid barcode rate for custom is about {valid_rate:.3f}%")
        else:
            valid_rate, fq_out = call_process(fq1_tmp, fq2_tmp, outdir_tmp, kwargs["chemistry"], **kwargs)
            _valid_rate[f"{_fq1},{_fq2}"] = valid_rate
            logger.info(f"{_fq1},{_fq2}\nValid barcode rate for --chemistry {kwargs['chemistry']} is about {valid_rate:.3f}%")

    for k, v in _valid_rate.items():
        if v < 30.0:
            for chem, options in CHEMISTRY.items():
                if chem.startswith("__"):
                    continue
                if chem == kwargs['chemistry']:
                    continue
                _kwargs = kwargs
                _kwargs.update(options)
                valid_rate, fq_out = call_process(fq1_tmp, fq2_tmp, outdir_tmp, chem, **_kwargs)
                logger.info(f"{k}\nValid barcode rate for --chemistry {chem} is about {valid_rate:.3f}%")

    return
