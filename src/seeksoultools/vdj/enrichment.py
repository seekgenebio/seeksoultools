from pathlib import Path
from collections import defaultdict
import pysam
import dnaio
import re
from ..utils.wrappers import cmd_execute
from ..utils.helper import logger

def STAR_build_wrapper(
    ref: str,
    outdir: str,
    core:int=4,
    star_path:str="STAR",
)->str:

    cmd = (f"{star_path} --runMode genomeGenerate --runThreadN {core}  " 
           f"--genomeFastaFiles {ref}  --genomeDir {outdir} ")
    cmd_execute(cmd, check=True)


def format_ref(ref: Path, outdir: Path) -> Path:
    """
    格式化参考序列的名称并写入到输出目录。

    该函数读取一个参考序列文件，格式化其中的序列名称，然后将这些序列写入到指定的输出目录中的一个新的FASTA文件里。
    
    参数:
    - ref: Path类型的参考序列文件路径。
    - outdir: Path类型的输出目录路径。
    
    返回:
    - Path: 新格式的参考序列文件的路径。
    """
    # 编译一个正则表达式，用于匹配并格式化序列名称
    regex = re.compile(r'((TR|IG)\S+)[\|\s:]')
    
    # 使用上下文管理器打开输入参考序列文件和输出目录中的新FASTA文件
    # 输入文件以适合的格式打开，输出文件也以写入模式打开

    ref_l = outdir / "raw_ref.fa"


    with dnaio.open(ref, fileformat="fasta") as fh_in, \
         dnaio.open(ref_l, mode="w", fileformat="fasta") as fh_out:
        # 初始化计数器，用于跟踪处理的记录数
        n = 0
        # 遍历输入文件中的每一条序列记录
        for record in fh_in:
            # 增加记录计数器
            n += 1
            # 搜索当前记录的名称中是否包含匹配正则表达式的部分
            m = regex.search(record.name)
            if m:
                # 如果找到匹配，格式化记录的名称并写入到输出文件
                record.name = f"{m.group(1)}_{n}"
                fh_out.write(record)
    return ref_l
def count_chain(bam: Path):
    """
    统计.bam文件中每个参考序列的读取次数。

    该函数处理一个.bam文件，统计每个参考序列（通过其名称的前三个字符识别）被映射为唯一读取或多重映射读取的次数。
    如果读取是唯一映射的，直接计数。如果是多重映射的，只计算其中的一个映射。

    参数:
    sam (Path): .sam文件的路径。

    返回:
    defaultdict: 一个字典，键为参考序列的前三个字符（'chain'），值为该参考序列的计数。
    """
    # 初始化默认字典d，用于存储最终的链计数
    d = defaultdict(int)
    # 初始化临时默认字典tmp，用于临时存储每个查询名称的链计数
    tmp = defaultdict(lambda: defaultdict(int))
    _default_verbosity = pysam.set_verbosity(0)
    # 打开.bam文件并处理其中的每一行读取
    with pysam.AlignmentFile(bam) as bamfh:
        for r in bamfh:
            # 提取读取的参考序列名称的前三个字符作为链标识
            chain = r.reference_name[:3]
            # 如果读取是唯一映射的（NH标签值为1），在d中增加该链的计数
            if r.get_tag("NH")==1:
                d[chain] += 1
            else:
                # 如果读取是多重映射的，在tmp中增加相应的计数
                tmp[r.qname][chain] += 1
    
    # 处理多重映射读取，决定如何计入tmp中的数据
    for k, v in tmp.items():
        # 如果一个读取只映射到一个链，将其计数加到d中
        if len(v)==1:
            chain, _count = v.popitem()
            d[chain] += 1
        else:
            # 如果有多个映射，选择最频繁的那个进行计数
            sorted_v = sorted(v.items(), key=lambda x: x[1], reverse=True)
            if sorted_v[0][1] > sorted_v[1][1]:
                chain, _count = sorted_v[0]
                d[chain] += 1
            else:
                # 如果最大计数相等，表示不确定，增加'unknown'链的计数
                d["unknown"] += 1
    
    # 返回最终的链计数字典
    return d
def enrichment_qc(
    fq: list,
    chain: str,
    fa: Path,
    samplename:str,
    outdir: Path,
    core: int=4,
    star_path: str="STAR"
):

    ref_dir = outdir / "fa_ref"
    ref_dir.mkdir(exist_ok=True, parents=True)
    ref_l = format_ref(fa, ref_dir)

    STAR_build_wrapper(ref_l, ref_dir, core, star_path)

    prefix = f"{outdir}/{samplename}_"
    from ..utils.wrappers import STAR_wrapper
    STAR_wrapper(
        fq=list(fq), genomeDir=ref_dir, prefix=prefix, scoremin=0.33,
        matchnmin=0.33, core = core, star_path=star_path        
    )
    bamfile = f"{prefix}Aligned.out.bam"
    d = count_chain(bamfile)
    with open(prefix + "chain_summary.txt", "w") as fh:
        for k, v in d.items():
            fh.write(f"{k}: {v}\n")
