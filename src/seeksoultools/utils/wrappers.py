
import os
import shutil
from subprocess import run
from .helper import logger

def cmd_execute(
    args, check:bool=False, text:bool=True,
    capture_output:bool=True,env:os.environ=None
):
    """execute cmd, output cmd in debug mode."""
    if isinstance(args, list):
        args = [str(_) for _ in args]
        logger.debug(" ".join(args))
        _call = run(args, check=check, text=text, capture_output=capture_output, env=env)
    elif isinstance(args, str):
        _call = run(args, shell=True, check=check, text=text, capture_output=capture_output, env=env)
    return _call

def STAR_wrapper(
    fq:str, genomeDir:str, prefix:str,
    core:int=4, star_path:str="STAR"
    ):
    """wrapper for STAR"""
    args = [
        star_path, '--runThreadN', core, '--limitOutSJcollapsed', 5000000,
        '--genomeDir', genomeDir, '--readFilesIn', fq, '--readFilesCommand',
        'zcat', '--outFileNamePrefix', prefix, '--outSAMtype', 'BAM',
        'Unsorted'
    ]
    _ = cmd_execute(args)
    return f'{prefix}Aligned.out.bam', f'{prefix}Log.final.out'


def samtools_sort_wrapper(
    inbam:str, outbam:str, byname:bool=False, clean:bool=True,
    core:int=4, samtools_path:str="samtools"
)->str:
    """wrapper for samtools sort"""
    args = [
        samtools_path, "sort", "-O", "BAM", "-@",  core, "-o", outbam, inbam
    ]
    if byname:
        args.insert(2, '-n')

    _call = cmd_execute(args, check=True)

    if _call.returncode == 0 & clean:
        os.remove(inbam)
    return outbam

def unzip_wrapper(gzfile:str, outfile:str):
    args = [
        "gzip", "-dc", gzfile, ">", outfile
    ]
    cmd_execute(args, check=True)
    return outfile

def qualimap_wrapper(
    bam:str, gtf:str, outdir:str, SC5P:bool, qualimap_path:str="qualimap"
    ):
    strand = {
        'f': 'strand-specific-forward',
        'r': 'strand-specific-reverse',
        'non': 'non-strand-specific'
    }

    # 3'
    s = "f"

    # 5'
    if SC5P:
        s="r"

    # gtf.gz?
    if gtf.endswith('.gz'):
        gtf_file = os.path.join(outdir, 'tmp.gtf')
        unzip_wrapper(gtf, gtf_file)
    else:
        gtf_file = gtf

    args = [
        qualimap_path, 'rnaseq', '-outformat', 'PDF', '-outdir', outdir, '-bam',
        bam, '-gtf', gtf, '-p', strand[s], '--java-mem-size=8G'
    ]
    my_env = os.environ.copy()
    if 'DISPLAY' in my_env:
        del my_env['DISPLAY']
    
    cmd_execute(args, check=False, env=my_env)
    return os.path.join(outdir, "rnaseq_qc_results.txt")

def featureCounts_wrapper(
    bam:str, gtf:str, samplename:str, outdir:str, region:str, SC5P:bool,
    core:int=4, featureCounts_path="featureCounts", **kwargs
    ):
    outcounts = os.path.join(outdir, 'counts.txt')

    # 3'
    s="1"

    # 5'
    if SC5P:
        s="2"

    args = [
        featureCounts_path, "-T", core, "-t", region , "-s", s, "-M", "-O", "-g", "gene_id",
        "--fracOverlap", 0.5, "-a", gtf, "-o", outcounts, "-R", "BAM", bam
    ]

    cmd_execute(args, check=True)
    return os.path.join(outdir, f"{samplename}_SortedByCoordinate.bam.featureCounts.bam")


def bowtie2_wrapper(
    fq:str, ref:str, bam:str, core:int=1, local_mode=True,
    bowtie2_path:str="bowtie2", samtools_path:str="samtools"
)->str:
    """使用bowtie2比对，并用samtools转bam"""
    local_option = ""
    if local_mode:
        local_option = "--local"
    args = (
        f"{bowtie2_path} -p {core} -x {ref} -U {fq} {local_option}|"
        f"samtools view -b > {bam}"
    )
    cmd_execute(args, check=True)
    return bam

def bowtie2_build_wrapper(
    ref:str,
    outdir:str,
    bowtie2_path:str="bowtie2",
)->str:
    os.makedirs(outdir, exist_ok=True)
    if os.path.dirname(os.path.abspath(ref)) != os.path.dirname(os.path.abspath(outdir)):
        ref_l = shutil.copy(ref, outdir)
    args = [f"{bowtie2_path}-build", ref_l, ref_l]
    cmd_execute(args, check=True)
    return ref_l

def igblast_wrapper(
    input:str,
    output:str,
    core:int=4,
    organism:str="",
    igblastn_path:str="igblastn",
    auxiliary_data:str="",
    internal_data:str="",
):
    outdir = os.path.dirname(output)
    cmd = (
        f"cd {outdir}; ln -s {internal_data} .; "
        f"{igblastn_path} -query {input} -auxiliary_data {auxiliary_data} -outfmt 19 "
        f"-germline_db_J internal_data/human/human_J -germline_db_D internal_data/human/human_D "
        f"-germline_db_V internal_data/human/human_V -c_region_db internal_data/human/human_C "
        f"-num_threads {core} -organism {organism} -out {output}"
    )
    cmd_execute(cmd)
    return output
