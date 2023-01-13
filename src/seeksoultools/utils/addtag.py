import pysam
import argparse
import re
from collections import defaultdict


def sep_umi(umi_cand):
    return umi_cand.split(":")


def get_corrected_umi(umi_f):
    #read corrected umi file into dict
    with open(umi_f, "r") as f:
        umi_corrected_dict = defaultdict(lambda: defaultdict(dict))
        for l in f:
            bc, ensemble_id, umi_pair_1, umi_pair_2 = l.strip().split("\t")
            umi_1, umi_count_1 = sep_umi(umi_pair_1)
            umi_2, umi_count_2 = sep_umi(umi_pair_2)
            umi_corrected_dict[bc][ensemble_id].update({umi_1:int(umi_count_1), umi_2:int(umi_count_2)})
    return umi_corrected_dict


def parse_alt(alt, bc, regex = re.compile(r'(\d+)([NACGT])')):
    #get raw barcode based on alteration signiture
    #alt: M means match
    #     1N3T means corrected barcode replaced N at 1st place and T at 3rd place(0 based)
    if alt == "M":
        bc_raw = bc
    else:
        groups = regex.findall(alt)
        for pos, base in groups:
            bc = bc[:int(pos)] + base + bc[int(pos)+1:]
        bc_raw = bc
    return bc_raw

def add_tag(inbam:str, outbam:str, umi_f:str=None):
    """bam添加BC、BR、UC、UR标签，用于细胞速率

    """

    samfile = pysam.AlignmentFile(inbam, "rb")
    samfile_mod = pysam.AlignmentFile(outbam, "wb", template=samfile)
    if umi_f:
        umi_correct_flag = True
        umi_corrected_dict = get_corrected_umi(umi_f)
    
    for read in samfile:
        #read samfile
        bc, umi_raw, alt, readname = read.query_name.split("_", 3)
        bc_raw = parse_alt(alt, bc)

        if umi_correct_flag and read.has_tag("XT"):
            # check if umi correction is needed
            ensemble_id = read.get_tag("XT")
            #try if ensemble_id is contained in dict
            #if contained, do umi correction
            try:
                umi_dict = umi_corrected_dict[bc][ensemble_id]
                umi = umi if umi_raw not in umi_dict.keys() else max(umi_dict, key=umi_dict.get)
            except:
                umi = umi_raw
        else:
            umi = umi_raw

        read.set_tags([("BC", bc), ("BR", bc_raw), ("UC", umi), ("UR", umi_raw)])
        #change read.query_name by removing barcode, umi and alteration
        read.query_name = readname
        samfile_mod.write(read)
        
    samfile_mod.close()
    samfile.close()    
    return 
