import os
import re
import gzip
import json
from collections import defaultdict
from itertools import groupby
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from ..utils.helper import logger, hamming_distance, get_file_path
from ..utils.wrappers import cmd_execute

__srcdir = get_file_path(__file__)

def read_gtf(gtf):
    gene_table = []
    mt_regex = re.compile("^(MT|mt|Mt)-")
    if gtf.endswith(".gz"):
        _open = gzip.open
    else:
        _open = open
    with _open(gtf, "rt") as fh:
        for line in fh:
            if not line: continue
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if tmp[2] == "gene":
                #print(tmp[-1])
                if "gene_id" in tmp[-1]:
                    gene_id=re.findall("gene_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\";",tmp[-1])[0]
                    gene_names=re.findall("gene_name \"([A-Za-z0-9_\.\-\:/\ ()]+)\"",tmp[-1])
                    if len(gene_names)==0:
                          gene_name = gene_id
                    else:
                          gene_name = gene_names[0]

                    if tmp[1] in ("chrM", "MT", "mt"):
                        if not mt_regex.match(gene_id):
                            gene_id = f"MT-{gene_id}"
                            gene_name = f"MT-{gene_name}"
                    gene_table.append([ gene_id, gene_name ])
    return gene_table

def umi_correct(geneid_umi_dict, bc, umi_correct_detail_fh):
    counts_dict = defaultdict(lambda: [0, 0])
    for gene_id, umi_dict in geneid_umi_dict.items():
        sorted_umis = sorted(umi_dict.keys(), key=lambda x: umi_dict[x], reverse=True)
        while len(sorted_umis) > 1:
            umi_s = sorted_umis.pop()
            for umi in sorted_umis:
                if hamming_distance(umi_s, umi) == 1:
                    umi_dict[umi] += umi_dict[umi_s]
                    umi_correct_detail_fh.write(f"{bc}\t{gene_id}\t{umi_s}:{umi_dict[umi_s]}\t{umi}:{umi_dict[umi]}\n")
                    del umi_dict[umi_s]
                    break
        counts_dict[gene_id][0] = len(umi_dict.keys())
        counts_dict[gene_id][1] = sum(umi_dict.values())
    return counts_dict, geneid_umi_dict

def umi_count(reads_group, umi_correct_detail_fh):
    assigned_dict = defaultdict(lambda: defaultdict(int))
    for barcode_umi, g in groupby(reads_group, key=lambda x: x.qname):
        _, umi = barcode_umi.split("_")[:2]
        if umi == umi[0]*len(umi):   # poly
            continue
        tmp_dict = defaultdict(int)
        n = 0
        for r in g:
            n += 1
            if r.has_tag("XT"):
                XT =r.get_tag("XT").split("XT:Z:")[-1]
                if "," in XT:
                    break
                gene_id = XT.split("XT:Z:")[-1]
                tmp_dict[gene_id] = r.mapping_quality
        if len(tmp_dict) == 1:
            if tmp_dict[gene_id] == 255 or n > 1:
                assigned_dict[gene_id][umi] += 1
    counts_dict, geneid_umi_dict = umi_correct(assigned_dict, _, umi_correct_detail_fh)
    return counts_dict, geneid_umi_dict

def bam2table(bam, detail_file, counts_file, umi_correct_detail):
    sam_file = pysam.AlignmentFile(bam, "rb")

    umi_correct_detail_fh = open(umi_correct_detail, "w")
    with open(detail_file, "w") as fh1, open(counts_file, "w") as fh2:
        fh1.write("\t".join(["cellID", "geneID", "UMI", "Num"]) + "\n")
        fh2.write("\t".join(["cellID", "geneID", "UMINum", "ReadsNum"]) + "\n")
        for barcode, g in groupby(sam_file, key=lambda x: x.qname.split("_", 1)[0]):
            counts_dict, geneid_umi_dict = umi_count(g, umi_correct_detail_fh)
            for gene_id in geneid_umi_dict:
                for umi in geneid_umi_dict[gene_id]:
                    raw_umi_count = geneid_umi_dict[gene_id][umi]
                    fh1.write(f"{barcode}\t{gene_id}\t{umi}\t{raw_umi_count}\n")
            for gene_id in counts_dict:
                umi_num, reads_num = counts_dict[gene_id]
                fh2.write(f"{barcode}\t{gene_id}\t{umi_num}\t{reads_num}\n")
    umi_correct_detail_fh.close()


def write_raw_matrix(counts_file, raw_matrix_dir, gtf):
    name_df = pd.DataFrame(read_gtf(gtf), columns=["geneID", "Symbol"])

    gene_dict = name_df.reset_index().set_index("geneID")["index"].to_dict()
    row, col, data = [], [], []
    barcodes = []
    n = 0
    with open(counts_file) as fh:
        fh.readline()
        for k, g in groupby(fh, lambda x: x.split("\t")[0]):
            barcodes.append(k)
            for _ in g:
                tmp = _.split("\t")
                data.append(int(tmp[2]))
                row.append(gene_dict[tmp[1]])
                col.append(n)
            n += 1
    mat = coo_matrix((data, (row, col)), shape=(name_df.shape[0], n))

    matrix_file = os.path.join(raw_matrix_dir, "matrix.mtx.gz")
    with gzip.open(matrix_file, "w") as fh:
        mmwrite(fh, mat)

    name_df["type"] = "Gene Expression"
    features_file = os.path.join(raw_matrix_dir, "features.tsv.gz")
    name_df.to_csv(features_file, sep="\t", index=False, header=False)

    barcodes_file = os.path.join(raw_matrix_dir, "barcodes.tsv.gz")
    with gzip.open(barcodes_file, "wt") as fh:
        fh.write("\n".join(barcodes))
        fh.write("\n")

def calculate_metrics(counts_file, detail_file, filterd_barcodes_file, filterd_features_file):
    summary = defaultdict()
    summary["Estimated Number of Cells"] = 0
    summary["Fraction Reads in Cells"] = 0
    summary["Mean Reads per Cell"] = 0
    summary["Median Genes per Cell"] = 0
    summary["Median UMI Counts per Cell"] = 0
    summary["Total Genes Detected"] = 0

    barcodes = pd.read_csv(filterd_barcodes_file, header=None, sep="\t")
    summary["Estimated Number of Cells"] = barcodes.shape[0]

    df0 = pd.read_csv(
        counts_file,
        dtype={
            "cellID": "category",
            "geneID": "category",
            "UMINum": "int32",
            "ReadsNum": "int32"
        },
        sep="\t"
    )
    summary["Sequencing Saturation"] = 1 - df0.UMINum.sum()/df0.ReadsNum.sum()
    
    df0 = df0.loc[df0["cellID"].isin(barcodes[0]), :].reset_index(drop=True)
    umi_median = int(df0.groupby(["cellID"], observed=True)["UMINum"].sum().median())
    summary["Median UMI Counts per Cell"] = umi_median
    gene_total = int(df0[["geneID"]].nunique())
    summary["Total Genes Detected"] = gene_total
    gene_median = int(df0.groupby(["cellID"], observed=True)["geneID"].nunique().median())
    summary["Median Genes per Cell"] = gene_median
    
    del df0

    df = pd.read_csv(
        detail_file,
        dtype={
            "cellID": "category",
            "geneID": "category",
            "UMI": "category",
            "Num": "int32"
        },
        sep="\t"
    )
    mapped_reads_total = df["Num"].sum()

    df = df.loc[df["cellID"].isin(barcodes[0]), :].reset_index(drop=True)
    cell_reads_total = df["Num"].sum()
    cell_reads_ratio = cell_reads_total / mapped_reads_total
    summary["Fraction Reads in Cells"] = cell_reads_ratio
    rep = df["Num"]
    df = df.drop(["Num"], axis=1)
    idx = df.index.repeat(rep)
    df = df.iloc[idx].reset_index(drop=True)
    del rep, idx

    # shuffle
    df = df.sample(frac=1.0).reset_index(drop=True)
    percentage_sampling = []
    saturation_sampling = []
    median_sampling = []

    # downsample
    for n, interval in enumerate(np.array_split(np.arange(df.shape[0]), 10)):
        idx = interval[-1]
        percentage = (n + 1) / 10
        percentage_sampling.append(percentage)
        sampled_df = df.iloc[:idx]
        UMIcounts = sampled_df.groupby(
            [sampled_df["cellID"], sampled_df["geneID"], sampled_df["UMI"]],
            observed=True)["UMI"].count().reset_index(drop=True)
        saturation = ((UMIcounts[UMIcounts > 1] - 1).sum() + 0.0) / idx * 100
        saturation_sampling.append(float(saturation))

        median = sampled_df.groupby(
            [sampled_df["cellID"]],
            observed=True)["geneID"].nunique().reset_index(drop=True).median()
        median_sampling.append(int(median))

    return summary, {
        "percentage": percentage_sampling,
        "saturation": saturation_sampling,
        "median": median_sampling
    }

def count(bam, outdir, samplename, gtf, expectNum=3000, **kwargs):
    basedir = os.path.join(outdir, "step3")
    os.makedirs(basedir, exist_ok=True)

    detail_file = os.path.join(basedir, "detail.xls")
    counts_file = os.path.join(basedir, "counts.xls")
    umi_file = os.path.join(basedir, "umi.xls")

    bam2table(bam=bam, detail_file=detail_file, counts_file=counts_file, umi_correct_detail=umi_file)

    raw_matrix_dir = os.path.join(basedir, "raw_feature_bc_matrix")
    os.makedirs(raw_matrix_dir, exist_ok=True)
    logger.info("write raw matrix started!")
    write_raw_matrix(counts_file, raw_matrix_dir, gtf)
    if kwargs["forceCell"] != None:
        Rapp = os.path.join(__srcdir, "force.cell_identify.R")
        args = ["Rscript", Rapp, "-i", raw_matrix_dir, "-f", kwargs["forceCell"]]
        logger.info("force cell.")
    else:
        Rapp = os.path.join(__srcdir, "cell_identify.R")
        args = ["Rscript", Rapp, "-i", raw_matrix_dir, "-e", expectNum, "-p", 0.01]
    logger.info("call cell started!")
    cmd_execute(args, check=False)

    filterd_barcodes_file = os.path.join(
        basedir, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    filterd_features_file = os.path.join(
        basedir, "filtered_feature_bc_matrix/features.tsv.gz")

    logger.info("calculate metrics started!")
    summary_tmp, downsample = calculate_metrics(
        counts_file,
        detail_file,
        filterd_barcodes_file,
        filterd_features_file
    )

    with open(os.path.join(outdir, samplename+"_summary.json"), "r") as fh:
        summary = json.load(fh)
        # Total = int(summary["Sequencing"]["Number of Reads"].replace(",", ""))
        Total = summary["stat"]["total"]
    with open(os.path.join(outdir, samplename+"_summary.json"), "w") as fh:
        estimated_cell_num = summary_tmp["Estimated Number of Cells"]
        mean_reads_per_cell = Total / estimated_cell_num
        summary_tmp["Mean Reads per Cell"] = int(mean_reads_per_cell)
        summary["cells"] = summary_tmp
        downsample["Reads"] = [
            int(summary_tmp["Mean Reads per Cell"] * p)
            for p in downsample["percentage"]
        ]
        summary["downsample"] = downsample
        json.dump(
            summary,
            fh,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )

    logger.info("count done!")

    return os.path.join(basedir, "filtered_feature_bc_matrix")
