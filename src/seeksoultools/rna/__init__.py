import os
import click
from ..utils.chemistry import CHEMISTRY, ADAPTERS
from ..utils.helper import logger, include_introns_callback


_steps = {
    "step1": [],
    "step2": ["STAR", "SortByPos", "FeatureCounts", "SortByName"],
    "step3": [],
    "step4": []
}

@click.group(help="quantifies singlecell gene expression")
@click.option("--steps", default=None, type=click.Path(), help="")
@click.pass_obj
def rna(obj, steps):
    if  steps:
        from yaml import load
        try:
            from yaml import CLoader as Loader
        except ImportError:
            from yaml import Loader
        with open(steps) as fh:
            obj["steps"] = load(fh, Loader=Loader)
    else:
        obj["steps"] = _steps


@rna.command(help="extract cell barcode and umi.")
@click.option("--fq1", "fq1", required=True, type=click.Path(), multiple=True, help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(), multiple=True, help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--shift", is_flag=True, default=False, show_default=True, help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, help="Barcode white list file, can specify multiple times.")
@click.option("--structure", help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, help="Linker white list file, can specify multiple times.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--chemistry", type=click.Choice(["DDV1", "DDV2", "DD5V1", "MM", "MM-D", "DD-Q"]), help="DDV1, DDV2, DD5V1, MM, MM-D, DD-Q.")
@click.pass_obj
def step1(obj, **kwargs):
    if kwargs["chemistry"]:
        kwargs.update(CHEMISTRY[kwargs["chemistry"]])

    # 设置默认adapters
    kwargs.update({"adapters": ADAPTERS})
    from ..utils.barcode import barcode_main
    from ..utils.helper import r1_structure_validation

    # 测试valid barcode
    r1_structure_validation(**kwargs)
    barcode_main(**kwargs)

@rna.command(help="align reads to genome.")
@click.option("--fq", required=True, help="Read2 fq file")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(), help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.option("--sc5p",is_flag=True,default=False,show_default=True,help="If set, the single cell data is considered as 5' data.")
@click.pass_obj
def step2(obj, **kwargs):
    from .step2 import align
    align(**kwargs)

@rna.command(help="quantifies.")
@click.option("--bam", required=True, help="Bam file which contain alignment info.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def step3(obj, **kwargs):
    from .step3 import count
    count(**kwargs)

@rna.command(help="seurat.")
@click.option("--matrix", type=click.Path(), help="Feature-barcode matrix.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--dims", default=15, show_default=True, help="Number of dimension used as input in dimensional reduction.")
@click.option("--minpct", default=0.1, show_default=True, help="Minimum percentage that a feature to be detected in either of the two groups of cells.")
@click.option("--logfc", default=0.25, show_default=True, help="Limit testing to genes to this number of fold difference(log_scale) between the two groups of cells")
@click.pass_obj
def step4(obj, **kwargs):
    from .step4 import do_seurat
    do_seurat(**kwargs)


@rna.command(help="report.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.pass_obj
def report(obj, **kwargs):
    from .report import report
    report(**kwargs)


@rna.command(help="run all steps.")
@click.pass_obj
@click.option("--fq1", "fq1", required=True, type=click.Path(), multiple=True, help="Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(), multiple=True, help="Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--shift", is_flag=True, default=False, help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, help="Barcode white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--structure", help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, help="Linker white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.option("--sc5p", is_flag=True, default=False, show_default=True, help="If set, the single cell data is considered as 5' data.")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(), help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--chemistry", type=click.Choice(["DDV1", "DDV2", "DD5V1", "MM", "MM-D", "DD-Q"]), help="DDV1, DDV2, DD5V1, MM, MM-D, DD-Q.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
def run(obj, **kwargs):
    if kwargs["chemistry"]:
        kwargs.update(CHEMISTRY[kwargs["chemistry"]])

    # 设置默认adapters
    kwargs.update({"adapters": ADAPTERS})

    kwargs["outdir"] = os.path.join(kwargs["outdir"], kwargs["samplename"])
    if "step1" in obj["steps"]:
        from ..utils.barcode import barcode_main
        from ..utils.helper import r1_structure_validation
        # test valid barcode
        r1_structure_validation(**kwargs)
        
        barcode_main(**kwargs)
    fq = os.path.join(kwargs["outdir"], "step1", f"{kwargs['samplename']}_2.fq.gz")
    kwargs["fq"] = fq

    if "step2" in obj["steps"]:
        from .step2 import align
        align(flags=_steps["step2"], **kwargs)
    bam = os.path.join(kwargs["outdir"], "step2", "featureCounts",  f"{kwargs['samplename']}_SortedByName.bam")
    kwargs["bam"] = bam

    if "step3" in obj["steps"]:
        from .step3 import count
        count(**kwargs)

    matrix = os.path.join(kwargs["outdir"], "step3", "filtered_feature_bc_matrix")
    kwargs["matrix"] = matrix
    if "step4" in obj["steps"]:
        from .step4 import do_seurat
        do_seurat(**kwargs)

        from .report import report
        report(**kwargs)
