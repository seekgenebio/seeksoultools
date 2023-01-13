import click

@click.group(help="utils.")
@click.pass_context
def utils(ctx):
    pass

@utils.command(help="gene type summary.")
@click.option('--gtf', type=click.Path(), required=True, help="gtf file.")
@click.option('--feature', default='gene', show_default=True, help="feature, e.g., gene, transcript")
@click.option('--key', default='gene_type', show_default=True, help='attribution key, e.g., gene_type')
def gtfstat(gtf, feature, key='gene_biotype'):
    from .mkref import gtfstat
    gtfstat(gtf, feature, key)

@utils.command(help="gtf filter.")
@click.option('--gtf', type=click.Path(), required=True, help="gtf file.")
@click.option('--biotype',  multiple=True, help="biotype to keep,  can specify multiple times.")
@click.option('--key', default='gene_type', show_default=True, help='attribution key, e.g., gene_type')
def gtffilter(gtf, biotype, key='gene_type'):
    from .mkref import gtffilter
    gtffilter(gtf, biotype, key)

@utils.command(help="make ref.")
@click.option('--fa', type=click.Path(), required=True, help="fa file.")
@click.option('--gtf', type=click.Path(), required=True, help="gtf file.")
@click.option('--genomeDir', 'genomeDir', type=click.Path(), required=True, help="genomeDir.")
@click.option('--runThreadN', 'runThreadN', default=8, help="runThreadN.")
@click.option('--star_path', default='STAR', help="STAR path.")
@click.option('--star_opt', default='', help="STAR opt.")
def mkref(fa, gtf, genomeDir, runThreadN=8, star_path='STAR', star_opt=''):
    from .mkref import mkref
    mkref(fa, gtf, genomeDir, runThreadN, star_path, star_opt)

@utils.command(help="addtag.")
@click.option("--inbam", type=click.Path(), required=True, help="inbam file, barcode and umi include in reads header.")
@click.option("--outbam", type=click.Path(), required=True, help="outbam file.")
@click.option("--umifile", type=click.Path(), help="fa file.")
def addtag(inbam, outbam, umifile=None):
    from .addtag import add_tag
    add_tag(inbam, outbam, umifile)
