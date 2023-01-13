import os
import importlib
import click
from seeksoultools.rna import rna
from seeksoultools.utils import utils
from seeksoultools.utils._version import __version__
from seeksoultools.utils import helper

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=__version__, message="%(version)s")
@click.option("--debug", is_flag=True, default=False, help="debug flag.")
@click.pass_context
def cli(ctx, debug):
    ctx.ensure_object(dict)
    ctx.obj["version"] = __version__
    if debug:
        os.environ["seeksoultools_debug"] = "True"
    else:
        os.environ["seeksoultools_debug"] = "False"   
    importlib.reload(helper)

cli.add_command(rna)
cli.add_command(utils)
