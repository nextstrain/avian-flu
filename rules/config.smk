"""
Functions and logic related to finding, parsing and interpreting configfiles
and other config-related stuff.
"""

class InvalidConfigError(Exception):
    pass

AVIAN_FLU_DIR = os.path.normpath(os.path.join(workflow.current_basedir, ".."))
# NOTE: `workflow.basedir` is the Snakemake entry point, i.e. the directory of the first encountered Snakefile

if workflow.basedir == AVIAN_FLU_DIR:
    raise InvalidConfigError("The avian-flu phylogenetic workflows can no longer be run using the top-level Snakefile as the entrypoint. "
    "(See README for updated instructions)")

if os.path.exists(os.path.join(workflow.basedir, 'config.yaml')):
    configfile: os.path.join(workflow.basedir, 'config.yaml')
