"""
Functions and logic related to finding, parsing and interpreting configfiles
and other config-related stuff.
"""

class InvalidConfigError(Exception):
    pass

AVIAN_FLU_DIR = os.path.normpath(os.path.join(workflow.current_basedir, ".."))
# NOTE: `workflow.basedir` is the Snakemake entry point, i.e. the directory of the first encountered Snakefile

# load the default config which must exist.
# NOTE: despite any `--config` or `--configfile` having been already loaded into the `config` variable by this stage,
# the following `configfile: ...` expression doesn't override this, rather it recomputes the entire `config` structure
# (by merging) such that the precedence of "configfile: < --configfile < --config" is maintained.
if not os.path.exists(os.path.join(workflow.basedir, 'config.yaml')):
    raise InvalidConfigError("No default config - this is the result of an error/bug in the underlying workflow.")
configfile: os.path.join(workflow.basedir, 'config.yaml')

# load a config.yaml file if it exists in the current working directory
# (current working directory is most often the directory where you ran `snakemake` from, but can be changed via
# Snakemake's --directory arg)
if os.path.exists("config.yaml"):
    configfile: "config.yaml"


def resolve_config_path(path):
    """
    Resolve a relative *path* (string) given in a configuration value. Returns a
    function which takes a single argument *wildcards* and returns the resolved
    path with any '{x}' substrings replaced by their corresponding wildcards
    (using Snakemake's `expand`).

    Search order (first match returned):
        1. Relative to the analysis directory
        2. Relative to the `avian-flu` directory

    The returned function will raise an `InvalidConfigError` if a match is not found
    """
    if not isinstance(path, str):
        raise InvalidConfigError(f"Config path provided to resolve_config_path must be a string. Provided value: {str(path)}")

    def resolve(wildcards):
        try:
            path_expanded = expand(path, **wildcards)[0]
        except snakemake.exceptions.WildcardError as e:
            # str(e) looks like "No values given for wildcard 'subtypes'."
            raise InvalidConfigError(f"resolve_config_path called with path {path!r} however {str(e)}")

        # check if the path exists relative to the working analysis directory
        if os.path.exists(path_expanded):
            return path_expanded

        # Check if the path exists relative to the avian-flu directory
        if os.path.exists(p:=os.path.join(AVIAN_FLU_DIR, path_expanded)):
            return p

        raise InvalidConfigError(f"Unable to resolve the config-provided path {path!r}, expanded to {path_expanded!r} after filling in wildcards. "
            f"The following directories were searched:\n"
            f"\t1. {os.path.abspath(os.curdir)} (current working directory)\n"
            f"\t2. {AVIAN_FLU_DIR} (the avian-flu repo)\n")

    return resolve

def script(path):
    """
    Resolve a provided script *path* (string)
    
    Search order (first match returned):
        1. Relative to the 'scripts' directory in the avian-flu repo (`AVIAN_FLU_DIR`)
        2. Relative to the avian-flu repo (`AVIAN_FLU_DIR`)

    An `InvalidConfigError` is raised if a match is not found
    """

    if os.path.exists(p:=os.path.join(AVIAN_FLU_DIR, "scripts", path)):
        return p
    
    if os.path.exists(p:=os.path.join(AVIAN_FLU_DIR, path)):
        return p

    raise InvalidConfigError(f"Unable to resolve the provided script {path!r}. "
        f"The following directories were searched:\n"
        f"\t1. {os.path.join(AVIAN_FLU_DIR, 'scripts')}\n"
        f"\t2. {AVIAN_FLU_DIR}\n")
