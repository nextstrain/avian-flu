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

def resolve_path(path, wildcards):
    """
    Resolves a relative *path* (string) by first expanding wildcards ('{x}'
    substrings) using Snakemake's `expand` functionality and then returning
    the appropriate path to the file by searching various locations.

    Search order (first match returned):
        1. Relative to the analysis directory
        2. Relative to the `avian-flu` directory

    Will raise an `InvalidConfigError` if a suitable file is not found

    NOTE: Consider using the more ergonomic `resolve_config_path` instead of this
    function
    """
    if not isinstance(path, str):
        raise InvalidConfigError(f"Config path provided to resolve_path must be a string. Provided value: {str(path)}")

    try:
        path_expanded = expand(path, **wildcards)[0]
    except snakemake.exceptions.WildcardError as e:
        # str(e) looks like "No values given for wildcard 'subtypes'."
        raise InvalidConfigError(f"resolve_path called with path {path!r} however {str(e)}")

    # check if the path exists relative to the working analysis directory
    if os.path.exists(path_expanded):
        return path_expanded

    # Note for developers: For some workflows it may be better to now search for
    # the path within the workflow directory itself (`workflow.current_basedir`)
    # (in addition to the AVIAN_FLU_DIR)

    # Check if the path exists relative to the avian-flu directory
    if os.path.exists(p:=os.path.join(AVIAN_FLU_DIR, path_expanded)):
        return p

    raise InvalidConfigError(f"Unable to resolve the config-provided path {path!r}, expanded to {path_expanded!r} after filling in wildcards. "
        f"The following directories were searched:\n"
    f"\t1. {os.path.abspath(os.curdir)} (current working directory)\n"
    f"\t2. {AVIAN_FLU_DIR} (the avian-flu repo)\n")

def is_scalar(x):
    return any([isinstance(x, t) for t in [float, int, bool, str]])

def resolve_config_value(*rule_parts, sep="/"):
    """
    Resolve a config value defined by the *rule_parts* of the config, e.g.
    rule_parts = ['filter', 'min_length'] then we expect a scalar or a
    dictionary to be present at config['filter']['min_length']. If the config
    value is to be used as a path please use the `resolve_config_path` function
    instead.

    If a scalar then that value is returned, i.e. it's always the same no matter
    what the wildcards are. Note that this scalar may have wildcards in it which will
    be filled in by Snakemake when the relevant rule is evaluated.

    If a dictionary we search it for the relevant value by finding the closest matching
    key once wildcards have been considered. For instance in a three-tiered wildcard pipeline
    such as this with subtype, segment and time we interpret these as ordered in specificity.
    For instance, if only one wildcard value is specified in the config (the others are '*')
    then matching subtype is more specific than segment. Given example 
    wildcard values of {subtype=h5nx, segment=pb2, time=2y} then we have a search order of:
    - 'h5nx/pb2/2y'   ─ all 3 wildcard values specified
    - 'h5nx/pb2/*'    ┐
    - 'h5nx/*/2y'     ├ 2/3 wildcard values specified
    - '*/pb2/2y'      ┘
    - 'h5nx/*/*'      ┐
    - '*/pb2/*'       ├ 1/3 wildcard values specified
    - '*/*/2y'        ┘
    - '*/*/*'         ─ default / fall-back
    and the first key present in the config is used.

    Examples:
        params:
            correction = resolve_config_value(['traits', 'sampling_bias_correction'])
    """
    try:
        config_lookup = config
        for i,rule_key in enumerate(rule_parts): # or use functools.reduce etc
            config_lookup = config_lookup[rule_key]
    except KeyError:
        raise InvalidConfigError('Config missing entire entry for config'+''.join(['["'+rule_parts[j]+'"]' for j in range(0,i+1)]))

    def resolve(wildcards):
        if is_scalar(config_lookup):
            return config_lookup

        if not isinstance(config_lookup, dict):
            raise InvalidConfigError(f"ERROR: config under {'.'.join(rule_parts)} must be a scalar value or a dictionary")

        wild_keys = ['subtype', 'segment', 'time'] # workflow specific
        search_keys = [                            # workflow independent, as long as there are 3 or fewer wildcard categories
            sep.join([wildcards[k] for k in wild_keys]),
            *([sep.join(['*' if i==k else wildcards[key] for k,key in enumerate(wild_keys)])
                for i in range(len(wild_keys)-1, -1, -1)] if len(wild_keys)>=2 else []),
            *([sep.join(['*' if i!=k else wildcards[key] for k,key in enumerate(wild_keys)])
                for i in range(0, len(wild_keys))] if len(wild_keys)==3 else []),
            sep.join(['*']*len(wild_keys))
        ]

        for key in search_keys:
            if key in config_lookup:
                return config_lookup[key]
        msg  =  'Config structure incorrect or incomplete for config'+''.join(['["'+rule_parts[j]+'"]' for j in range(0,i+1)])
        msg += f'\n\tThe dictionary is missing a matching key for the current target of {search_keys[0]!r}, or a fuzzy match (i.e. using "*" placeholders)'
        msg +=  '\n\tP.S. If you want to use a single value across all builds then set a scalar value (number, string, boolean)'
        raise InvalidConfigError(msg)
    return resolve

def resolve_config_path(*fields):
    """
    A helper function intended to be used as directly as a value within
    snakemake input/params blocks.

    Given an array of config *fields* (keys), we return a function with a single
    argument *wildcards* which resolves to a appropriate local file path
    (string). These are accomplished via calls to helper functions
    `resolve_config_value` and `resolve_path`, respectively; see the docstrings
    of those functions for more details.

    Examples:
        input:
            colors = resolve_config_path(['colors', 'hardcoded']),
            lat_longs = resolve_config_path(["lat_longs"]),
    """
    assert all([isinstance(f,str) for f in fields]), \
        "Arguments to `resolve_config_path` must be strings"

    def resolve(wildcards):
        raw_value = resolve_config_value(*fields)(wildcards)
        if not raw_value: # falsey -> don't resolve to a path!
            return ""
        return resolve_path(raw_value, wildcards)

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


def as_list(x):
    return x if isinstance(x, list) else [x]


def expand_target_patterns():
    """
    iteratively create workflow targets from config.builds for the `all` rule
    you can over-ride this by specifying targets (filenames) on the command line
    """
    targets = []

    if not isinstance(config.get('target_patterns', None), list) and not isinstance(config.get('target_patterns', None), str):
        raise InvalidConfigError('config["target_patterns"] must be defined (either as a list or a single string)')
    target_patterns = as_list(config['target_patterns'])

    if 'builds' not in config:
        raise InvalidConfigError('config["builds"] is not defined!')
    if not isinstance(config['builds'], list):
        raise InvalidConfigError('config["builds"] must be a list')

    if 'segments' in config:
        raise InvalidConfigError('config["segments"] is no longer used. Please remove it and encode the target segments within config["builds"]')

    for i,subconfig in enumerate(config['builds']):
        required_keys = ['subtype', 'segment']
        optional_keys = ['time']
        if not isinstance(subconfig, dict):
            raise InvalidConfigError(f'config["builds"][{i}] must be a dictionary!')
        if not all([k in subconfig for k in required_keys]):
            raise InvalidConfigError(f'config["builds"][{i}] must have {", ".join(required_keys)} keys')
        if not all([isinstance(v, list) or is_scalar(v) for v in subconfig.values()]):
            raise InvalidConfigError(f'config["builds"][{i}] values must all be scalars or lists')

        for subtype in as_list(subconfig['subtype']):
            for segment in as_list(subconfig['segment']):
                # Some builds (GISAID) have a time component, some (cattle-outbreak) don't
                for time in (as_list(subconfig['time']) if 'time' in subconfig else [None]):
                    for target_pattern in target_patterns:
                        if time is None and '{time}' in target_pattern:
                            raise InvalidConfigError(f'target pattern {target_pattern!r} specifies time, but config["builds"][{i}] doesn\'t!')
                        target = target_pattern.format(subtype=subtype, segment=segment, time=time)
                        if target not in targets:
                            targets.append(target)

    return targets


# REMOVE BELOW LINE # ------------------------------------------------------------------------------------------------


def old_resolve_config_value_old(*rule_parts, fallback="FALLBACK"):
    """
    Note that the underlying algorithm for finding the config value,
    and the config syntax itself, is going to be significantly
    modified in <https://github.com/nextstrain/avian-flu/pull/104>

    Given an array of config *fields* (keys), we return a function with a single
    argument *wildcards* which resolves to the config value (any type). If the
    config value is to be used as a path please use the `resolve_config_path`
    function instead.

    Examples:
        params:
            correction = resolve_config_value(['traits', 'sampling_bias_correction'])
    """
    def resolve(wildcards):
        rule_name = rule_parts[0]
        assert rule_name in config, f"Config missing top-level {rule_name} key"
        if len(rule_parts)==1:
            # Ths form of config syntax cannot use the dictionary-based wildcard-dependent syntax
            # As per the docstring, this will be redone in <https://github.com/nextstrain/avian-flu/pull/104>
            assert isinstance(config[rule_name], str), f"Invalid config spec for config.{rule_name}"
            return config[rule_name]
        rule_key = rule_parts[1]
        assert rule_key in config[rule_name], f"Config missing entry for {rule_name}.{rule_key}"
        try:
            return config[rule_name][rule_key][wildcards.subtype][wildcards.time]
        except KeyError:
            assert fallback in config[rule_name][rule_key], f"config.{rule_name!r}.{rule_key!r} either needs " \
                f"an entry for {wildcards.subtype!r}.{wildcards.time!r} added or a (default) {fallback!r} key."
            return config[rule_name][rule_key][fallback]
    return resolve

def old_resolve_config_path_old(*fields):
    """
    A helper function intended to be used as directly as a value within
    snakemake input/params blocks.

    Given an array of config *fields* (keys), we return a function with a single
    argument *wildcards* which resolves to a appropriate local file path
    (string). These are accomplished via calls to helper functions
    `resolve_config_value` and `resolve_path`, respectively; see the docstrings
    of those functions for more details.

    Examples:
        input:
            colors = old_resolve_config_path_old(['colors', 'hardcoded']),
            lat_longs = old_resolve_config_path_old(["lat_longs"]),
    """
    assert all([isinstance(f,str) for f in fields]), \
        "Arguments to `resolve_config_path` must be strings"

    def resolve(wildcards):
        raw_value = old_resolve_config_value_old(*fields)(wildcards)
        if not raw_value: # falsey -> don't resolve to a path!
            return ""
        return resolve_path(raw_value, wildcards)

    return resolve