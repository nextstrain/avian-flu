"""
Functions and logic related to finding, parsing and interpreting configfiles
and other config-related stuff.
"""

include: "../shared/vendored/snakemake/config.smk"

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

def is_scalar(x):
    return any([isinstance(x, t) for t in [float, int, bool, str]])

def resolve_config_value(*rule_parts, sep="/"):
    """
    A helper function intended to be used as directly as a Snakemake Input
    function to resolve the appropriate config value.

    Given an array of config *rule_parts* (keys), we return a function with a
    single argument *wildcards* which resolves to the appropriate config value.
    For instance the config value(s) for `config['filter']['min_length']` would
    be accessed via

        # within a Snakemake rule:
        params:
            min_length = resolve_config_value('filter', 'min_length')

        # within a python function:
        min_length = resolve_config_value('filter', 'min_length')(wildcards)

    The underlying config value may be structured in two ways:

    1. A scalar, e.g. `config['filter']['min_length'] = some_scalar`. In
       this case the scalar value is simply returned, i.e. it's always the same no
       matter what the wildcards are.

    2. A dictionary, e.g. `config['filter']['min_length'] = {...}`), which
       allows the resolved value to vary according to the wildcards. We search
       the dictionary for the relevant value by finding the closest matching key
       once wildcards have been considered. For instance in a three-tiered
       wildcard pipeline such as this with subtype, segment and time we
       interpret these as ordered in specificity. For instance, if only one
       wildcard value is specified in the config (the others are '*') then
       matching subtype is more specific than segment. Given example wildcard
       values of {subtype=h5nx, segment=pb2, time=2y} then we have a search
       order of:
            - 'h5nx/pb2/2y'   ─ all 3 wildcard values specified
            - 'h5nx/pb2/*'    ┐
            - 'h5nx/*/2y'     ├ 2/3 wildcard values specified
            - '*/pb2/2y'      ┘
            - 'h5nx/*/*'      ┐
            - '*/pb2/*'       ├ 1/3 wildcard values specified
            - '*/*/2y'        ┘
            - '*/*/*'         ─ default / fall-back
       and the first key present in the config is used.

    Note that in both cases, the resolved value may be a string with wildcard
    placeholders in it which may (separately) be filled in by Snakemake or a
    call to `format` or `expand`.
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

def resolve_config_fields_path(*fields):
    """
    A helper function intended to be used as directly as a Snakemake Input
    function to resolve the appropriate config path.

    Given an array of config *fields* (keys), we return a function with a single
    argument *wildcards* which resolves to a appropriate local file path
    (string). These are accomplished via calls to helper functions
    `resolve_config_value` and `resolve_filepath` (from `../shared/vendored/snakemake/config.smk`),
    respectively; see the docstrings of those functions for more details.

    Examples:
        # within a Snakemake rule
        input:
            colors = resolve_config_fields_path('colors', 'hardcoded'),
            lat_longs = resolve_config_fields_path("lat_longs"),

        # within a python function
        colors = resolve_config_fields_path('colors', 'hardcoded')(wildcards)
    """
    assert all([isinstance(f,str) for f in fields]), \
        "Arguments to `resolve_config_fields_path` must be strings"

    def resolve(wildcards):
        raw_value = resolve_config_value(*fields)(wildcards)
        if not raw_value: # falsey -> don't resolve to a path!
            return ""
        return resolve_filepath(raw_value)(wildcards)

    return resolve


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
