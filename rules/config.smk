"""
Functions and logic related to finding, parsing and interpreting configfiles
and other config-related stuff.
"""
import copy
from typing import Any, Literal, TypedDict

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

DATASET_LEVELS = ['subtype', 'segment', 'time']

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

    The first *rule_part* names a top-level config section; any remaining parts
    name keys within the resolved value. The section may be structured in three
    ways:

    1. A scalar, e.g. `config['reference'] = some_scalar`. The scalar is
       returned as-is, i.e. it's always the same no matter what the wildcards
       are.

    2. A "dataset-first" mapping of patterns to config layers, e.g.
       `config['filter'] = {'*/*/*': {...}, '*/pb2/*': {...}}`. Every pattern
       matching the current dataset (built from the wildcards in the order of
       `DATASET_LEVELS`) contributes a layer; the layers are deep-merged in the
       order they appear in the config (later, more specific patterns win). The
       remaining *rule_parts* then index into the merged mapping.

    3. A mapping of patterns directly to (non-dict) values, e.g.
       `config['subtype_query'] = {'h5nx/*/*': '...', ...}`. The value from the
       last matching pattern is returned.

    Patterns are slash-delimited with one part per `DATASET_LEVELS` entry. Each
    part may be a literal value, '*' (matches anything), or a whole-part
    multivalue like '(h5nx|h5n1)'.

    Note that in all cases, the resolved value may be a string with wildcard
    placeholders in it which may (separately) be filled in by Snakemake or a
    call to `format` or `expand`.
    """
    section, *param_path = rule_parts
    try:
        section_config = config[section]
    except KeyError:
        raise InvalidConfigError(f'Config missing entire entry for config["{section}"]')

    def resolve(wildcards):
        if is_scalar(section_config):
            value = section_config
        else:
            if not isinstance(section_config, dict):
                raise InvalidConfigError(f"config[{section!r}] must be a scalar value or a dictionary")

            dataset = tuple(wildcards[k] for k in DATASET_LEVELS)
            matches = matching_pattern_values(section_config, dataset)
            if not matches:
                raise InvalidConfigError(
                    f'Config structure incorrect or incomplete for config[{section!r}]\n'
                    f'\tThe dictionary is missing a matching pattern for the current '
                    f'target of {sep.join(dataset)!r} (literal or "*" placeholders).')

            # Dataset-first config layers are deep-merged; a plain pattern→value
            # mapping resolves to the last (most specific) matching value.
            if all(isinstance(m, dict) for m in matches):
                value = merge_layers(matches)
            else:
                value = matches[-1]

        for key in param_path:
            try:
                value = value[key]
            except (KeyError, TypeError):
                raise InvalidConfigError(
                    f'Config missing entry {key!r} for config'
                    + ''.join(f'[{p!r}]' for p in rule_parts))
        return value

    return resolve

def resolve_config_fields_path(*fields):
    """
    A helper function intended to be used as directly as a Snakemake Input
    function to resolve the appropriate config path.

    Given an array of config *fields* (keys), we return a function with a single
    argument *wildcards* which resolves to a appropriate local file path
    (string). These are accomplished via calls to helper functions
    `resolve_config_value` and `resolve_config_path` (from `../shared/vendored/snakemake/config.smk`),
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
        return resolve_config_path(raw_value, AVIAN_FLU_DIR)(wildcards)

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

# FIXME: everything below is independent of VALID_DATASET_LEVELS/DATASET_LEVELS_TO_RUN and can be moved to shared/vendored/config.smk or Augur

ExactDataset = tuple[str, ...]
"""Exact dataset values, ordered to match VALID_DATASET_LEVELS."""

class DatasetPatternPart(TypedDict):
    type: Literal["wildcard", "literal", "multivalue"]
    matches: tuple[str, ...] | None


class DatasetLevel(TypedDict):
    name: str
    values: list[str]


def validate_rule_config(
    rule_name: str,
    rule_config: dict[str, Any],
    dataset_levels: list[DatasetLevel],
) -> None:
    if not isinstance(rule_config, dict):
        raise InvalidConfigError(f"'{rule_name}' must be a mapping of dataset patterns to config layers.")

    for pattern, config_layer in rule_config.items():
        if not isinstance(pattern, str):
            raise InvalidConfigError(f"{rule_name} pattern {pattern!r} must be a string.")
        _validate_dataset_pattern(pattern, dataset_levels, rule_name)

        if not isinstance(config_layer, dict):
            raise InvalidConfigError(f"{rule_name} config for pattern {pattern!r} must be a mapping.")


def _validate_dataset_pattern(
    pattern: str,
    dataset_levels: list[DatasetLevel],
    context: str,
) -> None:
    pattern_parts = parse_dataset_pattern(pattern)
    if len(pattern_parts) != len(dataset_levels):
        raise InvalidConfigError(dedent(f"""\
            Invalid {context} dataset pattern {pattern!r}.
            Expected {len(dataset_levels)} slash-separated parts matching:
                {'/'.join(level['name'] for level in dataset_levels)}"""))

    for pattern_part, level in zip(pattern_parts, dataset_levels):
        if pattern_part["type"] == "wildcard":
            continue

        invalid_values = sorted(set(pattern_part["matches"]) - set(level["values"]))
        if invalid_values:
            raise InvalidConfigError(dedent(f"""\
                Invalid {context} dataset value(s) {invalid_values!r} in pattern {pattern!r}.
                Expected {level['name']} values from: {level['values']}"""))


def get_datasets(levels: list[DatasetLevel]) -> list[ExactDataset]:
    """
    Return all datasets requested by config, in the given levels order.
    """
    return product(*(level["values"] for level in levels))


def write_command_configs(
    config: dict[str, Any],
    rule_name: str,
    datasets: list[ExactDataset],
    output_dir: str,
) -> None:
    """
    Write a per-dataset Augur command config, one file per dataset.
    """
    for dataset in datasets:
        out = get_rule_config(config, rule_name, dataset)
        path = dataset_config_path(output_dir, dataset, rule_name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            print(f"# {'/'.join(dataset)}", file=f)
            yaml.dump(out, f, sort_keys=False, Dumper=NoAliasDumper)
        print(f"Saved {rule_name} config to {path!r}.", file=sys.stderr)


def dataset_config_path(output_dir: str, dataset: ExactDataset, rule_name: str) -> str:
    """
    Path of the augur config written for a dataset and rule.

    The path is '<output_dir>/<dataset>/<rule_name>_config.yaml', where
    '<dataset>' is the dataset's slash-joined values.
    """
    return f"{output_dir}/{'/'.join(dataset)}/{rule_name}_config.yaml"


def matching_pattern_values(
    patterned_config: dict[str, Any],
    dataset: ExactDataset,
) -> list[Any]:
    """
    Return the values from a {pattern: value} mapping whose pattern matches the
    dataset, preserving insertion order.
    """
    return [
        value
        for pattern, value in patterned_config.items()
        if pattern_matches_dataset(pattern, dataset)
    ]


def get_rule_config(
    config: dict[str, Any],
    rule_name: str,
    dataset: ExactDataset,
) -> dict[str, Any]:
    """
    Build the config for a rule and dataset.

    A matching 'custom_<rule>' replaces '<rule>' entirely: if any matching layer
    defines it, the '<rule>' layers are discarded and the config is built from
    the 'custom_<rule>' layers alone (e.g. 'custom_subsample' replaces the
    default 'subsample'). Otherwise the config is built from the '<rule>' layers.

    Within either, layers are merged top-to-bottom with later values overriding
    earlier ones.
    """
    if custom_layers := get_rule_layers(config, rule_override_name(rule_name), dataset):
        return merge_layers(custom_layers)

    return merge_layers(get_rule_layers(config, rule_name, dataset))


def rule_override_name(rule_name: str) -> str:
    return f"custom_{rule_name}"


def get_rule_layers(
    config: dict[str, Any],
    rule_name: str,
    dataset: ExactDataset,
) -> list[dict[str, Any]]:
    """
    Config layers for a rule and dataset, lowest priority first.

    The top-level '<rule>' key maps dataset patterns to config layers.
    Matching '<rule>.<pattern>' configs apply in pattern order.
    """
    layers = []
    if rule_name in config:
        layers += matching_pattern_values(config[rule_name], dataset)
    return layers


def merge_layers(layers: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Deep-merge config layers in order, with later layers winning.
    """
    merged: dict[str, Any] = {}
    for layer in layers:
        deep_merge(merged, layer)
    return merged


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> None:
    """
    Recursively merge 'override' into 'base' (mutating 'base').

    Mappings are merged recursively; scalars and lists overwrite. Dicts brought
    in from 'override' are deep-copied so the source config is never mutated.
    """
    for key, value in override.items():
        # FIXME: allow null-deletion?
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_merge(base[key], value)
        elif isinstance(value, dict):
            base[key] = copy.deepcopy(value)
        else:
            base[key] = value


def pattern_matches_dataset(
    pattern: str,
    dataset: ExactDataset,
) -> bool:
    """
    Return whether a dataset pattern matches an exact dataset.
    """
    pattern_parts = parse_dataset_pattern(pattern)
    if len(pattern_parts) != len(dataset):
        return False

    for pattern_part, dataset_value in zip(pattern_parts, dataset):
        if pattern_part["type"] == "wildcard":
            continue

        if dataset_value not in pattern_part["matches"]:
            return False

    return True


def parse_dataset_pattern(pattern: str) -> tuple[DatasetPatternPart, ...]:
    """
    Parse a slash-delimited dataset pattern.
    """
    return tuple(parse_dataset_pattern_part(part) for part in pattern.split("/"))


def parse_dataset_pattern_part(part: str) -> DatasetPatternPart:
    """
    Parse one part of a dataset pattern.

    Supported syntax:
    1. A literal value : 6y
    2. Multiple values : (6y|3y)
    3. All values      : *
    """
    if part == "*":
        return {"type": "wildcard", "matches": None}

    if part.startswith("(") and part.endswith(")"):
        values = tuple(part[1:-1].split("|"))
        if not values or any(not value for value in values):
            raise InvalidConfigError(f"Invalid multivalue dataset part {part!r}.")
        return {"type": "multivalue", "matches": values}

    if any(char in part for char in "()|"):
        raise InvalidConfigError(dedent(f"""\
            Invalid subsample dataset part {part!r}.
            Use '*', a literal value, or a whole-part multivalue like '(genome|N450)'."""))

    return {"type": "literal", "matches": (part,)}
