"""
Shared functions to be used within a Snakemake workflow for parsing
workflow configs.
"""
import os.path
from collections.abc import Callable
from snakemake.io import Wildcards
from typing import Optional
from textwrap import dedent, indent


class InvalidConfigError(Exception):
    pass


def resolve_config_path(path: str) -> Callable[[Wildcards], str]:
    """
    Resolve a relative *path* given in a configuration value. Will always try to
    resolve *path* after expanding wildcards with Snakemake's `expand` functionality.

    Returns the path for the first existing file found relative to any directory
    in the AUGUR_SEARCH_PATHS environment variable (searched in order).
    """
    global workflow

    def _resolve_config_path(wildcards):
        try:
            expanded_path = expand(path, **wildcards)[0]
        except snakemake.exceptions.WildcardError as e:
            available_wildcards = "\n".join(f"  - {wildcard}" for wildcard in wildcards)
            raise snakemake.exceptions.WildcardError(indent(dedent(f"""\
                {str(e)}

                However, resolve_config_path({{path}}) requires the wildcard.

                Wildcards available for this path are:

                {{available_wildcards}}

                Hint: Check that the config path value does not misspell the wildcard name
                and that the rule actually uses the wildcard name.
                """.lstrip("\n").rstrip()).format(path=repr(path), available_wildcards=available_wildcards), " " * 4))

        search_paths = os.environ.get('AUGUR_SEARCH_PATHS', '').split(':')

        for search_path in search_paths:
            candidate_path = os.path.join(search_path, expanded_path)
            if os.path.exists(candidate_path):
                return candidate_path

        raise InvalidConfigError(indent(dedent(f"""\
            Unable to resolve the config-provided path {path!r},
            expanded to {expanded_path!r} after filling in wildcards.

            Searched in the following locations from AUGUR_SEARCH_PATHS:
            {search_paths}
            """), " " * 4))

    return _resolve_config_path
