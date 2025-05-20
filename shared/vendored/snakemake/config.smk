"""
Shared functions to be used within a Snakemake workflow for parsing
workflow configs.
"""
import os.path
from collections.abc import Callable
from snakemake.io import Wildcards
from typing import Optional, List


class InvalidConfigError(Exception):
    pass


def resolve_config_path(path: str, other_prefixes: Optional[List[str]] = None) -> Callable[[Wildcards], str]:
    """
    Resolve a relative *path* given in a configuration value.

    Resolves *path* as relative to the workflow's ``defaults/`` directory (i.e.
    ``os.path.join(workflow.basedir, "defaults", path)``) if it doesn't exist
    in the workflow's analysis directory (i.e. the current working
    directory, or workdir, usually given by ``--directory`` (``-d``)).

    This behaviour allows a default configuration value to point to a default
    auxiliary file while also letting the file used be overridden either by
    setting an alternate file path in the configuration or by creating a file
    with the conventional name in the workflow's analysis directory.

    Will always try to resolve *path* or the default path after expanding
    wildcards with Snakemake's `expand` functionality.

    If *other_prefixes* are provided, then will also try to resolve *path*
    relative to the provided prefixes if the *path* and default path do not exist.
    """
    global workflow

    def _resolve_config_path(wildcards):
        expanded_path = expand(path, **wildcards)[0]
        if os.path.exists(expanded_path):
            return expanded_path

        # Special-case defaults/… for backwards compatibility with older
        # configs.  We could achieve the same behaviour with a symlink
        # (defaults/defaults → .) but that seems less clear.
        if path.startswith("defaults/"):
            defaults_path = os.path.join(workflow.basedir, expanded_path)
        else:
            defaults_path = os.path.join(workflow.basedir, "defaults", expanded_path)

        if os.path.exists(defaults_path):
            return defaults_path

        checked_paths = [expanded_path, defaults_path]
        if other_prefixes:
            for prefix in other_prefixes:
                prefixed_path = os.path.join(prefix, expanded_path)

                if os.path.exists(prefixed_path):
                    return prefixed_path

                checked_paths.append(prefixed_path)

        raise InvalidConfigError(
            "Unable to resolve config provided path. Checked for the following files:\n" + \
            "\n".join("\t" + f"{index + 1}. {path!r}" for index, path in enumerate(checked_paths)))

    return _resolve_config_path
