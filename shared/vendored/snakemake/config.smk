"""
Shared functions to be used within a Snakemake workflow for parsing
workflow configs.
"""
import os
import sys
from collections.abc import Callable
from snakemake.io import Wildcards
from textwrap import dedent, indent

# FIXME: move to augur.io.resolve_filepath ?
from augur.subsample import _resolve_filepath


# Set search paths for Augur
if "AUGUR_SEARCH_PATHS" in os.environ:
    print(dedent(f"""\
        Using existing search paths in AUGUR_SEARCH_PATHS:

            {os.environ["AUGUR_SEARCH_PATHS"]!r}
        """), file=sys.stderr)
else:
    search_paths = [
        # User analysis directory
        Path.cwd(),

        # Workflow defaults folder
        Path(workflow.basedir) / "defaults",

        # Workflow root (contains Snakefile)
        Path(workflow.basedir),
    ]

    repo_root = Path(workflow.basedir) / ".."
    if (repo_root / "nextstrain-pathogen.yaml").is_file():
        search_paths.extend([
            # Pathogen repo root
            repo_root,
        ])

    search_paths = [path.resolve() for path in search_paths if path.is_dir()]

    os.environ["AUGUR_SEARCH_PATHS"] = ":".join(map(str, search_paths))


class InvalidConfigError(Exception):
    pass


def resolve_filepath(path: str) -> Callable[[Wildcards], str]:
    """
    Resolve a relative *path*. Will always try to
    resolve *path* after expanding wildcards with Snakemake's `expand` functionality.

    Returns the path for the first existing file found relative to any directory
    in the AUGUR_SEARCH_PATHS environment variable (searched in order).
    """
    search_paths = [Path(p) for p in os.environ["AUGUR_SEARCH_PATHS"].split(":")]

    def callable(wildcards):
        try:
            expanded_path = expand(path, **wildcards)[0]
        except snakemake.exceptions.WildcardError as e:
            available_wildcards = "\n".join(f"  - {wildcard}" for wildcard in wildcards)
            raise snakemake.exceptions.WildcardError(indent(dedent(f"""\
                {str(e)}

                However, resolve_filepath({{path}}) requires the wildcard.

                Wildcards available for this path are:

                {{available_wildcards}}

                Hint: Check that the path value does not misspell the wildcard name
                and that the rule actually uses the wildcard name.
                """.lstrip("\n").rstrip()).format(path=repr(path), available_wildcards=available_wildcards), " " * 4))

        return str(_resolve_filepath(Path(expanded_path), search_paths))

    return callable
