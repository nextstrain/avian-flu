# ------------- helper functions to collect, merge & download input files ------------------- #

def _parse_config_input(input):
    """
    Parses information from an individual config-defined input, i.e. an element within `config.inputs` or `config.additional_inputs`
    and returns information snakemake rules can use to obtain the underlying data.

    The structure of `input` is a dictionary with keys:
    - name:string (required)
    - metadata:string (optional) - a s3 URI or a local file path
    - sequences:string|dict[string,string] (optional) - either a s3 URI or a local file path, in which case
      it must include a '{segment}' wildcard substring, or a dict of segment → s3 URI or local file path,
      in which case it must not include the wildcard substring.

    Returns a dictionary with optional keys:
    - metadata:string - path or url to the metadata file.
    - sequences:function. Takes in wildcards and returns path or url to the sequences FASTA for the provided
      segment wildcard, or returns `None` if this input doesn't define sequences for the provided segment.

    Raises InvalidConfigError
    """
    info = {
        "name": input["name"],
        "metadata": path_or_url(input["metadata"]) if input.get("metadata") else None,
        "sequences": None,
    }

    if location:=input.get('sequences', False):
        if isinstance(location, dict):
            info['sequences'] = lambda w: path_or_url(location[w.segment]) \
                if w.segment in location \
                else None
        elif isinstance(location, str):
            info['sequences'] = lambda w: path_or_url(location)
        else:
            raise InvalidConfigError(f"Config input for {info['name']} specifies sequences in an unknown format; must be dict or string")

    return info

def _gather_inputs():
    all_inputs = [*config['inputs'], *config.get('additional_inputs', [])]

    if len(all_inputs)==0:
        raise InvalidConfigError("Config must define at least one element in config.inputs or config.additional_inputs lists")
    if not all([isinstance(i, dict) for i in all_inputs]):
        raise InvalidConfigError("All of the elements in config.inputs and config.additional_inputs lists must be dictionaries. "
            "If you've used a command line '--config' double check your quoting.")
    if len({i['name'] for i in all_inputs})!=len(all_inputs):
        raise InvalidConfigError("Names of inputs (config.inputs and config.additional_inputs) must be unique")
    if not all(['name' in i and ('sequences' in i or 'metadata' in i) for i in all_inputs]):
        raise InvalidConfigError("Each input (config.inputs and config.additional_inputs) must have a 'name' and 'metadata' and/or 'sequences'")
    if not any(['metadata' in i for i in all_inputs]):
        raise InvalidConfigError("At least one input must have 'metadata'")
    if not any (['sequences' in i for i in all_inputs]):
        raise InvalidConfigError("At least one input must have 'sequences'")

    available_keys = set(['name', 'metadata', 'sequences'])
    if any([len(set(el.keys())-available_keys)>0 for el in all_inputs]):
        raise InvalidConfigError(f"Each input (config.inputs and config.additional_inputs) can only include keys of {', '.join(available_keys)}")

    return {i['name']: _parse_config_input(i) for i in all_inputs}

input_sources = _gather_inputs()

def input_metadata(wildcards):
    inputs = [info['metadata'] for info in input_sources.values() if info.get('metadata', None)]
    return inputs[0] if len(inputs)==1 else "results/metadata_merged.tsv"

def input_sequences(wildcards):
    inputs = list(filter(None, [info['sequences'](wildcards) for info in input_sources.values() if info.get('sequences', None)]))
    return inputs[0] if len(inputs)==1 else "results/sequences_merged_{segment}.fasta"

rule merge_metadata:
    """
    This rule should only be invoked if there are multiple defined metadata inputs
    (config.inputs + config.additional_inputs)
    """
    input:
        **{name: info['metadata'] for name,info in input_sources.items() if info.get('metadata', None)}
    params:
        metadata = lambda w, input: list(map("=".join, input.items()))
    output:
        metadata = "results/metadata_merged.tsv"
    shell:
        r"""
        augur merge \
            --metadata {params.metadata:q} \
            --source-columns 'input_{{NAME}}' \
            --output-metadata {output.metadata}
        """

rule merge_sequences:
    """
    This rule should only be invoked if there are multiple defined metadata inputs
    (config.inputs + config.additional_inputs) for this particular segment
    """
    input:
        lambda w: list(filter(None, [info['sequences'](w) for info in input_sources.values()]))
    output:
        sequences = "results/sequences_merged_{segment}.fasta"
    shell:
        r"""
        augur merge \
            --sequences {input:q} \
            --output-sequences {output.sequences:q}
        """

# -------------------------------------------------------------------------------------------- #
