
class ValidateError(Exception):
    pass

def validate(config):
    # prototype based on <https://github.com/nextstrain/augur/blob/master/augur/validate.py>
    from importlib import metadata
    # <https://github.com/nextstrain/augur/issues/1358>
    assert str(metadata.version("jsonschema")).startswith('3.'), "jsonschema must be version 3"
    import jsonschema
    import jsonschema.exceptions
    import yte
    from os import path

    with open(path.join(path.dirname(path.realpath(__file__)), "../config/schema.yaml"), encoding='utf-8') as f:
        schema = yte.process_yaml(f, require_use_yte=True)

    Validator = jsonschema.validators.validator_for(schema)

    try:
        Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise ValidateError(f"Internal error: config schema is not a valid JSON Schema ({Validator.META_SCHEMA['$schema']}). Error: {err}")

    # Here we're validating the merged schema. We could also validate the user config on its own by making (all?) properties optional?
    from augur.validate import validate_json
    validate_json(config, Validator(schema), "config")

    # There are more checks we can do by using code, as desired.
