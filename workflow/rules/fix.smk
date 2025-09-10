# FIX for setting default values (https://github.com/snakemake/snakemake/issues/3614)
# FROM
# https://python-jsonschema.readthedocs.io/en/stable/faq/#why-doesn-t-my-schema-s-default-property-set-the-default-on-my-instance

from jsonschema import Draft202012Validator, validators
import yaml


def extend_with_default(validator_class):
  validate_properties = validator_class.VALIDATORS["properties"]

  def set_defaults(validator, properties, instance, schema):
    for property, subschema in properties.items():
      if "default" in subschema:
        instance.setdefault(property,subschema["default"])

    for error in validate_properties(
            validator,properties,instance,schema,
    ):
      yield error

  return validators.extend(
    validator_class,{"properties": set_defaults},
  )
DefaultValidatingValidator = extend_with_default(Draft202012Validator)


# validate schema
with open(workflow.basedir + "/schemas/config.yaml") as stream:
  config_schema = yaml.safe_load(stream)
  DefaultValidatingValidator(config_schema).validate(config)

# FIXME validate pep