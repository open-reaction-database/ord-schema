# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Functions for creating Datasets by enumerating a template with a spreadsheet.

The templating code has specific expectations for how the reaction pbtxt and
spreadsheet are defined, namely that placeholder values in the pbtxt begin and
end with a "$" (dollar sign) and that these match a unique column header in the
spreadsheet file.
"""

import os
import re
from collections.abc import Mapping
from typing import BinaryIO, Optional, Union

import pandas as pd
from google.protobuf import text_format  # pytype: disable=import-error

import ord_schema
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


def read_spreadsheet(file_name_or_buffer: Union[str, BinaryIO], suffix: Optional[str] = None) -> pd.DataFrame:
    """Reads a {csv, xls, xlsx} spreadsheet file.

    Args:
        file_name_or_buffer: Filename or buffer. Note that a buffer is only
            allowed if suffix is not None.
        suffix: Filename suffix, used to determine the data encoding.

    Returns:
        DataFrame containing the reaction spreadsheet data.
    """
    if suffix is None:
        _, suffix = os.path.splitext(file_name_or_buffer)
    if suffix in [".xls", ".xlsx"]:
        return pd.read_excel(file_name_or_buffer)
    return pd.read_csv(file_name_or_buffer)


def _is_null(value: Union[float, str]) -> bool:
    """Returns whether a value is null."""
    return pd.isnull(value) or (isinstance(value, str) and (value == "nan" or not value.strip()))


def _fill_template(string: str, substitutions: Mapping[str, ord_schema.ScalarType]) -> reaction_pb2.Reaction:
    """Performs substring substitutions according to a dictionary.

    If any pattern has a null replacement value (i.e. this is an empty cell in
    the spreadsheet), apply the following edits in the given order:

        1. If after templating, any identifier does not have a defined value,
           remove the identifier.
        2. If a compound doesn't have any identifiers, remove that compound.
        3. If a compound amount is NaN, remove that compound.
        4. If an input has no components, remove that input.

    These edits are logged for easier error checking.

    Note that these edits are focused on *inputs* and are not intended to be
    exhaustive (e.g. by covering all Compound messages anywhere in Reaction).

    Args:
        string: A string whose contents should be modified.
        substitutions: A dictionary where each (key, value) pair defines
            a substring to replace and what its replacement should be.

    Returns:
        Reaction message with substitutions filled in.

    Raises:
        ValueError: If the substituted reaction template cannot be parsed.
    """
    check_null = False
    for pattern, value in substitutions.items():
        if pd.isnull(value):
            check_null = True
        string = string.replace(pattern, repr(value).strip("'"))
    try:
        reaction = text_format.Parse(string, reaction_pb2.Reaction())
    except text_format.ParseError as error:
        raise ValueError(f"Failed to parse reaction pbtxt after templating: {error}") from error
    if check_null:
        for message in reaction.inputs.values():
            for component in message.components:
                for identifier in list(component.identifiers):
                    if _is_null(identifier.value):
                        component.identifiers.remove(identifier)
            for component in list(message.components):
                kind = component.amount.WhichOneof("kind")
                if _is_null(getattr(component.amount, kind).value) or not component.identifiers:
                    message.components.remove(component)
        for key in list(reaction.inputs.keys()):
            if not reaction.inputs[key].components:
                del reaction.inputs[key]
    return reaction


def generate_dataset(template_string: str, df: pd.DataFrame, validate: bool = True) -> dataset_pb2.Dataset:
    """Generates a Dataset by enumerating a template reaction.

    Args:
        template_string: The contents of a Reaction pbtxt where placeholder
            values to be replaced are defined between dollar signs. For example,
            a SMILES identifier value could be "$product_smiles$". PLaceholders
            may only use letters, numbers, and underscores.
        df: Pandas Dataframe where each row corresponds to one reaction and
            column names match placeholders in the template_string.
        validate: Optional Boolean controlling whether Reaction messages should
            be validated as they are defined. Defaults to True.

    Returns:
        A Dataset message.

    Raises:
        ValueError: If there is no match for a placeholder string in df.
        ValueError: If validate is True and there are validation errors when
            validating an enumerated Reaction message.

    """
    placeholders = set(re.findall(r"\$\w+\$", template_string))
    for placeholder in placeholders:
        if placeholder not in df.columns:
            # Allow "$my_placeholder$" to match "my_placeholder" in df.
            if placeholder[1:-1] not in df.columns:
                raise ValueError(f"Placeholder {placeholder} not found as a" " column in dataset spreadsheet")
            df.rename(columns={placeholder[1:-1]: placeholder}, inplace=True)

    reactions = []
    for _, substitutions in df[list(placeholders)].iterrows():
        reaction = _fill_template(template_string, substitutions)
        if validate:
            output = validations.validate_message(reaction, raise_on_error=False)
            if output.errors:
                raise ValueError(f"Enumerated Reaction is not valid: {output.errors}")
        reactions.append(reaction)

    return dataset_pb2.Dataset(reactions=reactions)
