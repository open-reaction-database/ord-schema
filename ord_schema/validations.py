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
"""Helpers validating specific Message types."""

import dataclasses
import datetime
import functools
import math
import pathlib
import re
from collections.abc import Callable, Mapping
from enum import IntEnum
from typing import Any

from dateutil import parser
from rdkit import (
    # (RDKIT_VERSION reads better than rdkit_version)
    __version__ as RDKIT_VERSION,  # noqa: N812
)

import ord_schema
from ord_schema import message_helpers, parquet
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2, reaction_pb2

logger = get_logger(__name__)

# Record timestamps repeat across a dataset's reactions, but the distinct set is
# small next to the identifier caches in message_helpers, so this stays modest.
_DATETIME_CACHE_SIZE = 10_000
# Fills the fields a partial DateTime leaves out. Fixed rather than the current date
# so parsing is a pure function of the string and safe to memoize; a leap year so
# "Feb 29" still parses, and naive to match what a full value produces.
_DATETIME_DEFAULT = datetime.datetime(2000, 1, 1)  # noqa: DTZ001


@dataclasses.dataclass
class ValidationOptions:
    """Options for message validation."""

    # Check that Dataset and Reaction IDs are well-formed.
    validate_ids: bool = False
    # Require ReactionProvenance for Reactions.
    require_provenance: bool = True
    # Allow reactions with valid reaction SMILES and nothing else.
    allow_reaction_smiles_only: bool = True


@dataclasses.dataclass
class ValidationOutput:
    """Validation output: errors and warnings."""

    errors: list[str] = dataclasses.field(default_factory=list)
    warnings: list[str] = dataclasses.field(default_factory=list)

    def extend(self, other: "ValidationOutput") -> None:
        """Appends the errors and warnings from another output to this one."""
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)


class Severity(IntEnum):
    """How serious a validation finding is.

    Ordered, so a caller can threshold on ``>= Severity.ERROR`` rather than
    enumerating members.
    """

    WARNING = 1
    ERROR = 2


class ValidationError(Exception):
    """Raised when validation is asked to fail on invalid data.

    A finding's severity is :class:`Severity`; this is only the exception that
    carries an error out to the caller, via ``raise_on_error`` or
    ``validate_datasets``.
    """


@dataclasses.dataclass
class ValidationContext:
    """Where a validator reports to, and the toggles it validates under.

    Passed to every validator rather than held in module state, so what a call
    can see and affect is visible in its signature. ``validate_message`` reads
    back one message's findings by slicing off whatever its validator appended.
    """

    options: ValidationOptions = dataclasses.field(default_factory=ValidationOptions)
    findings: list[tuple[str, Severity]] = dataclasses.field(default_factory=list)

    def error(self, message: str) -> None:
        """Records a finding that makes the message invalid."""
        self.findings.append((message, Severity.ERROR))

    def warn(self, message: str) -> None:
        """Records a finding worth surfacing that does not fail validation."""
        self.findings.append((message, Severity.WARNING))


def validate_datasets(
    datasets: Mapping[str, dataset_pb2.Dataset | parquet.DatasetView],
    write_errors: bool = False,
    options: ValidationOptions | None = None,
) -> None:
    """Runs validation for a set of datasets.

    Args:
        datasets: Dict mapping text filenames to Dataset protos.
        write_errors: If True, errors are written to disk.
        options: ValidationOptions.

    Raises:
        ValidationError: if any Dataset does not pass validation.
    """
    all_errors = []
    for filename, dataset in datasets.items():
        basename = pathlib.Path(filename).name
        errors = _validate_datasets(dataset, label=basename, options=options)
        if errors:
            all_errors.extend(f"{filename}: {error}" for error in errors)
            if write_errors:
                with pathlib.Path(f"{filename}.error").open("w") as f:
                    f.writelines(f"{error}\n" for error in errors)
    # NOTE(kearnes): We run validation for all datasets before exiting if there
    # are errors.
    if all_errors:
        error_string = "\n".join(all_errors)
        raise ValidationError(f"validation encountered errors:\n{error_string}")


def _validate_datasets(
    dataset: dataset_pb2.Dataset | parquet.DatasetView,
    label: str = "dataset",
    options: ValidationOptions | None = None,
) -> list[str]:
    """Validates Reaction messages and cross-references in a Dataset.

    ``dataset`` may be a ``dataset_pb2.Dataset`` or a
    ``parquet.DatasetView``; the view re-iterates ``.reactions``
    from disk on each access, so the two iterations below (per-Reaction +
    cross-reference) both stream.

    Args:
        dataset: Dataset message or compatible stand-in.
        label: string label for logging purposes only.
        options: ValidationOptions.

    Returns:
        List of validation error messages.
    """
    errors = []
    # Reaction-level validation.
    num_bad_reactions = 0
    for i, reaction in enumerate(dataset.reactions):
        reaction_output = validate_message(
            reaction, raise_on_error=False, options=options
        )
        if reaction_output.errors:
            num_bad_reactions += 1
        for error in reaction_output.errors:
            errors.append(error)
            logger.error(f"Validation error for {label}[{i}]: {error}")
    # Dataset-level validation of cross-references. Call ``_validate_dataset``
    # directly rather than going through ``validate_message``, which insists
    # on a proto type (``_VALIDATOR_SWITCH`` lookup, ``DESCRIPTOR`` access) and
    # would reject a non-proto stand-in like ``DatasetView``.
    context = ValidationContext(options=options or ValidationOptions())
    _validate_dataset(dataset, context)
    for text, severity in context.findings:
        if severity < Severity.ERROR:
            continue
        error = f"Dataset: {text}"
        errors.append(error)
        logger.error(f"Validation error for {label}: {error}")

    return errors


def validate_message(
    message: ord_schema.Message,
    recurse: bool = True,
    raise_on_error: bool = True,
    options: ValidationOptions | None = None,
    trace: tuple[str, ...] | None = None,
    context: ValidationContext | None = None,
) -> ValidationOutput:
    """Template function for validating custom messages in the reaction_pb2.

    Messages are not validated to check enum values, since these are enforced
    by the schema. Instead, we only check for validity of items that cannot be
    enforced in the schema (e.g., non-negativity of certain measurements,
    consistency of cross-referenced keys).

    The message may be modified in place with any unambiguous changes needed to
    ensure validity.

    Args:
        message: A message to validate.
        recurse: If True, also validate submessages, meaning fields that are
            themselves messages.
        raise_on_error: If True, raises a ValidationError exception when errors
            are encountered. If False, the user must manually check the return
            value to identify validation errors.
        options: Toggles for the checks that are not always applied; see
            ``ValidationOptions``.
        trace: Tuple containing a string "stack trace" to track the position of
            the current message relative to the recursion root.
        context: Where findings are reported, and the options they are validated
            under. Created for the root call and threaded through the recursion;
            callers do not normally pass one. Supplying both this and ``options``
            validates under ``context.options`` and ignores ``options``.

    Returns:
        Errors and warnings accumulated over the message and, when recursing,
        its submessages.

    Raises:
        ValidationError: If any fields are invalid.
    """
    if context is None:
        context = ValidationContext(options=options or ValidationOptions())
    if trace is None:
        root_desc = type(message).DESCRIPTOR
        assert root_desc is not None  # Type hint.
        trace = (root_desc.name,)
    output = ValidationOutput()
    if recurse:
        for field, value in message.ListFields():
            if field.type == field.TYPE_MESSAGE:
                _validate_message(
                    field=field,
                    value=value,
                    output=output,
                    raise_on_error=raise_on_error,
                    context=context,
                    trace=trace,
                )

    # Message-specific validation
    if not isinstance(message, tuple(_VALIDATOR_SWITCH.keys())):
        # NOTE(ccoley): I made the conscious decision to raise an error here,
        # rather than assume that the message is valid. If a message does not
        # require any message-level checks (not uncommon), it should still be
        # registered in ``_VALIDATOR_SWITCH`` (mapped to ``_skip_validation``).
        # This forces us to think about what is necessary if/when new messages
        # are added.
        raise NotImplementedError(f"Don't know how to validate {type(message)}")

    # Children have already taken theirs, so whatever the validator appends from
    # here is exactly this message's.
    start = len(context.findings)
    _VALIDATOR_SWITCH[type(message)](message, context)
    stack = ".".join(trace)
    mine = context.findings[start:]
    del context.findings[start:]
    for text, severity in mine:
        warning_text = f"{stack}: {text}"
        # Thresholded so a level added above ERROR would fail validation rather
        # than be filed as a warning.
        if severity >= Severity.ERROR:
            if raise_on_error:
                raise ValidationError(warning_text)
            output.errors.append(warning_text)
        else:
            output.warnings.append(warning_text)
    return output


def _validate_message(
    field: ord_schema.FieldDescriptor,
    value: Any,
    output: ValidationOutput,
    raise_on_error: bool,
    context: ValidationContext,
    trace: tuple[str, ...],
) -> None:
    """Validates a single message field and its children.

    Args:
        field: FieldDescriptor instance.
        value: The value of the current message field.
        output: ValidationOutput.
        raise_on_error: If True, raises a ValidationError exception when errors
            are encountered. If False, the user must manually check the return
            value to identify validation errors.
        context: Where findings are reported, and the options they are
            validated under.
        trace: Tuple containing a string "stack trace" to track the position of
            the current message relative to the recursion root.
    """
    if field.label == field.LABEL_REPEATED:
        if field.message_type.GetOptions().map_entry:
            if field.message_type.fields_by_name["value"].type == field.TYPE_MESSAGE:
                for key, submessage in value.items():
                    this_trace = (*trace, f'{field.name}["{key}"]')
                    this_output = validate_message(
                        submessage,
                        raise_on_error=raise_on_error,
                        context=context,
                        trace=this_trace,
                    )
                    output.extend(this_output)
            else:  # value is a primitive
                pass
        else:  # Just a repeated message
            for index, submessage in enumerate(value):
                this_trace = (*trace, f"{field.name}[{index}]")
                this_output = validate_message(
                    submessage,
                    raise_on_error=raise_on_error,
                    context=context,
                    trace=this_trace,
                )
                output.extend(this_output)
    else:  # no recursion needed
        this_trace = (*trace, field.name)
        this_output = validate_message(
            value, raise_on_error=raise_on_error, context=context, trace=this_trace
        )
        output.extend(this_output)


def is_empty(message: ord_schema.Message) -> bool:
    """Returns whether the given message is empty."""
    empty = type(message)().SerializeToString()
    return message.SerializeToString(deterministic=True) == empty


def _ensure_float_nonnegative(
    message: ord_schema.Message, field: str, context: ValidationContext
) -> None:
    """Warns if the given numeric field of the message is negative.

    Args:
        message: The message whose field is checked.
        field: Name of the numeric field to check.
        context: Where the finding is recorded.
    """
    if getattr(message, field) < 0:
        desc = type(message).DESCRIPTOR
        assert desc is not None  # Type hint.
        context.error(f"Field {field} of message {desc.name} must be non-negative")


def _ensure_float_range(
    message: ord_schema.Message,
    field: str,
    context: ValidationContext,
    min_value: float = -math.inf,
    max_value: float = math.inf,
) -> None:
    """Warns if a numeric field of the message is outside [min_value, max_value].

    Args:
        message: The message whose field is checked.
        field: Name of the numeric field to check.
        context: Where the finding is recorded.
        min_value: Inclusive lower bound for the field value.
        max_value: Inclusive upper bound for the field value.
    """
    if getattr(message, field) < min_value or getattr(message, field) > max_value:
        desc = type(message).DESCRIPTOR
        assert desc is not None  # Type hint.
        context.error(
            f"Field {field} of message {desc.name} must be between "
            f"{min_value} and {max_value}"
        )


def _check_value_and_units(
    message: ord_schema.UnitMessage, context: ValidationContext
) -> None:
    """Checks that value/units messages are complete."""
    if not message.HasField("value"):
        context.error(f"{type(message)} requires `value` to be set")
    if message.units == message.UNSPECIFIED:
        context.error(f"{type(message)} requires `units` to be set")


def _check_type_and_details(
    message: ord_schema.TypeDetailsMessage, context: ValidationContext
) -> None:
    """Checks that type/details messages are complete."""
    if is_empty(message):
        return
    if message.type == message.UNSPECIFIED:
        context.error(f"{type(message)} requires `type` to be set")
    if message.type == message.CUSTOM and not message.details:
        context.error(f"{type(message)} has type CUSTOM but details field is empty")


def _validate_unit(message: ord_schema.UnitMessage, context: ValidationContext) -> None:
    """Validates a value/units measurement with non-negative value and precision.

    Covers the unit message types that share this exact contract (Time, Mass,
    Volume, etc.); Temperature is validated separately because of its per-unit
    absolute-zero bounds.

    Args:
        message: A unit message to validate.
        context: Where findings are recorded.
    """
    _check_value_and_units(message, context)
    _ensure_float_nonnegative(message, "value", context)
    _ensure_float_nonnegative(message, "precision", context)


def _skip_validation(message: ord_schema.Message, context: ValidationContext) -> None:
    """No-op validator for message types that need no message-level checks.

    Registered explicitly in ``_VALIDATOR_SWITCH`` so that every message type has
    a deliberate entry; see the note in ``validate_message`` about forcing a
    decision whenever a new message type is added.

    Args:
        message: The message that requires no validation.
        context: Unused; present so every dispatched validator shares a signature.
    """
    del message, context  # No message-level checks for this type.


def reaction_has_internal_standard(message: reaction_pb2.Reaction) -> bool:
    """Whether any reaction component uses the internal standard role."""
    for reaction_input in message.inputs.values():
        for compound in reaction_input.components:
            if compound.reaction_role == reaction_pb2.ReactionRole.INTERNAL_STANDARD:
                return True
    for workup in message.workups:
        if workup.input:
            for compound in workup.input.components:
                if (
                    compound.reaction_role
                    == reaction_pb2.ReactionRole.INTERNAL_STANDARD
                ):
                    return True
    return False


def reaction_has_limiting_component(message: reaction_pb2.Reaction) -> bool:
    """Whether any reaction input compound is limiting."""
    for reaction_input in message.inputs.values():
        for compound in reaction_input.components:
            if compound.is_limiting:
                return True
    return False


def reaction_needs_internal_standard(message: reaction_pb2.Reaction) -> bool:
    """Whether any analysis uses an internal standard."""
    for outcome in message.outcomes:
        for product in outcome.products:
            for measurement in product.measurements:
                if measurement.uses_internal_standard:
                    return True
    return False


def get_referenced_reaction_ids(message: reaction_pb2.Reaction) -> set[str]:
    """Return the set of reaction IDs that are referenced in a Reaction."""
    referenced_ids = set()
    for reaction_input in message.inputs.values():
        for component in reaction_input.components:
            for preparation in component.preparations:
                if preparation.reaction_id:
                    referenced_ids.add(preparation.reaction_id)
        for crude_component in reaction_input.crude_components:
            referenced_ids.add(crude_component.reaction_id)
    return referenced_ids


def is_valid_reaction_id(reaction_id: str) -> bool:
    """Returns whether a reaction ID matches the ord-<32 hex digits> format."""
    match = re.fullmatch("^ord-[0-9a-f]{32}$", reaction_id)
    return bool(match)


def is_valid_dataset_id(dataset_id: str) -> bool:
    """Returns whether a dataset ID matches the ord_dataset-<32 hex digits> format."""
    match = re.fullmatch("^ord_dataset-[0-9a-f]{32}$", dataset_id)
    return bool(match)


def is_valid_orcid(orcid: str) -> bool:
    """Returns whether an ORCID is well-formed, including its checksum.

    The final character is an ISO 7064 MOD 11-2 check digit over the preceding
    15 digits; see https://support.orcid.org/hc/en-us/articles/360006897674.

    Args:
        orcid: ORCID string, expected as 0000-0000-0000-0000.

    Returns:
        True if ``orcid`` is well-formed and the checksum is correct.
    """
    if not re.fullmatch(r"[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]", orcid):
        return False
    digits = orcid.replace("-", "")
    total = 0
    for digit in digits[:-1]:
        total = (total + int(digit)) * 2
    expected = (12 - total % 11) % 11
    check = "X" if expected == 10 else str(expected)
    return check == digits[-1]


def is_url(value: str) -> bool:
    """Returns whether a string looks like an http(s) URL with a host."""
    return bool(re.fullmatch(r"https?://[^\s/]+(?:/.*)?", value, flags=re.IGNORECASE))


def has_atom_mapping(smiles: str) -> bool:
    """Returns whether a SMILES string contains atom-map numbers."""
    # Atom maps in SMILES are always written inside brackets as e.g. [CH3:1].
    return bool(re.search(r"\[[^][]*:[0-9]+]", smiles))


@dataclasses.dataclass
class DatasetCrossRefState:
    """Aggregated cross-reference observations for a Dataset.

    A worker validating a slice of reactions feeds each one into ``observe`` and returns
    the resulting state. The master process merges the per-slice states with ``merge``
    and then ``report`` records a finding per duplicate occurrence, per self-reference,
    and one summary finding if any referenced reaction_ids are undefined.
    This keeps the streaming path behaviorally equivalent to the in-memory path.
    """

    defined_ids: set[str] = dataclasses.field(default_factory=set)
    referenced_ids: set[str] = dataclasses.field(default_factory=set)
    duplicate_count: int = 0
    self_reference_count: int = 0

    def observe(self, reaction: reaction_pb2.Reaction) -> None:
        """Records one reaction's defined ID, referenced IDs, and self-references."""
        if reaction.reaction_id:
            if reaction.reaction_id in self.defined_ids:
                self.duplicate_count += 1
            self.defined_ids.add(reaction.reaction_id)
        referenced = get_referenced_reaction_ids(reaction)
        if any(_id == reaction.reaction_id for _id in referenced):
            self.self_reference_count += 1
        self.referenced_ids |= referenced

    def merge(self, other: "DatasetCrossRefState") -> None:
        """Merges another state into this one, counting cross-slice duplicate IDs."""
        self.duplicate_count += other.duplicate_count + len(
            self.defined_ids & other.defined_ids
        )
        self.defined_ids |= other.defined_ids
        self.referenced_ids |= other.referenced_ids
        self.self_reference_count += other.self_reference_count

    def report(self, context: "ValidationContext") -> None:
        """Reports duplicate IDs, self-references, and undefined references.

        Args:
            context: Where the findings are recorded.
        """
        for _ in range(self.duplicate_count):
            context.error("Multiple Reactions should never have the same IDs")
        for _ in range(self.self_reference_count):
            context.error("A Reaction should not reference its own ID")
        undefined = self.referenced_ids - self.defined_ids
        if undefined:
            context.error(
                f"Reactions in the Dataset refer to undefined reaction_ids {undefined}"
            )


def _validate_dataset_scalars(
    *,
    context: ValidationContext,
    name: str,
    description: str,
    dataset_id: str,
    reaction_ids: list[str],
    has_reactions: bool,
) -> None:
    """Issues the dataset-level scalar/structural warnings.

    Pure function over the scalar fields so the streaming path can reuse it without re-
    iterating reactions.
    """
    if not name:
        context.error("Dataset name is required")
    if not description:
        context.error("Dataset description is required")
    if not has_reactions and not reaction_ids:
        context.error("Dataset requires reactions or reaction_ids")
    elif has_reactions and reaction_ids:
        context.error("Dataset requires reactions or reaction_ids, not both")
    if reaction_ids:
        for reaction_id in reaction_ids:
            if not is_valid_reaction_id(reaction_id):
                context.error("Reaction ID is malformed")
    if context.options.validate_ids:
        # The dataset_id is a 32-character uuid4 hex string.
        if not is_valid_dataset_id(dataset_id):
            context.error("Dataset ID is malformed")


def _validate_dataset(
    message: dataset_pb2.Dataset | parquet.DatasetView, context: ValidationContext
) -> None:
    """Validates a Dataset's scalar fields, reactions, and cross-references."""
    _validate_dataset_scalars(
        context=context,
        name=message.name,
        description=message.description,
        dataset_id=message.dataset_id,
        reaction_ids=list(message.reaction_ids),
        has_reactions=bool(message.reactions),
    )
    state = DatasetCrossRefState()
    for reaction in message.reactions:
        state.observe(reaction)
    state.report(context)


def validate_dataset_streaming(
    *,
    context: ValidationContext,
    name: str,
    description: str,
    dataset_id: str,
    reaction_ids: list[str],
    has_reactions: bool,
    state: DatasetCrossRefState,
) -> None:
    """Dataset-level validation for callers that have already streamed reactions.

    Equivalent to ``_validate_dataset`` for a Dataset whose reactions have been iterated
    in slices (e.g., per Parquet row group) by upstream workers, with each worker
    contributing a ``DatasetCrossRefState`` that the caller has merged.
    ``has_reactions`` should reflect the source's row count (e.g.,
    ``len(parquet.DatasetView(...).reactions)`` for parquet); inferring it from
    ``state`` would misclassify reactions without reaction_ids or references as empty.
    Pass ``reaction_ids=[]`` for the typical streaming case (parquet does not persist
    Dataset.reaction_ids).

    Validates under ``context.options``; there is no separate options argument, so
    this and ``_validate_dataset`` cannot disagree about which toggles are in force.
    """
    _validate_dataset_scalars(
        context=context,
        name=name,
        description=description,
        dataset_id=dataset_id,
        reaction_ids=reaction_ids,
        has_reactions=has_reactions,
    )
    state.report(context)


def _validate_dataset_example(
    message: dataset_pb2.DatasetExample, context: ValidationContext
) -> None:
    """Validates that a DatasetExample has description, url, and created set."""
    if not message.description:
        context.error("DatasetExample.description is required")
    if not message.url:
        context.error("DatasetExample.url is required")
    if not message.HasField("created"):
        context.error("DatasetExample.created is required")


def _validate_reaction(
    message: reaction_pb2.Reaction, context: ValidationContext
) -> None:
    """Validates a Reaction's inputs, outcomes, identifiers, and provenance."""
    options = context.options
    # A reaction-SMILES-only record is allowed to omit inputs and outcomes.
    smiles_only = (
        options.allow_reaction_smiles_only
        and message_helpers.get_reaction_smiles(message)
        and not message.inputs
        and not message.outcomes
    )
    if not smiles_only:
        if not message.inputs:
            context.error("Reactions should have at least 1 reaction input")
        if not message.outcomes:
            context.error("Reactions should have at least 1 reaction outcome")
    for input_ in message.inputs:
        for component in message.inputs[input_].components:
            if not component.amount.WhichOneof("kind"):
                context.error("All reaction input components require an amount")
    if reaction_needs_internal_standard(message) and not reaction_has_internal_standard(
        message
    ):
        context.error(
            "Reaction analysis uses an internal standard, but no "
            "component (as reaction input or workup) uses the "
            "reaction role INTERNAL_STANDARD"
        )
    if any(
        outcome.HasField("conversion") for outcome in message.outcomes
    ) and not reaction_has_limiting_component(message):
        context.error(
            "If reaction conversion is specified, at least one reaction input "
            "component must be labeled is_limiting"
        )
    if options.validate_ids:
        if not is_valid_reaction_id(message.reaction_id):
            context.error("Reaction ID is malformed")
    if options.require_provenance:
        if not message.HasField("provenance"):
            context.error("Reaction requires provenance")


# Block formats whose leading and trailing newlines are part of the format.
_WHITESPACE_EXEMPT_TYPES = frozenset(
    {
        reaction_pb2.CompoundIdentifier.MOLBLOCK,
        reaction_pb2.CompoundIdentifier.XYZ,
        reaction_pb2.ReactionIdentifier.RDFILE,
    }
)


def _check_surrounding_whitespace(
    message: reaction_pb2.CompoundIdentifier | reaction_pb2.ReactionIdentifier,
    context: ValidationContext,
) -> None:
    """Warns when an identifier value is padded, which defeats exact-match lookups."""
    if message.type in _WHITESPACE_EXEMPT_TYPES:
        return
    if message.value and message.value != message.value.strip():
        context.warn("value has leading or trailing whitespace")


def _check_cxsmiles_type(
    value: str, cxsmiles_type: str, context: ValidationContext
) -> None:
    """Warns when a plain-SMILES identifier carries a CXSMILES extension block.

    RDKit parses the block, so one recorded under the wrong type is otherwise silent.
    The value is still validated as recorded, since a malformed block is an error.
    """
    if message_helpers.split_cxsmiles_extension(value)[1] is not None:
        context.warn(f"value carries a CXSMILES extension block; use {cxsmiles_type}")


def _validate_reaction_identifier(
    message: reaction_pb2.ReactionIdentifier, context: ValidationContext
) -> None:
    """Validates a ReactionIdentifier's SMILES and atom-mapping consistency."""
    _check_type_and_details(message, context)
    _check_surrounding_whitespace(message, context)
    if message.type in [message.REACTION_SMILES, message.REACTION_CXSMILES]:
        # Parsed as recorded under either type: RDKit reads CXSMILES, and a malformed
        # extension block is an error rather than something to strip and ignore.
        if message.type == message.REACTION_SMILES:
            _check_cxsmiles_type(message.value, "REACTION_CXSMILES", context)
        try:
            message_helpers.validate_reaction_smiles(message.value)
        except ValueError as error:
            context.error(str(error))
        # Atom maps live in the SMILES half, never in the extension block.
        bare, _ = message_helpers.split_cxsmiles_extension(message.value)
        has_mapping = has_atom_mapping(bare)
        if message.is_mapped and not has_mapping:
            context.warn(
                "ReactionIdentifier is marked is_mapped but the SMILES "
                "contains no atom maps"
            )
        elif has_mapping and not message.is_mapped:
            context.warn(
                "ReactionIdentifier SMILES contains atom maps but is_mapped "
                "is not set to True"
            )
    if not message.value:
        context.error("value must be set")


class _StateOfMatter(IntEnum):
    GAS = 1
    LIQUID = 2
    SOLID = 3


# Coarse state of matter implied by each Texture type, for input/component
# consistency checks; UNSPECIFIED/CUSTOM map to None (unknown).
_TEXTURE_TO_STATE = {
    reaction_pb2.Texture.UNSPECIFIED: None,
    reaction_pb2.Texture.CUSTOM: None,
    reaction_pb2.Texture.GAS: _StateOfMatter.GAS,
    reaction_pb2.Texture.OIL: _StateOfMatter.LIQUID,
    reaction_pb2.Texture.FOAM: _StateOfMatter.LIQUID,
    reaction_pb2.Texture.LIQUID: _StateOfMatter.LIQUID,
    reaction_pb2.Texture.POWDER: _StateOfMatter.SOLID,
    reaction_pb2.Texture.CRYSTAL: _StateOfMatter.SOLID,
    reaction_pb2.Texture.WAX: _StateOfMatter.SOLID,
    reaction_pb2.Texture.AMORPHOUS_SOLID: _StateOfMatter.SOLID,
    reaction_pb2.Texture.SEMI_SOLID: _StateOfMatter.SOLID,
    reaction_pb2.Texture.SOLID: _StateOfMatter.SOLID,
}


def _validate_reaction_input(
    message: reaction_pb2.ReactionInput, context: ValidationContext
) -> None:
    """Validates ReactionInput component counts and texture consistency."""
    if len(message.components) + len(message.crude_components) == 0:
        context.error("Reaction inputs must have at least one component")
    elif len(message.components) + len(message.crude_components) == 1:
        for component in message.components:
            if (
                component.amount.WhichOneof("kind") == "unmeasured"
                and component.amount.unmeasured.type
                == reaction_pb2.UnmeasuredAmount.SATURATED
            ):
                context.warn(
                    "SATURATED compound amounts should only be used "
                    "for solutes when another component (solvent) is present"
                )

    input_state_code = _TEXTURE_TO_STATE[message.texture.type]
    if input_state_code is not None:
        components = [*message.components, *message.crude_components]
        component_state_codes = [_TEXTURE_TO_STATE[c.texture.type] for c in components]
        if (
            component_state_codes
            and None not in component_state_codes
            and max(component_state_codes) < input_state_code
        ):
            context.warn(
                f"the ReactionInput has texture type of: {message.texture.type}, "
                f"but its components are: {[c.texture.type for c in components]}; "
                "this seems unlikely"
            )


def _validate_amount(message: reaction_pb2.Amount, context: ValidationContext) -> None:
    """Validates that volume_includes_solutes is only set for volume amounts."""
    if (
        message.HasField("volume_includes_solutes")
        and message.WhichOneof("kind") != "volume"
    ):
        context.error("volume_includes_solutes should only be set for volume amounts")


def _validate_crude_component(
    message: reaction_pb2.CrudeComponent, context: ValidationContext
) -> None:
    """Validates that a CrudeComponent has a reaction_id and a consistent amount."""
    if not message.reaction_id:
        context.error("CrudeComponents must specify a reaction_id")
    if message.has_derived_amount and message.amount.HasField("kind"):
        context.error(
            "CrudeComponents with derived amounts cannot have their mass or "
            "volume specified explicitly"
        )
    if (
        not message.HasField("has_derived_amount") or not message.has_derived_amount
    ) and not message.amount.HasField("kind"):
        context.error(
            "Crude components should either have a derived amount or a "
            "specified mass or volume"
        )
    if message.amount.WhichOneof("kind") not in [None, "mass", "volume"]:
        context.error("Crude component amounts must be specified by mass or volume")
    if message.amount.HasField("volume_includes_solutes"):
        context.error("volume_includes_solutes should only be used for input Compounds")


def _validate_compound_identifiers(
    message: reaction_pb2.Compound | reaction_pb2.ProductCompound,
    context: ValidationContext,
) -> None:
    """Warns if a compound lacks identifiers, has only NAME, or inconsistent ones."""
    if len(message.identifiers) == 0:
        context.error("Compounds must have at least one identifier")
    if all(identifier.type == identifier.NAME for identifier in message.identifiers):
        context.warn(
            "Compounds should have more specific identifiers than NAME "
            "whenever possible"
        )
    try:
        message_helpers.check_compound_identifiers(message)
    except ValueError as error:
        context.warn(str(error))


def _validate_compound(
    message: reaction_pb2.Compound, context: ValidationContext
) -> None:
    """Validates that a Compound has usable identifiers."""
    _validate_compound_identifiers(message, context)


def _validate_compound_preparation(
    message: reaction_pb2.CompoundPreparation, context: ValidationContext
) -> None:
    """Validates CompoundPreparation type/details and reaction_id usage."""
    _check_type_and_details(message, context)
    if message.reaction_id and message.type != message.SYNTHESIZED:
        context.error(
            "Reaction IDs should only be specified in compound preparations "
            "when SYNTHESIZED"
        )


def _validate_compound_identifier(
    message: reaction_pb2.CompoundIdentifier, context: ValidationContext
) -> None:
    """Validates a CompoundIdentifier's value and type-specific format."""
    _check_type_and_details(message, context)
    _check_surrounding_whitespace(message, context)
    if not message.value:
        context.error("value must be set")
    if message.type in (
        message.SMILES,
        message.CXSMILES,
        message.INCHI,
        message.MOLBLOCK,
    ):
        # MolFromSmiles reads CXSMILES, so both types validate as recorded.
        identifier_type = {
            message.SMILES: "SMILES",
            message.CXSMILES: "CXSMILES",
            message.INCHI: "InChI",
            message.MOLBLOCK: "MolBlock",
        }[message.type]
        if message.type == message.SMILES:
            _check_cxsmiles_type(message.value, "CXSMILES", context)
        # Memoized, and shared with the consistency check in
        # check_compound_identifiers, which otherwise parses the same string again.
        if (
            message_helpers.canonical_smiles_for_identifier(message.type, message.value)
            is None
        ):
            if not message_helpers.identifier_parses_unsanitized(
                message.type, message.value
            ):
                context.error(
                    f"RDKit {RDKIT_VERSION} could not validate {identifier_type} "
                    f"identifier {message.value}"
                )
            else:
                context.warn(
                    f"RDKit {RDKIT_VERSION} could not sanitize {identifier_type} "
                    f"identifier {message.value}"
                )
    elif message.type == message.CAS_NUMBER:
        # CAS numbers are 2-7 digits, 2 digits, and a single check digit,
        # separated by hyphens (e.g., 64-17-5 for ethanol).
        if not re.fullmatch(r"\d{2,7}-\d{2}-\d", message.value):
            context.warn(
                f"CAS number {message.value} is malformed (expected e.g. 64-17-5)"
            )
    elif message.type == message.INCHI_KEY:
        if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", message.value):
            context.warn(f"InChIKey {message.value} is malformed")
    elif message.type in (message.PUBCHEM_CID, message.CHEMSPIDER_ID):
        identifier_type = (
            "PubChem CID" if message.type == message.PUBCHEM_CID else "ChemSpider ID"
        )
        if not message.value.isdecimal():
            context.warn(f"{identifier_type} {message.value} should be an integer")


def _validate_reaction_conditions(
    message: reaction_pb2.ReactionConditions, context: ValidationContext
) -> None:
    """Validates ReactionConditions dynamic-details pairing and pH range."""
    if message.conditions_are_dynamic and not message.details:
        context.error(
            "Reaction conditions are dynamic, but no details"
            " provided to explain how procedure deviates from"
            " normal single-step reaction conditions."
        )
    if message.details and not message.conditions_are_dynamic:
        context.warn(
            "Reaction condition details provided but field "
            "conditions_are_dynamic is False. If the conditions "
            "cannot be fully captured by the schema, set to True."
        )
    if message.HasField("ph") and not 0 <= message.ph <= 14:
        context.warn(f"Reaction pH ({message.ph}) is outside the expected range (0-14)")


def _validate_stirring_rate(
    message: reaction_pb2.StirringConditions.StirringRate,
    context: ValidationContext,
) -> None:
    """Validates that the stirring rate (rpm) is non-negative."""
    _ensure_float_nonnegative(message, "rpm", context)


def _validate_illumination_conditions(
    message: reaction_pb2.IlluminationConditions,
    context: ValidationContext,
) -> None:
    """Validates IlluminationConditions type and peak_wavelength usage."""
    _check_type_and_details(message, context)
    if message.type in (
        reaction_pb2.IlluminationConditions.DARK,
        reaction_pb2.IlluminationConditions.AMBIENT,
    ) and message.HasField("peak_wavelength"):
        context.warn(
            "peak_wavelength should not be specified for DARK or AMBIENT illumination"
        )


def _validate_electrochemistry_conditions(
    message: reaction_pb2.ElectrochemistryConditions,
    context: ValidationContext,
) -> None:
    """Validates ElectrochemistryConditions type/field consistency."""
    _check_type_and_details(message, context)
    if (
        message.type == reaction_pb2.ElectrochemistryConditions.CONSTANT_CURRENT
        and not message.HasField("current")
    ):
        context.warn(
            "CONSTANT_CURRENT electrochemistry conditions should specify the current"
        )
    if (
        message.type == reaction_pb2.ElectrochemistryConditions.CONSTANT_VOLTAGE
        and not message.HasField("voltage")
    ):
        context.warn(
            "CONSTANT_VOLTAGE electrochemistry conditions should specify the voltage"
        )


def _validate_reaction_workup(
    message: reaction_pb2.ReactionWorkup, context: ValidationContext
) -> None:
    """Validates a ReactionWorkup's type-specific required fields and pH range."""
    _check_type_and_details(message, context)
    if message.type == reaction_pb2.ReactionWorkup.WAIT and not message.duration.value:
        context.warn("WAIT workup steps should have a defined duration")
    if message.type == reaction_pb2.ReactionWorkup.TEMPERATURE and not message.HasField(
        "temperature"
    ):
        context.warn(
            "TEMPERATURE workup steps should have defined temperature conditions"
        )
    if (
        message.type
        in (
            reaction_pb2.ReactionWorkup.EXTRACTION,
            reaction_pb2.ReactionWorkup.FILTRATION,
        )
        and not message.keep_phase
    ):
        context.warn(
            "Workup step EXTRACTION or FILTRATION missing a recommended field "
            "keep_phase"
        )
    if (
        message.type
        in (
            reaction_pb2.ReactionWorkup.ADDITION,
            reaction_pb2.ReactionWorkup.WASH,
            reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL,
            reaction_pb2.ReactionWorkup.SCAVENGING,
            reaction_pb2.ReactionWorkup.DISSOLUTION,
            reaction_pb2.ReactionWorkup.PH_ADJUST,
        )
        and not message.input.components
    ):
        context.warn("Workup step missing recommended inputs definition")
    if message.type == reaction_pb2.ReactionWorkup.STIRRING and not message.HasField(
        "stirring"
    ):
        context.warn("Stirring workup step missing stirring definition")
    if message.type == reaction_pb2.ReactionWorkup.PH_ADJUST and not message.HasField(
        "target_ph"
    ):
        context.warn("pH adjustment workup missing target pH")
    if message.HasField("target_ph") and not 0 <= message.target_ph <= 14:
        context.warn(
            f"Workup target pH ({message.target_ph}) is outside the expected "
            "range (0-14)"
        )
    if message.type == reaction_pb2.ReactionWorkup.ALIQUOT:
        if message.amount.WhichOneof("kind") is None:
            context.warn("Aliquot workup step missing volume/mass amount")
        elif message.amount.WhichOneof("kind") not in ["mass", "volume"]:
            context.warn("Aliquot amounts should be specified by mass or volume")
        if message.amount.HasField("volume_includes_solutes"):
            context.warn(
                "volume_includes_solutes should only be used for input Compounds"
            )
    # Question: Are there other reaction workup types with specifiable amounts?
    if message.amount.WhichOneof("kind") is not None and message.type not in (
        reaction_pb2.ReactionWorkup.ALIQUOT,
        reaction_pb2.ReactionWorkup.CUSTOM,
    ):
        context.warn(
            "Workup amount should only be specified if workup type is ALIQUOT or CUSTOM"
        )


def _validate_reaction_outcome(
    message: reaction_pb2.ReactionOutcome, context: ValidationContext
) -> None:
    """Validates ReactionOutcome products, analysis keys, and conversion."""
    # *Usually* there should be at most one PRODUCT & is_desired_product
    ndp = sum(
        product.is_desired_product
        for product in message.products
        if product.reaction_role == reaction_pb2.ReactionRole.ReactionRoleType.PRODUCT
    )
    if ndp > 1:
        context.warn(
            f"Usually at most one (reaction_role == PRODUCT & is_desired_product) "
            f"product, but we have: {ndp}"
        )

    # Check key values for product analyses
    # NOTE(ccoley): Could use any(), but using expanded loops for clarity
    analysis_keys = list(message.analyses.keys())
    for product in message.products:
        for measurement in product.measurements:
            if (
                measurement.analysis_key
                and measurement.analysis_key not in analysis_keys
            ):
                context.error(
                    f"analysis key {measurement.analysis_key} does not match "
                    f"any known analysis ({analysis_keys})"
                )
    # TODO(ccoley): While we do not currently check whether the parent Reaction
    # is *actually* used in a multistep reaction within a Dataset (i.e., in a
    # CrudeComponent); this is an additional check that could be added to the
    # submission pipeline.
    if not message.products and not message.HasField("conversion"):
        context.warn(
            "No products or conversion are specified for reaction; this is "
            "permissible only for multistep reactions"
        )


def _validate_product_compound(
    message: reaction_pb2.ProductCompound, context: ValidationContext
) -> None:
    """Validates a ProductCompound's identifiers and desired-product role."""
    _validate_compound_identifiers(message, context)
    if message.is_desired_product:
        if (
            message.reaction_role
            == reaction_pb2.ReactionRole.ReactionRoleType.SIDE_PRODUCT
        ):
            context.error("a product cannot be (SIDE_PRODUCT & is_desired_product)")


def _validate_product_measurement(
    message: reaction_pb2.ProductMeasurement, context: ValidationContext
) -> None:
    """Validates a ProductMeasurement's type-specific value fields."""
    _check_type_and_details(message, context)
    if not message.analysis_key:
        context.warn(
            "Product measurements should be associated with an analysis through "
            "its analysis_key"
        )
    if message.type == reaction_pb2.ProductMeasurement.IDENTITY:
        if message.WhichOneof("value"):
            context.error(
                "Product measurements to confirm IDENTITY should not have any "
                "values defined"
            )
    elif message.type == reaction_pb2.ProductMeasurement.YIELD:
        if message.WhichOneof("value") != "percentage":
            context.warn(
                "YIELD measurements should be defined as percentage values if possible"
            )
    elif message.type == reaction_pb2.ProductMeasurement.PURITY:
        if message.WhichOneof("value") != "percentage":
            context.warn(
                "PURITY measurements should be defined as percentage values if possible"
            )
    elif message.type in (
        reaction_pb2.ProductMeasurement.AREA,
        reaction_pb2.ProductMeasurement.COUNTS,
        reaction_pb2.ProductMeasurement.INTENSITY,
    ):
        if message.WhichOneof("value") not in ("percentage", "float_value"):
            context.error(
                "Product measurements of type AREA, COUNTS, or "
                "INTENSITY must use numeric values (percentage or float_value)"
            )
    if message.HasField("selectivity") and (
        message.type != reaction_pb2.ProductMeasurement.SELECTIVITY
    ):
        context.error(
            "The selectivity_type field should only be used for a product "
            "measurement with type SELECTIVITY"
        )


def _validate_mass_spec_measurement_type(
    message: reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails,
    context: ValidationContext,
) -> None:
    """Validates mass spec m/z ranges and EIC/TIC mass usage."""
    _check_type_and_details(message, context)
    _ensure_float_nonnegative(message, "tic_minimum_mz", context)
    _ensure_float_nonnegative(message, "tic_maximum_mz", context)
    if (
        message.HasField("tic_minimum_mz")
        and message.HasField("tic_maximum_mz")
        and message.tic_minimum_mz > message.tic_maximum_mz
    ):
        context.error(
            f"tic_minimum_mz ({message.tic_minimum_mz}) must not exceed "
            f"tic_maximum_mz ({message.tic_maximum_mz})"
        )
    if any(mass < 0 for mass in message.eic_masses):
        context.error("eic_masses must be non-negative")
    if (
        message.type == reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails.EIC
        and not message.eic_masses
    ):
        context.warn("EIC mass spec measurements should specify eic_masses")
    if (
        message.type
        in (
            reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails.TIC,
            reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails.TIC_POSITIVE,
            reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails.TIC_NEGATIVE,
        )
        and message.eic_masses
    ):
        context.warn(
            "eic_masses should only be specified for EIC mass spec measurements"
        )


@functools.lru_cache(maxsize=_DATETIME_CACHE_SIZE)
def _parse_date_time(value: str) -> datetime.datetime | None:
    """Parses a DateTime value, returning None if it is unparseable.

    Memoized: a DateTime string is parsed once by this module's per-message walk
    and again by :func:`_validate_reaction_provenance` to order the timestamps,
    and record timestamps repeat across the reactions of a dataset. Caching is
    safe here because ``datetime`` is immutable.

    A partial value names only some fields ("10:30" carries no date); the rest
    come from ``_DATETIME_DEFAULT`` rather than today, which is what makes this
    safe to memoize. The filled-in date is arbitrary, and provenance compares
    timestamps only against each other.

    Args:
        value: The DateTime value.

    Returns:
        The parsed datetime, or None if ``value`` could not be parsed.
    """
    try:
        return parser.parse(value, default=_DATETIME_DEFAULT)
    except parser.ParserError:
        return None


def _parse_date_time_or_raise(value: str) -> datetime.datetime:
    """Parses a DateTime value through the cache, raising as ``parser.parse`` does.

    Lets callers that short-circuit a sequence of parses on the first bad value
    keep doing so while still going through the cache.

    Args:
        value: The DateTime value.

    Returns:
        The parsed datetime.

    Raises:
        parser.ParserError: If ``value`` could not be parsed.
    """
    parsed = _parse_date_time(value)
    if parsed is None:
        raise parser.ParserError("Unknown string format: %s", value)
    return parsed


def _validate_date_time(
    message: reaction_pb2.DateTime, context: ValidationContext
) -> None:
    """Validates that a DateTime value is parseable."""
    if message.value and _parse_date_time(message.value) is None:
        context.error(f"Could not parse DateTime string {message.value}")


def _validate_analysis(
    message: reaction_pb2.Analysis, context: ValidationContext
) -> None:
    """Validates an Analysis message's type and details."""
    # TODO(ccoley): Will be lots to expand here if we add structured data.
    _check_type_and_details(message, context)


def _validate_reaction_provenance(
    message: reaction_pb2.ReactionProvenance, context: ValidationContext
) -> None:
    """Validates ReactionProvenance timestamps, emails, DOI, and URL."""
    # Prepare datetimes
    if not message.HasField("record_created"):
        context.error("Reactions must have record_created defined")
    experiment_start = None
    record_created = None
    record_modified = None
    try:
        if message.experiment_start.value:
            experiment_start = _parse_date_time_or_raise(message.experiment_start.value)
        if message.record_created.time.value:
            record_created = _parse_date_time_or_raise(
                message.record_created.time.value
            )
        for record in message.record_modified:
            # Use the last record as the most recent modification time.
            record_modified = _parse_date_time_or_raise(record.time.value)
    except parser.ParserError:
        context.warn("Failed to parse DateTime string(s)")
    # Check signs of time differences
    if experiment_start and record_created:
        if (record_created - experiment_start).total_seconds() < 0:
            context.error("Record creation time should be after experiment")
    if record_modified and record_created:
        if (record_modified - record_created).total_seconds() < 0:
            context.error("Record modified time should be after creation")
    if not message.record_created.person.email:
        context.error("User email is required for record_created")
    for record in message.record_modified:
        if not record.person.email:
            context.error("User email is required for record_modified")
    if message.doi:
        try:
            parsed_doi = message_helpers.parse_doi(message.doi)
            if message.doi != parsed_doi:
                context.error(f"DOI should be trimmed ({message.doi} -> {parsed_doi})")
        except ValueError as error:
            context.error(str(error))
    if message.publication_url and not is_url(message.publication_url):
        context.warn(
            f"publication_url does not look like a valid URL: {message.publication_url}"
        )


def _validate_record_event(
    message: reaction_pb2.RecordEvent, context: ValidationContext
) -> None:
    """Validates that a RecordEvent has a time and an identifiable person."""
    if not message.time.value:
        context.error("RecordEvent must have `time` specified")
    person = message.person
    if not (person.username or person.name or person.orcid):
        context.error(
            "Person must have at least one of `username`, `name`, or `orcid` specified"
        )
    if not person.email:
        context.error("Person must have `email` specified")


def _validate_person(message: reaction_pb2.Person, context: ValidationContext) -> None:
    """Validates a Person's ORCID and email formats."""
    if message.orcid:
        if not is_valid_orcid(message.orcid):
            context.error("Invalid ORCID: Enter as 0000-0000-0000-0000")
    if message.email:
        # Based on https://www.regular-expressions.info/email.html.
        # Added optional "[bot]" suffix to the username for GitHub actions.
        if not re.fullmatch(
            r"[a-zA-Z0-9._+-]+(?:\[bot\])?@[a-zA-Z0-9.-]+\.[a-z]{2,}", message.email
        ):
            context.error(f"Invalid email address: {message.email}")


def _validate_temperature(
    message: reaction_pb2.Temperature, context: ValidationContext
) -> None:
    """Validates a Temperature, enforcing absolute-zero lower bounds per unit."""
    _check_value_and_units(message, context)
    if message.units == message.CELSIUS:
        _ensure_float_range(message, "value", context, min_value=-273.15)
    elif message.units == message.FAHRENHEIT:
        _ensure_float_range(message, "value", context, min_value=-459)
    elif message.units == message.KELVIN:
        _ensure_float_range(message, "value", context, min_value=0)
    _ensure_float_nonnegative(message, "precision", context)


def _validate_percentage(
    message: reaction_pb2.Percentage, context: ValidationContext
) -> None:
    """Validates that a Percentage value is within 0-100 (and not a fraction)."""
    if not message.HasField("value"):
        context.error(f"{type(message)} requires `value` to be set")
    if 0 < message.value < 1:
        context.warn(
            f"Percentage values are 0-100, not fractions ({message.value} used)"
        )
    if message.value < 0 or message.value > 100:
        context.warn(
            f"Percentage value ({message.value}) is outside the expected range (0-100)"
        )
    _ensure_float_nonnegative(message, "precision", context)


def _validate_float_value(
    message: reaction_pb2.FloatValue, context: ValidationContext
) -> None:
    """Validates that a FloatValue's precision is non-negative."""
    _ensure_float_nonnegative(message, "precision", context)


def _validate_data(message: reaction_pb2.Data, context: ValidationContext) -> None:
    """Validates that a Data message has a value and a valid URL/format."""
    if not message.WhichOneof("kind"):
        context.error("Data requires one of {value, bytes_value, url}")
    if message.bytes_value and not message.format:
        context.error("Data format is required for bytes_data")
    if message.WhichOneof("kind") == "url" and not is_url(message.url):
        context.warn(f"Data URL does not look like a valid URL: {message.url}")


_VALIDATOR_SWITCH: dict[type, Callable[..., None]] = {
    dataset_pb2.Dataset: _validate_dataset,
    dataset_pb2.DatasetExample: _validate_dataset_example,
    reaction_pb2.Reaction: _validate_reaction,
    # Basics
    reaction_pb2.ReactionIdentifier: _validate_reaction_identifier,
    reaction_pb2.ReactionInput: _validate_reaction_input,
    reaction_pb2.ReactionInput.AdditionDevice: _check_type_and_details,
    reaction_pb2.ReactionInput.AdditionSpeed: _skip_validation,
    # Compounds
    reaction_pb2.Amount: _validate_amount,
    reaction_pb2.UnmeasuredAmount: _check_type_and_details,
    reaction_pb2.Texture: _check_type_and_details,
    reaction_pb2.CrudeComponent: _validate_crude_component,
    reaction_pb2.Compound: _validate_compound,
    reaction_pb2.CompoundPreparation: _validate_compound_preparation,
    reaction_pb2.CompoundIdentifier: _validate_compound_identifier,
    reaction_pb2.Compound.Source: _skip_validation,
    # Setup
    reaction_pb2.Vessel: _check_type_and_details,
    reaction_pb2.VesselMaterial: _check_type_and_details,
    reaction_pb2.VesselAttachment: _check_type_and_details,
    reaction_pb2.VesselPreparation: _check_type_and_details,
    reaction_pb2.ReactionSetup: _skip_validation,
    reaction_pb2.ReactionSetup.ReactionEnvironment: _check_type_and_details,
    # Conditions
    reaction_pb2.ReactionConditions: _validate_reaction_conditions,
    reaction_pb2.TemperatureConditions: _skip_validation,
    reaction_pb2.TemperatureConditions.TemperatureControl: _check_type_and_details,
    reaction_pb2.TemperatureConditions.TemperatureMeasurement: _check_type_and_details,
    reaction_pb2.PressureConditions: _skip_validation,
    reaction_pb2.PressureConditions.PressureControl: _check_type_and_details,
    reaction_pb2.PressureConditions.Atmosphere: _check_type_and_details,
    reaction_pb2.PressureConditions.PressureMeasurement: _check_type_and_details,
    reaction_pb2.StirringConditions: _check_type_and_details,
    reaction_pb2.StirringConditions.StirringRate: _validate_stirring_rate,
    reaction_pb2.IlluminationConditions: _validate_illumination_conditions,
    reaction_pb2.ElectrochemistryConditions: _validate_electrochemistry_conditions,
    reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell: _check_type_and_details,  # noqa: E501
    reaction_pb2.ElectrochemistryConditions.ElectrochemistryMeasurement: _skip_validation,  # noqa: E501
    reaction_pb2.FlowConditions: _check_type_and_details,
    reaction_pb2.FlowConditions.Tubing: _check_type_and_details,
    # Annotations
    reaction_pb2.ReactionNotes: _skip_validation,
    reaction_pb2.ReactionObservation: _skip_validation,
    # Outcome
    reaction_pb2.ReactionWorkup: _validate_reaction_workup,
    reaction_pb2.ReactionOutcome: _validate_reaction_outcome,
    reaction_pb2.ProductCompound: _validate_product_compound,
    reaction_pb2.ProductMeasurement: _validate_product_measurement,
    reaction_pb2.ProductMeasurement.Selectivity: _check_type_and_details,
    reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails: _validate_mass_spec_measurement_type,  # noqa: E501
    reaction_pb2.DateTime: _validate_date_time,
    reaction_pb2.Analysis: _validate_analysis,
    # Metadata
    reaction_pb2.ReactionProvenance: _validate_reaction_provenance,
    reaction_pb2.RecordEvent: _validate_record_event,
    reaction_pb2.Person: _validate_person,
    # Units
    reaction_pb2.Time: _validate_unit,
    reaction_pb2.Mass: _validate_unit,
    reaction_pb2.Moles: _validate_unit,
    reaction_pb2.Volume: _validate_unit,
    reaction_pb2.Concentration: _validate_unit,
    reaction_pb2.Pressure: _validate_unit,
    reaction_pb2.Temperature: _validate_temperature,
    reaction_pb2.Current: _validate_unit,
    reaction_pb2.Voltage: _validate_unit,
    reaction_pb2.Length: _validate_unit,
    reaction_pb2.Wavelength: _validate_unit,
    reaction_pb2.FlowRate: _validate_unit,
    reaction_pb2.Percentage: _validate_percentage,
    reaction_pb2.FloatValue: _validate_float_value,
    reaction_pb2.Data: _validate_data,
}
