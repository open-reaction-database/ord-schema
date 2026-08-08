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
"""Tests for ord_schema.units."""

import typing

import pytest
from google.protobuf import descriptor

import ord_schema
from ord_schema import units
from ord_schema.proto import reaction_pb2


@pytest.fixture
def resolver() -> units.UnitResolver:
    return units.UnitResolver()


@pytest.mark.parametrize(
    ("string", "expected"),
    [
        (
            "15.0 ML",
            reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER),
        ),
        ("24 H", reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR)),
        ("32.1g", reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM)),
        (
            "   32.1      \t   g  ",
            reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM),
        ),
        ("10 min.", reaction_pb2.Time(value=10.0, units=reaction_pb2.Time.MINUTE)),
        (
            "-78°C",
            reaction_pb2.Temperature(
                value=-78.0, units=reaction_pb2.Temperature.CELSIUS
            ),
        ),
        (
            "10±5g",
            reaction_pb2.Mass(value=10, precision=5, units=reaction_pb2.Mass.GRAM),
        ),
        (" 10 meter", reaction_pb2.Length(value=10, units=reaction_pb2.Length.METER)),
        ("1.2e-3g", reaction_pb2.Mass(value=0.0012, units=reaction_pb2.Mass.GRAM)),
        (
            "0.12 nL",
            reaction_pb2.Volume(value=0.12, units=reaction_pb2.Volume.NANOLITER),
        ),
    ],
)
def test_resolve(resolver, string, expected):
    assert resolver.resolve(string) == expected


@pytest.mark.parametrize(
    ("message", "new_units", "expected"),
    [
        (
            reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER),
            "L",
            reaction_pb2.Volume(value=0.015, units=reaction_pb2.Volume.LITER),
        ),
        (
            reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR),
            "s",
            reaction_pb2.Time(value=86400, units=reaction_pb2.Time.SECOND),
        ),
        (
            reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR),
            reaction_pb2.Time.SECOND,
            reaction_pb2.Time(value=86400, units=reaction_pb2.Time.SECOND),
        ),
        (
            reaction_pb2.Temperature(
                value=-78.0, units=reaction_pb2.Temperature.CELSIUS
            ),
            "K",
            reaction_pb2.Temperature(
                value=195.15, units=reaction_pb2.Temperature.KELVIN
            ),
        ),
        (
            reaction_pb2.Temperature(
                value=-78.0, units=reaction_pb2.Temperature.CELSIUS, precision=5.2
            ),
            "K",
            reaction_pb2.Temperature(
                value=195.15, units=reaction_pb2.Temperature.KELVIN, precision=5.2
            ),
        ),
        (
            reaction_pb2.Temperature(
                value=25, units=reaction_pb2.Temperature.CELSIUS, precision=5
            ),
            "F",
            reaction_pb2.Temperature(
                value=77, units=reaction_pb2.Temperature.FAHRENHEIT, precision=9
            ),
        ),
        (
            reaction_pb2.Mass(value=0.0012, units=reaction_pb2.Mass.GRAM),
            "mg",
            reaction_pb2.Mass(value=1.2, units=reaction_pb2.Mass.MILLIGRAM),
        ),
        (
            reaction_pb2.Wavelength(value=500, units=reaction_pb2.Wavelength.NANOMETER),
            "cm^-1",
            reaction_pb2.Wavelength(
                value=20000, units=reaction_pb2.Wavelength.WAVENUMBER
            ),
        ),
        (
            reaction_pb2.Moles(value=0.0003, units=reaction_pb2.Moles.MOLE),
            "umol",
            reaction_pb2.Moles(value=300, units=reaction_pb2.Moles.MICROMOLE),
        ),
        (
            reaction_pb2.Pressure(
                value=1, units=reaction_pb2.Pressure.ATMOSPHERE, precision=0.1
            ),
            "bar",
            reaction_pb2.Pressure(
                value=1.01325, units=reaction_pb2.Pressure.BAR, precision=0.101325
            ),
        ),
        (
            reaction_pb2.FlowRate(
                value=2, units=reaction_pb2.FlowRate.MILLILITER_PER_SECOND
            ),
            "ul/s",
            reaction_pb2.FlowRate(
                value=2000, units=reaction_pb2.FlowRate.MICROLITER_PER_SECOND
            ),
        ),
        (
            reaction_pb2.Voltage(
                value=50, units=reaction_pb2.Voltage.VOLT, precision=2
            ),
            "mv",
            reaction_pb2.Voltage(
                value=50000, units=reaction_pb2.Voltage.MILLIVOLT, precision=2000
            ),
        ),
    ],
)
def test_convert(resolver, message, new_units, expected):
    assert resolver.convert(message, new_units) == expected


@pytest.mark.parametrize(
    ("message", "new_units", "expected"),
    [
        (
            reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER),
            "smoot",
            "unrecognized units",
        ),
        (
            reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER),
            "gram",
            "different types",
        ),
        (
            reaction_pb2.Temperature(value=5, units=reaction_pb2.Temperature.CELSIUS),
            "meter",
            "different types",
        ),
    ],
)
def test_convert_should_fail(resolver, message, new_units, expected):
    with pytest.raises((KeyError, ValueError), match=expected):
        resolver.convert(message, new_units)


@pytest.mark.parametrize(
    ("volume", "concentration", "expected"),
    [
        ("1L", "1 molar", "1 mol"),
        ("3mL", "0.1 molar", "300 micromoles"),
        ("100mL", "0.1 molar", "10 millimoles"),
    ],
)
def test_compute_solute_quantity(resolver, volume, concentration, expected):
    conc_resolver = units.UnitResolver(unit_synonyms=units.CONCENTRATION_UNIT_SYNONYMS)
    volume_pb = resolver.resolve(volume)
    concentration_pb = conc_resolver.resolve(concentration)
    assert isinstance(volume_pb, reaction_pb2.Volume)  # Type hint.
    assert isinstance(concentration_pb, reaction_pb2.Concentration)  # Type hint.
    assert units.compute_solute_quantity(
        volume=volume_pb,
        concentration=concentration_pb,
    ) == reaction_pb2.Amount(moles=resolver.resolve(expected))


@pytest.mark.parametrize(
    ("string", "expected"),
    [
        (
            "1-2 h",
            reaction_pb2.Time(value=1.5, precision=0.5, units=reaction_pb2.Time.HOUR),
        )
    ],
)
def test_resolve_allow_range(resolver, string, expected):
    assert resolver.resolve(string, allow_range=True) == expected


@pytest.mark.parametrize(
    ("string", "expected"),
    [
        ("1.21 GW", "unrecognized units"),
        ("15.0 ML 20.0 L", "string does not contain a value with units"),
        ("15.0. ML", "string does not contain a value with units"),
        ("5.2 m", "ambiguous"),
    ],
)
def test_resolve_should_fail(resolver, string, expected):
    with pytest.raises((KeyError, ValueError), match=expected):
        resolver.resolve(string)


def test_resolve_range_without_allow_range(resolver):
    with pytest.raises(ValueError, match="appears to contain a range"):
        resolver.resolve("1-2 h")


def test_resolver_init_rejects_duplicated_unit():
    duplicate_synonyms: dict[
        type[ord_schema.UnitMessage], dict[units.ProtoEnumMember, list[str]]
    ] = {
        reaction_pb2.Mass: {
            reaction_pb2.Mass.GRAM: ["g", "gram"],
            reaction_pb2.Mass.MILLIGRAM: ["g"],  # duplicate of "g".
        },
    }
    with pytest.raises(KeyError, match="duplicated unit"):
        units.UnitResolver(unit_synonyms=duplicate_synonyms, forbidden_units={})


@pytest.mark.parametrize(
    ("message", "new_units", "expected"),
    [
        (
            reaction_pb2.Temperature(
                value=77, units=reaction_pb2.Temperature.FAHRENHEIT, precision=9
            ),
            "C",
            reaction_pb2.Temperature(
                value=25, units=reaction_pb2.Temperature.CELSIUS, precision=5
            ),
        ),
        (
            reaction_pb2.Temperature(
                value=300, units=reaction_pb2.Temperature.KELVIN, precision=2
            ),
            "C",
            reaction_pb2.Temperature(
                value=26.85, units=reaction_pb2.Temperature.CELSIUS, precision=2
            ),
        ),
        (
            reaction_pb2.Wavelength(
                value=500, units=reaction_pb2.Wavelength.NANOMETER, precision=10
            ),
            "nm",
            reaction_pb2.Wavelength(
                value=500, units=reaction_pb2.Wavelength.NANOMETER, precision=10
            ),
        ),
    ],
)
def test_convert_precision(resolver, message, new_units, expected):
    actual = resolver.convert(message, new_units)
    assert actual.units == expected.units
    assert actual.value == pytest.approx(expected.value)
    assert actual.precision == pytest.approx(expected.precision)


def test_convert_wavelength_with_precision(resolver):
    # Expected per the implementation's documented formula:
    # (documented formula)
    # 10000000 / 2 * (1 / (value - precision) + 1 / (value + precision))  # noqa: ERA001
    # = 5e6 * (1/490 + 1/510) ≈ 20008.0032.
    converted = resolver.convert(
        reaction_pb2.Wavelength(
            value=500, units=reaction_pb2.Wavelength.NANOMETER, precision=10
        ),
        "cm^-1",
    )
    assert converted.units == reaction_pb2.Wavelength.WAVENUMBER
    assert converted.value == pytest.approx(20000)
    assert converted.precision == pytest.approx(20008.003201, rel=1e-6)


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        (reaction_pb2.Mass(value=1.5, units=reaction_pb2.Mass.GRAM), "1.5 g"),
        (
            reaction_pb2.Mass(value=1.5, units=reaction_pb2.Mass.GRAM, precision=0.1),
            "1.5 (± 0.1) g",
        ),
        (reaction_pb2.Time(value=10, units=reaction_pb2.Time.MINUTE), "10 min"),
    ],
)
def test_format_message(message, expected):
    assert units.format_message(message) == expected


def test_format_message_unspecified_returns_none():
    # Default-constructed message has UNSPECIFIED units.
    assert units.format_message(reaction_pb2.Mass()) is None


@pytest.mark.parametrize(
    ("string", "expected"),
    [
        (
            "760 torr",
            reaction_pb2.Pressure(value=760, units=reaction_pb2.Pressure.TORR),
        ),
        (
            "760 Torr",
            reaction_pb2.Pressure(value=760, units=reaction_pb2.Pressure.TORR),
        ),
        (
            "760 mmHg",
            reaction_pb2.Pressure(value=760, units=reaction_pb2.Pressure.MM_HG),
        ),
        (
            "760 mm Hg",
            reaction_pb2.Pressure(value=760, units=reaction_pb2.Pressure.MM_HG),
        ),
    ],
)
def test_resolve_mercury_pressure_units(resolver, string, expected):
    assert resolver.resolve(string) == expected


@pytest.mark.parametrize(
    "units_enum", [reaction_pb2.Pressure.TORR, reaction_pb2.Pressure.MM_HG]
)
def test_convert_mercury_pressure_units(resolver, units_enum):
    # 1 atm is 760 torr and, to the precision reaction data carries, 760 mmHg.
    message = reaction_pb2.Pressure(value=760, units=units_enum)
    assert resolver.convert(message, "atm").value == pytest.approx(1.0)
    assert resolver.convert(message, "kPa").value == pytest.approx(101.325, rel=1e-6)


def _unit_enum_values(message_type):
    """Returns {enum value: name} for a united message, minus UNSPECIFIED."""
    descriptor = message_type.DESCRIPTOR
    return {
        value.number: value.name
        for enum_type in descriptor.enum_types
        for value in enum_type.values
        if value.name != "UNSPECIFIED"
    }


def _united_message_types():
    """Returns every message in the schema carrying a value and a units enum.

    Read from the proto descriptor rather than from any hand-maintained list, so that a
    united message added to the schema is covered without a second edit here.
    """
    types = []
    for message_descriptor in reaction_pb2.DESCRIPTOR.message_types_by_name.values():
        fields = {field.name: field for field in message_descriptor.fields}
        units_field = fields.get("units")
        if (
            units_field is not None
            and units_field.type == descriptor.FieldDescriptor.TYPE_ENUM
            and "value" in fields
        ):
            types.append(getattr(reaction_pb2, message_descriptor.name))
    return types


# The checks below are parametrized from the schema, not from the tables under test:
# driving them from _UNIT_CONVERSIONS or _UNIT_SYNONYMS would generate no case at all
# for a message type missing from those, which is the failure they exist to catch. They
# assert on behavior, so it also does not matter which mechanism implements a
# conversion -- Temperature goes through Celsius and Wavelength through a wavenumber
# reciprocal, neither of which appears in _UNIT_CONVERSIONS.
_UNITED_MESSAGE_TYPES = _united_message_types()


def test_unit_message_alias_covers_the_schema():
    """``ord_schema.UnitMessage`` lists exactly the schema's united messages."""
    assert set(typing.get_args(ord_schema.UnitMessage)) == set(_UNITED_MESSAGE_TYPES)


def _synonyms_for(message_type):
    """Returns the synonym table covering a united message type.

    Concentration spellings ("M", "molar") are ambiguous with mass and mole spellings,
    so they live in a separate table used by a separate resolver.
    """
    if message_type in units.CONCENTRATION_UNIT_SYNONYMS:
        return units.CONCENTRATION_UNIT_SYNONYMS
    return units._UNIT_SYNONYMS


@pytest.mark.parametrize(
    "message_type", _UNITED_MESSAGE_TYPES, ids=lambda t: t.__name__
)
def test_every_unit_resolves_from_its_spellings(message_type):
    """Every unit in the schema has at least one usable spelling that maps back to it.

    Checked through ``resolve_unit`` rather than ``resolve``: a few spellings ("ul/s",
    "cm^-1") carry characters the value-plus-unit pattern does not accept, so they are
    reachable only when the unit is supplied on its own.
    """
    synonyms = _synonyms_for(message_type)
    resolver = units.UnitResolver(unit_synonyms=synonyms)
    for value, name in _unit_enum_values(message_type).items():
        spellings = synonyms.get(message_type, {}).get(value)
        assert spellings, f"{message_type.__name__}.{name} has no spelling"
        # "M" for molar collides with meter/minute and is forbidden outright.
        usable = [s for s in spellings if s.lower() not in units._FORBIDDEN_UNITS]
        assert usable, f"{message_type.__name__}.{name} has no usable spelling"
        for spelling in usable:
            assert resolver.resolve_unit(spelling) == (message_type, value)


@pytest.mark.parametrize(
    "message_type", _UNITED_MESSAGE_TYPES, ids=lambda t: t.__name__
)
def test_every_unit_converts_to_every_other(resolver, message_type):
    """Every pair of units within a message type converts."""
    unit_values = list(_unit_enum_values(message_type))
    for source in unit_values:
        message = message_type(value=1.0, units=source)
        for target in unit_values:
            assert resolver.convert(message, target).units == target


@pytest.mark.parametrize(
    ("message", "target", "expected"),
    [
        # An offset scale: the precision converts by the ratio alone, so the same
        # multiplier applied to the value would be wrong.
        (
            reaction_pb2.Temperature(
                value=77, units=reaction_pb2.Temperature.FAHRENHEIT, precision=1.8
            ),
            "K",
            1.0,
        ),
        (
            reaction_pb2.Temperature(
                value=25, units=reaction_pb2.Temperature.CELSIUS, precision=0.5
            ),
            "K",
            0.5,
        ),
        (
            reaction_pb2.Mass(
                value=500, units=reaction_pb2.Mass.MILLIGRAM, precision=10
            ),
            "g",
            0.01,
        ),
        (
            reaction_pb2.Time(value=90, units=reaction_pb2.Time.MINUTE, precision=1.5),
            "s",
            90.0,
        ),
        # Recorded zero is a precision the source stated, not a missing one.
        (
            reaction_pb2.Mass(value=1, units=reaction_pb2.Mass.GRAM, precision=0.0),
            "g",
            0.0,
        ),
    ],
)
def test_canonical_precision(message, target, expected):
    assert units.canonical_precision(message, target) == pytest.approx(expected)


@pytest.mark.parametrize(
    "message",
    [
        # No precision recorded.
        reaction_pb2.Mass(value=1, units=reaction_pb2.Mass.GRAM),
        # Precision recorded, but nothing says what unit it is in.
        reaction_pb2.Mass(value=1, precision=0.1),
        # Precision recorded against no value: an uncertainty published beside a null
        # value reads as a measurement nobody took but everybody bounded.
        reaction_pb2.Mass(precision=0.1, units=reaction_pb2.Mass.GRAM),
    ],
)
def test_canonical_precision_is_null(message):
    assert units.canonical_precision(message, "g") is None


@pytest.mark.parametrize(
    "message",
    [
        reaction_pb2.Temperature(precision=0.5, units=reaction_pb2.Temperature.CELSIUS),
        reaction_pb2.Time(precision=1.5, units=reaction_pb2.Time.MINUTE),
        reaction_pb2.Wavelength(precision=10, units=reaction_pb2.Wavelength.NANOMETER),
    ],
)
def test_precision_without_a_value_is_null_for_every_conversion_kind(message):
    # Offset scales, multiplier scales, and the inverting one all reach the value
    # through a different branch of convert(), so each is checked.
    target = {"Temperature": "K", "Time": "s", "Wavelength": "nm"}[
        type(message).__name__
    ]
    assert units.canonical_value(message, target) is None
    assert units.canonical_precision(message, target) is None
