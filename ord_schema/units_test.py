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

import pytest

from ord_schema import units
from ord_schema.proto import reaction_pb2


@pytest.fixture
def resolver() -> units.UnitResolver:
    return units.UnitResolver()


@pytest.mark.parametrize(
    "string,expected",
    (
        ("15.0 ML", reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER)),
        ("24 H", reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR)),
        ("32.1g", reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM)),
        ("   32.1      \t   g  ", reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM)),
        ("10 min.", reaction_pb2.Time(value=10.0, units=reaction_pb2.Time.MINUTE)),
        ("-78°C", reaction_pb2.Temperature(value=-78.0, units=reaction_pb2.Temperature.CELSIUS)),
        ("10±5g", reaction_pb2.Mass(value=10, precision=5, units=reaction_pb2.Mass.GRAM)),
        (" 10 meter", reaction_pb2.Length(value=10, units=reaction_pb2.Length.METER)),
        ("1.2e-3g", reaction_pb2.Mass(value=0.0012, units=reaction_pb2.Mass.GRAM)),
        ("0.12 nL", reaction_pb2.Volume(value=0.12, units=reaction_pb2.Volume.NANOLITER)),
    ),
)
def test_resolve(resolver, string, expected):
    assert resolver.resolve(string) == expected


@pytest.mark.parametrize(
    "volume,concentration,expected",
    (("1L", "1 molar", "1 mol"), ("3mL", "0.1 molar", "300 micromoles"), ("100mL", "0.1 molar", "10 millimoles")),
)
def test_compute_solute_quantity(resolver, volume, concentration, expected):
    conc_resolver = units.UnitResolver(unit_synonyms=units.CONCENTRATION_UNIT_SYNONYMS)
    assert units.compute_solute_quantity(
        volume=resolver.resolve(volume),
        concentration=conc_resolver.resolve(concentration),
    ) == reaction_pb2.Amount(moles=resolver.resolve(expected))


@pytest.mark.parametrize(
    "string,expected",
    (("1-2 h", reaction_pb2.Time(value=1.5, precision=0.5, units=reaction_pb2.Time.HOUR)),),
)
def test_resolve_allow_range(resolver, string, expected):
    assert resolver.resolve(string, allow_range=True) == expected


@pytest.mark.parametrize(
    "string,expected",
    (
        ("1.21 GW", "unrecognized units"),
        ("15.0 ML 20.0 L", "string does not contain a value with units"),
        ("15.0. ML", "string does not contain a value with units"),
        ("5.2 m", "ambiguous"),
    ),
)
def test_resolve_should_fail(resolver, string, expected):
    with pytest.raises((KeyError, ValueError), match=expected):
        resolver.resolve(string)
