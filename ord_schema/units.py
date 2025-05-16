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
"""Helpers for translating strings with units."""

import re
from typing import Optional, Tuple, Type, Union

import numpy as np

import ord_schema
from ord_schema.proto import reaction_pb2

# Accepted synonyms for units. Note that all values will be converted to
# lowercase.
_UNIT_SYNONYMS = {
    reaction_pb2.Time: {
        reaction_pb2.Time.DAY: ["d", "day", "days"],
        reaction_pb2.Time.HOUR: ["h", "hour", "hours", "hr", "hrs"],
        reaction_pb2.Time.MINUTE: ["min", "mins", "minute", "minutes"],
        reaction_pb2.Time.SECOND: ["s", "sec", "secs", "second", "seconds"],
    },
    reaction_pb2.Mass: {
        reaction_pb2.Mass.GRAM: ["g", "gram", "grams", "gs", "gm", "gms", "gr"],
        reaction_pb2.Mass.MILLIGRAM: ["mg", "mgs", "milligrams", "milligram"],
        reaction_pb2.Mass.MICROGRAM: [
            "μg",
            "ug",
            "ugs",
            "micg",
            "micgs",
            "micrograms",
            "microgram",
        ],
        reaction_pb2.Mass.KILOGRAM: ["kg", "kgs", "kilogram", "kilograms"],
    },
    reaction_pb2.Moles: {
        reaction_pb2.Moles.MOLE: ["mol", "mols", "mole", "moles"],
        reaction_pb2.Moles.MILLIMOLE: [
            "mmol",
            "millimoles",
            "mmols",
            "mmole",
            "mmoles",
        ],
        reaction_pb2.Moles.MICROMOLE: ["μmol", "umol", "umols", "micromoles"],
        reaction_pb2.Moles.NANOMOLE: ["nmol", "nanomoles"],
    },
    reaction_pb2.Volume: {
        reaction_pb2.Volume.MILLILITER: [
            "mL",
            "milliliter",
            "milliliters",
            "cc",
            "cm3",
            "mls",
        ],
        reaction_pb2.Volume.MICROLITER: [
            "μL",
            "uL",
            "micl",
            "microliter",
            "microliters",
        ],
        reaction_pb2.Volume.LITER: ["L", "liter", "liters", "litres"],
        reaction_pb2.Volume.NANOLITER: ["nL", "nanoliter", "nanoliters"],
    },
    reaction_pb2.Length: {
        reaction_pb2.Length.CENTIMETER: ["cm", "centimeter"],
        reaction_pb2.Length.MILLIMETER: ["millimeter", "millimeters"],
        reaction_pb2.Length.METER: ["meter", "meters"],
        reaction_pb2.Length.INCH: ["in", "inch", "inches"],
        reaction_pb2.Length.FOOT: ["ft", "foot", "feet"],
    },
    reaction_pb2.Pressure: {
        reaction_pb2.Pressure.BAR: ["bar", "barg", "bars"],
        reaction_pb2.Pressure.ATMOSPHERE: ["atm", "atmosphere", "atmospheres"],
        reaction_pb2.Pressure.PSI: ["psi"],
        reaction_pb2.Pressure.KPSI: ["kpsi"],
        reaction_pb2.Pressure.PASCAL: ["Pa", "pascal", "pascals", "pas"],
        reaction_pb2.Pressure.KILOPASCAL: ["kPa", "kilopascals", "kPas"],
    },
    reaction_pb2.Temperature: {
        reaction_pb2.Temperature.CELSIUS: [
            "°C",
            "C",
            "degC",
            "°celsius",
            "celsius",
            "degrees C",
        ],
        reaction_pb2.Temperature.FAHRENHEIT: ["°F", "F", "degF", "fahrenheit"],
        reaction_pb2.Temperature.KELVIN: ["K", "degK", "Kelvin"],
    },
    reaction_pb2.Current: {
        reaction_pb2.Current.AMPERE: ["A", "ampere", "amps", "amp"],
        reaction_pb2.Current.MILLIAMPERE: [
            "mA",
            "milliampere",
            "milliamp",
            "milliamps",
        ],
    },
    reaction_pb2.Voltage: {
        reaction_pb2.Voltage.VOLT: ["V", "volt", "volts"],
        reaction_pb2.Voltage.MILLIVOLT: ["mV", "millivolt", "millivolts"],
    },
    reaction_pb2.Wavelength: {
        reaction_pb2.Wavelength.NANOMETER: ["nm", "nanometer", "nanometers"],
        reaction_pb2.Wavelength.WAVENUMBER: [
            "cm⁻¹",
            "cm^-1",
            "cm-1",
            "wavenumber",
            "1/cm",
        ],
    },
    reaction_pb2.FlowRate: {
        reaction_pb2.FlowRate.MICROLITER_PER_MINUTE: ["μL/min", "uL/min"],
        reaction_pb2.FlowRate.MICROLITER_PER_SECOND: ["μL/s", "uL/s"],
        reaction_pb2.FlowRate.MILLILITER_PER_MINUTE: ["mL/min"],
        reaction_pb2.FlowRate.MILLILITER_PER_SECOND: ["mL/s"],
        reaction_pb2.FlowRate.MICROLITER_PER_HOUR: ["μL/h", "uL/h"],
    },
}

_FORBIDDEN_UNITS = {
    "m": "ambiguous between meter and minute",
}

_UNIT_CONVERSIONS = {
    reaction_pb2.Time: {
        reaction_pb2.Time.DAY: 24,
        reaction_pb2.Time.HOUR: 1,
        reaction_pb2.Time.MINUTE: 1 / 60,
        reaction_pb2.Time.SECOND: 1 / 3600,
    },
    reaction_pb2.Mass: {
        reaction_pb2.Mass.GRAM: 1000,
        reaction_pb2.Mass.MILLIGRAM: 1,
        reaction_pb2.Mass.MICROGRAM: 0.001,
        reaction_pb2.Mass.KILOGRAM: 1000000,
    },
    reaction_pb2.Moles: {
        reaction_pb2.Moles.MOLE: 1000,
        reaction_pb2.Moles.MILLIMOLE: 1,
        reaction_pb2.Moles.MICROMOLE: 0.001,
        reaction_pb2.Moles.NANOMOLE: 0.000001,
    },
    reaction_pb2.Volume: {
        reaction_pb2.Volume.MILLILITER: 1,
        reaction_pb2.Volume.MICROLITER: 0.001,
        reaction_pb2.Volume.LITER: 1000,
        reaction_pb2.Volume.NANOLITER: 0.000001,
    },
    reaction_pb2.Length: {
        reaction_pb2.Length.CENTIMETER: 1,
        reaction_pb2.Length.MILLIMETER: 0.1,
        reaction_pb2.Length.METER: 100,
        reaction_pb2.Length.INCH: 2.54,
        reaction_pb2.Length.FOOT: 30.48,
    },
    reaction_pb2.Pressure: {
        reaction_pb2.Pressure.BAR: 1 / 1.01325,
        reaction_pb2.Pressure.ATMOSPHERE: 1,
        reaction_pb2.Pressure.PSI: 1 / 14.69595,
        reaction_pb2.Pressure.KPSI: 1000 / 14.69595,
        reaction_pb2.Pressure.PASCAL: 1 / 101325,
        reaction_pb2.Pressure.KILOPASCAL: 1000 / 101325,
    },
    reaction_pb2.Current: {
        reaction_pb2.Current.AMPERE: 1,
        reaction_pb2.Current.MILLIAMPERE: 0.001,
    },
    reaction_pb2.Voltage: {
        reaction_pb2.Voltage.VOLT: 1,
        reaction_pb2.Voltage.MILLIVOLT: 0.001,
    },
    reaction_pb2.FlowRate: {
        reaction_pb2.FlowRate.MICROLITER_PER_MINUTE: 1,
        reaction_pb2.FlowRate.MICROLITER_PER_SECOND: 60,
        reaction_pb2.FlowRate.MILLILITER_PER_MINUTE: 1000,
        reaction_pb2.FlowRate.MILLILITER_PER_SECOND: 60000,
        reaction_pb2.FlowRate.MICROLITER_PER_HOUR: 1 / 60,
    },
}

# Concentration units are defined separately since they are not needed for any
# native fields in the reaction schema.
CONCENTRATION_UNIT_SYNONYMS = {
    reaction_pb2.Concentration: {
        reaction_pb2.Concentration.MOLAR: ["M", "molar"],
        reaction_pb2.Concentration.MILLIMOLAR: ["mM", "millimolar"],
        reaction_pb2.Concentration.MICROMOLAR: ["uM", "micromolar"],
    },
}

CONCENTRATION_M_PER_UNIT = {
    reaction_pb2.Concentration.ConcentrationUnit.MOLAR: 1,
    reaction_pb2.Concentration.ConcentrationUnit.MILLIMOLAR: 1e-3,
    reaction_pb2.Concentration.ConcentrationUnit.MICROMOLAR: 1e-6,
}


class UnitResolver:
    """Resolver class for translating value+unit strings into messages."""

    def __init__(
        self,
        unit_synonyms: Optional[dict[Type[ord_schema.UnitMessage], dict[ord_schema.Message, list[str]]]] = None,
        forbidden_units: Optional[dict[str, str]] = None,
    ):
        """Initializes a UnitResolver.

        Args:
            unit_synonyms: A dictionary of dictionaries that defines, for each
                message type (first key) and for each unit option (second key),
                a list of strings that defines how that unit can be written.
                Defaults to None. If None, uses default _UNIT_SYNONYMS dict.
            forbidden_units: A dictionary where each key is a string that is a
                plausible way of writing a unit and a value explaining why
                the UnitResolver cannot resolve that unit. The prototypical
                case is one of ambiguity (e.g., "m" can mean meter or minute).
                Defaults to None. If None, uses default _FORBIDDEN_UNITS dict.
                If no units are forbidden, an empty dictionary should be used.
        """
        if unit_synonyms is None:
            unit_synonyms = _UNIT_SYNONYMS
        if forbidden_units is None:
            forbidden_units = _FORBIDDEN_UNITS
        self._forbidden_units = forbidden_units
        self._resolver = {}
        for message in unit_synonyms:
            for unit in unit_synonyms[message]:
                for string_unit in unit_synonyms[message][unit]:
                    string_unit = string_unit.lower()
                    if string_unit in self._resolver:
                        raise KeyError(f"duplicated unit: {string_unit}")
                    self._resolver[string_unit] = (message, unit)
        # Values must have zero or one decimal point. Whitespace between the
        # value and the unit is optional.
        self._pattern = re.compile(r"(-?\d+\.?\d*(?:[eE]-?\d+)?)(?:[-±](-?\d+\.?\d*))?\s*([\w\sμ°]+)\.?")

    def resolve(self, string: str, allow_range: bool = False) -> ord_schema.UnitMessage:
        """Resolves a string into a message containing a value with units.

        Args:
            string: The string to parse; must contain a numeric value and a
                string unit. For example: "1.25 h".
            allow_range: If True, ranges like "1-2 h" can be provided and the
                average value will be reported along with the standard
                deviation.

        Returns:
            Message containing a numeric value with units listed in the schema.

        Raises:
            ValueError: if string does not contain a value with units, or if
                the value is invalid.
        """
        # NOTE(kearnes): Use fullmatch() to catch cases with multiple matches.
        match = self._pattern.fullmatch(string.strip().replace("−", "-"))
        if not match:
            raise ValueError(f"string does not contain a value with units: {string}")
        value, range_value, string_unit = match.groups()
        precision = None
        if range_value is not None:
            if "±" in string:
                value = float(value)
                precision = float(range_value)
            elif not allow_range:
                raise ValueError("string appears to contain a range of values " f"but allow_range is False: {string}")
            else:
                values = np.asarray([value, range_value], dtype=float)
                value = values.mean()
                precision = values.std()
        else:
            value = float(value)
        assert string_unit is not None  # Type hint.
        message, unit = self.resolve_unit(string_unit)
        if value < 0.0 and message != reaction_pb2.Temperature:
            raise ValueError(f"negative values are only allowed for temperature: {string}")
        if precision:
            return message(value=value, precision=precision, units=unit)
        return message(value=value, units=unit)

    def resolve_unit(self, string_unit: str) -> Tuple[ord_schema.UnitMessage, ord_schema.Message]:
        """Resolves a unit string into its message type and unit ENUM value.

        Args:
            string_unit: The string unit to parse; for example: "gram".

        Returns:
            Tuple containing the message type and unit ENUM value.

        Raises:
            KeyError: if string unit cannot be parsed.
        """
        string_unit = string_unit.lower()
        if string_unit in self._forbidden_units:
            raise KeyError(f"forbidden units: {string_unit}: " f"({self._forbidden_units[string_unit]})")
        if string_unit not in self._resolver:
            raise KeyError(f"unrecognized units: {string_unit}")
        message, unit = self._resolver[string_unit]
        return (message, unit)

    def convert(self, message: ord_schema.UnitMessage, new_units: Union[str, int]) -> ord_schema.UnitMessage:
        """Converts a united message into another united message of the same
        type, but with different units.

        Args:
            message: a message with units, e.g., Mass, Length.
            new_units: the desired units of the new message, expressed
                either as a string or an integer (ENUM value). Use of a
                string is recommended due to the ambiguity of using
                ENUM values; for example, Mass.GRAM == Time.MINUTE.

        Returns:
            A new message with units, e.g., Mass, Length.
        """
        # pylint: disable=too-many-branches
        message_type = type(message)
        if isinstance(new_units, str):
            (new_message_type, new_units) = self.resolve_unit(new_units)
            if new_message_type != message_type:
                raise ValueError("Cannot convert units between messages of different types")
        new_message = message_type(units=new_units)

        if message_type == reaction_pb2.Temperature:
            temperature_celsius = {
                reaction_pb2.Temperature.CELSIUS: lambda T: T,
                reaction_pb2.Temperature.FAHRENHEIT: lambda T: (T - 32) * 5 / 9,
                reaction_pb2.Temperature.KELVIN: lambda T: T - 273.15,
            }[message.units](message.value)
            new_message.value = {
                reaction_pb2.Temperature.CELSIUS: lambda T: T,
                reaction_pb2.Temperature.FAHRENHEIT: lambda T: T * 9 / 5 + 32,
                reaction_pb2.Temperature.KELVIN: lambda T: T + 273.15,
            }[new_units](temperature_celsius)
            if message.precision:
                if message.units == reaction_pb2.Temperature.FAHRENHEIT and new_units in (
                    reaction_pb2.Temperature.CELSIUS,
                    reaction_pb2.Temperature.KELVIN,
                ):
                    new_message.precision = message.precision * 5 / 9
                elif (
                    message.units in (reaction_pb2.Temperature.CELSIUS, reaction_pb2.Temperature.KELVIN)
                    and new_units == reaction_pb2.Temperature.FAHRENHEIT
                ):
                    new_message.precision = message.precision * 9 / 5
                else:
                    new_message.precision = message.precision

        elif message_type == reaction_pb2.Wavelength:
            if message.units != new_units:
                new_message.value = 10000000 / message.value
                # Approximate precision as stdev of previous value \pm precision
                if message.precision:
                    new_message.precision = (
                        10000000
                        / 2
                        * (1 / (message.value - message.precision) + 1 / (message.value + message.precision))
                    )
            else:
                new_message.value = message.value
                if message.precision:
                    new_message.precision = message.precision

        else:
            # All other cases can be handled by a multiplier conversion.
            new_message.value = (
                message.value
                * _UNIT_CONVERSIONS[message_type][message.units]
                / _UNIT_CONVERSIONS[message_type][new_units]
            )
            if message.precision:
                new_message.precision = (
                    message.precision
                    * _UNIT_CONVERSIONS[message_type][message.units]
                    / _UNIT_CONVERSIONS[message_type][new_units]
                )
        return new_message


def format_message(message: ord_schema.UnitMessage) -> Optional[str]:
    """Formats a united message into a string.

    Args:
        message: a message with units, e.g., Mass, Length.

    Returns:
        A string describing the value, e.g., "5.0 (p/m 0.1) mL" using the
            first unit synonym listed in _UNIT_SYNONYMS.
    """
    if message.units == getattr(type(message)(), "UNSPECIFIED"):
        return None
    txt = f"{message.value:.7g} "
    if message.precision:
        txt += f"(± {message.precision:.7g}) "
    txt += _UNIT_SYNONYMS[type(message)][message.units][0]
    return txt


def compute_solute_quantity(
    volume: reaction_pb2.Volume, concentration: reaction_pb2.Concentration
) -> reaction_pb2.Amount:
    """Computes the quantity of a solute, given volume and concentration."""
    volume_conversion = _UNIT_CONVERSIONS[reaction_pb2.Volume][volume.units]
    volume_liter = volume.value * volume_conversion / 1000
    concentration_conversion = CONCENTRATION_M_PER_UNIT[concentration.units]
    concentration_molar = concentration.value * concentration_conversion

    moles = volume_liter * concentration_molar
    if moles < 1e-6:
        value = moles * 1e9
        unit = reaction_pb2.Moles.NANOMOLE
    elif moles < 1e-3:
        value = moles * 1e6
        unit = reaction_pb2.Moles.MICROMOLE
    elif moles < 1:
        value = moles * 1e3
        unit = reaction_pb2.Moles.MILLIMOLE
    else:
        value = moles
        unit = reaction_pb2.Moles.MOLE
    return reaction_pb2.Amount(moles=reaction_pb2.Moles(value=value, units=unit))
