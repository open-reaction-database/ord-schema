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

from ord_schema.proto import reaction_pb2

# Accepted synonyms for units. Note that all values will be converted to
# lowercase.
_UNIT_SYNONYMS = {
    reaction_pb2.Time: {
        reaction_pb2.Time.HOUR: ['h', 'hour', 'hours', 'hr', 'hrs'],
        reaction_pb2.Time.MINUTE: ['min', 'mins', 'minute', 'minutes'],
        reaction_pb2.Time.SECOND: ['s', 'sec', 'secs', 'second', 'seconds'],
    },
    reaction_pb2.Mass: {
        reaction_pb2.Mass.GRAM: ['g', 'gram', 'grams', 'gs', 'gm', 'gms'],
        reaction_pb2.Mass.MILLIGRAM: ['mg', 'mgs', 'milligrams', 'milligram'],
        reaction_pb2.Mass.MICROGRAM: [
            'ug', 'ugs', 'micg', 'micgs', 'micrograms', 'microgram'
        ],
        reaction_pb2.Mass.KILOGRAM: ['kg', 'kgs', 'kilogram', 'kilograms'],
    },
    reaction_pb2.Moles: {
        reaction_pb2.Moles.MOLE: ['mol', 'mols', 'mole', 'moles'],
        reaction_pb2.Moles.MILLIMOLE: ['mmol', 'millimoles', 'mmols'],
        reaction_pb2.Moles.MICROMOLE: ['umol', 'umols', 'micromoles'],
        reaction_pb2.Moles.NANOMOLE: ['nmol', 'nanomoles'],
    },
    reaction_pb2.Volume: {
        reaction_pb2.Volume.MILLILITER: ['mL', 'milliliters'],
        reaction_pb2.Volume.MICROLITER: ['uL', 'micl', 'microliters'],
        reaction_pb2.Volume.LITER: ['L', 'liters', 'litres'],
    },
    reaction_pb2.Length: {
        reaction_pb2.Length.CENTIMETER: ['cm', 'centimeter'],
        reaction_pb2.Length.MILLIMETER: ['millimeter', 'millimeters'],
        reaction_pb2.Length.METER: ['meter', 'meters'],
        reaction_pb2.Length.INCH: ['in', 'inch', 'inches'],
        reaction_pb2.Length.FOOT: ['ft', 'foot', 'feet'],
    },
    reaction_pb2.Pressure: {
        reaction_pb2.Pressure.BAR: ['bar', 'barg', 'bars'],
        reaction_pb2.Pressure.ATMOSPHERE: ['atm', 'atmosphere', 'atmospheres'],
        reaction_pb2.Pressure.PSI: ['psi'],
        reaction_pb2.Pressure.KPSI: ['kpsi'],
        reaction_pb2.Pressure.PASCAL: ['Pa', 'pascal', 'pascals', 'pas'],
        reaction_pb2.Pressure.KILOPASCAL: ['kPa', 'kilopascals', 'kPas'],
    },
    reaction_pb2.Temperature: {
        reaction_pb2.Temperature.CELSIUS: ['C', 'degC', 'celsius'],
        reaction_pb2.Temperature.FAHRENHEIT: ['F', 'degF', 'fahrenheit'],
        reaction_pb2.Temperature.KELVIN: ['K', 'degK', 'Kelvin'],
    },
    reaction_pb2.Current: {
        reaction_pb2.Current.AMPERE: ['A', 'ampere', 'amps', 'amp'],
        reaction_pb2.Current.MILLIAMPERE: [
            'mA', 'milliampere', 'milliamp', 'milliamps'
        ],
    },
    reaction_pb2.Voltage: {
        reaction_pb2.Voltage.VOLT: ['V', 'volt', 'volts'],
        reaction_pb2.Voltage.MILLIVOLT: ['mV', 'millivolt', 'millivolts'],
    },
    reaction_pb2.Wavelength: {
        reaction_pb2.Wavelength.NANOMETER: ['nm', 'nanometer', 'nanometers'],
        reaction_pb2.Wavelength.WAVENUMBER: ['cm-1', 'wavenumber', '1/cm'],
    },
    reaction_pb2.FlowRate: {
        reaction_pb2.FlowRate.MICROLITER_PER_MINUTE: ['uL/min'],
        reaction_pb2.FlowRate.MICROLITER_PER_SECOND: ['uL/s'],
        reaction_pb2.FlowRate.MILLILITER_PER_MINUTE: ['mL/min'],
        reaction_pb2.FlowRate.MILLILITER_PER_SECOND: ['mL/s'],
        reaction_pb2.FlowRate.MICROLITER_PER_HOUR: ['uL/h'],
    },
}

_FORBIDDEN_UNITS = {
    'm': 'ambiguous between meter and minute',
}

# Concentration units are defined separately since they are not needed for any
# native fields in the reaction schema.
CONCENTRATION_UNIT_SYNONYMS = {
    reaction_pb2.Concentration: {
        reaction_pb2.Concentration.MOLAR: ['M', 'molar'],
        reaction_pb2.Concentration.MILLIMOLAR: ['mM', 'millimolar'],
        reaction_pb2.Concentration.MICROMOLAR: ['uM', 'micromolar'],
    },
}


class UnitResolver:
    """Resolver class for translating value+unit strings into messages."""
    def __init__(self, unit_synonyms=None, forbidden_units=None):
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

        Returns:
            None
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
                        raise KeyError(f'duplicated unit: {string_unit}')
                    self._resolver[string_unit] = (message, unit)
        # Values must have zero or one decimal point. Whitespace between the
        # value and the unit is optional.
        self._pattern = re.compile(r'(\d+.?\d*)\s*(\w+)')

    def resolve(self, string):
        """Resolves a string into a message containing a value with units.

        Args:
            string: The string to parse; must contain a numeric value and a
                string unit. For example: "1.25 h".

        Returns:
            Message containing a numeric value with units listed in the schema.

        Raises:
            ValueError: If string does not contain a value with units.
        """
        # NOTE(kearnes): Use fullmatch() to catch cases with multiple matches.
        match = self._pattern.fullmatch(string.strip())
        if not match:
            raise ValueError(
                f'string does not contain a value with units: {string}')
        value, string_unit = match.groups()
        string_unit = string_unit.lower()
        if string_unit in self._forbidden_units:
            raise KeyError(f'forbidden units: {string_unit}: '
                           f'({self._forbidden_units[string_unit]})')
        if string_unit not in self._resolver:
            raise KeyError(f'unrecognized units: {string_unit}')
        message, unit = self._resolver[string_unit]
        return message(value=float(value), units=unit)


def format_message(message):
    """Formats a united message into a string.

    Args:
        message: a message with units, e.g., Mass, Length.

    Returns:
        A string describing the value, e.g., "5.0 (p/m 0.1) mL" using the
            first unit synonym listed in _UNIT_SYNONYMS.
    """
    if message.units == getattr(type(message)(), 'UNSPECIFIED'):
        return None
    txt = f'{message.value:.4g} '
    if message.precision:
        txt += f'(p/m {message.precision}) '
    txt += _UNIT_SYNONYMS[type(message)][message.units][0]
    return txt
