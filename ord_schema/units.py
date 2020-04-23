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
        reaction_pb2.Mass.MICROGRAM: ['ug', 'ugs', 'micg', 'micgs',
                                      'micrograms', 'microgram'],
        reaction_pb2.Mass.KILOGRAM: ['kg', 'kgs', 'kilogram', 'kilograms'],
    },
    reaction_pb2.Moles: {
        reaction_pb2.Moles.MOLES: ['mol', 'mols', 'mole', 'moles'],
        reaction_pb2.Moles.MILLIMOLES: ['mmol', 'millimoles', 'mmols'],
        reaction_pb2.Moles.MICROMOLES: ['umol', 'umols', 'micromoles'],
        reaction_pb2.Moles.NANOMOLES: ['nmol', 'nanomoles'],
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
    reaction_pb2.Concentration: {
        reaction_pb2.Concentration.MOLAR: ['molar'],
        reaction_pb2.Concentration.MILLIMOLAR: ['millimolar'],
        reaction_pb2.Concentration.MICROMOLAR: ['micromolar'],
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
        reaction_pb2.Current.MILLIAMPERE: ['mA', 'milliampere', 'milliamp',
                                           'milliamps'],
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
    'mm': 'ambiguous between millimeter and millimolar',
}


class UnitResolver:
    """Resolver class for translating value+unit strings into messages."""

    def __init__(self):
        self._resolver = {}
        for message in _UNIT_SYNONYMS:
            for unit in _UNIT_SYNONYMS[message]:
                for string_unit in _UNIT_SYNONYMS[message][unit]:
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
        if string_unit in _FORBIDDEN_UNITS:
            raise KeyError(f'forbidden units: {string_unit}: '
                           f'({_FORBIDDEN_UNITS[string_unit]})')
        if string_unit not in self._resolver:
            raise KeyError(f'unrecognized units: {string_unit}')
        message, unit = self._resolver[string_unit]
        return message(value=float(value), units=unit)
