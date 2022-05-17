# Copyright 2020 The Open Reaction Database Authors
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
"""Wrappers and utilities for handling protocol buffers in Python."""

import collections
import dataclasses
from typing import Union

import ord_schema

_MESSAGE_TYPES = (
    collections.abc.MutableMapping,  # Proto map.
    ord_schema.Message,  # Generic submessage.
)

# pytype: disable=attribute-error


@dataclasses.dataclass(eq=True, frozen=True)
class FrozenMessage(collections.abc.Mapping):
    """Container for a protocol buffer that does not allow edits.

    Notes:
        * For standard scalar values, it is not possible to distinguish between
          default values and explicitly set values that match the default. If
          the default is a valid value, add the `optional` label to the field.
          See https://github.com/Open-Reaction-Database/ord-schema/pull/174.
        * For `optional` scalar values and all submessage fields, exceptions
          are raised if the user attempts to access an undefined attribute
          (AttributeError), access an undefined map key (KeyError), or set
          any attribute or map value (dataclasses.FrozenInstanceError).
        * I considered adding a `raise_on_error` option that would return None
          instead of raising AttributeError or KeyError when requesting unset
          values. However, this breaks the guarantee that `hasattr` returns
          False for unset `optional` scalar values and submessages.
    """

    _message: Union[_MESSAGE_TYPES]

    def __getattr__(self, name: str):
        """Fetches a message attribute, if it exists.

        Notes:
            * If `name` is a submessage, a new FrozenMessage wrapping the
              submessage is returned.
            * If `name` has not been set explicitly, an AttributeError is
              raised. This allows for distinguishing between default values and
              values that were set explicitly (to the default).

        Args:
            name: Text attribute name.

        Returns:
            The requested attribute. If the attribute is a submessage, a new
            FrozenMessage is returned.

        Raises:
            AttributeError: if `name` has not been set explicitly.
        """
        try:
            if not self._message.HasField(name):
                raise AttributeError(f'attribute "{name}" has not been set')
        except ValueError:
            pass  # The requested attribute is not a submessage.
        value = getattr(self._message, name)
        if isinstance(value, _MESSAGE_TYPES):
            return FrozenMessage(value)
        if hasattr(value, "append"):
            # Make repeated fields (and their elements) immutable.
            repeated_values = []
            for element in value:
                if isinstance(element, _MESSAGE_TYPES):
                    repeated_values.append(FrozenMessage(element))
                else:
                    repeated_values.append(element)
            return tuple(repeated_values)
        return value

    def __iter__(self):
        return self._message.__iter__()

    def __len__(self):
        return len(self._message)

    def __getitem__(self, key):
        """Fetches a message key, if it exists.

        Notes:
            * If `key` is a submessage, a new FrozenMessage wrapping the
              submessage is returned.
            * If `key` has not been set explicitly, an AttributeError is
              raised. This allows for distinguishing between default values and
              values that were set explicitly (to the default).

        Args:
            key: Text key into the underlying message.

        Returns:
            The requested value. If the value is a submessage, a new
            FrozenMessage is returned.

        Raises:
            AttributeError: if `key` has not been set explicitly.
        """
        if key not in self._message:
            raise KeyError(f'key "{key}" has not been set')
        value = self._message[key]
        if isinstance(value, _MESSAGE_TYPES):
            return FrozenMessage(value)
        return value
