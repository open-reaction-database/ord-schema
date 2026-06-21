# Copyright 2026 Open Reaction Database Project Authors
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
"""Deprecated alias module for :mod:`ord_schema.parquet`.

``ord_schema.parquet_dataset`` was renamed to ``ord_schema.parquet``. This
module forwards attribute access to the new module for backwards compatibility
and will be removed in a future minor release.
"""

import warnings
from typing import Any

from ord_schema import parquet


def __getattr__(name: str) -> Any:
    warnings.warn(
        "ord_schema.parquet_dataset was renamed to ord_schema.parquet; update your imports.",
        DeprecationWarning,
        stacklevel=2,
    )
    return getattr(parquet, name)
