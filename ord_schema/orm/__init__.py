# Copyright 2022 Open Reaction Database Project Authors
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

"""Base ORM objects.

Deprecated: query the derived Parquet artifacts with ``ord_schema.search`` instead; see
``ord_schema/orm/README.md``.
"""

import warnings

from sqlalchemy.orm import declarative_base

from ord_schema.logging import silence_rdkit_logs

# stacklevel=2 attributes the warning to the importing module: the warnings module skips
# importlib's frames, so Python's default filters show it to code in __main__.
warnings.warn(
    "ord_schema.orm is deprecated and may be removed in a future release; query the "
    "derived Parquet artifacts with ord_schema.search instead. See "
    "https://github.com/open-reaction-database/ord-schema/blob/main/ord_schema/orm/README.md",
    DeprecationWarning,
    stacklevel=2,
)

silence_rdkit_logs()
Base = declarative_base()
