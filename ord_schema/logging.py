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
"""Logging utilities."""

import logging

_initialized = False


def silence_rdkit_logs(pattern: str = "rdApp.*") -> None:
    """Disables noisy RDKit logs (use getattr; RDKit stubs omit ``DisableLog``)."""
    from rdkit import RDLogger

    getattr(RDLogger, "DisableLog")(pattern)


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Creates a Logger."""
    global _initialized
    if not _initialized:
        logging.basicConfig(format="%(levelname)s %(asctime)s %(filename)s:%(lineno)d: %(message)s")
        _initialized = True
    logger = logging.getLogger(name)
    logger.setLevel(level)
    return logger
