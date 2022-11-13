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

"""Tests for ord_schema.orm.query."""
import pytest

from ord_schema.orm.mappers import Percentage, ProductCompound, ProductMeasurement, Reaction, ReactionOutcome
from ord_schema.orm.query import get_join_path


@pytest.mark.parametrize(
    "source,target,expected",
    (
        (
            Reaction,
            ProductMeasurement.percentage,
            [Reaction, ReactionOutcome, ProductCompound, ProductMeasurement, Percentage],
        ),
        (
            Reaction,
            Percentage.value,
            [Reaction, ReactionOutcome, Percentage],  # Returns the first path it finds!
        ),
    ),
)
def test_get_join_path(source, target, expected):
    assert get_join_path(source, target) == expected


@pytest.mark.parametrize(
    "source,target",
    ((ProductMeasurement, ProductCompound.measurements),),
)
def test_get_join_path_fails(source, target):
    with pytest.raises(ValueError, match="could not find a path"):
        get_join_path(source, target)
