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

"""Hash-partition predicate shared by the sharded derived passes.

Kept in its own module so ``database`` and ``reaction_class`` can both use it without an import
cycle (``database`` imports ``reaction_class``).
"""


def shard_predicate(
    column: str, shard: tuple[int, int] | None
) -> tuple[str, dict[str, int]]:
    """Returns an ``AND`` predicate keeping 1/num_shards of rows by ``column``, or ``('', {})``.

    Partitions a dataset's rows deterministically and disjointly by a hash of ``column``, so
    independent workers can each process their shard without coordinating. ``shard`` is
    ``(index, num_shards)``; ``None`` disables sharding (whole dataset).

    ``column`` is interpolated into the SQL, so it must be a trusted, hardcoded column reference
    (e.g. ``"ord.reaction.id"``) -- never a user-supplied string. The shard bounds are passed as
    bind parameters.
    """
    if shard is None:
        return "", {}
    shard_index, num_shards = shard
    # ((h % n) + n) % n keeps the bucket in [0, num_shards) without abs() overflow.
    predicate = (
        f"AND ((hashtextextended({column}::text, 0) % :num_shards) + :num_shards) "
        "% :num_shards = :shard_index"
    )
    return predicate, {"num_shards": num_shards, "shard_index": shard_index}
