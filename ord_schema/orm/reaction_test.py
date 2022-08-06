"""Tests for ord_schema.orm.reaction."""
import os

from ord_schema.message_helpers import load_message
from ord_schema.orm.reaction import Reaction, from_proto, to_proto
from ord_schema.proto.dataset_pb2 import Dataset


def test_orm():
    Reaction()


def test_round_trip():
    dataset = load_message(os.path.join(os.path.dirname(__file__), "testdata", "full.pbtxt"), Dataset)
    for reaction in dataset.reactions:
        assert reaction == to_proto(from_proto(reaction, Reaction))


if __name__ == "__main__":
    test_orm()
