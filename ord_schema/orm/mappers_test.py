"""Tests for ord_schema.orm."""
import os
import pytest

from ord_schema.message_helpers import load_message
from ord_schema.orm.mappers import from_proto, to_proto
from ord_schema.proto.dataset_pb2 import Dataset


@pytest.mark.parametrize(
    "filename",
    (
        os.path.join(os.path.dirname(__file__), "testdata", "empty.pbtxt"),
        os.path.join(os.path.dirname(__file__), "testdata", "full.pbtxt"),
        os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"),
    ),
)
def test_round_trip(filename):
    dataset = load_message(filename, Dataset)
    assert dataset == to_proto(from_proto(dataset))
