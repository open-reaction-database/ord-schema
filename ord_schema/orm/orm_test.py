"""Tests for ord_schema.orm."""
import os

from ord_schema.message_helpers import load_message
from ord_schema.orm import from_proto, to_proto
from ord_schema.proto.dataset_pb2 import Dataset


def test_round_trip():
    dataset = load_message(os.path.join(os.path.dirname(__file__), "testdata", "full.pbtxt"), Dataset)
    assert dataset == to_proto(from_proto(dataset))
