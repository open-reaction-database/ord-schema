"""Tests for ord_schema.orm.reaction."""
from ord_schema.orm import reaction
from ord_schema.proto import reaction_pb2


def test_orm():
    reaction.ProductMeasurementPercentage(value=23.4)


def test_float_value():
    message = reaction_pb2.FloatValue(value=1.2, precision=3.4)
    assert message == reaction.FloatValue.from_proto(message).to_proto()


if __name__ == "__main__":
    test_orm()
