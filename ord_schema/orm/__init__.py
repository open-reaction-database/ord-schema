"""ord-schema ORM."""
from inspect import getmro
from typing import Optional, Type

from google.protobuf.descriptor import FieldDescriptor
from google.protobuf.pyext._message import Message, MessageMapContainer

from ord_schema.orm import dataset
from ord_schema.orm import reaction
from ord_schema.orm.reaction import Base
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

MAPPERS: dict[Type[Message], Type[Base]] = {
    dataset_pb2.Dataset: dataset.Dataset,
    dataset_pb2.DatasetExample: dataset.DatasetExample,
    reaction_pb2.Reaction: reaction.Reaction,
    reaction_pb2.ReactionIdentifier: reaction.ReactionIdentifier,
    reaction_pb2.ReactionInput: reaction.ReactionInput,
    reaction_pb2.ReactionInput.AdditionSpeed: reaction.AdditionSpeed,
    reaction_pb2.ReactionInput.AdditionDevice: reaction.AdditionDevice,
    reaction_pb2.Amount: reaction.Amount,
    reaction_pb2.UnmeasuredAmount: reaction.UnmeasuredAmount,
    reaction_pb2.CrudeComponent: reaction.CrudeComponent,
    reaction_pb2.Compound: reaction.Compound,
    reaction_pb2.Compound.Source: reaction.Source,
    reaction_pb2.CompoundPreparation: reaction.CompoundPreparation,
    reaction_pb2.CompoundIdentifier: reaction.CompoundIdentifier,
    reaction_pb2.Vessel: reaction.Vessel,
    reaction_pb2.VesselMaterial: reaction.VesselMaterial,
    reaction_pb2.VesselAttachment: reaction.VesselAttachment,
    reaction_pb2.VesselPreparation: reaction.VesselPreparation,
    reaction_pb2.ReactionSetup: reaction.ReactionSetup,
    reaction_pb2.ReactionSetup.ReactionEnvironment: reaction.ReactionEnvironment,
    reaction_pb2.ReactionConditions: reaction.ReactionConditions,
    reaction_pb2.TemperatureConditions: reaction.TemperatureConditions,
    reaction_pb2.TemperatureConditions.TemperatureControl: reaction.TemperatureControl,
    reaction_pb2.TemperatureConditions.Measurement: reaction.TemperatureConditionsMeasurement,
    reaction_pb2.PressureConditions: reaction.PressureConditions,
    reaction_pb2.PressureConditions.PressureControl: reaction.PressureControl,
    reaction_pb2.PressureConditions.Atmosphere: reaction.Atmosphere,
    reaction_pb2.PressureConditions.Measurement: reaction.PressureConditionsMeasurement,
    reaction_pb2.StirringConditions: reaction.StirringConditions,
    reaction_pb2.StirringConditions.StirringRate: reaction.StirringRate,
    reaction_pb2.IlluminationConditions: reaction.IlluminationConditions,
    reaction_pb2.ElectrochemistryConditions: reaction.ElectrochemistryConditions,
    reaction_pb2.ElectrochemistryConditions.Measurement: reaction.ElectrochemistryConditionsMeasurement,
    reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell: reaction.ElectrochemistryCell,
    reaction_pb2.FlowConditions: reaction.FlowConditions,
    reaction_pb2.FlowConditions.Tubing: reaction.Tubing,
    reaction_pb2.ReactionNotes: reaction.ReactionNotes,
    reaction_pb2.ReactionObservation: reaction.ReactionObservation,
    reaction_pb2.ReactionWorkup: reaction.ReactionWorkup,
    reaction_pb2.ReactionOutcome: reaction.ReactionOutcome,
    reaction_pb2.ProductCompound: reaction.ProductCompound,
    reaction_pb2.ProductCompound.Texture: reaction.Texture,
    reaction_pb2.ProductMeasurement: reaction.ProductMeasurement,
    reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails: reaction.MassSpecMeasurementDetails,
    reaction_pb2.ProductMeasurement.Selectivity: reaction.Selectivity,
    reaction_pb2.DateTime: reaction.DateTime,
    reaction_pb2.Analysis: reaction.Analysis,
    reaction_pb2.ReactionProvenance: reaction.ReactionProvenance,
    reaction_pb2.Person: reaction.Person,
    reaction_pb2.RecordEvent: reaction.RecordEvent,
    reaction_pb2.Time: reaction.Time,
    reaction_pb2.Mass: reaction.Mass,
    reaction_pb2.Moles: reaction.Moles,
    reaction_pb2.Volume: reaction.Volume,
    reaction_pb2.Concentration: reaction.Concentration,
    reaction_pb2.Pressure: reaction.Pressure,
    reaction_pb2.Temperature: reaction.Temperature,
    reaction_pb2.Current: reaction.Current,
    reaction_pb2.Voltage: reaction.Voltage,
    reaction_pb2.Length: reaction.Length,
    reaction_pb2.Wavelength: reaction.Wavelength,
    reaction_pb2.FlowRate: reaction.FlowRate,
    reaction_pb2.Percentage: reaction.Percentage,
    reaction_pb2.FloatValue: reaction.FloatValue,
    reaction_pb2.Data: reaction.Data,
}
PROTOS: dict[Type[Base], Type[Message]] = {value: key for key, value in MAPPERS.items()}

MAPPER_RENAMES: dict[tuple[Type[Message], str], str] = {
    (reaction_pb2.Compound.Source, "id"): "vendor_id",
}
PROTO_RENAMES: dict[tuple[Type[Base], str], str] = {
    (MAPPERS[key[0]], value): key[1] for key, value in MAPPER_RENAMES.items()
}


def from_proto(message: Message, mapper: Optional[Type[Base]] = None, key: Optional[str] = None) -> Base:
    if mapper is None:
        mapper = MAPPERS[type(message)]
    kwargs = {}
    if key is not None:
        kwargs["key"] = key
    if mapper == reaction.Reaction:
        kwargs["proto"] = message.SerializeToString()
    for field, value in message.ListFields():
        field_name = MAPPER_RENAMES.get((type(message), field.name), field.name)
        if field_name == "eic_masses":
            # Convert repeated float to repeated FloatValue.
            kwargs[field_name] = [reaction.MassSpecMeasurementDetailsEicMasses(value=v) for v in value]
        elif field_name == "reaction_ids":
            # Convert repeated string to repeated ReactionId.
            kwargs[field_name] = [dataset.ReactionId(reaction_id=v) for v in value]
        elif field.type == FieldDescriptor.TYPE_MESSAGE:
            field_mapper = getattr(mapper, field_name).mapper.class_
            if isinstance(value, MessageMapContainer):
                kwargs[field_name] = [from_proto(v, mapper=field_mapper, key=k) for k, v in value.items()]
            elif field.label == FieldDescriptor.LABEL_REPEATED:
                kwargs[field_name] = [from_proto(v, mapper=field_mapper) for v in value]
            else:
                kwargs[field_name] = from_proto(value, mapper=field_mapper)
        elif field.type == FieldDescriptor.TYPE_ENUM:
            kwargs[field_name] = field.enum_type.values_by_number[value].name
        else:
            kwargs[field_name] = value
    return mapper(**kwargs)


def to_proto(base: reaction.Base) -> Message:
    kwargs = {}
    proto = None
    for mapper in getmro(type(base)):
        if not issubclass(mapper, reaction.Child):
            proto = PROTOS[mapper]
            break
    assert issubclass(proto, Message)
    for field in proto.DESCRIPTOR.fields:
        mapper_field_name = MAPPER_RENAMES.get((proto, field.name), field.name)
        value = getattr(base, mapper_field_name)
        field_name = PROTO_RENAMES.get((type(base), mapper_field_name), mapper_field_name)
        if field_name == "eic_masses":
            # Convert repeated FloatValue to repeated float.
            kwargs[field_name] = [v.value for v in value]
        elif field_name == "reaction_ids":
            # Convert repeated ReactionId to repeated string.
            kwargs[field_name] = [v.reaction_id for v in value]
        elif isinstance(value, list):
            if not len(value):
                continue
            if hasattr(value[0], "key"):
                kwargs[field_name] = {v.key: to_proto(v) for v in value}
            else:
                kwargs[field_name] = [to_proto(v) for v in value]
        elif isinstance(value, Base):
            kwargs[field_name] = to_proto(value)
        else:
            kwargs[field_name] = value
    return proto(**kwargs)
