from ord_schema.proto import reaction_pb2 as _reaction_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class Dataset(_message.Message):
    __slots__ = ["name", "description", "reactions", "reaction_ids", "dataset_id"]
    NAME_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    REACTIONS_FIELD_NUMBER: _ClassVar[int]
    REACTION_IDS_FIELD_NUMBER: _ClassVar[int]
    DATASET_ID_FIELD_NUMBER: _ClassVar[int]
    name: str
    description: str
    reactions: _containers.RepeatedCompositeFieldContainer[_reaction_pb2.Reaction]
    reaction_ids: _containers.RepeatedScalarFieldContainer[str]
    dataset_id: str
    def __init__(self, name: _Optional[str] = ..., description: _Optional[str] = ..., reactions: _Optional[_Iterable[_Union[_reaction_pb2.Reaction, _Mapping]]] = ..., reaction_ids: _Optional[_Iterable[str]] = ..., dataset_id: _Optional[str] = ...) -> None: ...

class DatasetExample(_message.Message):
    __slots__ = ["dataset_id", "description", "url", "created"]
    DATASET_ID_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    URL_FIELD_NUMBER: _ClassVar[int]
    CREATED_FIELD_NUMBER: _ClassVar[int]
    dataset_id: str
    description: str
    url: str
    created: _reaction_pb2.RecordEvent
    def __init__(self, dataset_id: _Optional[str] = ..., description: _Optional[str] = ..., url: _Optional[str] = ..., created: _Optional[_Union[_reaction_pb2.RecordEvent, _Mapping]] = ...) -> None: ...
