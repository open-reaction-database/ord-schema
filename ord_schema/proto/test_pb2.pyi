from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class Scalar(_message.Message):
    __slots__ = ["int32_value", "int64_value", "float_value", "string_value", "bytes_value", "bool_value"]
    INT32_VALUE_FIELD_NUMBER: _ClassVar[int]
    INT64_VALUE_FIELD_NUMBER: _ClassVar[int]
    FLOAT_VALUE_FIELD_NUMBER: _ClassVar[int]
    STRING_VALUE_FIELD_NUMBER: _ClassVar[int]
    BYTES_VALUE_FIELD_NUMBER: _ClassVar[int]
    BOOL_VALUE_FIELD_NUMBER: _ClassVar[int]
    int32_value: int
    int64_value: int
    float_value: float
    string_value: str
    bytes_value: bytes
    bool_value: bool
    def __init__(self, int32_value: _Optional[int] = ..., int64_value: _Optional[int] = ..., float_value: _Optional[float] = ..., string_value: _Optional[str] = ..., bytes_value: _Optional[bytes] = ..., bool_value: bool = ...) -> None: ...

class RepeatedScalar(_message.Message):
    __slots__ = ["values"]
    VALUES_FIELD_NUMBER: _ClassVar[int]
    values: _containers.RepeatedScalarFieldContainer[float]
    def __init__(self, values: _Optional[_Iterable[float]] = ...) -> None: ...

class Enum(_message.Message):
    __slots__ = ["value"]
    class EnumValues(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Enum.EnumValues]
        FIRST: _ClassVar[Enum.EnumValues]
        SECOND: _ClassVar[Enum.EnumValues]
    UNSPECIFIED: Enum.EnumValues
    FIRST: Enum.EnumValues
    SECOND: Enum.EnumValues
    VALUE_FIELD_NUMBER: _ClassVar[int]
    value: Enum.EnumValues
    def __init__(self, value: _Optional[_Union[Enum.EnumValues, str]] = ...) -> None: ...

class RepeatedEnum(_message.Message):
    __slots__ = ["values"]
    class EnumValues(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[RepeatedEnum.EnumValues]
        FIRST: _ClassVar[RepeatedEnum.EnumValues]
        SECOND: _ClassVar[RepeatedEnum.EnumValues]
    UNSPECIFIED: RepeatedEnum.EnumValues
    FIRST: RepeatedEnum.EnumValues
    SECOND: RepeatedEnum.EnumValues
    VALUES_FIELD_NUMBER: _ClassVar[int]
    values: _containers.RepeatedScalarFieldContainer[RepeatedEnum.EnumValues]
    def __init__(self, values: _Optional[_Iterable[_Union[RepeatedEnum.EnumValues, str]]] = ...) -> None: ...

class Nested(_message.Message):
    __slots__ = ["child"]
    class Child(_message.Message):
        __slots__ = ["value"]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        value: float
        def __init__(self, value: _Optional[float] = ...) -> None: ...
    CHILD_FIELD_NUMBER: _ClassVar[int]
    child: Nested.Child
    def __init__(self, child: _Optional[_Union[Nested.Child, _Mapping]] = ...) -> None: ...

class RepeatedNested(_message.Message):
    __slots__ = ["children"]
    class Child(_message.Message):
        __slots__ = ["value"]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        value: float
        def __init__(self, value: _Optional[float] = ...) -> None: ...
    CHILDREN_FIELD_NUMBER: _ClassVar[int]
    children: _containers.RepeatedCompositeFieldContainer[RepeatedNested.Child]
    def __init__(self, children: _Optional[_Iterable[_Union[RepeatedNested.Child, _Mapping]]] = ...) -> None: ...

class Map(_message.Message):
    __slots__ = ["values"]
    class ValuesEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: float
        def __init__(self, key: _Optional[str] = ..., value: _Optional[float] = ...) -> None: ...
    VALUES_FIELD_NUMBER: _ClassVar[int]
    values: _containers.ScalarMap[str, float]
    def __init__(self, values: _Optional[_Mapping[str, float]] = ...) -> None: ...

class MapNested(_message.Message):
    __slots__ = ["children"]
    class Child(_message.Message):
        __slots__ = ["value"]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        value: float
        def __init__(self, value: _Optional[float] = ...) -> None: ...
    class ChildrenEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: MapNested.Child
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[MapNested.Child, _Mapping]] = ...) -> None: ...
    CHILDREN_FIELD_NUMBER: _ClassVar[int]
    children: _containers.MessageMap[str, MapNested.Child]
    def __init__(self, children: _Optional[_Mapping[str, MapNested.Child]] = ...) -> None: ...
