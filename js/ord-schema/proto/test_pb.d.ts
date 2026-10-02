/**
 * Copyright 2024 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// package: ord_test
// file: ord-schema/proto/test.proto

import * as jspb from "google-protobuf";

export class Scalar extends jspb.Message {
  getInt32Value(): number;
  setInt32Value(value: number): void;

  getInt64Value(): number;
  setInt64Value(value: number): void;

  hasFloatValue(): boolean;
  clearFloatValue(): void;
  getFloatValue(): number;
  setFloatValue(value: number): void;

  getStringValue(): string;
  setStringValue(value: string): void;

  getBytesValue(): Uint8Array | string;
  getBytesValue_asU8(): Uint8Array;
  getBytesValue_asB64(): string;
  setBytesValue(value: Uint8Array | string): void;

  hasBoolValue(): boolean;
  clearBoolValue(): void;
  getBoolValue(): boolean;
  setBoolValue(value: boolean): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Scalar.AsObject;
  static toObject(includeInstance: boolean, msg: Scalar): Scalar.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Scalar, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Scalar;
  static deserializeBinaryFromReader(message: Scalar, reader: jspb.BinaryReader): Scalar;
}

export namespace Scalar {
  export type AsObject = {
    int32Value: number,
    int64Value: number,
    floatValue: number,
    stringValue: string,
    bytesValue: Uint8Array | string,
    boolValue: boolean,
  }
}

export class RepeatedScalar extends jspb.Message {
  clearValuesList(): void;
  getValuesList(): Array<number>;
  setValuesList(value: Array<number>): void;
  addValues(value: number, index?: number): number;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): RepeatedScalar.AsObject;
  static toObject(includeInstance: boolean, msg: RepeatedScalar): RepeatedScalar.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: RepeatedScalar, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): RepeatedScalar;
  static deserializeBinaryFromReader(message: RepeatedScalar, reader: jspb.BinaryReader): RepeatedScalar;
}

export namespace RepeatedScalar {
  export type AsObject = {
    valuesList: Array<number>,
  }
}

export class Enum extends jspb.Message {
  getValue(): Enum.EnumValuesMap[keyof Enum.EnumValuesMap];
  setValue(value: Enum.EnumValuesMap[keyof Enum.EnumValuesMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Enum.AsObject;
  static toObject(includeInstance: boolean, msg: Enum): Enum.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Enum, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Enum;
  static deserializeBinaryFromReader(message: Enum, reader: jspb.BinaryReader): Enum;
}

export namespace Enum {
  export type AsObject = {
    value: Enum.EnumValuesMap[keyof Enum.EnumValuesMap],
  }

  export interface EnumValuesMap {
    UNSPECIFIED: 0;
    FIRST: 1;
    SECOND: 2;
  }

  export const EnumValues: EnumValuesMap;
}

export class RepeatedEnum extends jspb.Message {
  clearValuesList(): void;
  getValuesList(): Array<RepeatedEnum.EnumValuesMap[keyof RepeatedEnum.EnumValuesMap]>;
  setValuesList(value: Array<RepeatedEnum.EnumValuesMap[keyof RepeatedEnum.EnumValuesMap]>): void;
  addValues(value: RepeatedEnum.EnumValuesMap[keyof RepeatedEnum.EnumValuesMap], index?: number): RepeatedEnum.EnumValuesMap[keyof RepeatedEnum.EnumValuesMap];

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): RepeatedEnum.AsObject;
  static toObject(includeInstance: boolean, msg: RepeatedEnum): RepeatedEnum.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: RepeatedEnum, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): RepeatedEnum;
  static deserializeBinaryFromReader(message: RepeatedEnum, reader: jspb.BinaryReader): RepeatedEnum;
}

export namespace RepeatedEnum {
  export type AsObject = {
    valuesList: Array<RepeatedEnum.EnumValuesMap[keyof RepeatedEnum.EnumValuesMap]>,
  }

  export interface EnumValuesMap {
    UNSPECIFIED: 0;
    FIRST: 1;
    SECOND: 2;
  }

  export const EnumValues: EnumValuesMap;
}

export class Nested extends jspb.Message {
  hasChild(): boolean;
  clearChild(): void;
  getChild(): Nested.Child | undefined;
  setChild(value?: Nested.Child): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Nested.AsObject;
  static toObject(includeInstance: boolean, msg: Nested): Nested.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Nested, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Nested;
  static deserializeBinaryFromReader(message: Nested, reader: jspb.BinaryReader): Nested;
}

export namespace Nested {
  export type AsObject = {
    child?: Nested.Child.AsObject,
  }

  export class Child extends jspb.Message {
    hasValue(): boolean;
    clearValue(): void;
    getValue(): number;
    setValue(value: number): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Child.AsObject;
    static toObject(includeInstance: boolean, msg: Child): Child.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Child, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Child;
    static deserializeBinaryFromReader(message: Child, reader: jspb.BinaryReader): Child;
  }

  export namespace Child {
    export type AsObject = {
      value: number,
    }
  }
}

export class RepeatedNested extends jspb.Message {
  clearChildrenList(): void;
  getChildrenList(): Array<RepeatedNested.Child>;
  setChildrenList(value: Array<RepeatedNested.Child>): void;
  addChildren(value?: RepeatedNested.Child, index?: number): RepeatedNested.Child;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): RepeatedNested.AsObject;
  static toObject(includeInstance: boolean, msg: RepeatedNested): RepeatedNested.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: RepeatedNested, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): RepeatedNested;
  static deserializeBinaryFromReader(message: RepeatedNested, reader: jspb.BinaryReader): RepeatedNested;
}

export namespace RepeatedNested {
  export type AsObject = {
    childrenList: Array<RepeatedNested.Child.AsObject>,
  }

  export class Child extends jspb.Message {
    hasValue(): boolean;
    clearValue(): void;
    getValue(): number;
    setValue(value: number): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Child.AsObject;
    static toObject(includeInstance: boolean, msg: Child): Child.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Child, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Child;
    static deserializeBinaryFromReader(message: Child, reader: jspb.BinaryReader): Child;
  }

  export namespace Child {
    export type AsObject = {
      value: number,
    }
  }
}

export class Map extends jspb.Message {
  getValuesMap(): jspb.Map<string, number>;
  clearValuesMap(): void;
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Map.AsObject;
  static toObject(includeInstance: boolean, msg: Map): Map.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Map, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Map;
  static deserializeBinaryFromReader(message: Map, reader: jspb.BinaryReader): Map;
}

export namespace Map {
  export type AsObject = {
    valuesMap: Array<[string, number]>,
  }
}

export class MapNested extends jspb.Message {
  getChildrenMap(): jspb.Map<string, MapNested.Child>;
  clearChildrenMap(): void;
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): MapNested.AsObject;
  static toObject(includeInstance: boolean, msg: MapNested): MapNested.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: MapNested, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): MapNested;
  static deserializeBinaryFromReader(message: MapNested, reader: jspb.BinaryReader): MapNested;
}

export namespace MapNested {
  export type AsObject = {
    childrenMap: Array<[string, MapNested.Child.AsObject]>,
  }

  export class Child extends jspb.Message {
    hasValue(): boolean;
    clearValue(): void;
    getValue(): number;
    setValue(value: number): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Child.AsObject;
    static toObject(includeInstance: boolean, msg: Child): Child.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Child, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Child;
    static deserializeBinaryFromReader(message: Child, reader: jspb.BinaryReader): Child;
  }

  export namespace Child {
    export type AsObject = {
      value: number,
    }
  }
}

