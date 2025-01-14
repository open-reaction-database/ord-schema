/**
 * Copyright 2025 Open Reaction Database Project Authors
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

// package: ord
// file: ord-schema/proto/dataset.proto

import * as jspb from "google-protobuf";
import * as ord_schema_proto_reaction_pb from "../../ord-schema/proto/reaction_pb";

export class Dataset extends jspb.Message {
  getName(): string;
  setName(value: string): void;

  getDescription(): string;
  setDescription(value: string): void;

  clearReactionsList(): void;
  getReactionsList(): Array<ord_schema_proto_reaction_pb.Reaction>;
  setReactionsList(value: Array<ord_schema_proto_reaction_pb.Reaction>): void;
  addReactions(value?: ord_schema_proto_reaction_pb.Reaction, index?: number): ord_schema_proto_reaction_pb.Reaction;

  clearReactionIdsList(): void;
  getReactionIdsList(): Array<string>;
  setReactionIdsList(value: Array<string>): void;
  addReactionIds(value: string, index?: number): string;

  getDatasetId(): string;
  setDatasetId(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Dataset.AsObject;
  static toObject(includeInstance: boolean, msg: Dataset): Dataset.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Dataset, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Dataset;
  static deserializeBinaryFromReader(message: Dataset, reader: jspb.BinaryReader): Dataset;
}

export namespace Dataset {
  export type AsObject = {
    name: string,
    description: string,
    reactionsList: Array<ord_schema_proto_reaction_pb.Reaction.AsObject>,
    reactionIdsList: Array<string>,
    datasetId: string,
  }
}

export class DatasetExample extends jspb.Message {
  getDatasetId(): string;
  setDatasetId(value: string): void;

  getDescription(): string;
  setDescription(value: string): void;

  getUrl(): string;
  setUrl(value: string): void;

  hasCreated(): boolean;
  clearCreated(): void;
  getCreated(): ord_schema_proto_reaction_pb.RecordEvent | undefined;
  setCreated(value?: ord_schema_proto_reaction_pb.RecordEvent): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): DatasetExample.AsObject;
  static toObject(includeInstance: boolean, msg: DatasetExample): DatasetExample.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: DatasetExample, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): DatasetExample;
  static deserializeBinaryFromReader(message: DatasetExample, reader: jspb.BinaryReader): DatasetExample;
}

export namespace DatasetExample {
  export type AsObject = {
    datasetId: string,
    description: string,
    url: string,
    created?: ord_schema_proto_reaction_pb.RecordEvent.AsObject,
  }
}

