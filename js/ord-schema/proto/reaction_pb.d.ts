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

// package: ord
// file: ord-schema/proto/reaction.proto

import * as jspb from "google-protobuf";

export class Reaction extends jspb.Message {
  clearIdentifiersList(): void;
  getIdentifiersList(): Array<ReactionIdentifier>;
  setIdentifiersList(value: Array<ReactionIdentifier>): void;
  addIdentifiers(value?: ReactionIdentifier, index?: number): ReactionIdentifier;

  getInputsMap(): jspb.Map<string, ReactionInput>;
  clearInputsMap(): void;
  hasSetup(): boolean;
  clearSetup(): void;
  getSetup(): ReactionSetup | undefined;
  setSetup(value?: ReactionSetup): void;

  hasConditions(): boolean;
  clearConditions(): void;
  getConditions(): ReactionConditions | undefined;
  setConditions(value?: ReactionConditions): void;

  hasNotes(): boolean;
  clearNotes(): void;
  getNotes(): ReactionNotes | undefined;
  setNotes(value?: ReactionNotes): void;

  clearObservationsList(): void;
  getObservationsList(): Array<ReactionObservation>;
  setObservationsList(value: Array<ReactionObservation>): void;
  addObservations(value?: ReactionObservation, index?: number): ReactionObservation;

  clearWorkupsList(): void;
  getWorkupsList(): Array<ReactionWorkup>;
  setWorkupsList(value: Array<ReactionWorkup>): void;
  addWorkups(value?: ReactionWorkup, index?: number): ReactionWorkup;

  clearOutcomesList(): void;
  getOutcomesList(): Array<ReactionOutcome>;
  setOutcomesList(value: Array<ReactionOutcome>): void;
  addOutcomes(value?: ReactionOutcome, index?: number): ReactionOutcome;

  hasProvenance(): boolean;
  clearProvenance(): void;
  getProvenance(): ReactionProvenance | undefined;
  setProvenance(value?: ReactionProvenance): void;

  getReactionId(): string;
  setReactionId(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Reaction.AsObject;
  static toObject(includeInstance: boolean, msg: Reaction): Reaction.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Reaction, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Reaction;
  static deserializeBinaryFromReader(message: Reaction, reader: jspb.BinaryReader): Reaction;
}

export namespace Reaction {
  export type AsObject = {
    identifiersList: Array<ReactionIdentifier.AsObject>,
    inputsMap: Array<[string, ReactionInput.AsObject]>,
    setup?: ReactionSetup.AsObject,
    conditions?: ReactionConditions.AsObject,
    notes?: ReactionNotes.AsObject,
    observationsList: Array<ReactionObservation.AsObject>,
    workupsList: Array<ReactionWorkup.AsObject>,
    outcomesList: Array<ReactionOutcome.AsObject>,
    provenance?: ReactionProvenance.AsObject,
    reactionId: string,
  }
}

export class ReactionIdentifier extends jspb.Message {
  getType(): ReactionIdentifier.ReactionIdentifierTypeMap[keyof ReactionIdentifier.ReactionIdentifierTypeMap];
  setType(value: ReactionIdentifier.ReactionIdentifierTypeMap[keyof ReactionIdentifier.ReactionIdentifierTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  getValue(): string;
  setValue(value: string): void;

  hasIsMapped(): boolean;
  clearIsMapped(): void;
  getIsMapped(): boolean;
  setIsMapped(value: boolean): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionIdentifier.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionIdentifier): ReactionIdentifier.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionIdentifier, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionIdentifier;
  static deserializeBinaryFromReader(message: ReactionIdentifier, reader: jspb.BinaryReader): ReactionIdentifier;
}

export namespace ReactionIdentifier {
  export type AsObject = {
    type: ReactionIdentifier.ReactionIdentifierTypeMap[keyof ReactionIdentifier.ReactionIdentifierTypeMap],
    details: string,
    value: string,
    isMapped: boolean,
  }

  export interface ReactionIdentifierTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    REACTION_SMILES: 2;
    REACTION_CXSMILES: 6;
    RDFILE: 3;
    RINCHI: 4;
    REACTION_TYPE: 5;
  }

  export const ReactionIdentifierType: ReactionIdentifierTypeMap;
}

export class ReactionInput extends jspb.Message {
  clearComponentsList(): void;
  getComponentsList(): Array<Compound>;
  setComponentsList(value: Array<Compound>): void;
  addComponents(value?: Compound, index?: number): Compound;

  clearCrudeComponentsList(): void;
  getCrudeComponentsList(): Array<CrudeComponent>;
  setCrudeComponentsList(value: Array<CrudeComponent>): void;
  addCrudeComponents(value?: CrudeComponent, index?: number): CrudeComponent;

  getAdditionOrder(): number;
  setAdditionOrder(value: number): void;

  hasAdditionTime(): boolean;
  clearAdditionTime(): void;
  getAdditionTime(): Time | undefined;
  setAdditionTime(value?: Time): void;

  hasAdditionSpeed(): boolean;
  clearAdditionSpeed(): void;
  getAdditionSpeed(): ReactionInput.AdditionSpeed | undefined;
  setAdditionSpeed(value?: ReactionInput.AdditionSpeed): void;

  hasAdditionDuration(): boolean;
  clearAdditionDuration(): void;
  getAdditionDuration(): Time | undefined;
  setAdditionDuration(value?: Time): void;

  hasFlowRate(): boolean;
  clearFlowRate(): void;
  getFlowRate(): FlowRate | undefined;
  setFlowRate(value?: FlowRate): void;

  hasAdditionDevice(): boolean;
  clearAdditionDevice(): void;
  getAdditionDevice(): ReactionInput.AdditionDevice | undefined;
  setAdditionDevice(value?: ReactionInput.AdditionDevice): void;

  hasAdditionTemperature(): boolean;
  clearAdditionTemperature(): void;
  getAdditionTemperature(): Temperature | undefined;
  setAdditionTemperature(value?: Temperature): void;

  hasTexture(): boolean;
  clearTexture(): void;
  getTexture(): Texture | undefined;
  setTexture(value?: Texture): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionInput.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionInput): ReactionInput.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionInput, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionInput;
  static deserializeBinaryFromReader(message: ReactionInput, reader: jspb.BinaryReader): ReactionInput;
}

export namespace ReactionInput {
  export type AsObject = {
    componentsList: Array<Compound.AsObject>,
    crudeComponentsList: Array<CrudeComponent.AsObject>,
    additionOrder: number,
    additionTime?: Time.AsObject,
    additionSpeed?: ReactionInput.AdditionSpeed.AsObject,
    additionDuration?: Time.AsObject,
    flowRate?: FlowRate.AsObject,
    additionDevice?: ReactionInput.AdditionDevice.AsObject,
    additionTemperature?: Temperature.AsObject,
    texture?: Texture.AsObject,
  }

  export class AdditionSpeed extends jspb.Message {
    getType(): ReactionInput.AdditionSpeed.AdditionSpeedTypeMap[keyof ReactionInput.AdditionSpeed.AdditionSpeedTypeMap];
    setType(value: ReactionInput.AdditionSpeed.AdditionSpeedTypeMap[keyof ReactionInput.AdditionSpeed.AdditionSpeedTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): AdditionSpeed.AsObject;
    static toObject(includeInstance: boolean, msg: AdditionSpeed): AdditionSpeed.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: AdditionSpeed, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): AdditionSpeed;
    static deserializeBinaryFromReader(message: AdditionSpeed, reader: jspb.BinaryReader): AdditionSpeed;
  }

  export namespace AdditionSpeed {
    export type AsObject = {
      type: ReactionInput.AdditionSpeed.AdditionSpeedTypeMap[keyof ReactionInput.AdditionSpeed.AdditionSpeedTypeMap],
      details: string,
    }

    export interface AdditionSpeedTypeMap {
      UNSPECIFIED: 0;
      ALL_AT_ONCE: 1;
      FAST: 2;
      SLOW: 3;
      DROPWISE: 4;
      CONTINUOUS: 5;
      PORTIONWISE: 6;
    }

    export const AdditionSpeedType: AdditionSpeedTypeMap;
  }

  export class AdditionDevice extends jspb.Message {
    getType(): ReactionInput.AdditionDevice.AdditionDeviceTypeMap[keyof ReactionInput.AdditionDevice.AdditionDeviceTypeMap];
    setType(value: ReactionInput.AdditionDevice.AdditionDeviceTypeMap[keyof ReactionInput.AdditionDevice.AdditionDeviceTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): AdditionDevice.AsObject;
    static toObject(includeInstance: boolean, msg: AdditionDevice): AdditionDevice.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: AdditionDevice, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): AdditionDevice;
    static deserializeBinaryFromReader(message: AdditionDevice, reader: jspb.BinaryReader): AdditionDevice;
  }

  export namespace AdditionDevice {
    export type AsObject = {
      type: ReactionInput.AdditionDevice.AdditionDeviceTypeMap[keyof ReactionInput.AdditionDevice.AdditionDeviceTypeMap],
      details: string,
    }

    export interface AdditionDeviceTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      NONE: 2;
      SYRINGE: 3;
      CANNULA: 4;
      ADDITION_FUNNEL: 5;
      PIPETTE: 6;
      POSITIVE_DISPLACEMENT_PIPETTE: 7;
      PISTON_PUMP: 8;
      SYRINGE_PUMP: 9;
      PERISTALTIC_PUMP: 10;
    }

    export const AdditionDeviceType: AdditionDeviceTypeMap;
  }
}

export class Amount extends jspb.Message {
  hasMass(): boolean;
  clearMass(): void;
  getMass(): Mass | undefined;
  setMass(value?: Mass): void;

  hasMoles(): boolean;
  clearMoles(): void;
  getMoles(): Moles | undefined;
  setMoles(value?: Moles): void;

  hasVolume(): boolean;
  clearVolume(): void;
  getVolume(): Volume | undefined;
  setVolume(value?: Volume): void;

  hasUnmeasured(): boolean;
  clearUnmeasured(): void;
  getUnmeasured(): UnmeasuredAmount | undefined;
  setUnmeasured(value?: UnmeasuredAmount): void;

  hasVolumeIncludesSolutes(): boolean;
  clearVolumeIncludesSolutes(): void;
  getVolumeIncludesSolutes(): boolean;
  setVolumeIncludesSolutes(value: boolean): void;

  getKindCase(): Amount.KindCase;
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Amount.AsObject;
  static toObject(includeInstance: boolean, msg: Amount): Amount.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Amount, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Amount;
  static deserializeBinaryFromReader(message: Amount, reader: jspb.BinaryReader): Amount;
}

export namespace Amount {
  export type AsObject = {
    mass?: Mass.AsObject,
    moles?: Moles.AsObject,
    volume?: Volume.AsObject,
    unmeasured?: UnmeasuredAmount.AsObject,
    volumeIncludesSolutes: boolean,
  }

  export enum KindCase {
    KIND_NOT_SET = 0,
    MASS = 1,
    MOLES = 2,
    VOLUME = 3,
    UNMEASURED = 5,
  }
}

export class UnmeasuredAmount extends jspb.Message {
  getType(): UnmeasuredAmount.UnmeasuredAmountTypeMap[keyof UnmeasuredAmount.UnmeasuredAmountTypeMap];
  setType(value: UnmeasuredAmount.UnmeasuredAmountTypeMap[keyof UnmeasuredAmount.UnmeasuredAmountTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): UnmeasuredAmount.AsObject;
  static toObject(includeInstance: boolean, msg: UnmeasuredAmount): UnmeasuredAmount.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: UnmeasuredAmount, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): UnmeasuredAmount;
  static deserializeBinaryFromReader(message: UnmeasuredAmount, reader: jspb.BinaryReader): UnmeasuredAmount;
}

export namespace UnmeasuredAmount {
  export type AsObject = {
    type: UnmeasuredAmount.UnmeasuredAmountTypeMap[keyof UnmeasuredAmount.UnmeasuredAmountTypeMap],
    details: string,
  }

  export interface UnmeasuredAmountTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    SATURATED: 2;
    CATALYTIC: 3;
    TITRATED: 4;
  }

  export const UnmeasuredAmountType: UnmeasuredAmountTypeMap;
}

export class Texture extends jspb.Message {
  getType(): Texture.TextureTypeMap[keyof Texture.TextureTypeMap];
  setType(value: Texture.TextureTypeMap[keyof Texture.TextureTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Texture.AsObject;
  static toObject(includeInstance: boolean, msg: Texture): Texture.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Texture, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Texture;
  static deserializeBinaryFromReader(message: Texture, reader: jspb.BinaryReader): Texture;
}

export namespace Texture {
  export type AsObject = {
    type: Texture.TextureTypeMap[keyof Texture.TextureTypeMap],
    details: string,
  }

  export interface TextureTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    POWDER: 2;
    CRYSTAL: 3;
    OIL: 4;
    AMORPHOUS_SOLID: 5;
    FOAM: 6;
    WAX: 7;
    SEMI_SOLID: 8;
    SOLID: 9;
    LIQUID: 10;
    GAS: 11;
  }

  export const TextureType: TextureTypeMap;
}

export class CrudeComponent extends jspb.Message {
  getReactionId(): string;
  setReactionId(value: string): void;

  hasIncludesWorkup(): boolean;
  clearIncludesWorkup(): void;
  getIncludesWorkup(): boolean;
  setIncludesWorkup(value: boolean): void;

  hasHasDerivedAmount(): boolean;
  clearHasDerivedAmount(): void;
  getHasDerivedAmount(): boolean;
  setHasDerivedAmount(value: boolean): void;

  hasAmount(): boolean;
  clearAmount(): void;
  getAmount(): Amount | undefined;
  setAmount(value?: Amount): void;

  hasTexture(): boolean;
  clearTexture(): void;
  getTexture(): Texture | undefined;
  setTexture(value?: Texture): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): CrudeComponent.AsObject;
  static toObject(includeInstance: boolean, msg: CrudeComponent): CrudeComponent.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: CrudeComponent, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): CrudeComponent;
  static deserializeBinaryFromReader(message: CrudeComponent, reader: jspb.BinaryReader): CrudeComponent;
}

export namespace CrudeComponent {
  export type AsObject = {
    reactionId: string,
    includesWorkup: boolean,
    hasDerivedAmount: boolean,
    amount?: Amount.AsObject,
    texture?: Texture.AsObject,
  }
}

export class Compound extends jspb.Message {
  clearIdentifiersList(): void;
  getIdentifiersList(): Array<CompoundIdentifier>;
  setIdentifiersList(value: Array<CompoundIdentifier>): void;
  addIdentifiers(value?: CompoundIdentifier, index?: number): CompoundIdentifier;

  hasAmount(): boolean;
  clearAmount(): void;
  getAmount(): Amount | undefined;
  setAmount(value?: Amount): void;

  getReactionRole(): ReactionRole.ReactionRoleTypeMap[keyof ReactionRole.ReactionRoleTypeMap];
  setReactionRole(value: ReactionRole.ReactionRoleTypeMap[keyof ReactionRole.ReactionRoleTypeMap]): void;

  hasIsLimiting(): boolean;
  clearIsLimiting(): void;
  getIsLimiting(): boolean;
  setIsLimiting(value: boolean): void;

  clearPreparationsList(): void;
  getPreparationsList(): Array<CompoundPreparation>;
  setPreparationsList(value: Array<CompoundPreparation>): void;
  addPreparations(value?: CompoundPreparation, index?: number): CompoundPreparation;

  hasSource(): boolean;
  clearSource(): void;
  getSource(): Compound.Source | undefined;
  setSource(value?: Compound.Source): void;

  getFeaturesMap(): jspb.Map<string, Data>;
  clearFeaturesMap(): void;
  getAnalysesMap(): jspb.Map<string, Analysis>;
  clearAnalysesMap(): void;
  hasTexture(): boolean;
  clearTexture(): void;
  getTexture(): Texture | undefined;
  setTexture(value?: Texture): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Compound.AsObject;
  static toObject(includeInstance: boolean, msg: Compound): Compound.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Compound, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Compound;
  static deserializeBinaryFromReader(message: Compound, reader: jspb.BinaryReader): Compound;
}

export namespace Compound {
  export type AsObject = {
    identifiersList: Array<CompoundIdentifier.AsObject>,
    amount?: Amount.AsObject,
    reactionRole: ReactionRole.ReactionRoleTypeMap[keyof ReactionRole.ReactionRoleTypeMap],
    isLimiting: boolean,
    preparationsList: Array<CompoundPreparation.AsObject>,
    source?: Compound.Source.AsObject,
    featuresMap: Array<[string, Data.AsObject]>,
    analysesMap: Array<[string, Analysis.AsObject]>,
    texture?: Texture.AsObject,
  }

  export class Source extends jspb.Message {
    getVendor(): string;
    setVendor(value: string): void;

    getCatalogId(): string;
    setCatalogId(value: string): void;

    getLot(): string;
    setLot(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Source.AsObject;
    static toObject(includeInstance: boolean, msg: Source): Source.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Source, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Source;
    static deserializeBinaryFromReader(message: Source, reader: jspb.BinaryReader): Source;
  }

  export namespace Source {
    export type AsObject = {
      vendor: string,
      catalogId: string,
      lot: string,
    }
  }
}

export class ReactionRole extends jspb.Message {
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionRole.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionRole): ReactionRole.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionRole, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionRole;
  static deserializeBinaryFromReader(message: ReactionRole, reader: jspb.BinaryReader): ReactionRole;
}

export namespace ReactionRole {
  export type AsObject = {
  }

  export interface ReactionRoleTypeMap {
    UNSPECIFIED: 0;
    REACTANT: 1;
    REAGENT: 2;
    SOLVENT: 3;
    CATALYST: 4;
    WORKUP: 5;
    INTERNAL_STANDARD: 6;
    AUTHENTIC_STANDARD: 7;
    PRODUCT: 8;
    BYPRODUCT: 9;
    SIDE_PRODUCT: 10;
  }

  export const ReactionRoleType: ReactionRoleTypeMap;
}

export class CompoundPreparation extends jspb.Message {
  getType(): CompoundPreparation.CompoundPreparationTypeMap[keyof CompoundPreparation.CompoundPreparationTypeMap];
  setType(value: CompoundPreparation.CompoundPreparationTypeMap[keyof CompoundPreparation.CompoundPreparationTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  getReactionId(): string;
  setReactionId(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): CompoundPreparation.AsObject;
  static toObject(includeInstance: boolean, msg: CompoundPreparation): CompoundPreparation.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: CompoundPreparation, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): CompoundPreparation;
  static deserializeBinaryFromReader(message: CompoundPreparation, reader: jspb.BinaryReader): CompoundPreparation;
}

export namespace CompoundPreparation {
  export type AsObject = {
    type: CompoundPreparation.CompoundPreparationTypeMap[keyof CompoundPreparation.CompoundPreparationTypeMap],
    details: string,
    reactionId: string,
  }

  export interface CompoundPreparationTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    NONE: 2;
    REPURIFIED: 3;
    SPARGED: 4;
    DRIED: 5;
    SYNTHESIZED: 6;
  }

  export const CompoundPreparationType: CompoundPreparationTypeMap;
}

export class CompoundIdentifier extends jspb.Message {
  getType(): CompoundIdentifier.CompoundIdentifierTypeMap[keyof CompoundIdentifier.CompoundIdentifierTypeMap];
  setType(value: CompoundIdentifier.CompoundIdentifierTypeMap[keyof CompoundIdentifier.CompoundIdentifierTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  getValue(): string;
  setValue(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): CompoundIdentifier.AsObject;
  static toObject(includeInstance: boolean, msg: CompoundIdentifier): CompoundIdentifier.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: CompoundIdentifier, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): CompoundIdentifier;
  static deserializeBinaryFromReader(message: CompoundIdentifier, reader: jspb.BinaryReader): CompoundIdentifier;
}

export namespace CompoundIdentifier {
  export type AsObject = {
    type: CompoundIdentifier.CompoundIdentifierTypeMap[keyof CompoundIdentifier.CompoundIdentifierTypeMap],
    details: string,
    value: string,
  }

  export interface CompoundIdentifierTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    SMILES: 2;
    INCHI: 3;
    MOLBLOCK: 4;
    IUPAC_NAME: 5;
    NAME: 6;
    CAS_NUMBER: 7;
    PUBCHEM_CID: 8;
    CHEMSPIDER_ID: 9;
    CXSMILES: 10;
    INCHI_KEY: 11;
    XYZ: 12;
    UNIPROT_ID: 13;
    PDB_ID: 14;
    AMINO_ACID_SEQUENCE: 15;
    HELM: 16;
    MDL: 17;
  }

  export const CompoundIdentifierType: CompoundIdentifierTypeMap;
}

export class Vessel extends jspb.Message {
  getType(): Vessel.VesselTypeMap[keyof Vessel.VesselTypeMap];
  setType(value: Vessel.VesselTypeMap[keyof Vessel.VesselTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  hasMaterial(): boolean;
  clearMaterial(): void;
  getMaterial(): VesselMaterial | undefined;
  setMaterial(value?: VesselMaterial): void;

  clearPreparationsList(): void;
  getPreparationsList(): Array<VesselPreparation>;
  setPreparationsList(value: Array<VesselPreparation>): void;
  addPreparations(value?: VesselPreparation, index?: number): VesselPreparation;

  clearAttachmentsList(): void;
  getAttachmentsList(): Array<VesselAttachment>;
  setAttachmentsList(value: Array<VesselAttachment>): void;
  addAttachments(value?: VesselAttachment, index?: number): VesselAttachment;

  hasVolume(): boolean;
  clearVolume(): void;
  getVolume(): Volume | undefined;
  setVolume(value?: Volume): void;

  getVesselId(): string;
  setVesselId(value: string): void;

  getPosition(): string;
  setPosition(value: string): void;

  getRow(): string;
  setRow(value: string): void;

  getCol(): string;
  setCol(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Vessel.AsObject;
  static toObject(includeInstance: boolean, msg: Vessel): Vessel.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Vessel, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Vessel;
  static deserializeBinaryFromReader(message: Vessel, reader: jspb.BinaryReader): Vessel;
}

export namespace Vessel {
  export type AsObject = {
    type: Vessel.VesselTypeMap[keyof Vessel.VesselTypeMap],
    details: string,
    material?: VesselMaterial.AsObject,
    preparationsList: Array<VesselPreparation.AsObject>,
    attachmentsList: Array<VesselAttachment.AsObject>,
    volume?: Volume.AsObject,
    vesselId: string,
    position: string,
    row: string,
    col: string,
  }

  export interface VesselTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    ROUND_BOTTOM_FLASK: 2;
    VIAL: 3;
    WELL_PLATE: 4;
    MICROWAVE_VIAL: 5;
    TUBE: 6;
    CONTINUOUS_STIRRED_TANK_REACTOR: 7;
    PACKED_BED_REACTOR: 8;
    NMR_TUBE: 9;
    PRESSURE_FLASK: 10;
    PRESSURE_REACTOR: 11;
    ELECTROCHEMICAL_CELL: 12;
  }

  export const VesselType: VesselTypeMap;
}

export class VesselMaterial extends jspb.Message {
  getType(): VesselMaterial.VesselMaterialTypeMap[keyof VesselMaterial.VesselMaterialTypeMap];
  setType(value: VesselMaterial.VesselMaterialTypeMap[keyof VesselMaterial.VesselMaterialTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): VesselMaterial.AsObject;
  static toObject(includeInstance: boolean, msg: VesselMaterial): VesselMaterial.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: VesselMaterial, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): VesselMaterial;
  static deserializeBinaryFromReader(message: VesselMaterial, reader: jspb.BinaryReader): VesselMaterial;
}

export namespace VesselMaterial {
  export type AsObject = {
    type: VesselMaterial.VesselMaterialTypeMap[keyof VesselMaterial.VesselMaterialTypeMap],
    details: string,
  }

  export interface VesselMaterialTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    GLASS: 2;
    POLYPROPYLENE: 3;
    PLASTIC: 4;
    METAL: 5;
    QUARTZ: 6;
  }

  export const VesselMaterialType: VesselMaterialTypeMap;
}

export class VesselAttachment extends jspb.Message {
  getType(): VesselAttachment.VesselAttachmentTypeMap[keyof VesselAttachment.VesselAttachmentTypeMap];
  setType(value: VesselAttachment.VesselAttachmentTypeMap[keyof VesselAttachment.VesselAttachmentTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): VesselAttachment.AsObject;
  static toObject(includeInstance: boolean, msg: VesselAttachment): VesselAttachment.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: VesselAttachment, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): VesselAttachment;
  static deserializeBinaryFromReader(message: VesselAttachment, reader: jspb.BinaryReader): VesselAttachment;
}

export namespace VesselAttachment {
  export type AsObject = {
    type: VesselAttachment.VesselAttachmentTypeMap[keyof VesselAttachment.VesselAttachmentTypeMap],
    details: string,
  }

  export interface VesselAttachmentTypeMap {
    UNSPECIFIED: 0;
    NONE: 1;
    CUSTOM: 2;
    SEPTUM: 3;
    CAP: 4;
    MAT: 5;
    REFLUX_CONDENSER: 6;
    VENT_NEEDLE: 7;
    DEAN_STARK: 8;
    VACUUM_TUBE: 9;
    ADDITION_FUNNEL: 10;
    DRYING_TUBE: 11;
    ALUMINUM_FOIL: 12;
    THERMOCOUPLE: 13;
    BALLOON: 14;
    GAS_ADAPTER: 15;
    PRESSURE_REGULATOR: 16;
    RELEASE_VALVE: 17;
  }

  export const VesselAttachmentType: VesselAttachmentTypeMap;
}

export class VesselPreparation extends jspb.Message {
  getType(): VesselPreparation.VesselPreparationTypeMap[keyof VesselPreparation.VesselPreparationTypeMap];
  setType(value: VesselPreparation.VesselPreparationTypeMap[keyof VesselPreparation.VesselPreparationTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): VesselPreparation.AsObject;
  static toObject(includeInstance: boolean, msg: VesselPreparation): VesselPreparation.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: VesselPreparation, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): VesselPreparation;
  static deserializeBinaryFromReader(message: VesselPreparation, reader: jspb.BinaryReader): VesselPreparation;
}

export namespace VesselPreparation {
  export type AsObject = {
    type: VesselPreparation.VesselPreparationTypeMap[keyof VesselPreparation.VesselPreparationTypeMap],
    details: string,
  }

  export interface VesselPreparationTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    NONE: 2;
    OVEN_DRIED: 3;
    FLAME_DRIED: 4;
    EVACUATED_BACKFILLED: 5;
    PURGED: 6;
  }

  export const VesselPreparationType: VesselPreparationTypeMap;
}

export class ReactionSetup extends jspb.Message {
  hasVessel(): boolean;
  clearVessel(): void;
  getVessel(): Vessel | undefined;
  setVessel(value?: Vessel): void;

  hasIsAutomated(): boolean;
  clearIsAutomated(): void;
  getIsAutomated(): boolean;
  setIsAutomated(value: boolean): void;

  getAutomationPlatform(): string;
  setAutomationPlatform(value: string): void;

  getAutomationCodeMap(): jspb.Map<string, Data>;
  clearAutomationCodeMap(): void;
  hasEnvironment(): boolean;
  clearEnvironment(): void;
  getEnvironment(): ReactionSetup.ReactionEnvironment | undefined;
  setEnvironment(value?: ReactionSetup.ReactionEnvironment): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionSetup.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionSetup): ReactionSetup.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionSetup, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionSetup;
  static deserializeBinaryFromReader(message: ReactionSetup, reader: jspb.BinaryReader): ReactionSetup;
}

export namespace ReactionSetup {
  export type AsObject = {
    vessel?: Vessel.AsObject,
    isAutomated: boolean,
    automationPlatform: string,
    automationCodeMap: Array<[string, Data.AsObject]>,
    environment?: ReactionSetup.ReactionEnvironment.AsObject,
  }

  export class ReactionEnvironment extends jspb.Message {
    getType(): ReactionSetup.ReactionEnvironment.ReactionEnvironmentTypeMap[keyof ReactionSetup.ReactionEnvironment.ReactionEnvironmentTypeMap];
    setType(value: ReactionSetup.ReactionEnvironment.ReactionEnvironmentTypeMap[keyof ReactionSetup.ReactionEnvironment.ReactionEnvironmentTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): ReactionEnvironment.AsObject;
    static toObject(includeInstance: boolean, msg: ReactionEnvironment): ReactionEnvironment.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: ReactionEnvironment, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): ReactionEnvironment;
    static deserializeBinaryFromReader(message: ReactionEnvironment, reader: jspb.BinaryReader): ReactionEnvironment;
  }

  export namespace ReactionEnvironment {
    export type AsObject = {
      type: ReactionSetup.ReactionEnvironment.ReactionEnvironmentTypeMap[keyof ReactionSetup.ReactionEnvironment.ReactionEnvironmentTypeMap],
      details: string,
    }

    export interface ReactionEnvironmentTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      FUME_HOOD: 2;
      BENCH_TOP: 3;
      GLOVE_BOX: 4;
      GLOVE_BAG: 5;
    }

    export const ReactionEnvironmentType: ReactionEnvironmentTypeMap;
  }
}

export class ReactionConditions extends jspb.Message {
  hasTemperature(): boolean;
  clearTemperature(): void;
  getTemperature(): TemperatureConditions | undefined;
  setTemperature(value?: TemperatureConditions): void;

  hasPressure(): boolean;
  clearPressure(): void;
  getPressure(): PressureConditions | undefined;
  setPressure(value?: PressureConditions): void;

  hasStirring(): boolean;
  clearStirring(): void;
  getStirring(): StirringConditions | undefined;
  setStirring(value?: StirringConditions): void;

  hasIllumination(): boolean;
  clearIllumination(): void;
  getIllumination(): IlluminationConditions | undefined;
  setIllumination(value?: IlluminationConditions): void;

  hasElectrochemistry(): boolean;
  clearElectrochemistry(): void;
  getElectrochemistry(): ElectrochemistryConditions | undefined;
  setElectrochemistry(value?: ElectrochemistryConditions): void;

  hasFlow(): boolean;
  clearFlow(): void;
  getFlow(): FlowConditions | undefined;
  setFlow(value?: FlowConditions): void;

  hasReflux(): boolean;
  clearReflux(): void;
  getReflux(): boolean;
  setReflux(value: boolean): void;

  hasPh(): boolean;
  clearPh(): void;
  getPh(): number;
  setPh(value: number): void;

  hasConditionsAreDynamic(): boolean;
  clearConditionsAreDynamic(): void;
  getConditionsAreDynamic(): boolean;
  setConditionsAreDynamic(value: boolean): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionConditions.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionConditions): ReactionConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionConditions;
  static deserializeBinaryFromReader(message: ReactionConditions, reader: jspb.BinaryReader): ReactionConditions;
}

export namespace ReactionConditions {
  export type AsObject = {
    temperature?: TemperatureConditions.AsObject,
    pressure?: PressureConditions.AsObject,
    stirring?: StirringConditions.AsObject,
    illumination?: IlluminationConditions.AsObject,
    electrochemistry?: ElectrochemistryConditions.AsObject,
    flow?: FlowConditions.AsObject,
    reflux: boolean,
    ph: number,
    conditionsAreDynamic: boolean,
    details: string,
  }
}

export class TemperatureConditions extends jspb.Message {
  hasControl(): boolean;
  clearControl(): void;
  getControl(): TemperatureConditions.TemperatureControl | undefined;
  setControl(value?: TemperatureConditions.TemperatureControl): void;

  hasSetpoint(): boolean;
  clearSetpoint(): void;
  getSetpoint(): Temperature | undefined;
  setSetpoint(value?: Temperature): void;

  clearMeasurementsList(): void;
  getMeasurementsList(): Array<TemperatureConditions.TemperatureMeasurement>;
  setMeasurementsList(value: Array<TemperatureConditions.TemperatureMeasurement>): void;
  addMeasurements(value?: TemperatureConditions.TemperatureMeasurement, index?: number): TemperatureConditions.TemperatureMeasurement;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): TemperatureConditions.AsObject;
  static toObject(includeInstance: boolean, msg: TemperatureConditions): TemperatureConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: TemperatureConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): TemperatureConditions;
  static deserializeBinaryFromReader(message: TemperatureConditions, reader: jspb.BinaryReader): TemperatureConditions;
}

export namespace TemperatureConditions {
  export type AsObject = {
    control?: TemperatureConditions.TemperatureControl.AsObject,
    setpoint?: Temperature.AsObject,
    measurementsList: Array<TemperatureConditions.TemperatureMeasurement.AsObject>,
  }

  export class TemperatureControl extends jspb.Message {
    getType(): TemperatureConditions.TemperatureControl.TemperatureControlTypeMap[keyof TemperatureConditions.TemperatureControl.TemperatureControlTypeMap];
    setType(value: TemperatureConditions.TemperatureControl.TemperatureControlTypeMap[keyof TemperatureConditions.TemperatureControl.TemperatureControlTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): TemperatureControl.AsObject;
    static toObject(includeInstance: boolean, msg: TemperatureControl): TemperatureControl.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: TemperatureControl, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): TemperatureControl;
    static deserializeBinaryFromReader(message: TemperatureControl, reader: jspb.BinaryReader): TemperatureControl;
  }

  export namespace TemperatureControl {
    export type AsObject = {
      type: TemperatureConditions.TemperatureControl.TemperatureControlTypeMap[keyof TemperatureConditions.TemperatureControl.TemperatureControlTypeMap],
      details: string,
    }

    export interface TemperatureControlTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      AMBIENT: 2;
      OIL_BATH: 3;
      WATER_BATH: 4;
      SAND_BATH: 5;
      ICE_BATH: 6;
      DRY_ALUMINUM_PLATE: 7;
      MICROWAVE: 8;
      DRY_ICE_BATH: 9;
      AIR_FAN: 10;
      LIQUID_NITROGEN: 11;
    }

    export const TemperatureControlType: TemperatureControlTypeMap;
  }

  export class TemperatureMeasurement extends jspb.Message {
    getType(): TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementTypeMap[keyof TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementTypeMap];
    setType(value: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementTypeMap[keyof TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    hasTime(): boolean;
    clearTime(): void;
    getTime(): Time | undefined;
    setTime(value?: Time): void;

    hasTemperature(): boolean;
    clearTemperature(): void;
    getTemperature(): Temperature | undefined;
    setTemperature(value?: Temperature): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): TemperatureMeasurement.AsObject;
    static toObject(includeInstance: boolean, msg: TemperatureMeasurement): TemperatureMeasurement.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: TemperatureMeasurement, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): TemperatureMeasurement;
    static deserializeBinaryFromReader(message: TemperatureMeasurement, reader: jspb.BinaryReader): TemperatureMeasurement;
  }

  export namespace TemperatureMeasurement {
    export type AsObject = {
      type: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementTypeMap[keyof TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementTypeMap],
      details: string,
      time?: Time.AsObject,
      temperature?: Temperature.AsObject,
    }

    export interface TemperatureMeasurementTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      THERMOCOUPLE_INTERNAL: 2;
      THERMOCOUPLE_EXTERNAL: 3;
      INFRARED: 4;
    }

    export const TemperatureMeasurementType: TemperatureMeasurementTypeMap;
  }
}

export class PressureConditions extends jspb.Message {
  hasControl(): boolean;
  clearControl(): void;
  getControl(): PressureConditions.PressureControl | undefined;
  setControl(value?: PressureConditions.PressureControl): void;

  hasSetpoint(): boolean;
  clearSetpoint(): void;
  getSetpoint(): Pressure | undefined;
  setSetpoint(value?: Pressure): void;

  hasAtmosphere(): boolean;
  clearAtmosphere(): void;
  getAtmosphere(): PressureConditions.Atmosphere | undefined;
  setAtmosphere(value?: PressureConditions.Atmosphere): void;

  clearMeasurementsList(): void;
  getMeasurementsList(): Array<PressureConditions.PressureMeasurement>;
  setMeasurementsList(value: Array<PressureConditions.PressureMeasurement>): void;
  addMeasurements(value?: PressureConditions.PressureMeasurement, index?: number): PressureConditions.PressureMeasurement;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): PressureConditions.AsObject;
  static toObject(includeInstance: boolean, msg: PressureConditions): PressureConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: PressureConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): PressureConditions;
  static deserializeBinaryFromReader(message: PressureConditions, reader: jspb.BinaryReader): PressureConditions;
}

export namespace PressureConditions {
  export type AsObject = {
    control?: PressureConditions.PressureControl.AsObject,
    setpoint?: Pressure.AsObject,
    atmosphere?: PressureConditions.Atmosphere.AsObject,
    measurementsList: Array<PressureConditions.PressureMeasurement.AsObject>,
  }

  export class PressureControl extends jspb.Message {
    getType(): PressureConditions.PressureControl.PressureControlTypeMap[keyof PressureConditions.PressureControl.PressureControlTypeMap];
    setType(value: PressureConditions.PressureControl.PressureControlTypeMap[keyof PressureConditions.PressureControl.PressureControlTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): PressureControl.AsObject;
    static toObject(includeInstance: boolean, msg: PressureControl): PressureControl.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: PressureControl, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): PressureControl;
    static deserializeBinaryFromReader(message: PressureControl, reader: jspb.BinaryReader): PressureControl;
  }

  export namespace PressureControl {
    export type AsObject = {
      type: PressureConditions.PressureControl.PressureControlTypeMap[keyof PressureConditions.PressureControl.PressureControlTypeMap],
      details: string,
    }

    export interface PressureControlTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      AMBIENT: 2;
      SLIGHT_POSITIVE: 3;
      SEALED: 4;
      PRESSURIZED: 5;
    }

    export const PressureControlType: PressureControlTypeMap;
  }

  export class Atmosphere extends jspb.Message {
    getType(): PressureConditions.Atmosphere.AtmosphereTypeMap[keyof PressureConditions.Atmosphere.AtmosphereTypeMap];
    setType(value: PressureConditions.Atmosphere.AtmosphereTypeMap[keyof PressureConditions.Atmosphere.AtmosphereTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Atmosphere.AsObject;
    static toObject(includeInstance: boolean, msg: Atmosphere): Atmosphere.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Atmosphere, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Atmosphere;
    static deserializeBinaryFromReader(message: Atmosphere, reader: jspb.BinaryReader): Atmosphere;
  }

  export namespace Atmosphere {
    export type AsObject = {
      type: PressureConditions.Atmosphere.AtmosphereTypeMap[keyof PressureConditions.Atmosphere.AtmosphereTypeMap],
      details: string,
    }

    export interface AtmosphereTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      AIR: 2;
      NITROGEN: 3;
      ARGON: 4;
      OXYGEN: 5;
      HYDROGEN: 6;
      CARBON_MONOXIDE: 7;
      CARBON_DIOXIDE: 8;
      METHANE: 9;
      AMMONIA: 10;
      OZONE: 11;
      ETHYLENE: 12;
      ACETYLENE: 13;
    }

    export const AtmosphereType: AtmosphereTypeMap;
  }

  export class PressureMeasurement extends jspb.Message {
    getType(): PressureConditions.PressureMeasurement.PressureMeasurementTypeMap[keyof PressureConditions.PressureMeasurement.PressureMeasurementTypeMap];
    setType(value: PressureConditions.PressureMeasurement.PressureMeasurementTypeMap[keyof PressureConditions.PressureMeasurement.PressureMeasurementTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    hasTime(): boolean;
    clearTime(): void;
    getTime(): Time | undefined;
    setTime(value?: Time): void;

    hasPressure(): boolean;
    clearPressure(): void;
    getPressure(): Pressure | undefined;
    setPressure(value?: Pressure): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): PressureMeasurement.AsObject;
    static toObject(includeInstance: boolean, msg: PressureMeasurement): PressureMeasurement.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: PressureMeasurement, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): PressureMeasurement;
    static deserializeBinaryFromReader(message: PressureMeasurement, reader: jspb.BinaryReader): PressureMeasurement;
  }

  export namespace PressureMeasurement {
    export type AsObject = {
      type: PressureConditions.PressureMeasurement.PressureMeasurementTypeMap[keyof PressureConditions.PressureMeasurement.PressureMeasurementTypeMap],
      details: string,
      time?: Time.AsObject,
      pressure?: Pressure.AsObject,
    }

    export interface PressureMeasurementTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      PRESSURE_TRANSDUCER: 2;
    }

    export const PressureMeasurementType: PressureMeasurementTypeMap;
  }
}

export class StirringConditions extends jspb.Message {
  getType(): StirringConditions.StirringMethodTypeMap[keyof StirringConditions.StirringMethodTypeMap];
  setType(value: StirringConditions.StirringMethodTypeMap[keyof StirringConditions.StirringMethodTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  hasRate(): boolean;
  clearRate(): void;
  getRate(): StirringConditions.StirringRate | undefined;
  setRate(value?: StirringConditions.StirringRate): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): StirringConditions.AsObject;
  static toObject(includeInstance: boolean, msg: StirringConditions): StirringConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: StirringConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): StirringConditions;
  static deserializeBinaryFromReader(message: StirringConditions, reader: jspb.BinaryReader): StirringConditions;
}

export namespace StirringConditions {
  export type AsObject = {
    type: StirringConditions.StirringMethodTypeMap[keyof StirringConditions.StirringMethodTypeMap],
    details: string,
    rate?: StirringConditions.StirringRate.AsObject,
  }

  export class StirringRate extends jspb.Message {
    getType(): StirringConditions.StirringRate.StirringRateTypeMap[keyof StirringConditions.StirringRate.StirringRateTypeMap];
    setType(value: StirringConditions.StirringRate.StirringRateTypeMap[keyof StirringConditions.StirringRate.StirringRateTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    getRpm(): number;
    setRpm(value: number): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): StirringRate.AsObject;
    static toObject(includeInstance: boolean, msg: StirringRate): StirringRate.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: StirringRate, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): StirringRate;
    static deserializeBinaryFromReader(message: StirringRate, reader: jspb.BinaryReader): StirringRate;
  }

  export namespace StirringRate {
    export type AsObject = {
      type: StirringConditions.StirringRate.StirringRateTypeMap[keyof StirringConditions.StirringRate.StirringRateTypeMap],
      details: string,
      rpm: number,
    }

    export interface StirringRateTypeMap {
      UNSPECIFIED: 0;
      HIGH: 1;
      MEDIUM: 2;
      LOW: 3;
    }

    export const StirringRateType: StirringRateTypeMap;
  }

  export interface StirringMethodTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    NONE: 2;
    STIR_BAR: 3;
    OVERHEAD_MIXER: 4;
    AGITATION: 5;
    BALL_MILLING: 6;
    SONICATION: 7;
  }

  export const StirringMethodType: StirringMethodTypeMap;
}

export class IlluminationConditions extends jspb.Message {
  getType(): IlluminationConditions.IlluminationTypeMap[keyof IlluminationConditions.IlluminationTypeMap];
  setType(value: IlluminationConditions.IlluminationTypeMap[keyof IlluminationConditions.IlluminationTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  hasPeakWavelength(): boolean;
  clearPeakWavelength(): void;
  getPeakWavelength(): Wavelength | undefined;
  setPeakWavelength(value?: Wavelength): void;

  getColor(): string;
  setColor(value: string): void;

  hasDistanceToVessel(): boolean;
  clearDistanceToVessel(): void;
  getDistanceToVessel(): Length | undefined;
  setDistanceToVessel(value?: Length): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): IlluminationConditions.AsObject;
  static toObject(includeInstance: boolean, msg: IlluminationConditions): IlluminationConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: IlluminationConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): IlluminationConditions;
  static deserializeBinaryFromReader(message: IlluminationConditions, reader: jspb.BinaryReader): IlluminationConditions;
}

export namespace IlluminationConditions {
  export type AsObject = {
    type: IlluminationConditions.IlluminationTypeMap[keyof IlluminationConditions.IlluminationTypeMap],
    details: string,
    peakWavelength?: Wavelength.AsObject,
    color: string,
    distanceToVessel?: Length.AsObject,
  }

  export interface IlluminationTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    AMBIENT: 2;
    DARK: 3;
    LED: 4;
    HALOGEN_LAMP: 5;
    DEUTERIUM_LAMP: 6;
    SOLAR_SIMULATOR: 7;
    BROAD_SPECTRUM: 8;
  }

  export const IlluminationType: IlluminationTypeMap;
}

export class ElectrochemistryConditions extends jspb.Message {
  getType(): ElectrochemistryConditions.ElectrochemistryTypeMap[keyof ElectrochemistryConditions.ElectrochemistryTypeMap];
  setType(value: ElectrochemistryConditions.ElectrochemistryTypeMap[keyof ElectrochemistryConditions.ElectrochemistryTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  hasCurrent(): boolean;
  clearCurrent(): void;
  getCurrent(): Current | undefined;
  setCurrent(value?: Current): void;

  hasVoltage(): boolean;
  clearVoltage(): void;
  getVoltage(): Voltage | undefined;
  setVoltage(value?: Voltage): void;

  getAnodeMaterial(): string;
  setAnodeMaterial(value: string): void;

  getCathodeMaterial(): string;
  setCathodeMaterial(value: string): void;

  hasElectrodeSeparation(): boolean;
  clearElectrodeSeparation(): void;
  getElectrodeSeparation(): Length | undefined;
  setElectrodeSeparation(value?: Length): void;

  clearMeasurementsList(): void;
  getMeasurementsList(): Array<ElectrochemistryConditions.ElectrochemistryMeasurement>;
  setMeasurementsList(value: Array<ElectrochemistryConditions.ElectrochemistryMeasurement>): void;
  addMeasurements(value?: ElectrochemistryConditions.ElectrochemistryMeasurement, index?: number): ElectrochemistryConditions.ElectrochemistryMeasurement;

  hasCell(): boolean;
  clearCell(): void;
  getCell(): ElectrochemistryConditions.ElectrochemistryCell | undefined;
  setCell(value?: ElectrochemistryConditions.ElectrochemistryCell): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ElectrochemistryConditions.AsObject;
  static toObject(includeInstance: boolean, msg: ElectrochemistryConditions): ElectrochemistryConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ElectrochemistryConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ElectrochemistryConditions;
  static deserializeBinaryFromReader(message: ElectrochemistryConditions, reader: jspb.BinaryReader): ElectrochemistryConditions;
}

export namespace ElectrochemistryConditions {
  export type AsObject = {
    type: ElectrochemistryConditions.ElectrochemistryTypeMap[keyof ElectrochemistryConditions.ElectrochemistryTypeMap],
    details: string,
    current?: Current.AsObject,
    voltage?: Voltage.AsObject,
    anodeMaterial: string,
    cathodeMaterial: string,
    electrodeSeparation?: Length.AsObject,
    measurementsList: Array<ElectrochemistryConditions.ElectrochemistryMeasurement.AsObject>,
    cell?: ElectrochemistryConditions.ElectrochemistryCell.AsObject,
  }

  export class ElectrochemistryMeasurement extends jspb.Message {
    hasTime(): boolean;
    clearTime(): void;
    getTime(): Time | undefined;
    setTime(value?: Time): void;

    hasCurrent(): boolean;
    clearCurrent(): void;
    getCurrent(): Current | undefined;
    setCurrent(value?: Current): void;

    hasVoltage(): boolean;
    clearVoltage(): void;
    getVoltage(): Voltage | undefined;
    setVoltage(value?: Voltage): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): ElectrochemistryMeasurement.AsObject;
    static toObject(includeInstance: boolean, msg: ElectrochemistryMeasurement): ElectrochemistryMeasurement.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: ElectrochemistryMeasurement, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): ElectrochemistryMeasurement;
    static deserializeBinaryFromReader(message: ElectrochemistryMeasurement, reader: jspb.BinaryReader): ElectrochemistryMeasurement;
  }

  export namespace ElectrochemistryMeasurement {
    export type AsObject = {
      time?: Time.AsObject,
      current?: Current.AsObject,
      voltage?: Voltage.AsObject,
    }
  }

  export class ElectrochemistryCell extends jspb.Message {
    getType(): ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellTypeMap[keyof ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellTypeMap];
    setType(value: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellTypeMap[keyof ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): ElectrochemistryCell.AsObject;
    static toObject(includeInstance: boolean, msg: ElectrochemistryCell): ElectrochemistryCell.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: ElectrochemistryCell, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): ElectrochemistryCell;
    static deserializeBinaryFromReader(message: ElectrochemistryCell, reader: jspb.BinaryReader): ElectrochemistryCell;
  }

  export namespace ElectrochemistryCell {
    export type AsObject = {
      type: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellTypeMap[keyof ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellTypeMap],
      details: string,
    }

    export interface ElectrochemistryCellTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      DIVIDED_CELL: 2;
      UNDIVIDED_CELL: 3;
    }

    export const ElectrochemistryCellType: ElectrochemistryCellTypeMap;
  }

  export interface ElectrochemistryTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    CONSTANT_CURRENT: 2;
    CONSTANT_VOLTAGE: 3;
  }

  export const ElectrochemistryType: ElectrochemistryTypeMap;
}

export class FlowConditions extends jspb.Message {
  getType(): FlowConditions.FlowTypeMap[keyof FlowConditions.FlowTypeMap];
  setType(value: FlowConditions.FlowTypeMap[keyof FlowConditions.FlowTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  getPumpType(): string;
  setPumpType(value: string): void;

  hasTubing(): boolean;
  clearTubing(): void;
  getTubing(): FlowConditions.Tubing | undefined;
  setTubing(value?: FlowConditions.Tubing): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): FlowConditions.AsObject;
  static toObject(includeInstance: boolean, msg: FlowConditions): FlowConditions.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: FlowConditions, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): FlowConditions;
  static deserializeBinaryFromReader(message: FlowConditions, reader: jspb.BinaryReader): FlowConditions;
}

export namespace FlowConditions {
  export type AsObject = {
    type: FlowConditions.FlowTypeMap[keyof FlowConditions.FlowTypeMap],
    details: string,
    pumpType: string,
    tubing?: FlowConditions.Tubing.AsObject,
  }

  export class Tubing extends jspb.Message {
    getType(): FlowConditions.Tubing.TubingTypeMap[keyof FlowConditions.Tubing.TubingTypeMap];
    setType(value: FlowConditions.Tubing.TubingTypeMap[keyof FlowConditions.Tubing.TubingTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    hasDiameter(): boolean;
    clearDiameter(): void;
    getDiameter(): Length | undefined;
    setDiameter(value?: Length): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Tubing.AsObject;
    static toObject(includeInstance: boolean, msg: Tubing): Tubing.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Tubing, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Tubing;
    static deserializeBinaryFromReader(message: Tubing, reader: jspb.BinaryReader): Tubing;
  }

  export namespace Tubing {
    export type AsObject = {
      type: FlowConditions.Tubing.TubingTypeMap[keyof FlowConditions.Tubing.TubingTypeMap],
      details: string,
      diameter?: Length.AsObject,
    }

    export interface TubingTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      STEEL: 2;
      COPPER: 3;
      PFA: 4;
      FEP: 5;
      TEFLONAF: 6;
      PTFE: 7;
      GLASS: 8;
      QUARTZ: 9;
      SILICON: 10;
      PDMS: 11;
    }

    export const TubingType: TubingTypeMap;
  }

  export interface FlowTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    PLUG_FLOW_REACTOR: 2;
    CONTINUOUS_STIRRED_TANK_REACTOR: 3;
    PACKED_BED_REACTOR: 4;
  }

  export const FlowType: FlowTypeMap;
}

export class ReactionNotes extends jspb.Message {
  hasIsHeterogeneous(): boolean;
  clearIsHeterogeneous(): void;
  getIsHeterogeneous(): boolean;
  setIsHeterogeneous(value: boolean): void;

  hasFormsPrecipitate(): boolean;
  clearFormsPrecipitate(): void;
  getFormsPrecipitate(): boolean;
  setFormsPrecipitate(value: boolean): void;

  hasIsExothermic(): boolean;
  clearIsExothermic(): void;
  getIsExothermic(): boolean;
  setIsExothermic(value: boolean): void;

  hasOffgasses(): boolean;
  clearOffgasses(): void;
  getOffgasses(): boolean;
  setOffgasses(value: boolean): void;

  hasIsSensitiveToMoisture(): boolean;
  clearIsSensitiveToMoisture(): void;
  getIsSensitiveToMoisture(): boolean;
  setIsSensitiveToMoisture(value: boolean): void;

  hasIsSensitiveToOxygen(): boolean;
  clearIsSensitiveToOxygen(): void;
  getIsSensitiveToOxygen(): boolean;
  setIsSensitiveToOxygen(value: boolean): void;

  hasIsSensitiveToLight(): boolean;
  clearIsSensitiveToLight(): void;
  getIsSensitiveToLight(): boolean;
  setIsSensitiveToLight(value: boolean): void;

  getSafetyNotes(): string;
  setSafetyNotes(value: string): void;

  getProcedureDetails(): string;
  setProcedureDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionNotes.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionNotes): ReactionNotes.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionNotes, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionNotes;
  static deserializeBinaryFromReader(message: ReactionNotes, reader: jspb.BinaryReader): ReactionNotes;
}

export namespace ReactionNotes {
  export type AsObject = {
    isHeterogeneous: boolean,
    formsPrecipitate: boolean,
    isExothermic: boolean,
    offgasses: boolean,
    isSensitiveToMoisture: boolean,
    isSensitiveToOxygen: boolean,
    isSensitiveToLight: boolean,
    safetyNotes: string,
    procedureDetails: string,
  }
}

export class ReactionObservation extends jspb.Message {
  hasTime(): boolean;
  clearTime(): void;
  getTime(): Time | undefined;
  setTime(value?: Time): void;

  getComment(): string;
  setComment(value: string): void;

  hasImage(): boolean;
  clearImage(): void;
  getImage(): Data | undefined;
  setImage(value?: Data): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionObservation.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionObservation): ReactionObservation.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionObservation, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionObservation;
  static deserializeBinaryFromReader(message: ReactionObservation, reader: jspb.BinaryReader): ReactionObservation;
}

export namespace ReactionObservation {
  export type AsObject = {
    time?: Time.AsObject,
    comment: string,
    image?: Data.AsObject,
  }
}

export class ReactionWorkup extends jspb.Message {
  getType(): ReactionWorkup.ReactionWorkupTypeMap[keyof ReactionWorkup.ReactionWorkupTypeMap];
  setType(value: ReactionWorkup.ReactionWorkupTypeMap[keyof ReactionWorkup.ReactionWorkupTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  hasDuration(): boolean;
  clearDuration(): void;
  getDuration(): Time | undefined;
  setDuration(value?: Time): void;

  hasInput(): boolean;
  clearInput(): void;
  getInput(): ReactionInput | undefined;
  setInput(value?: ReactionInput): void;

  hasAmount(): boolean;
  clearAmount(): void;
  getAmount(): Amount | undefined;
  setAmount(value?: Amount): void;

  hasTemperature(): boolean;
  clearTemperature(): void;
  getTemperature(): TemperatureConditions | undefined;
  setTemperature(value?: TemperatureConditions): void;

  getKeepPhase(): string;
  setKeepPhase(value: string): void;

  hasStirring(): boolean;
  clearStirring(): void;
  getStirring(): StirringConditions | undefined;
  setStirring(value?: StirringConditions): void;

  hasTargetPh(): boolean;
  clearTargetPh(): void;
  getTargetPh(): number;
  setTargetPh(value: number): void;

  hasIsAutomated(): boolean;
  clearIsAutomated(): void;
  getIsAutomated(): boolean;
  setIsAutomated(value: boolean): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionWorkup.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionWorkup): ReactionWorkup.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionWorkup, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionWorkup;
  static deserializeBinaryFromReader(message: ReactionWorkup, reader: jspb.BinaryReader): ReactionWorkup;
}

export namespace ReactionWorkup {
  export type AsObject = {
    type: ReactionWorkup.ReactionWorkupTypeMap[keyof ReactionWorkup.ReactionWorkupTypeMap],
    details: string,
    duration?: Time.AsObject,
    input?: ReactionInput.AsObject,
    amount?: Amount.AsObject,
    temperature?: TemperatureConditions.AsObject,
    keepPhase: string,
    stirring?: StirringConditions.AsObject,
    targetPh: number,
    isAutomated: boolean,
  }

  export interface ReactionWorkupTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    ADDITION: 2;
    ALIQUOT: 3;
    TEMPERATURE: 4;
    CONCENTRATION: 5;
    EXTRACTION: 6;
    FILTRATION: 7;
    WASH: 8;
    DRY_IN_VACUUM: 9;
    DRY_WITH_MATERIAL: 10;
    FLASH_CHROMATOGRAPHY: 11;
    OTHER_CHROMATOGRAPHY: 12;
    SCAVENGING: 13;
    WAIT: 14;
    STIRRING: 15;
    PH_ADJUST: 16;
    DISSOLUTION: 17;
    DISTILLATION: 18;
  }

  export const ReactionWorkupType: ReactionWorkupTypeMap;
}

export class ReactionOutcome extends jspb.Message {
  hasReactionTime(): boolean;
  clearReactionTime(): void;
  getReactionTime(): Time | undefined;
  setReactionTime(value?: Time): void;

  hasConversion(): boolean;
  clearConversion(): void;
  getConversion(): Percentage | undefined;
  setConversion(value?: Percentage): void;

  clearProductsList(): void;
  getProductsList(): Array<ProductCompound>;
  setProductsList(value: Array<ProductCompound>): void;
  addProducts(value?: ProductCompound, index?: number): ProductCompound;

  getAnalysesMap(): jspb.Map<string, Analysis>;
  clearAnalysesMap(): void;
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionOutcome.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionOutcome): ReactionOutcome.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionOutcome, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionOutcome;
  static deserializeBinaryFromReader(message: ReactionOutcome, reader: jspb.BinaryReader): ReactionOutcome;
}

export namespace ReactionOutcome {
  export type AsObject = {
    reactionTime?: Time.AsObject,
    conversion?: Percentage.AsObject,
    productsList: Array<ProductCompound.AsObject>,
    analysesMap: Array<[string, Analysis.AsObject]>,
  }
}

export class ProductCompound extends jspb.Message {
  clearIdentifiersList(): void;
  getIdentifiersList(): Array<CompoundIdentifier>;
  setIdentifiersList(value: Array<CompoundIdentifier>): void;
  addIdentifiers(value?: CompoundIdentifier, index?: number): CompoundIdentifier;

  hasIsDesiredProduct(): boolean;
  clearIsDesiredProduct(): void;
  getIsDesiredProduct(): boolean;
  setIsDesiredProduct(value: boolean): void;

  clearMeasurementsList(): void;
  getMeasurementsList(): Array<ProductMeasurement>;
  setMeasurementsList(value: Array<ProductMeasurement>): void;
  addMeasurements(value?: ProductMeasurement, index?: number): ProductMeasurement;

  getIsolatedColor(): string;
  setIsolatedColor(value: string): void;

  hasTexture(): boolean;
  clearTexture(): void;
  getTexture(): Texture | undefined;
  setTexture(value?: Texture): void;

  getFeaturesMap(): jspb.Map<string, Data>;
  clearFeaturesMap(): void;
  getReactionRole(): ReactionRole.ReactionRoleTypeMap[keyof ReactionRole.ReactionRoleTypeMap];
  setReactionRole(value: ReactionRole.ReactionRoleTypeMap[keyof ReactionRole.ReactionRoleTypeMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ProductCompound.AsObject;
  static toObject(includeInstance: boolean, msg: ProductCompound): ProductCompound.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ProductCompound, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ProductCompound;
  static deserializeBinaryFromReader(message: ProductCompound, reader: jspb.BinaryReader): ProductCompound;
}

export namespace ProductCompound {
  export type AsObject = {
    identifiersList: Array<CompoundIdentifier.AsObject>,
    isDesiredProduct: boolean,
    measurementsList: Array<ProductMeasurement.AsObject>,
    isolatedColor: string,
    texture?: Texture.AsObject,
    featuresMap: Array<[string, Data.AsObject]>,
    reactionRole: ReactionRole.ReactionRoleTypeMap[keyof ReactionRole.ReactionRoleTypeMap],
  }
}

export class ProductMeasurement extends jspb.Message {
  getAnalysisKey(): string;
  setAnalysisKey(value: string): void;

  getType(): ProductMeasurement.ProductMeasurementTypeMap[keyof ProductMeasurement.ProductMeasurementTypeMap];
  setType(value: ProductMeasurement.ProductMeasurementTypeMap[keyof ProductMeasurement.ProductMeasurementTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  hasUsesInternalStandard(): boolean;
  clearUsesInternalStandard(): void;
  getUsesInternalStandard(): boolean;
  setUsesInternalStandard(value: boolean): void;

  hasIsNormalized(): boolean;
  clearIsNormalized(): void;
  getIsNormalized(): boolean;
  setIsNormalized(value: boolean): void;

  hasUsesAuthenticStandard(): boolean;
  clearUsesAuthenticStandard(): void;
  getUsesAuthenticStandard(): boolean;
  setUsesAuthenticStandard(value: boolean): void;

  hasAuthenticStandard(): boolean;
  clearAuthenticStandard(): void;
  getAuthenticStandard(): Compound | undefined;
  setAuthenticStandard(value?: Compound): void;

  hasPercentage(): boolean;
  clearPercentage(): void;
  getPercentage(): Percentage | undefined;
  setPercentage(value?: Percentage): void;

  hasFloatValue(): boolean;
  clearFloatValue(): void;
  getFloatValue(): FloatValue | undefined;
  setFloatValue(value?: FloatValue): void;

  hasStringValue(): boolean;
  clearStringValue(): void;
  getStringValue(): string;
  setStringValue(value: string): void;

  hasAmount(): boolean;
  clearAmount(): void;
  getAmount(): Amount | undefined;
  setAmount(value?: Amount): void;

  hasRetentionTime(): boolean;
  clearRetentionTime(): void;
  getRetentionTime(): Time | undefined;
  setRetentionTime(value?: Time): void;

  hasMassSpecDetails(): boolean;
  clearMassSpecDetails(): void;
  getMassSpecDetails(): ProductMeasurement.MassSpecMeasurementDetails | undefined;
  setMassSpecDetails(value?: ProductMeasurement.MassSpecMeasurementDetails): void;

  hasSelectivity(): boolean;
  clearSelectivity(): void;
  getSelectivity(): ProductMeasurement.Selectivity | undefined;
  setSelectivity(value?: ProductMeasurement.Selectivity): void;

  hasWavelength(): boolean;
  clearWavelength(): void;
  getWavelength(): Wavelength | undefined;
  setWavelength(value?: Wavelength): void;

  getValueCase(): ProductMeasurement.ValueCase;
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ProductMeasurement.AsObject;
  static toObject(includeInstance: boolean, msg: ProductMeasurement): ProductMeasurement.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ProductMeasurement, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ProductMeasurement;
  static deserializeBinaryFromReader(message: ProductMeasurement, reader: jspb.BinaryReader): ProductMeasurement;
}

export namespace ProductMeasurement {
  export type AsObject = {
    analysisKey: string,
    type: ProductMeasurement.ProductMeasurementTypeMap[keyof ProductMeasurement.ProductMeasurementTypeMap],
    details: string,
    usesInternalStandard: boolean,
    isNormalized: boolean,
    usesAuthenticStandard: boolean,
    authenticStandard?: Compound.AsObject,
    percentage?: Percentage.AsObject,
    floatValue?: FloatValue.AsObject,
    stringValue: string,
    amount?: Amount.AsObject,
    retentionTime?: Time.AsObject,
    massSpecDetails?: ProductMeasurement.MassSpecMeasurementDetails.AsObject,
    selectivity?: ProductMeasurement.Selectivity.AsObject,
    wavelength?: Wavelength.AsObject,
  }

  export class MassSpecMeasurementDetails extends jspb.Message {
    getType(): ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementTypeMap[keyof ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementTypeMap];
    setType(value: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementTypeMap[keyof ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    hasTicMinimumMz(): boolean;
    clearTicMinimumMz(): void;
    getTicMinimumMz(): number;
    setTicMinimumMz(value: number): void;

    hasTicMaximumMz(): boolean;
    clearTicMaximumMz(): void;
    getTicMaximumMz(): number;
    setTicMaximumMz(value: number): void;

    clearEicMassesList(): void;
    getEicMassesList(): Array<number>;
    setEicMassesList(value: Array<number>): void;
    addEicMasses(value: number, index?: number): number;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): MassSpecMeasurementDetails.AsObject;
    static toObject(includeInstance: boolean, msg: MassSpecMeasurementDetails): MassSpecMeasurementDetails.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: MassSpecMeasurementDetails, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): MassSpecMeasurementDetails;
    static deserializeBinaryFromReader(message: MassSpecMeasurementDetails, reader: jspb.BinaryReader): MassSpecMeasurementDetails;
  }

  export namespace MassSpecMeasurementDetails {
    export type AsObject = {
      type: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementTypeMap[keyof ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementTypeMap],
      details: string,
      ticMinimumMz: number,
      ticMaximumMz: number,
      eicMassesList: Array<number>,
    }

    export interface MassSpecMeasurementTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      TIC: 2;
      TIC_POSITIVE: 3;
      TIC_NEGATIVE: 4;
      EIC: 5;
    }

    export const MassSpecMeasurementType: MassSpecMeasurementTypeMap;
  }

  export class Selectivity extends jspb.Message {
    getType(): ProductMeasurement.Selectivity.SelectivityTypeMap[keyof ProductMeasurement.Selectivity.SelectivityTypeMap];
    setType(value: ProductMeasurement.Selectivity.SelectivityTypeMap[keyof ProductMeasurement.Selectivity.SelectivityTypeMap]): void;

    getDetails(): string;
    setDetails(value: string): void;

    serializeBinary(): Uint8Array;
    toObject(includeInstance?: boolean): Selectivity.AsObject;
    static toObject(includeInstance: boolean, msg: Selectivity): Selectivity.AsObject;
    static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
    static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
    static serializeBinaryToWriter(message: Selectivity, writer: jspb.BinaryWriter): void;
    static deserializeBinary(bytes: Uint8Array): Selectivity;
    static deserializeBinaryFromReader(message: Selectivity, reader: jspb.BinaryReader): Selectivity;
  }

  export namespace Selectivity {
    export type AsObject = {
      type: ProductMeasurement.Selectivity.SelectivityTypeMap[keyof ProductMeasurement.Selectivity.SelectivityTypeMap],
      details: string,
    }

    export interface SelectivityTypeMap {
      UNSPECIFIED: 0;
      CUSTOM: 1;
      EE: 2;
      ER: 3;
      DR: 4;
      EZ: 5;
      ZE: 6;
    }

    export const SelectivityType: SelectivityTypeMap;
  }

  export interface ProductMeasurementTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    IDENTITY: 2;
    YIELD: 3;
    SELECTIVITY: 4;
    PURITY: 5;
    AREA: 6;
    COUNTS: 7;
    INTENSITY: 8;
    AMOUNT: 9;
  }

  export const ProductMeasurementType: ProductMeasurementTypeMap;

  export enum ValueCase {
    VALUE_NOT_SET = 0,
    PERCENTAGE = 8,
    FLOAT_VALUE = 9,
    STRING_VALUE = 10,
    AMOUNT = 11,
  }
}

export class DateTime extends jspb.Message {
  getValue(): string;
  setValue(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): DateTime.AsObject;
  static toObject(includeInstance: boolean, msg: DateTime): DateTime.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: DateTime, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): DateTime;
  static deserializeBinaryFromReader(message: DateTime, reader: jspb.BinaryReader): DateTime;
}

export namespace DateTime {
  export type AsObject = {
    value: string,
  }
}

export class Analysis extends jspb.Message {
  getType(): Analysis.AnalysisTypeMap[keyof Analysis.AnalysisTypeMap];
  setType(value: Analysis.AnalysisTypeMap[keyof Analysis.AnalysisTypeMap]): void;

  getDetails(): string;
  setDetails(value: string): void;

  getChmoId(): number;
  setChmoId(value: number): void;

  hasIsOfIsolatedSpecies(): boolean;
  clearIsOfIsolatedSpecies(): void;
  getIsOfIsolatedSpecies(): boolean;
  setIsOfIsolatedSpecies(value: boolean): void;

  getDataMap(): jspb.Map<string, Data>;
  clearDataMap(): void;
  getInstrumentManufacturer(): string;
  setInstrumentManufacturer(value: string): void;

  hasInstrumentLastCalibrated(): boolean;
  clearInstrumentLastCalibrated(): void;
  getInstrumentLastCalibrated(): DateTime | undefined;
  setInstrumentLastCalibrated(value?: DateTime): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Analysis.AsObject;
  static toObject(includeInstance: boolean, msg: Analysis): Analysis.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Analysis, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Analysis;
  static deserializeBinaryFromReader(message: Analysis, reader: jspb.BinaryReader): Analysis;
}

export namespace Analysis {
  export type AsObject = {
    type: Analysis.AnalysisTypeMap[keyof Analysis.AnalysisTypeMap],
    details: string,
    chmoId: number,
    isOfIsolatedSpecies: boolean,
    dataMap: Array<[string, Data.AsObject]>,
    instrumentManufacturer: string,
    instrumentLastCalibrated?: DateTime.AsObject,
  }

  export interface AnalysisTypeMap {
    UNSPECIFIED: 0;
    CUSTOM: 1;
    LC: 2;
    GC: 3;
    IR: 4;
    NMR_1H: 5;
    NMR_13C: 6;
    NMR_OTHER: 7;
    MP: 8;
    UV: 9;
    TLC: 10;
    MS: 11;
    HRMS: 12;
    MSMS: 13;
    WEIGHT: 14;
    LCMS: 15;
    GCMS: 16;
    ELSD: 17;
    CD: 18;
    SFC: 19;
    EPR: 20;
    XRD: 21;
    RAMAN: 22;
    ED: 23;
    OPTICAL_ROTATION: 24;
    CAD: 25;
  }

  export const AnalysisType: AnalysisTypeMap;
}

export class ReactionProvenance extends jspb.Message {
  hasExperimenter(): boolean;
  clearExperimenter(): void;
  getExperimenter(): Person | undefined;
  setExperimenter(value?: Person): void;

  getCity(): string;
  setCity(value: string): void;

  hasExperimentStart(): boolean;
  clearExperimentStart(): void;
  getExperimentStart(): DateTime | undefined;
  setExperimentStart(value?: DateTime): void;

  getDoi(): string;
  setDoi(value: string): void;

  getPatent(): string;
  setPatent(value: string): void;

  getPublicationUrl(): string;
  setPublicationUrl(value: string): void;

  hasRecordCreated(): boolean;
  clearRecordCreated(): void;
  getRecordCreated(): RecordEvent | undefined;
  setRecordCreated(value?: RecordEvent): void;

  clearRecordModifiedList(): void;
  getRecordModifiedList(): Array<RecordEvent>;
  setRecordModifiedList(value: Array<RecordEvent>): void;
  addRecordModified(value?: RecordEvent, index?: number): RecordEvent;

  getReactionMetadataMap(): jspb.Map<string, Data>;
  clearReactionMetadataMap(): void;
  hasIsMined(): boolean;
  clearIsMined(): void;
  getIsMined(): boolean;
  setIsMined(value: boolean): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): ReactionProvenance.AsObject;
  static toObject(includeInstance: boolean, msg: ReactionProvenance): ReactionProvenance.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: ReactionProvenance, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): ReactionProvenance;
  static deserializeBinaryFromReader(message: ReactionProvenance, reader: jspb.BinaryReader): ReactionProvenance;
}

export namespace ReactionProvenance {
  export type AsObject = {
    experimenter?: Person.AsObject,
    city: string,
    experimentStart?: DateTime.AsObject,
    doi: string,
    patent: string,
    publicationUrl: string,
    recordCreated?: RecordEvent.AsObject,
    recordModifiedList: Array<RecordEvent.AsObject>,
    reactionMetadataMap: Array<[string, Data.AsObject]>,
    isMined: boolean,
  }
}

export class Person extends jspb.Message {
  getUsername(): string;
  setUsername(value: string): void;

  getName(): string;
  setName(value: string): void;

  getOrcid(): string;
  setOrcid(value: string): void;

  getOrganization(): string;
  setOrganization(value: string): void;

  getEmail(): string;
  setEmail(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Person.AsObject;
  static toObject(includeInstance: boolean, msg: Person): Person.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Person, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Person;
  static deserializeBinaryFromReader(message: Person, reader: jspb.BinaryReader): Person;
}

export namespace Person {
  export type AsObject = {
    username: string,
    name: string,
    orcid: string,
    organization: string,
    email: string,
  }
}

export class RecordEvent extends jspb.Message {
  hasTime(): boolean;
  clearTime(): void;
  getTime(): DateTime | undefined;
  setTime(value?: DateTime): void;

  hasPerson(): boolean;
  clearPerson(): void;
  getPerson(): Person | undefined;
  setPerson(value?: Person): void;

  getDetails(): string;
  setDetails(value: string): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): RecordEvent.AsObject;
  static toObject(includeInstance: boolean, msg: RecordEvent): RecordEvent.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: RecordEvent, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): RecordEvent;
  static deserializeBinaryFromReader(message: RecordEvent, reader: jspb.BinaryReader): RecordEvent;
}

export namespace RecordEvent {
  export type AsObject = {
    time?: DateTime.AsObject,
    person?: Person.AsObject,
    details: string,
  }
}

export class Time extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Time.TimeUnitMap[keyof Time.TimeUnitMap];
  setUnits(value: Time.TimeUnitMap[keyof Time.TimeUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Time.AsObject;
  static toObject(includeInstance: boolean, msg: Time): Time.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Time, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Time;
  static deserializeBinaryFromReader(message: Time, reader: jspb.BinaryReader): Time;
}

export namespace Time {
  export type AsObject = {
    value: number,
    precision: number,
    units: Time.TimeUnitMap[keyof Time.TimeUnitMap],
  }

  export interface TimeUnitMap {
    UNSPECIFIED: 0;
    DAY: 4;
    HOUR: 1;
    MINUTE: 2;
    SECOND: 3;
  }

  export const TimeUnit: TimeUnitMap;
}

export class Mass extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Mass.MassUnitMap[keyof Mass.MassUnitMap];
  setUnits(value: Mass.MassUnitMap[keyof Mass.MassUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Mass.AsObject;
  static toObject(includeInstance: boolean, msg: Mass): Mass.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Mass, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Mass;
  static deserializeBinaryFromReader(message: Mass, reader: jspb.BinaryReader): Mass;
}

export namespace Mass {
  export type AsObject = {
    value: number,
    precision: number,
    units: Mass.MassUnitMap[keyof Mass.MassUnitMap],
  }

  export interface MassUnitMap {
    UNSPECIFIED: 0;
    KILOGRAM: 1;
    GRAM: 2;
    MILLIGRAM: 3;
    MICROGRAM: 4;
  }

  export const MassUnit: MassUnitMap;
}

export class Moles extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Moles.MolesUnitMap[keyof Moles.MolesUnitMap];
  setUnits(value: Moles.MolesUnitMap[keyof Moles.MolesUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Moles.AsObject;
  static toObject(includeInstance: boolean, msg: Moles): Moles.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Moles, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Moles;
  static deserializeBinaryFromReader(message: Moles, reader: jspb.BinaryReader): Moles;
}

export namespace Moles {
  export type AsObject = {
    value: number,
    precision: number,
    units: Moles.MolesUnitMap[keyof Moles.MolesUnitMap],
  }

  export interface MolesUnitMap {
    UNSPECIFIED: 0;
    MOLE: 1;
    MILLIMOLE: 2;
    MICROMOLE: 3;
    NANOMOLE: 4;
  }

  export const MolesUnit: MolesUnitMap;
}

export class Volume extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Volume.VolumeUnitMap[keyof Volume.VolumeUnitMap];
  setUnits(value: Volume.VolumeUnitMap[keyof Volume.VolumeUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Volume.AsObject;
  static toObject(includeInstance: boolean, msg: Volume): Volume.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Volume, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Volume;
  static deserializeBinaryFromReader(message: Volume, reader: jspb.BinaryReader): Volume;
}

export namespace Volume {
  export type AsObject = {
    value: number,
    precision: number,
    units: Volume.VolumeUnitMap[keyof Volume.VolumeUnitMap],
  }

  export interface VolumeUnitMap {
    UNSPECIFIED: 0;
    LITER: 1;
    MILLILITER: 2;
    MICROLITER: 3;
    NANOLITER: 4;
  }

  export const VolumeUnit: VolumeUnitMap;
}

export class Concentration extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Concentration.ConcentrationUnitMap[keyof Concentration.ConcentrationUnitMap];
  setUnits(value: Concentration.ConcentrationUnitMap[keyof Concentration.ConcentrationUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Concentration.AsObject;
  static toObject(includeInstance: boolean, msg: Concentration): Concentration.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Concentration, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Concentration;
  static deserializeBinaryFromReader(message: Concentration, reader: jspb.BinaryReader): Concentration;
}

export namespace Concentration {
  export type AsObject = {
    value: number,
    precision: number,
    units: Concentration.ConcentrationUnitMap[keyof Concentration.ConcentrationUnitMap],
  }

  export interface ConcentrationUnitMap {
    UNSPECIFIED: 0;
    MOLAR: 1;
    MILLIMOLAR: 2;
    MICROMOLAR: 3;
  }

  export const ConcentrationUnit: ConcentrationUnitMap;
}

export class Pressure extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Pressure.PressureUnitMap[keyof Pressure.PressureUnitMap];
  setUnits(value: Pressure.PressureUnitMap[keyof Pressure.PressureUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Pressure.AsObject;
  static toObject(includeInstance: boolean, msg: Pressure): Pressure.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Pressure, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Pressure;
  static deserializeBinaryFromReader(message: Pressure, reader: jspb.BinaryReader): Pressure;
}

export namespace Pressure {
  export type AsObject = {
    value: number,
    precision: number,
    units: Pressure.PressureUnitMap[keyof Pressure.PressureUnitMap],
  }

  export interface PressureUnitMap {
    UNSPECIFIED: 0;
    BAR: 1;
    ATMOSPHERE: 2;
    PSI: 3;
    KPSI: 4;
    PASCAL: 5;
    KILOPASCAL: 6;
    TORR: 7;
    MM_HG: 8;
  }

  export const PressureUnit: PressureUnitMap;
}

export class Temperature extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Temperature.TemperatureUnitMap[keyof Temperature.TemperatureUnitMap];
  setUnits(value: Temperature.TemperatureUnitMap[keyof Temperature.TemperatureUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Temperature.AsObject;
  static toObject(includeInstance: boolean, msg: Temperature): Temperature.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Temperature, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Temperature;
  static deserializeBinaryFromReader(message: Temperature, reader: jspb.BinaryReader): Temperature;
}

export namespace Temperature {
  export type AsObject = {
    value: number,
    precision: number,
    units: Temperature.TemperatureUnitMap[keyof Temperature.TemperatureUnitMap],
  }

  export interface TemperatureUnitMap {
    UNSPECIFIED: 0;
    CELSIUS: 1;
    FAHRENHEIT: 2;
    KELVIN: 3;
  }

  export const TemperatureUnit: TemperatureUnitMap;
}

export class Current extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Current.CurrentUnitMap[keyof Current.CurrentUnitMap];
  setUnits(value: Current.CurrentUnitMap[keyof Current.CurrentUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Current.AsObject;
  static toObject(includeInstance: boolean, msg: Current): Current.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Current, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Current;
  static deserializeBinaryFromReader(message: Current, reader: jspb.BinaryReader): Current;
}

export namespace Current {
  export type AsObject = {
    value: number,
    precision: number,
    units: Current.CurrentUnitMap[keyof Current.CurrentUnitMap],
  }

  export interface CurrentUnitMap {
    UNSPECIFIED: 0;
    AMPERE: 1;
    MILLIAMPERE: 2;
  }

  export const CurrentUnit: CurrentUnitMap;
}

export class Voltage extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Voltage.VoltageUnitMap[keyof Voltage.VoltageUnitMap];
  setUnits(value: Voltage.VoltageUnitMap[keyof Voltage.VoltageUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Voltage.AsObject;
  static toObject(includeInstance: boolean, msg: Voltage): Voltage.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Voltage, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Voltage;
  static deserializeBinaryFromReader(message: Voltage, reader: jspb.BinaryReader): Voltage;
}

export namespace Voltage {
  export type AsObject = {
    value: number,
    precision: number,
    units: Voltage.VoltageUnitMap[keyof Voltage.VoltageUnitMap],
  }

  export interface VoltageUnitMap {
    UNSPECIFIED: 0;
    VOLT: 1;
    MILLIVOLT: 2;
  }

  export const VoltageUnit: VoltageUnitMap;
}

export class Length extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Length.LengthUnitMap[keyof Length.LengthUnitMap];
  setUnits(value: Length.LengthUnitMap[keyof Length.LengthUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Length.AsObject;
  static toObject(includeInstance: boolean, msg: Length): Length.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Length, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Length;
  static deserializeBinaryFromReader(message: Length, reader: jspb.BinaryReader): Length;
}

export namespace Length {
  export type AsObject = {
    value: number,
    precision: number,
    units: Length.LengthUnitMap[keyof Length.LengthUnitMap],
  }

  export interface LengthUnitMap {
    UNSPECIFIED: 0;
    CENTIMETER: 1;
    MILLIMETER: 2;
    METER: 3;
    INCH: 4;
    FOOT: 5;
  }

  export const LengthUnit: LengthUnitMap;
}

export class Wavelength extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): Wavelength.WavelengthUnitMap[keyof Wavelength.WavelengthUnitMap];
  setUnits(value: Wavelength.WavelengthUnitMap[keyof Wavelength.WavelengthUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Wavelength.AsObject;
  static toObject(includeInstance: boolean, msg: Wavelength): Wavelength.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Wavelength, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Wavelength;
  static deserializeBinaryFromReader(message: Wavelength, reader: jspb.BinaryReader): Wavelength;
}

export namespace Wavelength {
  export type AsObject = {
    value: number,
    precision: number,
    units: Wavelength.WavelengthUnitMap[keyof Wavelength.WavelengthUnitMap],
  }

  export interface WavelengthUnitMap {
    UNSPECIFIED: 0;
    NANOMETER: 1;
    WAVENUMBER: 2;
  }

  export const WavelengthUnit: WavelengthUnitMap;
}

export class FlowRate extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  getUnits(): FlowRate.FlowRateUnitMap[keyof FlowRate.FlowRateUnitMap];
  setUnits(value: FlowRate.FlowRateUnitMap[keyof FlowRate.FlowRateUnitMap]): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): FlowRate.AsObject;
  static toObject(includeInstance: boolean, msg: FlowRate): FlowRate.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: FlowRate, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): FlowRate;
  static deserializeBinaryFromReader(message: FlowRate, reader: jspb.BinaryReader): FlowRate;
}

export namespace FlowRate {
  export type AsObject = {
    value: number,
    precision: number,
    units: FlowRate.FlowRateUnitMap[keyof FlowRate.FlowRateUnitMap],
  }

  export interface FlowRateUnitMap {
    UNSPECIFIED: 0;
    MICROLITER_PER_MINUTE: 1;
    MICROLITER_PER_SECOND: 2;
    MILLILITER_PER_MINUTE: 3;
    MILLILITER_PER_SECOND: 4;
    MICROLITER_PER_HOUR: 5;
  }

  export const FlowRateUnit: FlowRateUnitMap;
}

export class Percentage extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Percentage.AsObject;
  static toObject(includeInstance: boolean, msg: Percentage): Percentage.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Percentage, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Percentage;
  static deserializeBinaryFromReader(message: Percentage, reader: jspb.BinaryReader): Percentage;
}

export namespace Percentage {
  export type AsObject = {
    value: number,
    precision: number,
  }
}

export class FloatValue extends jspb.Message {
  hasValue(): boolean;
  clearValue(): void;
  getValue(): number;
  setValue(value: number): void;

  hasPrecision(): boolean;
  clearPrecision(): void;
  getPrecision(): number;
  setPrecision(value: number): void;

  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): FloatValue.AsObject;
  static toObject(includeInstance: boolean, msg: FloatValue): FloatValue.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: FloatValue, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): FloatValue;
  static deserializeBinaryFromReader(message: FloatValue, reader: jspb.BinaryReader): FloatValue;
}

export namespace FloatValue {
  export type AsObject = {
    value: number,
    precision: number,
  }
}

export class Data extends jspb.Message {
  hasFloatValue(): boolean;
  clearFloatValue(): void;
  getFloatValue(): number;
  setFloatValue(value: number): void;

  hasIntegerValue(): boolean;
  clearIntegerValue(): void;
  getIntegerValue(): number;
  setIntegerValue(value: number): void;

  hasBytesValue(): boolean;
  clearBytesValue(): void;
  getBytesValue(): Uint8Array | string;
  getBytesValue_asU8(): Uint8Array;
  getBytesValue_asB64(): string;
  setBytesValue(value: Uint8Array | string): void;

  hasStringValue(): boolean;
  clearStringValue(): void;
  getStringValue(): string;
  setStringValue(value: string): void;

  hasUrl(): boolean;
  clearUrl(): void;
  getUrl(): string;
  setUrl(value: string): void;

  getDescription(): string;
  setDescription(value: string): void;

  getFormat(): string;
  setFormat(value: string): void;

  getKindCase(): Data.KindCase;
  serializeBinary(): Uint8Array;
  toObject(includeInstance?: boolean): Data.AsObject;
  static toObject(includeInstance: boolean, msg: Data): Data.AsObject;
  static extensions: {[key: number]: jspb.ExtensionFieldInfo<jspb.Message>};
  static extensionsBinary: {[key: number]: jspb.ExtensionFieldBinaryInfo<jspb.Message>};
  static serializeBinaryToWriter(message: Data, writer: jspb.BinaryWriter): void;
  static deserializeBinary(bytes: Uint8Array): Data;
  static deserializeBinaryFromReader(message: Data, reader: jspb.BinaryReader): Data;
}

export namespace Data {
  export type AsObject = {
    floatValue: number,
    integerValue: number,
    bytesValue: Uint8Array | string,
    stringValue: string,
    url: string,
    description: string,
    format: string,
  }

  export enum KindCase {
    KIND_NOT_SET = 0,
    FLOAT_VALUE = 1,
    INTEGER_VALUE = 2,
    BYTES_VALUE = 3,
    STRING_VALUE = 4,
    URL = 5,
  }
}

