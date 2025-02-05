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

import * as $protobuf from "protobufjs";
import Long = require("long");
/** Namespace ord. */
export namespace ord {

    /** Properties of a Reaction. */
    interface IReaction {

        /** Reaction identifiers */
        identifiers?: (ord.IReactionIdentifier[]|null);

        /** Reaction inputs */
        inputs?: ({ [k: string]: ord.IReactionInput }|null);

        /** Reaction setup */
        setup?: (ord.IReactionSetup|null);

        /** Reaction conditions */
        conditions?: (ord.IReactionConditions|null);

        /** Reaction notes */
        notes?: (ord.IReactionNotes|null);

        /** Reaction observations */
        observations?: (ord.IReactionObservation[]|null);

        /** Reaction workups */
        workups?: (ord.IReactionWorkup[]|null);

        /** Reaction outcomes */
        outcomes?: (ord.IReactionOutcome[]|null);

        /** Reaction provenance */
        provenance?: (ord.IReactionProvenance|null);

        /** Reaction reactionId */
        reactionId?: (string|null);
    }

    /**
     * Throughout this schema, we introduce enums to encourage consistency in
     * nomenclature and to avoid unnecessary downstream data processing that would
     * otherwise be required to consolidate equivalent entries. However, we do
     * not wish to restrict what users are able to specify if their synthesis
     * does not fit cleanly into a pre-existing enum field. For that reason, many
     * enums contain a CUSTOM field, which must be accompanied by setting the
     * 'details' field, where appropriate).
     *
     * NOTE(kearnes): In many places, we deliberately violate the style guide for
     * enums by nesting instead of prefixing; this is not done lightly. The primary
     * consideration is API consistency and the ability to use unqualified strings
     * as enum values. For instance, we want 'CUSTOM' to be a valid value for all
     * enums that support custom types.
     */
    class Reaction implements IReaction {

        /**
         * Constructs a new Reaction.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReaction);

        /** Reaction identifiers. */
        public identifiers: ord.IReactionIdentifier[];

        /** Reaction inputs. */
        public inputs: { [k: string]: ord.IReactionInput };

        /** Reaction setup. */
        public setup?: (ord.IReactionSetup|null);

        /** Reaction conditions. */
        public conditions?: (ord.IReactionConditions|null);

        /** Reaction notes. */
        public notes?: (ord.IReactionNotes|null);

        /** Reaction observations. */
        public observations: ord.IReactionObservation[];

        /** Reaction workups. */
        public workups: ord.IReactionWorkup[];

        /** Reaction outcomes. */
        public outcomes: ord.IReactionOutcome[];

        /** Reaction provenance. */
        public provenance?: (ord.IReactionProvenance|null);

        /** Reaction reactionId. */
        public reactionId: string;

        /**
         * Creates a new Reaction instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Reaction instance
         */
        public static create(properties?: ord.IReaction): ord.Reaction;

        /**
         * Encodes the specified Reaction message. Does not implicitly {@link ord.Reaction.verify|verify} messages.
         * @param message Reaction message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReaction, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Reaction message, length delimited. Does not implicitly {@link ord.Reaction.verify|verify} messages.
         * @param message Reaction message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReaction, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Reaction message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Reaction
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Reaction;

        /**
         * Decodes a Reaction message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Reaction
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Reaction;

        /**
         * Verifies a Reaction message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Reaction message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Reaction
         */
        public static fromObject(object: { [k: string]: any }): ord.Reaction;

        /**
         * Creates a plain object from a Reaction message. Also converts values to other types if specified.
         * @param message Reaction
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Reaction, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Reaction to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Reaction
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a ReactionIdentifier. */
    interface IReactionIdentifier {

        /** ReactionIdentifier type */
        type?: (ord.ReactionIdentifier.ReactionIdentifierType|null);

        /** ReactionIdentifier details */
        details?: (string|null);

        /** ReactionIdentifier value */
        value?: (string|null);

        /** ReactionIdentifier isMapped */
        isMapped?: (boolean|null);
    }

    /**
     * Reaction identifiers define descriptions of the overall reaction.
     * While we encourage the use of SMILES strings, these do not work well in
     * all cases. The <reaction_smiles> field should be able to be derived
     * from the information present in the ReactionInput and ReactionOutcome
     * fields of any Reaction message.
     */
    class ReactionIdentifier implements IReactionIdentifier {

        /**
         * Constructs a new ReactionIdentifier.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionIdentifier);

        /** ReactionIdentifier type. */
        public type: ord.ReactionIdentifier.ReactionIdentifierType;

        /** ReactionIdentifier details. */
        public details: string;

        /** ReactionIdentifier value. */
        public value: string;

        /** ReactionIdentifier isMapped. */
        public isMapped?: (boolean|null);

        /**
         * Creates a new ReactionIdentifier instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionIdentifier instance
         */
        public static create(properties?: ord.IReactionIdentifier): ord.ReactionIdentifier;

        /**
         * Encodes the specified ReactionIdentifier message. Does not implicitly {@link ord.ReactionIdentifier.verify|verify} messages.
         * @param message ReactionIdentifier message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionIdentifier, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionIdentifier message, length delimited. Does not implicitly {@link ord.ReactionIdentifier.verify|verify} messages.
         * @param message ReactionIdentifier message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionIdentifier, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionIdentifier message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionIdentifier
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionIdentifier;

        /**
         * Decodes a ReactionIdentifier message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionIdentifier
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionIdentifier;

        /**
         * Verifies a ReactionIdentifier message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionIdentifier message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionIdentifier
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionIdentifier;

        /**
         * Creates a plain object from a ReactionIdentifier message. Also converts values to other types if specified.
         * @param message ReactionIdentifier
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionIdentifier, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionIdentifier to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionIdentifier
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ReactionIdentifier {

        /** ReactionIdentifierType enum. */
        enum ReactionIdentifierType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            REACTION_SMILES = 2,
            REACTION_CXSMILES = 6,
            RDFILE = 3,
            RINCHI = 4,
            REACTION_TYPE = 5
        }
    }

    /** Properties of a ReactionInput. */
    interface IReactionInput {

        /** ReactionInput components */
        components?: (ord.ICompound[]|null);

        /** ReactionInput crudeComponents */
        crudeComponents?: (ord.ICrudeComponent[]|null);

        /** ReactionInput additionOrder */
        additionOrder?: (number|null);

        /** ReactionInput additionTime */
        additionTime?: (ord.ITime|null);

        /** ReactionInput additionSpeed */
        additionSpeed?: (ord.ReactionInput.IAdditionSpeed|null);

        /** ReactionInput additionDuration */
        additionDuration?: (ord.ITime|null);

        /** ReactionInput flowRate */
        flowRate?: (ord.IFlowRate|null);

        /** ReactionInput additionDevice */
        additionDevice?: (ord.ReactionInput.IAdditionDevice|null);

        /** ReactionInput additionTemperature */
        additionTemperature?: (ord.ITemperature|null);

        /** ReactionInput texture */
        texture?: (ord.ITexture|null);
    }

    /**
     * A reaction input is any pure substance, mixture, or solution that is
     * added to the reaction vessel.
     *
     * For example, suppose we are adding 3 mL of a 4 M solution of NaOH in water.
     * We would define one component for the solvent and one component for the
     * solute with the correct respective amounts.
     *
     * input {
     * components: {
     * identifiers: {type: IDENTIFIER_SMILES, value: "O"}
     * identifiers: {type: IDENTIFIER_NAME, value: "water"}
     * volume: {value: 3, units: MILLILITER}
     * volume_includes_solutes: true
     * }
     * components: {
     * identifiers: {type: IDENTIFIER_SMILES, value: "[Na+].[OH-]"}
     * identifiers: {type: IDENTIFIER_NAME, value: "sodium hydroxide"}
     * moles: {value: 12, units: MILLIMOLES}
     * }
     * }
     */
    class ReactionInput implements IReactionInput {

        /**
         * Constructs a new ReactionInput.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionInput);

        /** ReactionInput components. */
        public components: ord.ICompound[];

        /** ReactionInput crudeComponents. */
        public crudeComponents: ord.ICrudeComponent[];

        /** ReactionInput additionOrder. */
        public additionOrder: number;

        /** ReactionInput additionTime. */
        public additionTime?: (ord.ITime|null);

        /** ReactionInput additionSpeed. */
        public additionSpeed?: (ord.ReactionInput.IAdditionSpeed|null);

        /** ReactionInput additionDuration. */
        public additionDuration?: (ord.ITime|null);

        /** ReactionInput flowRate. */
        public flowRate?: (ord.IFlowRate|null);

        /** ReactionInput additionDevice. */
        public additionDevice?: (ord.ReactionInput.IAdditionDevice|null);

        /** ReactionInput additionTemperature. */
        public additionTemperature?: (ord.ITemperature|null);

        /** ReactionInput texture. */
        public texture?: (ord.ITexture|null);

        /**
         * Creates a new ReactionInput instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionInput instance
         */
        public static create(properties?: ord.IReactionInput): ord.ReactionInput;

        /**
         * Encodes the specified ReactionInput message. Does not implicitly {@link ord.ReactionInput.verify|verify} messages.
         * @param message ReactionInput message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionInput, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionInput message, length delimited. Does not implicitly {@link ord.ReactionInput.verify|verify} messages.
         * @param message ReactionInput message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionInput, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionInput message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionInput
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionInput;

        /**
         * Decodes a ReactionInput message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionInput
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionInput;

        /**
         * Verifies a ReactionInput message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionInput message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionInput
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionInput;

        /**
         * Creates a plain object from a ReactionInput message. Also converts values to other types if specified.
         * @param message ReactionInput
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionInput, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionInput to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionInput
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ReactionInput {

        /** Properties of an AdditionSpeed. */
        interface IAdditionSpeed {

            /** AdditionSpeed type */
            type?: (ord.ReactionInput.AdditionSpeed.AdditionSpeedType|null);

            /** AdditionSpeed details */
            details?: (string|null);
        }

        /** Represents an AdditionSpeed. */
        class AdditionSpeed implements IAdditionSpeed {

            /**
             * Constructs a new AdditionSpeed.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ReactionInput.IAdditionSpeed);

            /** AdditionSpeed type. */
            public type: ord.ReactionInput.AdditionSpeed.AdditionSpeedType;

            /** AdditionSpeed details. */
            public details: string;

            /**
             * Creates a new AdditionSpeed instance using the specified properties.
             * @param [properties] Properties to set
             * @returns AdditionSpeed instance
             */
            public static create(properties?: ord.ReactionInput.IAdditionSpeed): ord.ReactionInput.AdditionSpeed;

            /**
             * Encodes the specified AdditionSpeed message. Does not implicitly {@link ord.ReactionInput.AdditionSpeed.verify|verify} messages.
             * @param message AdditionSpeed message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ReactionInput.IAdditionSpeed, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified AdditionSpeed message, length delimited. Does not implicitly {@link ord.ReactionInput.AdditionSpeed.verify|verify} messages.
             * @param message AdditionSpeed message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ReactionInput.IAdditionSpeed, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes an AdditionSpeed message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns AdditionSpeed
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionInput.AdditionSpeed;

            /**
             * Decodes an AdditionSpeed message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns AdditionSpeed
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionInput.AdditionSpeed;

            /**
             * Verifies an AdditionSpeed message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates an AdditionSpeed message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns AdditionSpeed
             */
            public static fromObject(object: { [k: string]: any }): ord.ReactionInput.AdditionSpeed;

            /**
             * Creates a plain object from an AdditionSpeed message. Also converts values to other types if specified.
             * @param message AdditionSpeed
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ReactionInput.AdditionSpeed, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this AdditionSpeed to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for AdditionSpeed
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace AdditionSpeed {

            /** AdditionSpeedType enum. */
            enum AdditionSpeedType {
                UNSPECIFIED = 0,
                ALL_AT_ONCE = 1,
                FAST = 2,
                SLOW = 3,
                DROPWISE = 4,
                CONTINUOUS = 5,
                PORTIONWISE = 6
            }
        }

        /** Properties of an AdditionDevice. */
        interface IAdditionDevice {

            /** AdditionDevice type */
            type?: (ord.ReactionInput.AdditionDevice.AdditionDeviceType|null);

            /** AdditionDevice details */
            details?: (string|null);
        }

        /** Represents an AdditionDevice. */
        class AdditionDevice implements IAdditionDevice {

            /**
             * Constructs a new AdditionDevice.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ReactionInput.IAdditionDevice);

            /** AdditionDevice type. */
            public type: ord.ReactionInput.AdditionDevice.AdditionDeviceType;

            /** AdditionDevice details. */
            public details: string;

            /**
             * Creates a new AdditionDevice instance using the specified properties.
             * @param [properties] Properties to set
             * @returns AdditionDevice instance
             */
            public static create(properties?: ord.ReactionInput.IAdditionDevice): ord.ReactionInput.AdditionDevice;

            /**
             * Encodes the specified AdditionDevice message. Does not implicitly {@link ord.ReactionInput.AdditionDevice.verify|verify} messages.
             * @param message AdditionDevice message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ReactionInput.IAdditionDevice, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified AdditionDevice message, length delimited. Does not implicitly {@link ord.ReactionInput.AdditionDevice.verify|verify} messages.
             * @param message AdditionDevice message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ReactionInput.IAdditionDevice, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes an AdditionDevice message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns AdditionDevice
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionInput.AdditionDevice;

            /**
             * Decodes an AdditionDevice message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns AdditionDevice
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionInput.AdditionDevice;

            /**
             * Verifies an AdditionDevice message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates an AdditionDevice message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns AdditionDevice
             */
            public static fromObject(object: { [k: string]: any }): ord.ReactionInput.AdditionDevice;

            /**
             * Creates a plain object from an AdditionDevice message. Also converts values to other types if specified.
             * @param message AdditionDevice
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ReactionInput.AdditionDevice, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this AdditionDevice to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for AdditionDevice
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace AdditionDevice {

            /** AdditionDeviceType enum. */
            enum AdditionDeviceType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                NONE = 2,
                SYRINGE = 3,
                CANNULA = 4,
                ADDITION_FUNNEL = 5,
                PIPETTE = 6,
                POSITIVE_DISPLACEMENT_PIPETTE = 7,
                PISTON_PUMP = 8,
                SYRINGE_PUMP = 9,
                PERISTALTIC_PUMP = 10
            }
        }
    }

    /** Properties of an Amount. */
    interface IAmount {

        /** Amount mass */
        mass?: (ord.IMass|null);

        /** Amount moles */
        moles?: (ord.IMoles|null);

        /** Amount volume */
        volume?: (ord.IVolume|null);

        /** Amount unmeasured */
        unmeasured?: (ord.IUnmeasuredAmount|null);

        /** Amount volumeIncludesSolutes */
        volumeIncludesSolutes?: (boolean|null);
    }

    /**
     * The quantitative amount of a Compound used in a particular reaction.
     * Compounds added in their pure form should have their value defined by
     * mass, moles, or volume. Compounds prepared as solutions should be defined
     * in terms of their volume. Compounds prepared on solid supports should
     * define the total mass/volume including the support.
     */
    class Amount implements IAmount {

        /**
         * Constructs a new Amount.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IAmount);

        /** Amount mass. */
        public mass?: (ord.IMass|null);

        /** Amount moles. */
        public moles?: (ord.IMoles|null);

        /** Amount volume. */
        public volume?: (ord.IVolume|null);

        /** Amount unmeasured. */
        public unmeasured?: (ord.IUnmeasuredAmount|null);

        /** Amount volumeIncludesSolutes. */
        public volumeIncludesSolutes?: (boolean|null);

        /** Amount kind. */
        public kind?: ("mass"|"moles"|"volume"|"unmeasured");

        /**
         * Creates a new Amount instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Amount instance
         */
        public static create(properties?: ord.IAmount): ord.Amount;

        /**
         * Encodes the specified Amount message. Does not implicitly {@link ord.Amount.verify|verify} messages.
         * @param message Amount message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IAmount, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Amount message, length delimited. Does not implicitly {@link ord.Amount.verify|verify} messages.
         * @param message Amount message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IAmount, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes an Amount message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Amount
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Amount;

        /**
         * Decodes an Amount message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Amount
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Amount;

        /**
         * Verifies an Amount message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates an Amount message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Amount
         */
        public static fromObject(object: { [k: string]: any }): ord.Amount;

        /**
         * Creates a plain object from an Amount message. Also converts values to other types if specified.
         * @param message Amount
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Amount, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Amount to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Amount
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of an UnmeasuredAmount. */
    interface IUnmeasuredAmount {

        /** UnmeasuredAmount type */
        type?: (ord.UnmeasuredAmount.UnmeasuredAmountType|null);

        /** UnmeasuredAmount details */
        details?: (string|null);
    }

    /**
     * Compounds may be defined with qualitative amounts in situations where a
     * precise quantity is not measured or reported.
     */
    class UnmeasuredAmount implements IUnmeasuredAmount {

        /**
         * Constructs a new UnmeasuredAmount.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IUnmeasuredAmount);

        /** UnmeasuredAmount type. */
        public type: ord.UnmeasuredAmount.UnmeasuredAmountType;

        /** UnmeasuredAmount details. */
        public details: string;

        /**
         * Creates a new UnmeasuredAmount instance using the specified properties.
         * @param [properties] Properties to set
         * @returns UnmeasuredAmount instance
         */
        public static create(properties?: ord.IUnmeasuredAmount): ord.UnmeasuredAmount;

        /**
         * Encodes the specified UnmeasuredAmount message. Does not implicitly {@link ord.UnmeasuredAmount.verify|verify} messages.
         * @param message UnmeasuredAmount message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IUnmeasuredAmount, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified UnmeasuredAmount message, length delimited. Does not implicitly {@link ord.UnmeasuredAmount.verify|verify} messages.
         * @param message UnmeasuredAmount message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IUnmeasuredAmount, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes an UnmeasuredAmount message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns UnmeasuredAmount
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.UnmeasuredAmount;

        /**
         * Decodes an UnmeasuredAmount message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns UnmeasuredAmount
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.UnmeasuredAmount;

        /**
         * Verifies an UnmeasuredAmount message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates an UnmeasuredAmount message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns UnmeasuredAmount
         */
        public static fromObject(object: { [k: string]: any }): ord.UnmeasuredAmount;

        /**
         * Creates a plain object from an UnmeasuredAmount message. Also converts values to other types if specified.
         * @param message UnmeasuredAmount
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.UnmeasuredAmount, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this UnmeasuredAmount to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for UnmeasuredAmount
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace UnmeasuredAmount {

        /** UnmeasuredAmountType enum. */
        enum UnmeasuredAmountType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            SATURATED = 2,
            CATALYTIC = 3,
            TITRATED = 4
        }
    }

    /** Properties of a Texture. */
    interface ITexture {

        /** Texture type */
        type?: (ord.Texture.TextureType|null);

        /** Texture details */
        details?: (string|null);
    }

    /**
     * This qualitatively describes the apparent size and morphology of a Compound,
     * a ProductCompound, or a ReactionInput.
     */
    class Texture implements ITexture {

        /**
         * Constructs a new Texture.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ITexture);

        /** Texture type. */
        public type: ord.Texture.TextureType;

        /** Texture details. */
        public details: string;

        /**
         * Creates a new Texture instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Texture instance
         */
        public static create(properties?: ord.ITexture): ord.Texture;

        /**
         * Encodes the specified Texture message. Does not implicitly {@link ord.Texture.verify|verify} messages.
         * @param message Texture message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ITexture, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Texture message, length delimited. Does not implicitly {@link ord.Texture.verify|verify} messages.
         * @param message Texture message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ITexture, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Texture message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Texture
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Texture;

        /**
         * Decodes a Texture message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Texture
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Texture;

        /**
         * Verifies a Texture message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Texture message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Texture
         */
        public static fromObject(object: { [k: string]: any }): ord.Texture;

        /**
         * Creates a plain object from a Texture message. Also converts values to other types if specified.
         * @param message Texture
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Texture, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Texture to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Texture
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Texture {

        /** TextureType enum. */
        enum TextureType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            POWDER = 2,
            CRYSTAL = 3,
            OIL = 4,
            AMORPHOUS_SOLID = 5,
            FOAM = 6,
            WAX = 7,
            SEMI_SOLID = 8,
            SOLID = 9,
            LIQUID = 10,
            GAS = 11
        }
    }

    /** Properties of a CrudeComponent. */
    interface ICrudeComponent {

        /** CrudeComponent reactionId */
        reactionId?: (string|null);

        /** CrudeComponent includesWorkup */
        includesWorkup?: (boolean|null);

        /** CrudeComponent hasDerivedAmount */
        hasDerivedAmount?: (boolean|null);

        /** CrudeComponent amount */
        amount?: (ord.IAmount|null);

        /** CrudeComponent texture */
        texture?: (ord.ITexture|null);
    }

    /**
     * Crude components are used in multi-step or multi-stage reactions (no strong
     * distinction is made here) where one synthetic process must be described by
     * multiple "Reaction" messages. In these cases, we often carry the crude
     * product from one step/stage into the next. This message is only to be used
     * when there is not complete isolation of an intermediate molecule; if there
     * is complete isolation, then a regular Compound should be used with the
     * SYNTHESIED preparation type.
     */
    class CrudeComponent implements ICrudeComponent {

        /**
         * Constructs a new CrudeComponent.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ICrudeComponent);

        /** CrudeComponent reactionId. */
        public reactionId: string;

        /** CrudeComponent includesWorkup. */
        public includesWorkup?: (boolean|null);

        /** CrudeComponent hasDerivedAmount. */
        public hasDerivedAmount?: (boolean|null);

        /** CrudeComponent amount. */
        public amount?: (ord.IAmount|null);

        /** CrudeComponent texture. */
        public texture?: (ord.ITexture|null);

        /**
         * Creates a new CrudeComponent instance using the specified properties.
         * @param [properties] Properties to set
         * @returns CrudeComponent instance
         */
        public static create(properties?: ord.ICrudeComponent): ord.CrudeComponent;

        /**
         * Encodes the specified CrudeComponent message. Does not implicitly {@link ord.CrudeComponent.verify|verify} messages.
         * @param message CrudeComponent message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ICrudeComponent, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified CrudeComponent message, length delimited. Does not implicitly {@link ord.CrudeComponent.verify|verify} messages.
         * @param message CrudeComponent message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ICrudeComponent, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a CrudeComponent message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns CrudeComponent
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.CrudeComponent;

        /**
         * Decodes a CrudeComponent message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns CrudeComponent
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.CrudeComponent;

        /**
         * Verifies a CrudeComponent message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a CrudeComponent message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns CrudeComponent
         */
        public static fromObject(object: { [k: string]: any }): ord.CrudeComponent;

        /**
         * Creates a plain object from a CrudeComponent message. Also converts values to other types if specified.
         * @param message CrudeComponent
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.CrudeComponent, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this CrudeComponent to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for CrudeComponent
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a Compound. */
    interface ICompound {

        /** Compound identifiers */
        identifiers?: (ord.ICompoundIdentifier[]|null);

        /** Compound amount */
        amount?: (ord.IAmount|null);

        /** Compound reactionRole */
        reactionRole?: (ord.ReactionRole.ReactionRoleType|null);

        /** Compound isLimiting */
        isLimiting?: (boolean|null);

        /** Compound preparations */
        preparations?: (ord.ICompoundPreparation[]|null);

        /** Compound source */
        source?: (ord.Compound.ISource|null);

        /** Compound features */
        features?: ({ [k: string]: ord.IData }|null);

        /** Compound analyses */
        analyses?: ({ [k: string]: ord.IAnalysis }|null);

        /** Compound texture */
        texture?: (ord.ITexture|null);
    }

    /**
     * A Compound defines both the identity of a pure species and a quantitative
     * amount (mass, moles, volume). For compounds used in inputs, details can
     * be provided about how it was prepared and from where it was purchased.
     */
    class Compound implements ICompound {

        /**
         * Constructs a new Compound.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ICompound);

        /** Compound identifiers. */
        public identifiers: ord.ICompoundIdentifier[];

        /** Compound amount. */
        public amount?: (ord.IAmount|null);

        /** Compound reactionRole. */
        public reactionRole: ord.ReactionRole.ReactionRoleType;

        /** Compound isLimiting. */
        public isLimiting?: (boolean|null);

        /** Compound preparations. */
        public preparations: ord.ICompoundPreparation[];

        /** Compound source. */
        public source?: (ord.Compound.ISource|null);

        /** Compound features. */
        public features: { [k: string]: ord.IData };

        /** Compound analyses. */
        public analyses: { [k: string]: ord.IAnalysis };

        /** Compound texture. */
        public texture?: (ord.ITexture|null);

        /**
         * Creates a new Compound instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Compound instance
         */
        public static create(properties?: ord.ICompound): ord.Compound;

        /**
         * Encodes the specified Compound message. Does not implicitly {@link ord.Compound.verify|verify} messages.
         * @param message Compound message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ICompound, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Compound message, length delimited. Does not implicitly {@link ord.Compound.verify|verify} messages.
         * @param message Compound message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ICompound, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Compound message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Compound
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Compound;

        /**
         * Decodes a Compound message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Compound
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Compound;

        /**
         * Verifies a Compound message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Compound message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Compound
         */
        public static fromObject(object: { [k: string]: any }): ord.Compound;

        /**
         * Creates a plain object from a Compound message. Also converts values to other types if specified.
         * @param message Compound
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Compound, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Compound to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Compound
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Compound {

        /** Properties of a Source. */
        interface ISource {

            /** Source vendor */
            vendor?: (string|null);

            /** Source catalogId */
            catalogId?: (string|null);

            /** Source lot */
            lot?: (string|null);
        }

        /** Represents a Source. */
        class Source implements ISource {

            /**
             * Constructs a new Source.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.Compound.ISource);

            /** Source vendor. */
            public vendor: string;

            /** Source catalogId. */
            public catalogId: string;

            /** Source lot. */
            public lot: string;

            /**
             * Creates a new Source instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Source instance
             */
            public static create(properties?: ord.Compound.ISource): ord.Compound.Source;

            /**
             * Encodes the specified Source message. Does not implicitly {@link ord.Compound.Source.verify|verify} messages.
             * @param message Source message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.Compound.ISource, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Source message, length delimited. Does not implicitly {@link ord.Compound.Source.verify|verify} messages.
             * @param message Source message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.Compound.ISource, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a Source message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Source
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Compound.Source;

            /**
             * Decodes a Source message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Source
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Compound.Source;

            /**
             * Verifies a Source message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a Source message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Source
             */
            public static fromObject(object: { [k: string]: any }): ord.Compound.Source;

            /**
             * Creates a plain object from a Source message. Also converts values to other types if specified.
             * @param message Source
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.Compound.Source, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Source to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Source
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }
    }

    /** Properties of a ReactionRole. */
    interface IReactionRole {
    }

    /** Represents a ReactionRole. */
    class ReactionRole implements IReactionRole {

        /**
         * Constructs a new ReactionRole.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionRole);

        /**
         * Creates a new ReactionRole instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionRole instance
         */
        public static create(properties?: ord.IReactionRole): ord.ReactionRole;

        /**
         * Encodes the specified ReactionRole message. Does not implicitly {@link ord.ReactionRole.verify|verify} messages.
         * @param message ReactionRole message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionRole, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionRole message, length delimited. Does not implicitly {@link ord.ReactionRole.verify|verify} messages.
         * @param message ReactionRole message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionRole, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionRole message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionRole
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionRole;

        /**
         * Decodes a ReactionRole message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionRole
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionRole;

        /**
         * Verifies a ReactionRole message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionRole message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionRole
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionRole;

        /**
         * Creates a plain object from a ReactionRole message. Also converts values to other types if specified.
         * @param message ReactionRole
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionRole, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionRole to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionRole
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ReactionRole {

        /** ReactionRoleType enum. */
        enum ReactionRoleType {
            UNSPECIFIED = 0,
            REACTANT = 1,
            REAGENT = 2,
            SOLVENT = 3,
            CATALYST = 4,
            WORKUP = 5,
            INTERNAL_STANDARD = 6,
            AUTHENTIC_STANDARD = 7,
            PRODUCT = 8,
            BYPRODUCT = 9,
            SIDE_PRODUCT = 10
        }
    }

    /** Properties of a CompoundPreparation. */
    interface ICompoundPreparation {

        /** CompoundPreparation type */
        type?: (ord.CompoundPreparation.CompoundPreparationType|null);

        /** CompoundPreparation details */
        details?: (string|null);

        /** CompoundPreparation reactionId */
        reactionId?: (string|null);
    }

    /**
     * Compounds may undergo additional preparation before being used in a
     * reaction after being received from a supplier or vendor. We encourage
     * the use of the 'preparation' enum when possible, even if the description
     * is an oversimplification of the full procedure, which can be described
     * in the 'details' field.
     */
    class CompoundPreparation implements ICompoundPreparation {

        /**
         * Constructs a new CompoundPreparation.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ICompoundPreparation);

        /** CompoundPreparation type. */
        public type: ord.CompoundPreparation.CompoundPreparationType;

        /** CompoundPreparation details. */
        public details: string;

        /** CompoundPreparation reactionId. */
        public reactionId: string;

        /**
         * Creates a new CompoundPreparation instance using the specified properties.
         * @param [properties] Properties to set
         * @returns CompoundPreparation instance
         */
        public static create(properties?: ord.ICompoundPreparation): ord.CompoundPreparation;

        /**
         * Encodes the specified CompoundPreparation message. Does not implicitly {@link ord.CompoundPreparation.verify|verify} messages.
         * @param message CompoundPreparation message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ICompoundPreparation, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified CompoundPreparation message, length delimited. Does not implicitly {@link ord.CompoundPreparation.verify|verify} messages.
         * @param message CompoundPreparation message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ICompoundPreparation, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a CompoundPreparation message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns CompoundPreparation
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.CompoundPreparation;

        /**
         * Decodes a CompoundPreparation message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns CompoundPreparation
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.CompoundPreparation;

        /**
         * Verifies a CompoundPreparation message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a CompoundPreparation message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns CompoundPreparation
         */
        public static fromObject(object: { [k: string]: any }): ord.CompoundPreparation;

        /**
         * Creates a plain object from a CompoundPreparation message. Also converts values to other types if specified.
         * @param message CompoundPreparation
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.CompoundPreparation, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this CompoundPreparation to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for CompoundPreparation
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace CompoundPreparation {

        /** CompoundPreparationType enum. */
        enum CompoundPreparationType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            NONE = 2,
            REPURIFIED = 3,
            SPARGED = 4,
            DRIED = 5,
            SYNTHESIZED = 6
        }
    }

    /** Properties of a CompoundIdentifier. */
    interface ICompoundIdentifier {

        /** CompoundIdentifier type */
        type?: (ord.CompoundIdentifier.CompoundIdentifierType|null);

        /** CompoundIdentifier details */
        details?: (string|null);

        /** CompoundIdentifier value */
        value?: (string|null);
    }

    /**
     * Compound identifiers uniquely define a single (pure) chemical species.
     * While we encourage the use of SMILES strings, these do not work well in
     * all cases (e.g., handling tautomerism, axial chirality). Multiple
     * identifiers may be specified for a single compound to avoid ambiguity.
     * We discourage chemicals from being defined only by a name. For compounds
     * that are prepared or isolated as salts, the identifier should include
     * specification of which salt.
     */
    class CompoundIdentifier implements ICompoundIdentifier {

        /**
         * Constructs a new CompoundIdentifier.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ICompoundIdentifier);

        /** CompoundIdentifier type. */
        public type: ord.CompoundIdentifier.CompoundIdentifierType;

        /** CompoundIdentifier details. */
        public details: string;

        /** CompoundIdentifier value. */
        public value: string;

        /**
         * Creates a new CompoundIdentifier instance using the specified properties.
         * @param [properties] Properties to set
         * @returns CompoundIdentifier instance
         */
        public static create(properties?: ord.ICompoundIdentifier): ord.CompoundIdentifier;

        /**
         * Encodes the specified CompoundIdentifier message. Does not implicitly {@link ord.CompoundIdentifier.verify|verify} messages.
         * @param message CompoundIdentifier message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ICompoundIdentifier, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified CompoundIdentifier message, length delimited. Does not implicitly {@link ord.CompoundIdentifier.verify|verify} messages.
         * @param message CompoundIdentifier message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ICompoundIdentifier, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a CompoundIdentifier message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns CompoundIdentifier
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.CompoundIdentifier;

        /**
         * Decodes a CompoundIdentifier message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns CompoundIdentifier
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.CompoundIdentifier;

        /**
         * Verifies a CompoundIdentifier message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a CompoundIdentifier message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns CompoundIdentifier
         */
        public static fromObject(object: { [k: string]: any }): ord.CompoundIdentifier;

        /**
         * Creates a plain object from a CompoundIdentifier message. Also converts values to other types if specified.
         * @param message CompoundIdentifier
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.CompoundIdentifier, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this CompoundIdentifier to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for CompoundIdentifier
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace CompoundIdentifier {

        /** CompoundIdentifierType enum. */
        enum CompoundIdentifierType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            SMILES = 2,
            INCHI = 3,
            MOLBLOCK = 4,
            IUPAC_NAME = 5,
            NAME = 6,
            CAS_NUMBER = 7,
            PUBCHEM_CID = 8,
            CHEMSPIDER_ID = 9,
            CXSMILES = 10,
            INCHI_KEY = 11,
            XYZ = 12,
            UNIPROT_ID = 13,
            PDB_ID = 14,
            AMINO_ACID_SEQUENCE = 15,
            HELM = 16,
            MDL = 17
        }
    }

    /** Properties of a Vessel. */
    interface IVessel {

        /** Vessel type */
        type?: (ord.Vessel.VesselType|null);

        /** Vessel details */
        details?: (string|null);

        /** Vessel material */
        material?: (ord.IVesselMaterial|null);

        /** Vessel preparations */
        preparations?: (ord.IVesselPreparation[]|null);

        /** Vessel attachments */
        attachments?: (ord.IVesselAttachment[]|null);

        /** Vessel volume */
        volume?: (ord.IVolume|null);

        /** Vessel vesselId */
        vesselId?: (string|null);

        /** Vessel position */
        position?: (string|null);

        /** Vessel row */
        row?: (string|null);

        /** Vessel col */
        col?: (string|null);
    }

    /**
     * The Vessel defines the primary container within which the reaction was
     * performed, including the type of vessel, its primary material, any
     * preparation steps or vessel attachments, and its volume.
     */
    class Vessel implements IVessel {

        /**
         * Constructs a new Vessel.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IVessel);

        /** Vessel type. */
        public type: ord.Vessel.VesselType;

        /** Vessel details. */
        public details: string;

        /** Vessel material. */
        public material?: (ord.IVesselMaterial|null);

        /** Vessel preparations. */
        public preparations: ord.IVesselPreparation[];

        /** Vessel attachments. */
        public attachments: ord.IVesselAttachment[];

        /** Vessel volume. */
        public volume?: (ord.IVolume|null);

        /** Vessel vesselId. */
        public vesselId: string;

        /** Vessel position. */
        public position: string;

        /** Vessel row. */
        public row: string;

        /** Vessel col. */
        public col: string;

        /**
         * Creates a new Vessel instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Vessel instance
         */
        public static create(properties?: ord.IVessel): ord.Vessel;

        /**
         * Encodes the specified Vessel message. Does not implicitly {@link ord.Vessel.verify|verify} messages.
         * @param message Vessel message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IVessel, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Vessel message, length delimited. Does not implicitly {@link ord.Vessel.verify|verify} messages.
         * @param message Vessel message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IVessel, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Vessel message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Vessel
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Vessel;

        /**
         * Decodes a Vessel message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Vessel
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Vessel;

        /**
         * Verifies a Vessel message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Vessel message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Vessel
         */
        public static fromObject(object: { [k: string]: any }): ord.Vessel;

        /**
         * Creates a plain object from a Vessel message. Also converts values to other types if specified.
         * @param message Vessel
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Vessel, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Vessel to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Vessel
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Vessel {

        /** VesselType enum. */
        enum VesselType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            ROUND_BOTTOM_FLASK = 2,
            VIAL = 3,
            WELL_PLATE = 4,
            MICROWAVE_VIAL = 5,
            TUBE = 6,
            CONTINUOUS_STIRRED_TANK_REACTOR = 7,
            PACKED_BED_REACTOR = 8,
            NMR_TUBE = 9,
            PRESSURE_FLASK = 10,
            PRESSURE_REACTOR = 11,
            ELECTROCHEMICAL_CELL = 12
        }
    }

    /** Properties of a VesselMaterial. */
    interface IVesselMaterial {

        /** VesselMaterial type */
        type?: (ord.VesselMaterial.VesselMaterialType|null);

        /** VesselMaterial details */
        details?: (string|null);
    }

    /** Represents a VesselMaterial. */
    class VesselMaterial implements IVesselMaterial {

        /**
         * Constructs a new VesselMaterial.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IVesselMaterial);

        /** VesselMaterial type. */
        public type: ord.VesselMaterial.VesselMaterialType;

        /** VesselMaterial details. */
        public details: string;

        /**
         * Creates a new VesselMaterial instance using the specified properties.
         * @param [properties] Properties to set
         * @returns VesselMaterial instance
         */
        public static create(properties?: ord.IVesselMaterial): ord.VesselMaterial;

        /**
         * Encodes the specified VesselMaterial message. Does not implicitly {@link ord.VesselMaterial.verify|verify} messages.
         * @param message VesselMaterial message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IVesselMaterial, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified VesselMaterial message, length delimited. Does not implicitly {@link ord.VesselMaterial.verify|verify} messages.
         * @param message VesselMaterial message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IVesselMaterial, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a VesselMaterial message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns VesselMaterial
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.VesselMaterial;

        /**
         * Decodes a VesselMaterial message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns VesselMaterial
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.VesselMaterial;

        /**
         * Verifies a VesselMaterial message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a VesselMaterial message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns VesselMaterial
         */
        public static fromObject(object: { [k: string]: any }): ord.VesselMaterial;

        /**
         * Creates a plain object from a VesselMaterial message. Also converts values to other types if specified.
         * @param message VesselMaterial
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.VesselMaterial, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this VesselMaterial to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for VesselMaterial
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace VesselMaterial {

        /** VesselMaterialType enum. */
        enum VesselMaterialType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            GLASS = 2,
            POLYPROPYLENE = 3,
            PLASTIC = 4,
            METAL = 5,
            QUARTZ = 6
        }
    }

    /** Properties of a VesselAttachment. */
    interface IVesselAttachment {

        /** VesselAttachment type */
        type?: (ord.VesselAttachment.VesselAttachmentType|null);

        /** VesselAttachment details */
        details?: (string|null);
    }

    /** Represents a VesselAttachment. */
    class VesselAttachment implements IVesselAttachment {

        /**
         * Constructs a new VesselAttachment.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IVesselAttachment);

        /** VesselAttachment type. */
        public type: ord.VesselAttachment.VesselAttachmentType;

        /** VesselAttachment details. */
        public details: string;

        /**
         * Creates a new VesselAttachment instance using the specified properties.
         * @param [properties] Properties to set
         * @returns VesselAttachment instance
         */
        public static create(properties?: ord.IVesselAttachment): ord.VesselAttachment;

        /**
         * Encodes the specified VesselAttachment message. Does not implicitly {@link ord.VesselAttachment.verify|verify} messages.
         * @param message VesselAttachment message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IVesselAttachment, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified VesselAttachment message, length delimited. Does not implicitly {@link ord.VesselAttachment.verify|verify} messages.
         * @param message VesselAttachment message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IVesselAttachment, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a VesselAttachment message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns VesselAttachment
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.VesselAttachment;

        /**
         * Decodes a VesselAttachment message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns VesselAttachment
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.VesselAttachment;

        /**
         * Verifies a VesselAttachment message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a VesselAttachment message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns VesselAttachment
         */
        public static fromObject(object: { [k: string]: any }): ord.VesselAttachment;

        /**
         * Creates a plain object from a VesselAttachment message. Also converts values to other types if specified.
         * @param message VesselAttachment
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.VesselAttachment, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this VesselAttachment to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for VesselAttachment
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace VesselAttachment {

        /** VesselAttachmentType enum. */
        enum VesselAttachmentType {
            UNSPECIFIED = 0,
            NONE = 1,
            CUSTOM = 2,
            SEPTUM = 3,
            CAP = 4,
            MAT = 5,
            REFLUX_CONDENSER = 6,
            VENT_NEEDLE = 7,
            DEAN_STARK = 8,
            VACUUM_TUBE = 9,
            ADDITION_FUNNEL = 10,
            DRYING_TUBE = 11,
            ALUMINUM_FOIL = 12,
            THERMOCOUPLE = 13,
            BALLOON = 14,
            GAS_ADAPTER = 15,
            PRESSURE_REGULATOR = 16,
            RELEASE_VALVE = 17
        }
    }

    /** Properties of a VesselPreparation. */
    interface IVesselPreparation {

        /** VesselPreparation type */
        type?: (ord.VesselPreparation.VesselPreparationType|null);

        /** VesselPreparation details */
        details?: (string|null);
    }

    /** Represents a VesselPreparation. */
    class VesselPreparation implements IVesselPreparation {

        /**
         * Constructs a new VesselPreparation.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IVesselPreparation);

        /** VesselPreparation type. */
        public type: ord.VesselPreparation.VesselPreparationType;

        /** VesselPreparation details. */
        public details: string;

        /**
         * Creates a new VesselPreparation instance using the specified properties.
         * @param [properties] Properties to set
         * @returns VesselPreparation instance
         */
        public static create(properties?: ord.IVesselPreparation): ord.VesselPreparation;

        /**
         * Encodes the specified VesselPreparation message. Does not implicitly {@link ord.VesselPreparation.verify|verify} messages.
         * @param message VesselPreparation message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IVesselPreparation, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified VesselPreparation message, length delimited. Does not implicitly {@link ord.VesselPreparation.verify|verify} messages.
         * @param message VesselPreparation message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IVesselPreparation, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a VesselPreparation message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns VesselPreparation
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.VesselPreparation;

        /**
         * Decodes a VesselPreparation message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns VesselPreparation
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.VesselPreparation;

        /**
         * Verifies a VesselPreparation message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a VesselPreparation message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns VesselPreparation
         */
        public static fromObject(object: { [k: string]: any }): ord.VesselPreparation;

        /**
         * Creates a plain object from a VesselPreparation message. Also converts values to other types if specified.
         * @param message VesselPreparation
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.VesselPreparation, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this VesselPreparation to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for VesselPreparation
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace VesselPreparation {

        /** VesselPreparationType enum. */
        enum VesselPreparationType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            NONE = 2,
            OVEN_DRIED = 3,
            FLAME_DRIED = 4,
            EVACUATED_BACKFILLED = 5,
            PURGED = 6
        }
    }

    /** Properties of a ReactionSetup. */
    interface IReactionSetup {

        /** ReactionSetup vessel */
        vessel?: (ord.IVessel|null);

        /** ReactionSetup isAutomated */
        isAutomated?: (boolean|null);

        /** ReactionSetup automationPlatform */
        automationPlatform?: (string|null);

        /** ReactionSetup automationCode */
        automationCode?: ({ [k: string]: ord.IData }|null);

        /** ReactionSetup environment */
        environment?: (ord.ReactionSetup.IReactionEnvironment|null);
    }

    /** Represents a ReactionSetup. */
    class ReactionSetup implements IReactionSetup {

        /**
         * Constructs a new ReactionSetup.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionSetup);

        /** ReactionSetup vessel. */
        public vessel?: (ord.IVessel|null);

        /** ReactionSetup isAutomated. */
        public isAutomated?: (boolean|null);

        /** ReactionSetup automationPlatform. */
        public automationPlatform: string;

        /** ReactionSetup automationCode. */
        public automationCode: { [k: string]: ord.IData };

        /** ReactionSetup environment. */
        public environment?: (ord.ReactionSetup.IReactionEnvironment|null);

        /**
         * Creates a new ReactionSetup instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionSetup instance
         */
        public static create(properties?: ord.IReactionSetup): ord.ReactionSetup;

        /**
         * Encodes the specified ReactionSetup message. Does not implicitly {@link ord.ReactionSetup.verify|verify} messages.
         * @param message ReactionSetup message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionSetup, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionSetup message, length delimited. Does not implicitly {@link ord.ReactionSetup.verify|verify} messages.
         * @param message ReactionSetup message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionSetup, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionSetup message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionSetup
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionSetup;

        /**
         * Decodes a ReactionSetup message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionSetup
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionSetup;

        /**
         * Verifies a ReactionSetup message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionSetup message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionSetup
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionSetup;

        /**
         * Creates a plain object from a ReactionSetup message. Also converts values to other types if specified.
         * @param message ReactionSetup
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionSetup, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionSetup to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionSetup
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ReactionSetup {

        /** Properties of a ReactionEnvironment. */
        interface IReactionEnvironment {

            /** ReactionEnvironment type */
            type?: (ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType|null);

            /** ReactionEnvironment details */
            details?: (string|null);
        }

        /** Represents a ReactionEnvironment. */
        class ReactionEnvironment implements IReactionEnvironment {

            /**
             * Constructs a new ReactionEnvironment.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ReactionSetup.IReactionEnvironment);

            /** ReactionEnvironment type. */
            public type: ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType;

            /** ReactionEnvironment details. */
            public details: string;

            /**
             * Creates a new ReactionEnvironment instance using the specified properties.
             * @param [properties] Properties to set
             * @returns ReactionEnvironment instance
             */
            public static create(properties?: ord.ReactionSetup.IReactionEnvironment): ord.ReactionSetup.ReactionEnvironment;

            /**
             * Encodes the specified ReactionEnvironment message. Does not implicitly {@link ord.ReactionSetup.ReactionEnvironment.verify|verify} messages.
             * @param message ReactionEnvironment message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ReactionSetup.IReactionEnvironment, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified ReactionEnvironment message, length delimited. Does not implicitly {@link ord.ReactionSetup.ReactionEnvironment.verify|verify} messages.
             * @param message ReactionEnvironment message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ReactionSetup.IReactionEnvironment, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a ReactionEnvironment message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns ReactionEnvironment
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionSetup.ReactionEnvironment;

            /**
             * Decodes a ReactionEnvironment message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns ReactionEnvironment
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionSetup.ReactionEnvironment;

            /**
             * Verifies a ReactionEnvironment message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a ReactionEnvironment message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns ReactionEnvironment
             */
            public static fromObject(object: { [k: string]: any }): ord.ReactionSetup.ReactionEnvironment;

            /**
             * Creates a plain object from a ReactionEnvironment message. Also converts values to other types if specified.
             * @param message ReactionEnvironment
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ReactionSetup.ReactionEnvironment, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this ReactionEnvironment to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for ReactionEnvironment
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace ReactionEnvironment {

            /** ReactionEnvironmentType enum. */
            enum ReactionEnvironmentType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                FUME_HOOD = 2,
                BENCH_TOP = 3,
                GLOVE_BOX = 4,
                GLOVE_BAG = 5
            }
        }
    }

    /** Properties of a ReactionConditions. */
    interface IReactionConditions {

        /** ReactionConditions temperature */
        temperature?: (ord.ITemperatureConditions|null);

        /** ReactionConditions pressure */
        pressure?: (ord.IPressureConditions|null);

        /** ReactionConditions stirring */
        stirring?: (ord.IStirringConditions|null);

        /** ReactionConditions illumination */
        illumination?: (ord.IIlluminationConditions|null);

        /** ReactionConditions electrochemistry */
        electrochemistry?: (ord.IElectrochemistryConditions|null);

        /** ReactionConditions flow */
        flow?: (ord.IFlowConditions|null);

        /** ReactionConditions reflux */
        reflux?: (boolean|null);

        /** ReactionConditions ph */
        ph?: (number|null);

        /** ReactionConditions conditionsAreDynamic */
        conditionsAreDynamic?: (boolean|null);

        /** ReactionConditions details */
        details?: (string|null);
    }

    /** Represents a ReactionConditions. */
    class ReactionConditions implements IReactionConditions {

        /**
         * Constructs a new ReactionConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionConditions);

        /** ReactionConditions temperature. */
        public temperature?: (ord.ITemperatureConditions|null);

        /** ReactionConditions pressure. */
        public pressure?: (ord.IPressureConditions|null);

        /** ReactionConditions stirring. */
        public stirring?: (ord.IStirringConditions|null);

        /** ReactionConditions illumination. */
        public illumination?: (ord.IIlluminationConditions|null);

        /** ReactionConditions electrochemistry. */
        public electrochemistry?: (ord.IElectrochemistryConditions|null);

        /** ReactionConditions flow. */
        public flow?: (ord.IFlowConditions|null);

        /** ReactionConditions reflux. */
        public reflux?: (boolean|null);

        /** ReactionConditions ph. */
        public ph?: (number|null);

        /** ReactionConditions conditionsAreDynamic. */
        public conditionsAreDynamic?: (boolean|null);

        /** ReactionConditions details. */
        public details: string;

        /**
         * Creates a new ReactionConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionConditions instance
         */
        public static create(properties?: ord.IReactionConditions): ord.ReactionConditions;

        /**
         * Encodes the specified ReactionConditions message. Does not implicitly {@link ord.ReactionConditions.verify|verify} messages.
         * @param message ReactionConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionConditions message, length delimited. Does not implicitly {@link ord.ReactionConditions.verify|verify} messages.
         * @param message ReactionConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionConditions;

        /**
         * Decodes a ReactionConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionConditions;

        /**
         * Verifies a ReactionConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionConditions;

        /**
         * Creates a plain object from a ReactionConditions message. Also converts values to other types if specified.
         * @param message ReactionConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a TemperatureConditions. */
    interface ITemperatureConditions {

        /** TemperatureConditions control */
        control?: (ord.TemperatureConditions.ITemperatureControl|null);

        /** TemperatureConditions setpoint */
        setpoint?: (ord.ITemperature|null);

        /** TemperatureConditions measurements */
        measurements?: (ord.TemperatureConditions.ITemperatureMeasurement[]|null);
    }

    /** Represents a TemperatureConditions. */
    class TemperatureConditions implements ITemperatureConditions {

        /**
         * Constructs a new TemperatureConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ITemperatureConditions);

        /** TemperatureConditions control. */
        public control?: (ord.TemperatureConditions.ITemperatureControl|null);

        /** TemperatureConditions setpoint. */
        public setpoint?: (ord.ITemperature|null);

        /** TemperatureConditions measurements. */
        public measurements: ord.TemperatureConditions.ITemperatureMeasurement[];

        /**
         * Creates a new TemperatureConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns TemperatureConditions instance
         */
        public static create(properties?: ord.ITemperatureConditions): ord.TemperatureConditions;

        /**
         * Encodes the specified TemperatureConditions message. Does not implicitly {@link ord.TemperatureConditions.verify|verify} messages.
         * @param message TemperatureConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ITemperatureConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified TemperatureConditions message, length delimited. Does not implicitly {@link ord.TemperatureConditions.verify|verify} messages.
         * @param message TemperatureConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ITemperatureConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a TemperatureConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns TemperatureConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.TemperatureConditions;

        /**
         * Decodes a TemperatureConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns TemperatureConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.TemperatureConditions;

        /**
         * Verifies a TemperatureConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a TemperatureConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns TemperatureConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.TemperatureConditions;

        /**
         * Creates a plain object from a TemperatureConditions message. Also converts values to other types if specified.
         * @param message TemperatureConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.TemperatureConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this TemperatureConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for TemperatureConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace TemperatureConditions {

        /** Properties of a TemperatureControl. */
        interface ITemperatureControl {

            /** TemperatureControl type */
            type?: (ord.TemperatureConditions.TemperatureControl.TemperatureControlType|null);

            /** TemperatureControl details */
            details?: (string|null);
        }

        /** Represents a TemperatureControl. */
        class TemperatureControl implements ITemperatureControl {

            /**
             * Constructs a new TemperatureControl.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.TemperatureConditions.ITemperatureControl);

            /** TemperatureControl type. */
            public type: ord.TemperatureConditions.TemperatureControl.TemperatureControlType;

            /** TemperatureControl details. */
            public details: string;

            /**
             * Creates a new TemperatureControl instance using the specified properties.
             * @param [properties] Properties to set
             * @returns TemperatureControl instance
             */
            public static create(properties?: ord.TemperatureConditions.ITemperatureControl): ord.TemperatureConditions.TemperatureControl;

            /**
             * Encodes the specified TemperatureControl message. Does not implicitly {@link ord.TemperatureConditions.TemperatureControl.verify|verify} messages.
             * @param message TemperatureControl message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.TemperatureConditions.ITemperatureControl, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified TemperatureControl message, length delimited. Does not implicitly {@link ord.TemperatureConditions.TemperatureControl.verify|verify} messages.
             * @param message TemperatureControl message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.TemperatureConditions.ITemperatureControl, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a TemperatureControl message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns TemperatureControl
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.TemperatureConditions.TemperatureControl;

            /**
             * Decodes a TemperatureControl message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns TemperatureControl
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.TemperatureConditions.TemperatureControl;

            /**
             * Verifies a TemperatureControl message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a TemperatureControl message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns TemperatureControl
             */
            public static fromObject(object: { [k: string]: any }): ord.TemperatureConditions.TemperatureControl;

            /**
             * Creates a plain object from a TemperatureControl message. Also converts values to other types if specified.
             * @param message TemperatureControl
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.TemperatureConditions.TemperatureControl, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this TemperatureControl to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for TemperatureControl
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace TemperatureControl {

            /** TemperatureControlType enum. */
            enum TemperatureControlType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                AMBIENT = 2,
                OIL_BATH = 3,
                WATER_BATH = 4,
                SAND_BATH = 5,
                ICE_BATH = 6,
                DRY_ALUMINUM_PLATE = 7,
                MICROWAVE = 8,
                DRY_ICE_BATH = 9,
                AIR_FAN = 10,
                LIQUID_NITROGEN = 11
            }
        }

        /** Properties of a TemperatureMeasurement. */
        interface ITemperatureMeasurement {

            /** TemperatureMeasurement type */
            type?: (ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType|null);

            /** TemperatureMeasurement details */
            details?: (string|null);

            /** TemperatureMeasurement time */
            time?: (ord.ITime|null);

            /** TemperatureMeasurement temperature */
            temperature?: (ord.ITemperature|null);
        }

        /** Represents a TemperatureMeasurement. */
        class TemperatureMeasurement implements ITemperatureMeasurement {

            /**
             * Constructs a new TemperatureMeasurement.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.TemperatureConditions.ITemperatureMeasurement);

            /** TemperatureMeasurement type. */
            public type: ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType;

            /** TemperatureMeasurement details. */
            public details: string;

            /** TemperatureMeasurement time. */
            public time?: (ord.ITime|null);

            /** TemperatureMeasurement temperature. */
            public temperature?: (ord.ITemperature|null);

            /**
             * Creates a new TemperatureMeasurement instance using the specified properties.
             * @param [properties] Properties to set
             * @returns TemperatureMeasurement instance
             */
            public static create(properties?: ord.TemperatureConditions.ITemperatureMeasurement): ord.TemperatureConditions.TemperatureMeasurement;

            /**
             * Encodes the specified TemperatureMeasurement message. Does not implicitly {@link ord.TemperatureConditions.TemperatureMeasurement.verify|verify} messages.
             * @param message TemperatureMeasurement message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.TemperatureConditions.ITemperatureMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified TemperatureMeasurement message, length delimited. Does not implicitly {@link ord.TemperatureConditions.TemperatureMeasurement.verify|verify} messages.
             * @param message TemperatureMeasurement message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.TemperatureConditions.ITemperatureMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a TemperatureMeasurement message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns TemperatureMeasurement
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.TemperatureConditions.TemperatureMeasurement;

            /**
             * Decodes a TemperatureMeasurement message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns TemperatureMeasurement
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.TemperatureConditions.TemperatureMeasurement;

            /**
             * Verifies a TemperatureMeasurement message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a TemperatureMeasurement message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns TemperatureMeasurement
             */
            public static fromObject(object: { [k: string]: any }): ord.TemperatureConditions.TemperatureMeasurement;

            /**
             * Creates a plain object from a TemperatureMeasurement message. Also converts values to other types if specified.
             * @param message TemperatureMeasurement
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.TemperatureConditions.TemperatureMeasurement, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this TemperatureMeasurement to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for TemperatureMeasurement
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace TemperatureMeasurement {

            /** TemperatureMeasurementType enum. */
            enum TemperatureMeasurementType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                THERMOCOUPLE_INTERNAL = 2,
                THERMOCOUPLE_EXTERNAL = 3,
                INFRARED = 4
            }
        }
    }

    /** Properties of a PressureConditions. */
    interface IPressureConditions {

        /** PressureConditions control */
        control?: (ord.PressureConditions.IPressureControl|null);

        /** PressureConditions setpoint */
        setpoint?: (ord.IPressure|null);

        /** PressureConditions atmosphere */
        atmosphere?: (ord.PressureConditions.IAtmosphere|null);

        /** PressureConditions measurements */
        measurements?: (ord.PressureConditions.IPressureMeasurement[]|null);
    }

    /** Represents a PressureConditions. */
    class PressureConditions implements IPressureConditions {

        /**
         * Constructs a new PressureConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IPressureConditions);

        /** PressureConditions control. */
        public control?: (ord.PressureConditions.IPressureControl|null);

        /** PressureConditions setpoint. */
        public setpoint?: (ord.IPressure|null);

        /** PressureConditions atmosphere. */
        public atmosphere?: (ord.PressureConditions.IAtmosphere|null);

        /** PressureConditions measurements. */
        public measurements: ord.PressureConditions.IPressureMeasurement[];

        /**
         * Creates a new PressureConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns PressureConditions instance
         */
        public static create(properties?: ord.IPressureConditions): ord.PressureConditions;

        /**
         * Encodes the specified PressureConditions message. Does not implicitly {@link ord.PressureConditions.verify|verify} messages.
         * @param message PressureConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IPressureConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified PressureConditions message, length delimited. Does not implicitly {@link ord.PressureConditions.verify|verify} messages.
         * @param message PressureConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IPressureConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a PressureConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns PressureConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.PressureConditions;

        /**
         * Decodes a PressureConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns PressureConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.PressureConditions;

        /**
         * Verifies a PressureConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a PressureConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns PressureConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.PressureConditions;

        /**
         * Creates a plain object from a PressureConditions message. Also converts values to other types if specified.
         * @param message PressureConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.PressureConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this PressureConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for PressureConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace PressureConditions {

        /** Properties of a PressureControl. */
        interface IPressureControl {

            /** PressureControl type */
            type?: (ord.PressureConditions.PressureControl.PressureControlType|null);

            /** PressureControl details */
            details?: (string|null);
        }

        /** Represents a PressureControl. */
        class PressureControl implements IPressureControl {

            /**
             * Constructs a new PressureControl.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.PressureConditions.IPressureControl);

            /** PressureControl type. */
            public type: ord.PressureConditions.PressureControl.PressureControlType;

            /** PressureControl details. */
            public details: string;

            /**
             * Creates a new PressureControl instance using the specified properties.
             * @param [properties] Properties to set
             * @returns PressureControl instance
             */
            public static create(properties?: ord.PressureConditions.IPressureControl): ord.PressureConditions.PressureControl;

            /**
             * Encodes the specified PressureControl message. Does not implicitly {@link ord.PressureConditions.PressureControl.verify|verify} messages.
             * @param message PressureControl message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.PressureConditions.IPressureControl, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified PressureControl message, length delimited. Does not implicitly {@link ord.PressureConditions.PressureControl.verify|verify} messages.
             * @param message PressureControl message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.PressureConditions.IPressureControl, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a PressureControl message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns PressureControl
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.PressureConditions.PressureControl;

            /**
             * Decodes a PressureControl message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns PressureControl
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.PressureConditions.PressureControl;

            /**
             * Verifies a PressureControl message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a PressureControl message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns PressureControl
             */
            public static fromObject(object: { [k: string]: any }): ord.PressureConditions.PressureControl;

            /**
             * Creates a plain object from a PressureControl message. Also converts values to other types if specified.
             * @param message PressureControl
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.PressureConditions.PressureControl, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this PressureControl to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for PressureControl
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace PressureControl {

            /** PressureControlType enum. */
            enum PressureControlType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                AMBIENT = 2,
                SLIGHT_POSITIVE = 3,
                SEALED = 4,
                PRESSURIZED = 5
            }
        }

        /** Properties of an Atmosphere. */
        interface IAtmosphere {

            /** Atmosphere type */
            type?: (ord.PressureConditions.Atmosphere.AtmosphereType|null);

            /** Atmosphere details */
            details?: (string|null);
        }

        /** Represents an Atmosphere. */
        class Atmosphere implements IAtmosphere {

            /**
             * Constructs a new Atmosphere.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.PressureConditions.IAtmosphere);

            /** Atmosphere type. */
            public type: ord.PressureConditions.Atmosphere.AtmosphereType;

            /** Atmosphere details. */
            public details: string;

            /**
             * Creates a new Atmosphere instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Atmosphere instance
             */
            public static create(properties?: ord.PressureConditions.IAtmosphere): ord.PressureConditions.Atmosphere;

            /**
             * Encodes the specified Atmosphere message. Does not implicitly {@link ord.PressureConditions.Atmosphere.verify|verify} messages.
             * @param message Atmosphere message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.PressureConditions.IAtmosphere, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Atmosphere message, length delimited. Does not implicitly {@link ord.PressureConditions.Atmosphere.verify|verify} messages.
             * @param message Atmosphere message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.PressureConditions.IAtmosphere, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes an Atmosphere message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Atmosphere
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.PressureConditions.Atmosphere;

            /**
             * Decodes an Atmosphere message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Atmosphere
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.PressureConditions.Atmosphere;

            /**
             * Verifies an Atmosphere message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates an Atmosphere message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Atmosphere
             */
            public static fromObject(object: { [k: string]: any }): ord.PressureConditions.Atmosphere;

            /**
             * Creates a plain object from an Atmosphere message. Also converts values to other types if specified.
             * @param message Atmosphere
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.PressureConditions.Atmosphere, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Atmosphere to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Atmosphere
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace Atmosphere {

            /** AtmosphereType enum. */
            enum AtmosphereType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                AIR = 2,
                NITROGEN = 3,
                ARGON = 4,
                OXYGEN = 5,
                HYDROGEN = 6,
                CARBON_MONOXIDE = 7,
                CARBON_DIOXIDE = 8,
                METHANE = 9,
                AMMONIA = 10,
                OZONE = 11,
                ETHYLENE = 12,
                ACETYLENE = 13
            }
        }

        /** Properties of a PressureMeasurement. */
        interface IPressureMeasurement {

            /** PressureMeasurement type */
            type?: (ord.PressureConditions.PressureMeasurement.PressureMeasurementType|null);

            /** PressureMeasurement details */
            details?: (string|null);

            /** PressureMeasurement time */
            time?: (ord.ITime|null);

            /** PressureMeasurement pressure */
            pressure?: (ord.IPressure|null);
        }

        /** Represents a PressureMeasurement. */
        class PressureMeasurement implements IPressureMeasurement {

            /**
             * Constructs a new PressureMeasurement.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.PressureConditions.IPressureMeasurement);

            /** PressureMeasurement type. */
            public type: ord.PressureConditions.PressureMeasurement.PressureMeasurementType;

            /** PressureMeasurement details. */
            public details: string;

            /** PressureMeasurement time. */
            public time?: (ord.ITime|null);

            /** PressureMeasurement pressure. */
            public pressure?: (ord.IPressure|null);

            /**
             * Creates a new PressureMeasurement instance using the specified properties.
             * @param [properties] Properties to set
             * @returns PressureMeasurement instance
             */
            public static create(properties?: ord.PressureConditions.IPressureMeasurement): ord.PressureConditions.PressureMeasurement;

            /**
             * Encodes the specified PressureMeasurement message. Does not implicitly {@link ord.PressureConditions.PressureMeasurement.verify|verify} messages.
             * @param message PressureMeasurement message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.PressureConditions.IPressureMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified PressureMeasurement message, length delimited. Does not implicitly {@link ord.PressureConditions.PressureMeasurement.verify|verify} messages.
             * @param message PressureMeasurement message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.PressureConditions.IPressureMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a PressureMeasurement message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns PressureMeasurement
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.PressureConditions.PressureMeasurement;

            /**
             * Decodes a PressureMeasurement message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns PressureMeasurement
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.PressureConditions.PressureMeasurement;

            /**
             * Verifies a PressureMeasurement message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a PressureMeasurement message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns PressureMeasurement
             */
            public static fromObject(object: { [k: string]: any }): ord.PressureConditions.PressureMeasurement;

            /**
             * Creates a plain object from a PressureMeasurement message. Also converts values to other types if specified.
             * @param message PressureMeasurement
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.PressureConditions.PressureMeasurement, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this PressureMeasurement to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for PressureMeasurement
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace PressureMeasurement {

            /** PressureMeasurementType enum. */
            enum PressureMeasurementType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                PRESSURE_TRANSDUCER = 2
            }
        }
    }

    /** Properties of a StirringConditions. */
    interface IStirringConditions {

        /** StirringConditions type */
        type?: (ord.StirringConditions.StirringMethodType|null);

        /** StirringConditions details */
        details?: (string|null);

        /** StirringConditions rate */
        rate?: (ord.StirringConditions.IStirringRate|null);
    }

    /** Represents a StirringConditions. */
    class StirringConditions implements IStirringConditions {

        /**
         * Constructs a new StirringConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IStirringConditions);

        /** StirringConditions type. */
        public type: ord.StirringConditions.StirringMethodType;

        /** StirringConditions details. */
        public details: string;

        /** StirringConditions rate. */
        public rate?: (ord.StirringConditions.IStirringRate|null);

        /**
         * Creates a new StirringConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns StirringConditions instance
         */
        public static create(properties?: ord.IStirringConditions): ord.StirringConditions;

        /**
         * Encodes the specified StirringConditions message. Does not implicitly {@link ord.StirringConditions.verify|verify} messages.
         * @param message StirringConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IStirringConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified StirringConditions message, length delimited. Does not implicitly {@link ord.StirringConditions.verify|verify} messages.
         * @param message StirringConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IStirringConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a StirringConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns StirringConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.StirringConditions;

        /**
         * Decodes a StirringConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns StirringConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.StirringConditions;

        /**
         * Verifies a StirringConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a StirringConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns StirringConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.StirringConditions;

        /**
         * Creates a plain object from a StirringConditions message. Also converts values to other types if specified.
         * @param message StirringConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.StirringConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this StirringConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for StirringConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace StirringConditions {

        /** StirringMethodType enum. */
        enum StirringMethodType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            NONE = 2,
            STIR_BAR = 3,
            OVERHEAD_MIXER = 4,
            AGITATION = 5,
            BALL_MILLING = 6,
            SONICATION = 7
        }

        /** Properties of a StirringRate. */
        interface IStirringRate {

            /** StirringRate type */
            type?: (ord.StirringConditions.StirringRate.StirringRateType|null);

            /** StirringRate details */
            details?: (string|null);

            /** StirringRate rpm */
            rpm?: (number|null);
        }

        /** Represents a StirringRate. */
        class StirringRate implements IStirringRate {

            /**
             * Constructs a new StirringRate.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.StirringConditions.IStirringRate);

            /** StirringRate type. */
            public type: ord.StirringConditions.StirringRate.StirringRateType;

            /** StirringRate details. */
            public details: string;

            /** StirringRate rpm. */
            public rpm: number;

            /**
             * Creates a new StirringRate instance using the specified properties.
             * @param [properties] Properties to set
             * @returns StirringRate instance
             */
            public static create(properties?: ord.StirringConditions.IStirringRate): ord.StirringConditions.StirringRate;

            /**
             * Encodes the specified StirringRate message. Does not implicitly {@link ord.StirringConditions.StirringRate.verify|verify} messages.
             * @param message StirringRate message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.StirringConditions.IStirringRate, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified StirringRate message, length delimited. Does not implicitly {@link ord.StirringConditions.StirringRate.verify|verify} messages.
             * @param message StirringRate message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.StirringConditions.IStirringRate, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a StirringRate message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns StirringRate
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.StirringConditions.StirringRate;

            /**
             * Decodes a StirringRate message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns StirringRate
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.StirringConditions.StirringRate;

            /**
             * Verifies a StirringRate message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a StirringRate message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns StirringRate
             */
            public static fromObject(object: { [k: string]: any }): ord.StirringConditions.StirringRate;

            /**
             * Creates a plain object from a StirringRate message. Also converts values to other types if specified.
             * @param message StirringRate
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.StirringConditions.StirringRate, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this StirringRate to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for StirringRate
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace StirringRate {

            /** StirringRateType enum. */
            enum StirringRateType {
                UNSPECIFIED = 0,
                HIGH = 1,
                MEDIUM = 2,
                LOW = 3
            }
        }
    }

    /** Properties of an IlluminationConditions. */
    interface IIlluminationConditions {

        /** IlluminationConditions type */
        type?: (ord.IlluminationConditions.IlluminationType|null);

        /** IlluminationConditions details */
        details?: (string|null);

        /** IlluminationConditions peakWavelength */
        peakWavelength?: (ord.IWavelength|null);

        /** IlluminationConditions color */
        color?: (string|null);

        /** IlluminationConditions distanceToVessel */
        distanceToVessel?: (ord.ILength|null);
    }

    /** Represents an IlluminationConditions. */
    class IlluminationConditions implements IIlluminationConditions {

        /**
         * Constructs a new IlluminationConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IIlluminationConditions);

        /** IlluminationConditions type. */
        public type: ord.IlluminationConditions.IlluminationType;

        /** IlluminationConditions details. */
        public details: string;

        /** IlluminationConditions peakWavelength. */
        public peakWavelength?: (ord.IWavelength|null);

        /** IlluminationConditions color. */
        public color: string;

        /** IlluminationConditions distanceToVessel. */
        public distanceToVessel?: (ord.ILength|null);

        /**
         * Creates a new IlluminationConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns IlluminationConditions instance
         */
        public static create(properties?: ord.IIlluminationConditions): ord.IlluminationConditions;

        /**
         * Encodes the specified IlluminationConditions message. Does not implicitly {@link ord.IlluminationConditions.verify|verify} messages.
         * @param message IlluminationConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IIlluminationConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified IlluminationConditions message, length delimited. Does not implicitly {@link ord.IlluminationConditions.verify|verify} messages.
         * @param message IlluminationConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IIlluminationConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes an IlluminationConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns IlluminationConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.IlluminationConditions;

        /**
         * Decodes an IlluminationConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns IlluminationConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.IlluminationConditions;

        /**
         * Verifies an IlluminationConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates an IlluminationConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns IlluminationConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.IlluminationConditions;

        /**
         * Creates a plain object from an IlluminationConditions message. Also converts values to other types if specified.
         * @param message IlluminationConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.IlluminationConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this IlluminationConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for IlluminationConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace IlluminationConditions {

        /** IlluminationType enum. */
        enum IlluminationType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            AMBIENT = 2,
            DARK = 3,
            LED = 4,
            HALOGEN_LAMP = 5,
            DEUTERIUM_LAMP = 6,
            SOLAR_SIMULATOR = 7,
            BROAD_SPECTRUM = 8
        }
    }

    /** Properties of an ElectrochemistryConditions. */
    interface IElectrochemistryConditions {

        /** ElectrochemistryConditions type */
        type?: (ord.ElectrochemistryConditions.ElectrochemistryType|null);

        /** ElectrochemistryConditions details */
        details?: (string|null);

        /** ElectrochemistryConditions current */
        current?: (ord.ICurrent|null);

        /** ElectrochemistryConditions voltage */
        voltage?: (ord.IVoltage|null);

        /** ElectrochemistryConditions anodeMaterial */
        anodeMaterial?: (string|null);

        /** ElectrochemistryConditions cathodeMaterial */
        cathodeMaterial?: (string|null);

        /** ElectrochemistryConditions electrodeSeparation */
        electrodeSeparation?: (ord.ILength|null);

        /** ElectrochemistryConditions measurements */
        measurements?: (ord.ElectrochemistryConditions.IElectrochemistryMeasurement[]|null);

        /** ElectrochemistryConditions cell */
        cell?: (ord.ElectrochemistryConditions.IElectrochemistryCell|null);
    }

    /** Represents an ElectrochemistryConditions. */
    class ElectrochemistryConditions implements IElectrochemistryConditions {

        /**
         * Constructs a new ElectrochemistryConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IElectrochemistryConditions);

        /** ElectrochemistryConditions type. */
        public type: ord.ElectrochemistryConditions.ElectrochemistryType;

        /** ElectrochemistryConditions details. */
        public details: string;

        /** ElectrochemistryConditions current. */
        public current?: (ord.ICurrent|null);

        /** ElectrochemistryConditions voltage. */
        public voltage?: (ord.IVoltage|null);

        /** ElectrochemistryConditions anodeMaterial. */
        public anodeMaterial: string;

        /** ElectrochemistryConditions cathodeMaterial. */
        public cathodeMaterial: string;

        /** ElectrochemistryConditions electrodeSeparation. */
        public electrodeSeparation?: (ord.ILength|null);

        /** ElectrochemistryConditions measurements. */
        public measurements: ord.ElectrochemistryConditions.IElectrochemistryMeasurement[];

        /** ElectrochemistryConditions cell. */
        public cell?: (ord.ElectrochemistryConditions.IElectrochemistryCell|null);

        /**
         * Creates a new ElectrochemistryConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ElectrochemistryConditions instance
         */
        public static create(properties?: ord.IElectrochemistryConditions): ord.ElectrochemistryConditions;

        /**
         * Encodes the specified ElectrochemistryConditions message. Does not implicitly {@link ord.ElectrochemistryConditions.verify|verify} messages.
         * @param message ElectrochemistryConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IElectrochemistryConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ElectrochemistryConditions message, length delimited. Does not implicitly {@link ord.ElectrochemistryConditions.verify|verify} messages.
         * @param message ElectrochemistryConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IElectrochemistryConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes an ElectrochemistryConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ElectrochemistryConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ElectrochemistryConditions;

        /**
         * Decodes an ElectrochemistryConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ElectrochemistryConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ElectrochemistryConditions;

        /**
         * Verifies an ElectrochemistryConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates an ElectrochemistryConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ElectrochemistryConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.ElectrochemistryConditions;

        /**
         * Creates a plain object from an ElectrochemistryConditions message. Also converts values to other types if specified.
         * @param message ElectrochemistryConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ElectrochemistryConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ElectrochemistryConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ElectrochemistryConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ElectrochemistryConditions {

        /** ElectrochemistryType enum. */
        enum ElectrochemistryType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            CONSTANT_CURRENT = 2,
            CONSTANT_VOLTAGE = 3
        }

        /** Properties of an ElectrochemistryMeasurement. */
        interface IElectrochemistryMeasurement {

            /** ElectrochemistryMeasurement time */
            time?: (ord.ITime|null);

            /** ElectrochemistryMeasurement current */
            current?: (ord.ICurrent|null);

            /** ElectrochemistryMeasurement voltage */
            voltage?: (ord.IVoltage|null);
        }

        /** Represents an ElectrochemistryMeasurement. */
        class ElectrochemistryMeasurement implements IElectrochemistryMeasurement {

            /**
             * Constructs a new ElectrochemistryMeasurement.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ElectrochemistryConditions.IElectrochemistryMeasurement);

            /** ElectrochemistryMeasurement time. */
            public time?: (ord.ITime|null);

            /** ElectrochemistryMeasurement current. */
            public current?: (ord.ICurrent|null);

            /** ElectrochemistryMeasurement voltage. */
            public voltage?: (ord.IVoltage|null);

            /**
             * Creates a new ElectrochemistryMeasurement instance using the specified properties.
             * @param [properties] Properties to set
             * @returns ElectrochemistryMeasurement instance
             */
            public static create(properties?: ord.ElectrochemistryConditions.IElectrochemistryMeasurement): ord.ElectrochemistryConditions.ElectrochemistryMeasurement;

            /**
             * Encodes the specified ElectrochemistryMeasurement message. Does not implicitly {@link ord.ElectrochemistryConditions.ElectrochemistryMeasurement.verify|verify} messages.
             * @param message ElectrochemistryMeasurement message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ElectrochemistryConditions.IElectrochemistryMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified ElectrochemistryMeasurement message, length delimited. Does not implicitly {@link ord.ElectrochemistryConditions.ElectrochemistryMeasurement.verify|verify} messages.
             * @param message ElectrochemistryMeasurement message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ElectrochemistryConditions.IElectrochemistryMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes an ElectrochemistryMeasurement message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns ElectrochemistryMeasurement
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ElectrochemistryConditions.ElectrochemistryMeasurement;

            /**
             * Decodes an ElectrochemistryMeasurement message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns ElectrochemistryMeasurement
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ElectrochemistryConditions.ElectrochemistryMeasurement;

            /**
             * Verifies an ElectrochemistryMeasurement message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates an ElectrochemistryMeasurement message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns ElectrochemistryMeasurement
             */
            public static fromObject(object: { [k: string]: any }): ord.ElectrochemistryConditions.ElectrochemistryMeasurement;

            /**
             * Creates a plain object from an ElectrochemistryMeasurement message. Also converts values to other types if specified.
             * @param message ElectrochemistryMeasurement
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ElectrochemistryConditions.ElectrochemistryMeasurement, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this ElectrochemistryMeasurement to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for ElectrochemistryMeasurement
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        /** Properties of an ElectrochemistryCell. */
        interface IElectrochemistryCell {

            /** ElectrochemistryCell type */
            type?: (ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType|null);

            /** ElectrochemistryCell details */
            details?: (string|null);
        }

        /** Represents an ElectrochemistryCell. */
        class ElectrochemistryCell implements IElectrochemistryCell {

            /**
             * Constructs a new ElectrochemistryCell.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ElectrochemistryConditions.IElectrochemistryCell);

            /** ElectrochemistryCell type. */
            public type: ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType;

            /** ElectrochemistryCell details. */
            public details: string;

            /**
             * Creates a new ElectrochemistryCell instance using the specified properties.
             * @param [properties] Properties to set
             * @returns ElectrochemistryCell instance
             */
            public static create(properties?: ord.ElectrochemistryConditions.IElectrochemistryCell): ord.ElectrochemistryConditions.ElectrochemistryCell;

            /**
             * Encodes the specified ElectrochemistryCell message. Does not implicitly {@link ord.ElectrochemistryConditions.ElectrochemistryCell.verify|verify} messages.
             * @param message ElectrochemistryCell message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ElectrochemistryConditions.IElectrochemistryCell, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified ElectrochemistryCell message, length delimited. Does not implicitly {@link ord.ElectrochemistryConditions.ElectrochemistryCell.verify|verify} messages.
             * @param message ElectrochemistryCell message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ElectrochemistryConditions.IElectrochemistryCell, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes an ElectrochemistryCell message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns ElectrochemistryCell
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ElectrochemistryConditions.ElectrochemistryCell;

            /**
             * Decodes an ElectrochemistryCell message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns ElectrochemistryCell
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ElectrochemistryConditions.ElectrochemistryCell;

            /**
             * Verifies an ElectrochemistryCell message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates an ElectrochemistryCell message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns ElectrochemistryCell
             */
            public static fromObject(object: { [k: string]: any }): ord.ElectrochemistryConditions.ElectrochemistryCell;

            /**
             * Creates a plain object from an ElectrochemistryCell message. Also converts values to other types if specified.
             * @param message ElectrochemistryCell
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ElectrochemistryConditions.ElectrochemistryCell, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this ElectrochemistryCell to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for ElectrochemistryCell
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace ElectrochemistryCell {

            /** ElectrochemistryCellType enum. */
            enum ElectrochemistryCellType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                DIVIDED_CELL = 2,
                UNDIVIDED_CELL = 3
            }
        }
    }

    /** Properties of a FlowConditions. */
    interface IFlowConditions {

        /** FlowConditions type */
        type?: (ord.FlowConditions.FlowType|null);

        /** FlowConditions details */
        details?: (string|null);

        /** FlowConditions pumpType */
        pumpType?: (string|null);

        /** FlowConditions tubing */
        tubing?: (ord.FlowConditions.ITubing|null);
    }

    /** Represents a FlowConditions. */
    class FlowConditions implements IFlowConditions {

        /**
         * Constructs a new FlowConditions.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IFlowConditions);

        /** FlowConditions type. */
        public type: ord.FlowConditions.FlowType;

        /** FlowConditions details. */
        public details: string;

        /** FlowConditions pumpType. */
        public pumpType: string;

        /** FlowConditions tubing. */
        public tubing?: (ord.FlowConditions.ITubing|null);

        /**
         * Creates a new FlowConditions instance using the specified properties.
         * @param [properties] Properties to set
         * @returns FlowConditions instance
         */
        public static create(properties?: ord.IFlowConditions): ord.FlowConditions;

        /**
         * Encodes the specified FlowConditions message. Does not implicitly {@link ord.FlowConditions.verify|verify} messages.
         * @param message FlowConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IFlowConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified FlowConditions message, length delimited. Does not implicitly {@link ord.FlowConditions.verify|verify} messages.
         * @param message FlowConditions message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IFlowConditions, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a FlowConditions message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns FlowConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.FlowConditions;

        /**
         * Decodes a FlowConditions message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns FlowConditions
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.FlowConditions;

        /**
         * Verifies a FlowConditions message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a FlowConditions message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns FlowConditions
         */
        public static fromObject(object: { [k: string]: any }): ord.FlowConditions;

        /**
         * Creates a plain object from a FlowConditions message. Also converts values to other types if specified.
         * @param message FlowConditions
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.FlowConditions, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this FlowConditions to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for FlowConditions
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace FlowConditions {

        /** FlowType enum. */
        enum FlowType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            PLUG_FLOW_REACTOR = 2,
            CONTINUOUS_STIRRED_TANK_REACTOR = 3,
            PACKED_BED_REACTOR = 4
        }

        /** Properties of a Tubing. */
        interface ITubing {

            /** Tubing type */
            type?: (ord.FlowConditions.Tubing.TubingType|null);

            /** Tubing details */
            details?: (string|null);

            /** Tubing diameter */
            diameter?: (ord.ILength|null);
        }

        /** Represents a Tubing. */
        class Tubing implements ITubing {

            /**
             * Constructs a new Tubing.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.FlowConditions.ITubing);

            /** Tubing type. */
            public type: ord.FlowConditions.Tubing.TubingType;

            /** Tubing details. */
            public details: string;

            /** Tubing diameter. */
            public diameter?: (ord.ILength|null);

            /**
             * Creates a new Tubing instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Tubing instance
             */
            public static create(properties?: ord.FlowConditions.ITubing): ord.FlowConditions.Tubing;

            /**
             * Encodes the specified Tubing message. Does not implicitly {@link ord.FlowConditions.Tubing.verify|verify} messages.
             * @param message Tubing message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.FlowConditions.ITubing, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Tubing message, length delimited. Does not implicitly {@link ord.FlowConditions.Tubing.verify|verify} messages.
             * @param message Tubing message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.FlowConditions.ITubing, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a Tubing message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Tubing
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.FlowConditions.Tubing;

            /**
             * Decodes a Tubing message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Tubing
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.FlowConditions.Tubing;

            /**
             * Verifies a Tubing message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a Tubing message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Tubing
             */
            public static fromObject(object: { [k: string]: any }): ord.FlowConditions.Tubing;

            /**
             * Creates a plain object from a Tubing message. Also converts values to other types if specified.
             * @param message Tubing
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.FlowConditions.Tubing, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Tubing to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Tubing
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace Tubing {

            /** TubingType enum. */
            enum TubingType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                STEEL = 2,
                COPPER = 3,
                PFA = 4,
                FEP = 5,
                TEFLONAF = 6,
                PTFE = 7,
                GLASS = 8,
                QUARTZ = 9,
                SILICON = 10,
                PDMS = 11
            }
        }
    }

    /** Properties of a ReactionNotes. */
    interface IReactionNotes {

        /** ReactionNotes isHeterogeneous */
        isHeterogeneous?: (boolean|null);

        /** ReactionNotes formsPrecipitate */
        formsPrecipitate?: (boolean|null);

        /** ReactionNotes isExothermic */
        isExothermic?: (boolean|null);

        /** ReactionNotes offgasses */
        offgasses?: (boolean|null);

        /** ReactionNotes isSensitiveToMoisture */
        isSensitiveToMoisture?: (boolean|null);

        /** ReactionNotes isSensitiveToOxygen */
        isSensitiveToOxygen?: (boolean|null);

        /** ReactionNotes isSensitiveToLight */
        isSensitiveToLight?: (boolean|null);

        /** ReactionNotes safetyNotes */
        safetyNotes?: (string|null);

        /** ReactionNotes procedureDetails */
        procedureDetails?: (string|null);
    }

    /** Represents a ReactionNotes. */
    class ReactionNotes implements IReactionNotes {

        /**
         * Constructs a new ReactionNotes.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionNotes);

        /** ReactionNotes isHeterogeneous. */
        public isHeterogeneous?: (boolean|null);

        /** ReactionNotes formsPrecipitate. */
        public formsPrecipitate?: (boolean|null);

        /** ReactionNotes isExothermic. */
        public isExothermic?: (boolean|null);

        /** ReactionNotes offgasses. */
        public offgasses?: (boolean|null);

        /** ReactionNotes isSensitiveToMoisture. */
        public isSensitiveToMoisture?: (boolean|null);

        /** ReactionNotes isSensitiveToOxygen. */
        public isSensitiveToOxygen?: (boolean|null);

        /** ReactionNotes isSensitiveToLight. */
        public isSensitiveToLight?: (boolean|null);

        /** ReactionNotes safetyNotes. */
        public safetyNotes: string;

        /** ReactionNotes procedureDetails. */
        public procedureDetails: string;

        /**
         * Creates a new ReactionNotes instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionNotes instance
         */
        public static create(properties?: ord.IReactionNotes): ord.ReactionNotes;

        /**
         * Encodes the specified ReactionNotes message. Does not implicitly {@link ord.ReactionNotes.verify|verify} messages.
         * @param message ReactionNotes message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionNotes, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionNotes message, length delimited. Does not implicitly {@link ord.ReactionNotes.verify|verify} messages.
         * @param message ReactionNotes message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionNotes, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionNotes message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionNotes
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionNotes;

        /**
         * Decodes a ReactionNotes message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionNotes
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionNotes;

        /**
         * Verifies a ReactionNotes message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionNotes message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionNotes
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionNotes;

        /**
         * Creates a plain object from a ReactionNotes message. Also converts values to other types if specified.
         * @param message ReactionNotes
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionNotes, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionNotes to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionNotes
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a ReactionObservation. */
    interface IReactionObservation {

        /** ReactionObservation time */
        time?: (ord.ITime|null);

        /** ReactionObservation comment */
        comment?: (string|null);

        /** ReactionObservation image */
        image?: (ord.IData|null);
    }

    /** Represents a ReactionObservation. */
    class ReactionObservation implements IReactionObservation {

        /**
         * Constructs a new ReactionObservation.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionObservation);

        /** ReactionObservation time. */
        public time?: (ord.ITime|null);

        /** ReactionObservation comment. */
        public comment: string;

        /** ReactionObservation image. */
        public image?: (ord.IData|null);

        /**
         * Creates a new ReactionObservation instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionObservation instance
         */
        public static create(properties?: ord.IReactionObservation): ord.ReactionObservation;

        /**
         * Encodes the specified ReactionObservation message. Does not implicitly {@link ord.ReactionObservation.verify|verify} messages.
         * @param message ReactionObservation message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionObservation, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionObservation message, length delimited. Does not implicitly {@link ord.ReactionObservation.verify|verify} messages.
         * @param message ReactionObservation message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionObservation, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionObservation message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionObservation
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionObservation;

        /**
         * Decodes a ReactionObservation message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionObservation
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionObservation;

        /**
         * Verifies a ReactionObservation message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionObservation message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionObservation
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionObservation;

        /**
         * Creates a plain object from a ReactionObservation message. Also converts values to other types if specified.
         * @param message ReactionObservation
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionObservation, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionObservation to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionObservation
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a ReactionWorkup. */
    interface IReactionWorkup {

        /** ReactionWorkup type */
        type?: (ord.ReactionWorkup.ReactionWorkupType|null);

        /** ReactionWorkup details */
        details?: (string|null);

        /** ReactionWorkup duration */
        duration?: (ord.ITime|null);

        /** ReactionWorkup input */
        input?: (ord.IReactionInput|null);

        /** ReactionWorkup amount */
        amount?: (ord.IAmount|null);

        /** ReactionWorkup temperature */
        temperature?: (ord.ITemperatureConditions|null);

        /** ReactionWorkup keepPhase */
        keepPhase?: (string|null);

        /** ReactionWorkup stirring */
        stirring?: (ord.IStirringConditions|null);

        /** ReactionWorkup targetPh */
        targetPh?: (number|null);

        /** ReactionWorkup isAutomated */
        isAutomated?: (boolean|null);
    }

    /** Represents a ReactionWorkup. */
    class ReactionWorkup implements IReactionWorkup {

        /**
         * Constructs a new ReactionWorkup.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionWorkup);

        /** ReactionWorkup type. */
        public type: ord.ReactionWorkup.ReactionWorkupType;

        /** ReactionWorkup details. */
        public details: string;

        /** ReactionWorkup duration. */
        public duration?: (ord.ITime|null);

        /** ReactionWorkup input. */
        public input?: (ord.IReactionInput|null);

        /** ReactionWorkup amount. */
        public amount?: (ord.IAmount|null);

        /** ReactionWorkup temperature. */
        public temperature?: (ord.ITemperatureConditions|null);

        /** ReactionWorkup keepPhase. */
        public keepPhase: string;

        /** ReactionWorkup stirring. */
        public stirring?: (ord.IStirringConditions|null);

        /** ReactionWorkup targetPh. */
        public targetPh?: (number|null);

        /** ReactionWorkup isAutomated. */
        public isAutomated?: (boolean|null);

        /**
         * Creates a new ReactionWorkup instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionWorkup instance
         */
        public static create(properties?: ord.IReactionWorkup): ord.ReactionWorkup;

        /**
         * Encodes the specified ReactionWorkup message. Does not implicitly {@link ord.ReactionWorkup.verify|verify} messages.
         * @param message ReactionWorkup message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionWorkup, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionWorkup message, length delimited. Does not implicitly {@link ord.ReactionWorkup.verify|verify} messages.
         * @param message ReactionWorkup message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionWorkup, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionWorkup message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionWorkup
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionWorkup;

        /**
         * Decodes a ReactionWorkup message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionWorkup
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionWorkup;

        /**
         * Verifies a ReactionWorkup message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionWorkup message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionWorkup
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionWorkup;

        /**
         * Creates a plain object from a ReactionWorkup message. Also converts values to other types if specified.
         * @param message ReactionWorkup
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionWorkup, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionWorkup to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionWorkup
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ReactionWorkup {

        /** ReactionWorkupType enum. */
        enum ReactionWorkupType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            ADDITION = 2,
            ALIQUOT = 3,
            TEMPERATURE = 4,
            CONCENTRATION = 5,
            EXTRACTION = 6,
            FILTRATION = 7,
            WASH = 8,
            DRY_IN_VACUUM = 9,
            DRY_WITH_MATERIAL = 10,
            FLASH_CHROMATOGRAPHY = 11,
            OTHER_CHROMATOGRAPHY = 12,
            SCAVENGING = 13,
            WAIT = 14,
            STIRRING = 15,
            PH_ADJUST = 16,
            DISSOLUTION = 17,
            DISTILLATION = 18
        }
    }

    /** Properties of a ReactionOutcome. */
    interface IReactionOutcome {

        /** ReactionOutcome reactionTime */
        reactionTime?: (ord.ITime|null);

        /** ReactionOutcome conversion */
        conversion?: (ord.IPercentage|null);

        /** ReactionOutcome products */
        products?: (ord.IProductCompound[]|null);

        /** ReactionOutcome analyses */
        analyses?: ({ [k: string]: ord.IAnalysis }|null);
    }

    /**
     * The outcomes of a reaction describe the conversion, yield, and/or other
     * analyses of the resulting product mixture after workup step(s). Each
     * outcome is associated with a reaction/residence time. To allow for
     * one Reaction message to contain the results of a full kinetic profiling
     * experiment, this is a repeated field of the Reaction message.
     *
     * It is the parent message for product characterization and any analytical
     * data.
     */
    class ReactionOutcome implements IReactionOutcome {

        /**
         * Constructs a new ReactionOutcome.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionOutcome);

        /** ReactionOutcome reactionTime. */
        public reactionTime?: (ord.ITime|null);

        /** ReactionOutcome conversion. */
        public conversion?: (ord.IPercentage|null);

        /** ReactionOutcome products. */
        public products: ord.IProductCompound[];

        /** ReactionOutcome analyses. */
        public analyses: { [k: string]: ord.IAnalysis };

        /**
         * Creates a new ReactionOutcome instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionOutcome instance
         */
        public static create(properties?: ord.IReactionOutcome): ord.ReactionOutcome;

        /**
         * Encodes the specified ReactionOutcome message. Does not implicitly {@link ord.ReactionOutcome.verify|verify} messages.
         * @param message ReactionOutcome message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionOutcome, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionOutcome message, length delimited. Does not implicitly {@link ord.ReactionOutcome.verify|verify} messages.
         * @param message ReactionOutcome message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionOutcome, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionOutcome message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionOutcome
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionOutcome;

        /**
         * Decodes a ReactionOutcome message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionOutcome
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionOutcome;

        /**
         * Verifies a ReactionOutcome message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionOutcome message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionOutcome
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionOutcome;

        /**
         * Creates a plain object from a ReactionOutcome message. Also converts values to other types if specified.
         * @param message ReactionOutcome
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionOutcome, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionOutcome to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionOutcome
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a ProductCompound. */
    interface IProductCompound {

        /** ProductCompound identifiers */
        identifiers?: (ord.ICompoundIdentifier[]|null);

        /** ProductCompound isDesiredProduct */
        isDesiredProduct?: (boolean|null);

        /** ProductCompound measurements */
        measurements?: (ord.IProductMeasurement[]|null);

        /** ProductCompound isolatedColor */
        isolatedColor?: (string|null);

        /** ProductCompound texture */
        texture?: (ord.ITexture|null);

        /** ProductCompound features */
        features?: ({ [k: string]: ord.IData }|null);

        /** ProductCompound reactionRole */
        reactionRole?: (ord.ReactionRole.ReactionRoleType|null);
    }

    /** Represents a ProductCompound. */
    class ProductCompound implements IProductCompound {

        /**
         * Constructs a new ProductCompound.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IProductCompound);

        /** ProductCompound identifiers. */
        public identifiers: ord.ICompoundIdentifier[];

        /** ProductCompound isDesiredProduct. */
        public isDesiredProduct?: (boolean|null);

        /** ProductCompound measurements. */
        public measurements: ord.IProductMeasurement[];

        /** ProductCompound isolatedColor. */
        public isolatedColor: string;

        /** ProductCompound texture. */
        public texture?: (ord.ITexture|null);

        /** ProductCompound features. */
        public features: { [k: string]: ord.IData };

        /** ProductCompound reactionRole. */
        public reactionRole: ord.ReactionRole.ReactionRoleType;

        /**
         * Creates a new ProductCompound instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ProductCompound instance
         */
        public static create(properties?: ord.IProductCompound): ord.ProductCompound;

        /**
         * Encodes the specified ProductCompound message. Does not implicitly {@link ord.ProductCompound.verify|verify} messages.
         * @param message ProductCompound message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IProductCompound, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ProductCompound message, length delimited. Does not implicitly {@link ord.ProductCompound.verify|verify} messages.
         * @param message ProductCompound message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IProductCompound, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ProductCompound message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ProductCompound
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ProductCompound;

        /**
         * Decodes a ProductCompound message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ProductCompound
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ProductCompound;

        /**
         * Verifies a ProductCompound message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ProductCompound message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ProductCompound
         */
        public static fromObject(object: { [k: string]: any }): ord.ProductCompound;

        /**
         * Creates a plain object from a ProductCompound message. Also converts values to other types if specified.
         * @param message ProductCompound
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ProductCompound, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ProductCompound to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ProductCompound
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a ProductMeasurement. */
    interface IProductMeasurement {

        /** ProductMeasurement analysisKey */
        analysisKey?: (string|null);

        /** ProductMeasurement type */
        type?: (ord.ProductMeasurement.ProductMeasurementType|null);

        /** ProductMeasurement details */
        details?: (string|null);

        /** ProductMeasurement usesInternalStandard */
        usesInternalStandard?: (boolean|null);

        /** ProductMeasurement isNormalized */
        isNormalized?: (boolean|null);

        /** ProductMeasurement usesAuthenticStandard */
        usesAuthenticStandard?: (boolean|null);

        /** ProductMeasurement authenticStandard */
        authenticStandard?: (ord.ICompound|null);

        /** ProductMeasurement percentage */
        percentage?: (ord.IPercentage|null);

        /** ProductMeasurement floatValue */
        floatValue?: (ord.IFloatValue|null);

        /** ProductMeasurement stringValue */
        stringValue?: (string|null);

        /** ProductMeasurement amount */
        amount?: (ord.IAmount|null);

        /** ProductMeasurement retentionTime */
        retentionTime?: (ord.ITime|null);

        /** ProductMeasurement massSpecDetails */
        massSpecDetails?: (ord.ProductMeasurement.IMassSpecMeasurementDetails|null);

        /** ProductMeasurement selectivity */
        selectivity?: (ord.ProductMeasurement.ISelectivity|null);

        /** ProductMeasurement wavelength */
        wavelength?: (ord.IWavelength|null);
    }

    /** Represents a ProductMeasurement. */
    class ProductMeasurement implements IProductMeasurement {

        /**
         * Constructs a new ProductMeasurement.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IProductMeasurement);

        /** ProductMeasurement analysisKey. */
        public analysisKey: string;

        /** ProductMeasurement type. */
        public type: ord.ProductMeasurement.ProductMeasurementType;

        /** ProductMeasurement details. */
        public details: string;

        /** ProductMeasurement usesInternalStandard. */
        public usesInternalStandard?: (boolean|null);

        /** ProductMeasurement isNormalized. */
        public isNormalized?: (boolean|null);

        /** ProductMeasurement usesAuthenticStandard. */
        public usesAuthenticStandard?: (boolean|null);

        /** ProductMeasurement authenticStandard. */
        public authenticStandard?: (ord.ICompound|null);

        /** ProductMeasurement percentage. */
        public percentage?: (ord.IPercentage|null);

        /** ProductMeasurement floatValue. */
        public floatValue?: (ord.IFloatValue|null);

        /** ProductMeasurement stringValue. */
        public stringValue?: (string|null);

        /** ProductMeasurement amount. */
        public amount?: (ord.IAmount|null);

        /** ProductMeasurement retentionTime. */
        public retentionTime?: (ord.ITime|null);

        /** ProductMeasurement massSpecDetails. */
        public massSpecDetails?: (ord.ProductMeasurement.IMassSpecMeasurementDetails|null);

        /** ProductMeasurement selectivity. */
        public selectivity?: (ord.ProductMeasurement.ISelectivity|null);

        /** ProductMeasurement wavelength. */
        public wavelength?: (ord.IWavelength|null);

        /** ProductMeasurement value. */
        public value?: ("percentage"|"floatValue"|"stringValue"|"amount");

        /**
         * Creates a new ProductMeasurement instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ProductMeasurement instance
         */
        public static create(properties?: ord.IProductMeasurement): ord.ProductMeasurement;

        /**
         * Encodes the specified ProductMeasurement message. Does not implicitly {@link ord.ProductMeasurement.verify|verify} messages.
         * @param message ProductMeasurement message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IProductMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ProductMeasurement message, length delimited. Does not implicitly {@link ord.ProductMeasurement.verify|verify} messages.
         * @param message ProductMeasurement message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IProductMeasurement, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ProductMeasurement message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ProductMeasurement
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ProductMeasurement;

        /**
         * Decodes a ProductMeasurement message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ProductMeasurement
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ProductMeasurement;

        /**
         * Verifies a ProductMeasurement message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ProductMeasurement message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ProductMeasurement
         */
        public static fromObject(object: { [k: string]: any }): ord.ProductMeasurement;

        /**
         * Creates a plain object from a ProductMeasurement message. Also converts values to other types if specified.
         * @param message ProductMeasurement
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ProductMeasurement, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ProductMeasurement to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ProductMeasurement
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace ProductMeasurement {

        /** ProductMeasurementType enum. */
        enum ProductMeasurementType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            IDENTITY = 2,
            YIELD = 3,
            SELECTIVITY = 4,
            PURITY = 5,
            AREA = 6,
            COUNTS = 7,
            INTENSITY = 8,
            AMOUNT = 9
        }

        /** Properties of a MassSpecMeasurementDetails. */
        interface IMassSpecMeasurementDetails {

            /** MassSpecMeasurementDetails type */
            type?: (ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType|null);

            /** MassSpecMeasurementDetails details */
            details?: (string|null);

            /** MassSpecMeasurementDetails ticMinimumMz */
            ticMinimumMz?: (number|null);

            /** MassSpecMeasurementDetails ticMaximumMz */
            ticMaximumMz?: (number|null);

            /** MassSpecMeasurementDetails eicMasses */
            eicMasses?: (number[]|null);
        }

        /** Represents a MassSpecMeasurementDetails. */
        class MassSpecMeasurementDetails implements IMassSpecMeasurementDetails {

            /**
             * Constructs a new MassSpecMeasurementDetails.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ProductMeasurement.IMassSpecMeasurementDetails);

            /** MassSpecMeasurementDetails type. */
            public type: ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType;

            /** MassSpecMeasurementDetails details. */
            public details: string;

            /** MassSpecMeasurementDetails ticMinimumMz. */
            public ticMinimumMz?: (number|null);

            /** MassSpecMeasurementDetails ticMaximumMz. */
            public ticMaximumMz?: (number|null);

            /** MassSpecMeasurementDetails eicMasses. */
            public eicMasses: number[];

            /**
             * Creates a new MassSpecMeasurementDetails instance using the specified properties.
             * @param [properties] Properties to set
             * @returns MassSpecMeasurementDetails instance
             */
            public static create(properties?: ord.ProductMeasurement.IMassSpecMeasurementDetails): ord.ProductMeasurement.MassSpecMeasurementDetails;

            /**
             * Encodes the specified MassSpecMeasurementDetails message. Does not implicitly {@link ord.ProductMeasurement.MassSpecMeasurementDetails.verify|verify} messages.
             * @param message MassSpecMeasurementDetails message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ProductMeasurement.IMassSpecMeasurementDetails, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified MassSpecMeasurementDetails message, length delimited. Does not implicitly {@link ord.ProductMeasurement.MassSpecMeasurementDetails.verify|verify} messages.
             * @param message MassSpecMeasurementDetails message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ProductMeasurement.IMassSpecMeasurementDetails, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a MassSpecMeasurementDetails message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns MassSpecMeasurementDetails
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ProductMeasurement.MassSpecMeasurementDetails;

            /**
             * Decodes a MassSpecMeasurementDetails message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns MassSpecMeasurementDetails
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ProductMeasurement.MassSpecMeasurementDetails;

            /**
             * Verifies a MassSpecMeasurementDetails message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a MassSpecMeasurementDetails message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns MassSpecMeasurementDetails
             */
            public static fromObject(object: { [k: string]: any }): ord.ProductMeasurement.MassSpecMeasurementDetails;

            /**
             * Creates a plain object from a MassSpecMeasurementDetails message. Also converts values to other types if specified.
             * @param message MassSpecMeasurementDetails
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ProductMeasurement.MassSpecMeasurementDetails, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this MassSpecMeasurementDetails to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for MassSpecMeasurementDetails
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace MassSpecMeasurementDetails {

            /** MassSpecMeasurementType enum. */
            enum MassSpecMeasurementType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                TIC = 2,
                TIC_POSITIVE = 3,
                TIC_NEGATIVE = 4,
                EIC = 5
            }
        }

        /** Properties of a Selectivity. */
        interface ISelectivity {

            /** Selectivity type */
            type?: (ord.ProductMeasurement.Selectivity.SelectivityType|null);

            /** Selectivity details */
            details?: (string|null);
        }

        /** Represents a Selectivity. */
        class Selectivity implements ISelectivity {

            /**
             * Constructs a new Selectivity.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord.ProductMeasurement.ISelectivity);

            /** Selectivity type. */
            public type: ord.ProductMeasurement.Selectivity.SelectivityType;

            /** Selectivity details. */
            public details: string;

            /**
             * Creates a new Selectivity instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Selectivity instance
             */
            public static create(properties?: ord.ProductMeasurement.ISelectivity): ord.ProductMeasurement.Selectivity;

            /**
             * Encodes the specified Selectivity message. Does not implicitly {@link ord.ProductMeasurement.Selectivity.verify|verify} messages.
             * @param message Selectivity message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord.ProductMeasurement.ISelectivity, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Selectivity message, length delimited. Does not implicitly {@link ord.ProductMeasurement.Selectivity.verify|verify} messages.
             * @param message Selectivity message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord.ProductMeasurement.ISelectivity, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a Selectivity message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Selectivity
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ProductMeasurement.Selectivity;

            /**
             * Decodes a Selectivity message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Selectivity
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ProductMeasurement.Selectivity;

            /**
             * Verifies a Selectivity message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a Selectivity message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Selectivity
             */
            public static fromObject(object: { [k: string]: any }): ord.ProductMeasurement.Selectivity;

            /**
             * Creates a plain object from a Selectivity message. Also converts values to other types if specified.
             * @param message Selectivity
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord.ProductMeasurement.Selectivity, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Selectivity to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Selectivity
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }

        namespace Selectivity {

            /** SelectivityType enum. */
            enum SelectivityType {
                UNSPECIFIED = 0,
                CUSTOM = 1,
                EE = 2,
                ER = 3,
                DR = 4,
                EZ = 5,
                ZE = 6
            }
        }
    }

    /** Properties of a DateTime. */
    interface IDateTime {

        /** DateTime value */
        value?: (string|null);
    }

    /** Represents a DateTime. */
    class DateTime implements IDateTime {

        /**
         * Constructs a new DateTime.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IDateTime);

        /** DateTime value. */
        public value: string;

        /**
         * Creates a new DateTime instance using the specified properties.
         * @param [properties] Properties to set
         * @returns DateTime instance
         */
        public static create(properties?: ord.IDateTime): ord.DateTime;

        /**
         * Encodes the specified DateTime message. Does not implicitly {@link ord.DateTime.verify|verify} messages.
         * @param message DateTime message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IDateTime, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified DateTime message, length delimited. Does not implicitly {@link ord.DateTime.verify|verify} messages.
         * @param message DateTime message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IDateTime, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a DateTime message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns DateTime
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.DateTime;

        /**
         * Decodes a DateTime message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns DateTime
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.DateTime;

        /**
         * Verifies a DateTime message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a DateTime message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns DateTime
         */
        public static fromObject(object: { [k: string]: any }): ord.DateTime;

        /**
         * Creates a plain object from a DateTime message. Also converts values to other types if specified.
         * @param message DateTime
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.DateTime, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this DateTime to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for DateTime
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of an Analysis. */
    interface IAnalysis {

        /** Analysis type */
        type?: (ord.Analysis.AnalysisType|null);

        /** Analysis details */
        details?: (string|null);

        /** Analysis chmoId */
        chmoId?: (number|null);

        /** Analysis isOfIsolatedSpecies */
        isOfIsolatedSpecies?: (boolean|null);

        /** Analysis data */
        data?: ({ [k: string]: ord.IData }|null);

        /** Analysis instrumentManufacturer */
        instrumentManufacturer?: (string|null);

        /** Analysis instrumentLastCalibrated */
        instrumentLastCalibrated?: (ord.IDateTime|null);
    }

    /** Represents an Analysis. */
    class Analysis implements IAnalysis {

        /**
         * Constructs a new Analysis.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IAnalysis);

        /** Analysis type. */
        public type: ord.Analysis.AnalysisType;

        /** Analysis details. */
        public details: string;

        /** Analysis chmoId. */
        public chmoId: number;

        /** Analysis isOfIsolatedSpecies. */
        public isOfIsolatedSpecies?: (boolean|null);

        /** Analysis data. */
        public data: { [k: string]: ord.IData };

        /** Analysis instrumentManufacturer. */
        public instrumentManufacturer: string;

        /** Analysis instrumentLastCalibrated. */
        public instrumentLastCalibrated?: (ord.IDateTime|null);

        /**
         * Creates a new Analysis instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Analysis instance
         */
        public static create(properties?: ord.IAnalysis): ord.Analysis;

        /**
         * Encodes the specified Analysis message. Does not implicitly {@link ord.Analysis.verify|verify} messages.
         * @param message Analysis message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IAnalysis, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Analysis message, length delimited. Does not implicitly {@link ord.Analysis.verify|verify} messages.
         * @param message Analysis message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IAnalysis, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes an Analysis message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Analysis
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Analysis;

        /**
         * Decodes an Analysis message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Analysis
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Analysis;

        /**
         * Verifies an Analysis message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates an Analysis message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Analysis
         */
        public static fromObject(object: { [k: string]: any }): ord.Analysis;

        /**
         * Creates a plain object from an Analysis message. Also converts values to other types if specified.
         * @param message Analysis
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Analysis, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Analysis to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Analysis
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Analysis {

        /** AnalysisType enum. */
        enum AnalysisType {
            UNSPECIFIED = 0,
            CUSTOM = 1,
            LC = 2,
            GC = 3,
            IR = 4,
            NMR_1H = 5,
            NMR_13C = 6,
            NMR_OTHER = 7,
            MP = 8,
            UV = 9,
            TLC = 10,
            MS = 11,
            HRMS = 12,
            MSMS = 13,
            WEIGHT = 14,
            LCMS = 15,
            GCMS = 16,
            ELSD = 17,
            CD = 18,
            SFC = 19,
            EPR = 20,
            XRD = 21,
            RAMAN = 22,
            ED = 23,
            OPTICAL_ROTATION = 24,
            CAD = 25
        }
    }

    /** Properties of a ReactionProvenance. */
    interface IReactionProvenance {

        /** ReactionProvenance experimenter */
        experimenter?: (ord.IPerson|null);

        /** ReactionProvenance city */
        city?: (string|null);

        /** ReactionProvenance experimentStart */
        experimentStart?: (ord.IDateTime|null);

        /** ReactionProvenance doi */
        doi?: (string|null);

        /** ReactionProvenance patent */
        patent?: (string|null);

        /** ReactionProvenance publicationUrl */
        publicationUrl?: (string|null);

        /** ReactionProvenance recordCreated */
        recordCreated?: (ord.IRecordEvent|null);

        /** ReactionProvenance recordModified */
        recordModified?: (ord.IRecordEvent[]|null);

        /** ReactionProvenance reactionMetadata */
        reactionMetadata?: ({ [k: string]: ord.IData }|null);

        /** ReactionProvenance isMined */
        isMined?: (boolean|null);
    }

    /** Represents a ReactionProvenance. */
    class ReactionProvenance implements IReactionProvenance {

        /**
         * Constructs a new ReactionProvenance.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IReactionProvenance);

        /** ReactionProvenance experimenter. */
        public experimenter?: (ord.IPerson|null);

        /** ReactionProvenance city. */
        public city: string;

        /** ReactionProvenance experimentStart. */
        public experimentStart?: (ord.IDateTime|null);

        /** ReactionProvenance doi. */
        public doi: string;

        /** ReactionProvenance patent. */
        public patent: string;

        /** ReactionProvenance publicationUrl. */
        public publicationUrl: string;

        /** ReactionProvenance recordCreated. */
        public recordCreated?: (ord.IRecordEvent|null);

        /** ReactionProvenance recordModified. */
        public recordModified: ord.IRecordEvent[];

        /** ReactionProvenance reactionMetadata. */
        public reactionMetadata: { [k: string]: ord.IData };

        /** ReactionProvenance isMined. */
        public isMined?: (boolean|null);

        /**
         * Creates a new ReactionProvenance instance using the specified properties.
         * @param [properties] Properties to set
         * @returns ReactionProvenance instance
         */
        public static create(properties?: ord.IReactionProvenance): ord.ReactionProvenance;

        /**
         * Encodes the specified ReactionProvenance message. Does not implicitly {@link ord.ReactionProvenance.verify|verify} messages.
         * @param message ReactionProvenance message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IReactionProvenance, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified ReactionProvenance message, length delimited. Does not implicitly {@link ord.ReactionProvenance.verify|verify} messages.
         * @param message ReactionProvenance message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IReactionProvenance, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a ReactionProvenance message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns ReactionProvenance
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.ReactionProvenance;

        /**
         * Decodes a ReactionProvenance message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns ReactionProvenance
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.ReactionProvenance;

        /**
         * Verifies a ReactionProvenance message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a ReactionProvenance message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns ReactionProvenance
         */
        public static fromObject(object: { [k: string]: any }): ord.ReactionProvenance;

        /**
         * Creates a plain object from a ReactionProvenance message. Also converts values to other types if specified.
         * @param message ReactionProvenance
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.ReactionProvenance, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this ReactionProvenance to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for ReactionProvenance
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a Person. */
    interface IPerson {

        /** Person username */
        username?: (string|null);

        /** Person name */
        name?: (string|null);

        /** Person orcid */
        orcid?: (string|null);

        /** Person organization */
        organization?: (string|null);

        /** Person email */
        email?: (string|null);
    }

    /** Represents a Person. */
    class Person implements IPerson {

        /**
         * Constructs a new Person.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IPerson);

        /** Person username. */
        public username: string;

        /** Person name. */
        public name: string;

        /** Person orcid. */
        public orcid: string;

        /** Person organization. */
        public organization: string;

        /** Person email. */
        public email: string;

        /**
         * Creates a new Person instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Person instance
         */
        public static create(properties?: ord.IPerson): ord.Person;

        /**
         * Encodes the specified Person message. Does not implicitly {@link ord.Person.verify|verify} messages.
         * @param message Person message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IPerson, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Person message, length delimited. Does not implicitly {@link ord.Person.verify|verify} messages.
         * @param message Person message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IPerson, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Person message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Person
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Person;

        /**
         * Decodes a Person message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Person
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Person;

        /**
         * Verifies a Person message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Person message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Person
         */
        public static fromObject(object: { [k: string]: any }): ord.Person;

        /**
         * Creates a plain object from a Person message. Also converts values to other types if specified.
         * @param message Person
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Person, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Person to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Person
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a RecordEvent. */
    interface IRecordEvent {

        /** RecordEvent time */
        time?: (ord.IDateTime|null);

        /** RecordEvent person */
        person?: (ord.IPerson|null);

        /** RecordEvent details */
        details?: (string|null);
    }

    /** Represents a RecordEvent. */
    class RecordEvent implements IRecordEvent {

        /**
         * Constructs a new RecordEvent.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IRecordEvent);

        /** RecordEvent time. */
        public time?: (ord.IDateTime|null);

        /** RecordEvent person. */
        public person?: (ord.IPerson|null);

        /** RecordEvent details. */
        public details: string;

        /**
         * Creates a new RecordEvent instance using the specified properties.
         * @param [properties] Properties to set
         * @returns RecordEvent instance
         */
        public static create(properties?: ord.IRecordEvent): ord.RecordEvent;

        /**
         * Encodes the specified RecordEvent message. Does not implicitly {@link ord.RecordEvent.verify|verify} messages.
         * @param message RecordEvent message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IRecordEvent, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified RecordEvent message, length delimited. Does not implicitly {@link ord.RecordEvent.verify|verify} messages.
         * @param message RecordEvent message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IRecordEvent, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a RecordEvent message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns RecordEvent
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.RecordEvent;

        /**
         * Decodes a RecordEvent message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns RecordEvent
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.RecordEvent;

        /**
         * Verifies a RecordEvent message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a RecordEvent message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns RecordEvent
         */
        public static fromObject(object: { [k: string]: any }): ord.RecordEvent;

        /**
         * Creates a plain object from a RecordEvent message. Also converts values to other types if specified.
         * @param message RecordEvent
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.RecordEvent, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this RecordEvent to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for RecordEvent
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a Time. */
    interface ITime {

        /** Time value */
        value?: (number|null);

        /** Time precision */
        precision?: (number|null);

        /** Time units */
        units?: (ord.Time.TimeUnit|null);
    }

    /**
     * To allow users to describe synthetic processes in whatever units they find
     * most natural, we define a fixed list of allowable units for each measurement
     * type. Upon submission to a centralized database, or using a validation and
     * canonicalization script, we will convert all values to the default units
     * (the first nonzero item in each enum).
     *
     * Each message also contains a `precision` field, which specifies the precision
     * of the measurement in the same units as the measurement itself. Often the
     * precision will be the standard deviation from an instrument calibration.
     */
    class Time implements ITime {

        /**
         * Constructs a new Time.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ITime);

        /** Time value. */
        public value?: (number|null);

        /** Time precision. */
        public precision?: (number|null);

        /** Time units. */
        public units: ord.Time.TimeUnit;

        /**
         * Creates a new Time instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Time instance
         */
        public static create(properties?: ord.ITime): ord.Time;

        /**
         * Encodes the specified Time message. Does not implicitly {@link ord.Time.verify|verify} messages.
         * @param message Time message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ITime, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Time message, length delimited. Does not implicitly {@link ord.Time.verify|verify} messages.
         * @param message Time message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ITime, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Time message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Time
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Time;

        /**
         * Decodes a Time message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Time
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Time;

        /**
         * Verifies a Time message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Time message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Time
         */
        public static fromObject(object: { [k: string]: any }): ord.Time;

        /**
         * Creates a plain object from a Time message. Also converts values to other types if specified.
         * @param message Time
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Time, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Time to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Time
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Time {

        /** TimeUnit enum. */
        enum TimeUnit {
            UNSPECIFIED = 0,
            DAY = 4,
            HOUR = 1,
            MINUTE = 2,
            SECOND = 3
        }
    }

    /** Properties of a Mass. */
    interface IMass {

        /** Mass value */
        value?: (number|null);

        /** Mass precision */
        precision?: (number|null);

        /** Mass units */
        units?: (ord.Mass.MassUnit|null);
    }

    /** Represents a Mass. */
    class Mass implements IMass {

        /**
         * Constructs a new Mass.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IMass);

        /** Mass value. */
        public value?: (number|null);

        /** Mass precision. */
        public precision?: (number|null);

        /** Mass units. */
        public units: ord.Mass.MassUnit;

        /**
         * Creates a new Mass instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Mass instance
         */
        public static create(properties?: ord.IMass): ord.Mass;

        /**
         * Encodes the specified Mass message. Does not implicitly {@link ord.Mass.verify|verify} messages.
         * @param message Mass message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IMass, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Mass message, length delimited. Does not implicitly {@link ord.Mass.verify|verify} messages.
         * @param message Mass message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IMass, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Mass message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Mass
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Mass;

        /**
         * Decodes a Mass message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Mass
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Mass;

        /**
         * Verifies a Mass message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Mass message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Mass
         */
        public static fromObject(object: { [k: string]: any }): ord.Mass;

        /**
         * Creates a plain object from a Mass message. Also converts values to other types if specified.
         * @param message Mass
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Mass, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Mass to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Mass
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Mass {

        /** MassUnit enum. */
        enum MassUnit {
            UNSPECIFIED = 0,
            KILOGRAM = 1,
            GRAM = 2,
            MILLIGRAM = 3,
            MICROGRAM = 4
        }
    }

    /** Properties of a Moles. */
    interface IMoles {

        /** Moles value */
        value?: (number|null);

        /** Moles precision */
        precision?: (number|null);

        /** Moles units */
        units?: (ord.Moles.MolesUnit|null);
    }

    /** Represents a Moles. */
    class Moles implements IMoles {

        /**
         * Constructs a new Moles.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IMoles);

        /** Moles value. */
        public value?: (number|null);

        /** Moles precision. */
        public precision?: (number|null);

        /** Moles units. */
        public units: ord.Moles.MolesUnit;

        /**
         * Creates a new Moles instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Moles instance
         */
        public static create(properties?: ord.IMoles): ord.Moles;

        /**
         * Encodes the specified Moles message. Does not implicitly {@link ord.Moles.verify|verify} messages.
         * @param message Moles message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IMoles, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Moles message, length delimited. Does not implicitly {@link ord.Moles.verify|verify} messages.
         * @param message Moles message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IMoles, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Moles message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Moles
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Moles;

        /**
         * Decodes a Moles message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Moles
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Moles;

        /**
         * Verifies a Moles message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Moles message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Moles
         */
        public static fromObject(object: { [k: string]: any }): ord.Moles;

        /**
         * Creates a plain object from a Moles message. Also converts values to other types if specified.
         * @param message Moles
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Moles, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Moles to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Moles
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Moles {

        /** MolesUnit enum. */
        enum MolesUnit {
            UNSPECIFIED = 0,
            MOLE = 1,
            MILLIMOLE = 2,
            MICROMOLE = 3,
            NANOMOLE = 4
        }
    }

    /** Properties of a Volume. */
    interface IVolume {

        /** Volume value */
        value?: (number|null);

        /** Volume precision */
        precision?: (number|null);

        /** Volume units */
        units?: (ord.Volume.VolumeUnit|null);
    }

    /** Represents a Volume. */
    class Volume implements IVolume {

        /**
         * Constructs a new Volume.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IVolume);

        /** Volume value. */
        public value?: (number|null);

        /** Volume precision. */
        public precision?: (number|null);

        /** Volume units. */
        public units: ord.Volume.VolumeUnit;

        /**
         * Creates a new Volume instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Volume instance
         */
        public static create(properties?: ord.IVolume): ord.Volume;

        /**
         * Encodes the specified Volume message. Does not implicitly {@link ord.Volume.verify|verify} messages.
         * @param message Volume message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IVolume, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Volume message, length delimited. Does not implicitly {@link ord.Volume.verify|verify} messages.
         * @param message Volume message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IVolume, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Volume message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Volume
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Volume;

        /**
         * Decodes a Volume message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Volume
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Volume;

        /**
         * Verifies a Volume message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Volume message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Volume
         */
        public static fromObject(object: { [k: string]: any }): ord.Volume;

        /**
         * Creates a plain object from a Volume message. Also converts values to other types if specified.
         * @param message Volume
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Volume, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Volume to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Volume
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Volume {

        /** VolumeUnit enum. */
        enum VolumeUnit {
            UNSPECIFIED = 0,
            LITER = 1,
            MILLILITER = 2,
            MICROLITER = 3,
            NANOLITER = 4
        }
    }

    /** Properties of a Concentration. */
    interface IConcentration {

        /** Concentration value */
        value?: (number|null);

        /** Concentration precision */
        precision?: (number|null);

        /** Concentration units */
        units?: (ord.Concentration.ConcentrationUnit|null);
    }

    /** Represents a Concentration. */
    class Concentration implements IConcentration {

        /**
         * Constructs a new Concentration.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IConcentration);

        /** Concentration value. */
        public value?: (number|null);

        /** Concentration precision. */
        public precision?: (number|null);

        /** Concentration units. */
        public units: ord.Concentration.ConcentrationUnit;

        /**
         * Creates a new Concentration instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Concentration instance
         */
        public static create(properties?: ord.IConcentration): ord.Concentration;

        /**
         * Encodes the specified Concentration message. Does not implicitly {@link ord.Concentration.verify|verify} messages.
         * @param message Concentration message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IConcentration, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Concentration message, length delimited. Does not implicitly {@link ord.Concentration.verify|verify} messages.
         * @param message Concentration message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IConcentration, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Concentration message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Concentration
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Concentration;

        /**
         * Decodes a Concentration message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Concentration
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Concentration;

        /**
         * Verifies a Concentration message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Concentration message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Concentration
         */
        public static fromObject(object: { [k: string]: any }): ord.Concentration;

        /**
         * Creates a plain object from a Concentration message. Also converts values to other types if specified.
         * @param message Concentration
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Concentration, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Concentration to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Concentration
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Concentration {

        /** ConcentrationUnit enum. */
        enum ConcentrationUnit {
            UNSPECIFIED = 0,
            MOLAR = 1,
            MILLIMOLAR = 2,
            MICROMOLAR = 3
        }
    }

    /** Properties of a Pressure. */
    interface IPressure {

        /** Pressure value */
        value?: (number|null);

        /** Pressure precision */
        precision?: (number|null);

        /** Pressure units */
        units?: (ord.Pressure.PressureUnit|null);
    }

    /** Represents a Pressure. */
    class Pressure implements IPressure {

        /**
         * Constructs a new Pressure.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IPressure);

        /** Pressure value. */
        public value?: (number|null);

        /** Pressure precision. */
        public precision?: (number|null);

        /** Pressure units. */
        public units: ord.Pressure.PressureUnit;

        /**
         * Creates a new Pressure instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Pressure instance
         */
        public static create(properties?: ord.IPressure): ord.Pressure;

        /**
         * Encodes the specified Pressure message. Does not implicitly {@link ord.Pressure.verify|verify} messages.
         * @param message Pressure message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IPressure, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Pressure message, length delimited. Does not implicitly {@link ord.Pressure.verify|verify} messages.
         * @param message Pressure message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IPressure, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Pressure message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Pressure
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Pressure;

        /**
         * Decodes a Pressure message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Pressure
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Pressure;

        /**
         * Verifies a Pressure message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Pressure message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Pressure
         */
        public static fromObject(object: { [k: string]: any }): ord.Pressure;

        /**
         * Creates a plain object from a Pressure message. Also converts values to other types if specified.
         * @param message Pressure
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Pressure, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Pressure to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Pressure
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Pressure {

        /** PressureUnit enum. */
        enum PressureUnit {
            UNSPECIFIED = 0,
            BAR = 1,
            ATMOSPHERE = 2,
            PSI = 3,
            KPSI = 4,
            PASCAL = 5,
            KILOPASCAL = 6,
            TORR = 7,
            MM_HG = 8
        }
    }

    /** Properties of a Temperature. */
    interface ITemperature {

        /** Temperature value */
        value?: (number|null);

        /** Temperature precision */
        precision?: (number|null);

        /** Temperature units */
        units?: (ord.Temperature.TemperatureUnit|null);
    }

    /** Represents a Temperature. */
    class Temperature implements ITemperature {

        /**
         * Constructs a new Temperature.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ITemperature);

        /** Temperature value. */
        public value?: (number|null);

        /** Temperature precision. */
        public precision?: (number|null);

        /** Temperature units. */
        public units: ord.Temperature.TemperatureUnit;

        /**
         * Creates a new Temperature instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Temperature instance
         */
        public static create(properties?: ord.ITemperature): ord.Temperature;

        /**
         * Encodes the specified Temperature message. Does not implicitly {@link ord.Temperature.verify|verify} messages.
         * @param message Temperature message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ITemperature, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Temperature message, length delimited. Does not implicitly {@link ord.Temperature.verify|verify} messages.
         * @param message Temperature message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ITemperature, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Temperature message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Temperature
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Temperature;

        /**
         * Decodes a Temperature message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Temperature
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Temperature;

        /**
         * Verifies a Temperature message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Temperature message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Temperature
         */
        public static fromObject(object: { [k: string]: any }): ord.Temperature;

        /**
         * Creates a plain object from a Temperature message. Also converts values to other types if specified.
         * @param message Temperature
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Temperature, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Temperature to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Temperature
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Temperature {

        /** TemperatureUnit enum. */
        enum TemperatureUnit {
            UNSPECIFIED = 0,
            CELSIUS = 1,
            FAHRENHEIT = 2,
            KELVIN = 3
        }
    }

    /** Properties of a Current. */
    interface ICurrent {

        /** Current value */
        value?: (number|null);

        /** Current precision */
        precision?: (number|null);

        /** Current units */
        units?: (ord.Current.CurrentUnit|null);
    }

    /** Represents a Current. */
    class Current implements ICurrent {

        /**
         * Constructs a new Current.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ICurrent);

        /** Current value. */
        public value?: (number|null);

        /** Current precision. */
        public precision?: (number|null);

        /** Current units. */
        public units: ord.Current.CurrentUnit;

        /**
         * Creates a new Current instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Current instance
         */
        public static create(properties?: ord.ICurrent): ord.Current;

        /**
         * Encodes the specified Current message. Does not implicitly {@link ord.Current.verify|verify} messages.
         * @param message Current message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ICurrent, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Current message, length delimited. Does not implicitly {@link ord.Current.verify|verify} messages.
         * @param message Current message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ICurrent, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Current message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Current
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Current;

        /**
         * Decodes a Current message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Current
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Current;

        /**
         * Verifies a Current message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Current message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Current
         */
        public static fromObject(object: { [k: string]: any }): ord.Current;

        /**
         * Creates a plain object from a Current message. Also converts values to other types if specified.
         * @param message Current
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Current, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Current to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Current
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Current {

        /** CurrentUnit enum. */
        enum CurrentUnit {
            UNSPECIFIED = 0,
            AMPERE = 1,
            MILLIAMPERE = 2
        }
    }

    /** Properties of a Voltage. */
    interface IVoltage {

        /** Voltage value */
        value?: (number|null);

        /** Voltage precision */
        precision?: (number|null);

        /** Voltage units */
        units?: (ord.Voltage.VoltageUnit|null);
    }

    /** Represents a Voltage. */
    class Voltage implements IVoltage {

        /**
         * Constructs a new Voltage.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IVoltage);

        /** Voltage value. */
        public value?: (number|null);

        /** Voltage precision. */
        public precision?: (number|null);

        /** Voltage units. */
        public units: ord.Voltage.VoltageUnit;

        /**
         * Creates a new Voltage instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Voltage instance
         */
        public static create(properties?: ord.IVoltage): ord.Voltage;

        /**
         * Encodes the specified Voltage message. Does not implicitly {@link ord.Voltage.verify|verify} messages.
         * @param message Voltage message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IVoltage, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Voltage message, length delimited. Does not implicitly {@link ord.Voltage.verify|verify} messages.
         * @param message Voltage message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IVoltage, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Voltage message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Voltage
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Voltage;

        /**
         * Decodes a Voltage message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Voltage
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Voltage;

        /**
         * Verifies a Voltage message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Voltage message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Voltage
         */
        public static fromObject(object: { [k: string]: any }): ord.Voltage;

        /**
         * Creates a plain object from a Voltage message. Also converts values to other types if specified.
         * @param message Voltage
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Voltage, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Voltage to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Voltage
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Voltage {

        /** VoltageUnit enum. */
        enum VoltageUnit {
            UNSPECIFIED = 0,
            VOLT = 1,
            MILLIVOLT = 2
        }
    }

    /** Properties of a Length. */
    interface ILength {

        /** Length value */
        value?: (number|null);

        /** Length precision */
        precision?: (number|null);

        /** Length units */
        units?: (ord.Length.LengthUnit|null);
    }

    /** Represents a Length. */
    class Length implements ILength {

        /**
         * Constructs a new Length.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.ILength);

        /** Length value. */
        public value?: (number|null);

        /** Length precision. */
        public precision?: (number|null);

        /** Length units. */
        public units: ord.Length.LengthUnit;

        /**
         * Creates a new Length instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Length instance
         */
        public static create(properties?: ord.ILength): ord.Length;

        /**
         * Encodes the specified Length message. Does not implicitly {@link ord.Length.verify|verify} messages.
         * @param message Length message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.ILength, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Length message, length delimited. Does not implicitly {@link ord.Length.verify|verify} messages.
         * @param message Length message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.ILength, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Length message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Length
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Length;

        /**
         * Decodes a Length message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Length
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Length;

        /**
         * Verifies a Length message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Length message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Length
         */
        public static fromObject(object: { [k: string]: any }): ord.Length;

        /**
         * Creates a plain object from a Length message. Also converts values to other types if specified.
         * @param message Length
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Length, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Length to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Length
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Length {

        /** LengthUnit enum. */
        enum LengthUnit {
            UNSPECIFIED = 0,
            CENTIMETER = 1,
            MILLIMETER = 2,
            METER = 3,
            INCH = 4,
            FOOT = 5
        }
    }

    /** Properties of a Wavelength. */
    interface IWavelength {

        /** Wavelength value */
        value?: (number|null);

        /** Wavelength precision */
        precision?: (number|null);

        /** Wavelength units */
        units?: (ord.Wavelength.WavelengthUnit|null);
    }

    /** Represents a Wavelength. */
    class Wavelength implements IWavelength {

        /**
         * Constructs a new Wavelength.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IWavelength);

        /** Wavelength value. */
        public value?: (number|null);

        /** Wavelength precision. */
        public precision?: (number|null);

        /** Wavelength units. */
        public units: ord.Wavelength.WavelengthUnit;

        /**
         * Creates a new Wavelength instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Wavelength instance
         */
        public static create(properties?: ord.IWavelength): ord.Wavelength;

        /**
         * Encodes the specified Wavelength message. Does not implicitly {@link ord.Wavelength.verify|verify} messages.
         * @param message Wavelength message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IWavelength, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Wavelength message, length delimited. Does not implicitly {@link ord.Wavelength.verify|verify} messages.
         * @param message Wavelength message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IWavelength, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Wavelength message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Wavelength
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Wavelength;

        /**
         * Decodes a Wavelength message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Wavelength
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Wavelength;

        /**
         * Verifies a Wavelength message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Wavelength message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Wavelength
         */
        public static fromObject(object: { [k: string]: any }): ord.Wavelength;

        /**
         * Creates a plain object from a Wavelength message. Also converts values to other types if specified.
         * @param message Wavelength
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Wavelength, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Wavelength to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Wavelength
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Wavelength {

        /** WavelengthUnit enum. */
        enum WavelengthUnit {
            UNSPECIFIED = 0,
            NANOMETER = 1,
            WAVENUMBER = 2
        }
    }

    /** Properties of a FlowRate. */
    interface IFlowRate {

        /** FlowRate value */
        value?: (number|null);

        /** FlowRate precision */
        precision?: (number|null);

        /** FlowRate units */
        units?: (ord.FlowRate.FlowRateUnit|null);
    }

    /** Represents a FlowRate. */
    class FlowRate implements IFlowRate {

        /**
         * Constructs a new FlowRate.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IFlowRate);

        /** FlowRate value. */
        public value?: (number|null);

        /** FlowRate precision. */
        public precision?: (number|null);

        /** FlowRate units. */
        public units: ord.FlowRate.FlowRateUnit;

        /**
         * Creates a new FlowRate instance using the specified properties.
         * @param [properties] Properties to set
         * @returns FlowRate instance
         */
        public static create(properties?: ord.IFlowRate): ord.FlowRate;

        /**
         * Encodes the specified FlowRate message. Does not implicitly {@link ord.FlowRate.verify|verify} messages.
         * @param message FlowRate message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IFlowRate, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified FlowRate message, length delimited. Does not implicitly {@link ord.FlowRate.verify|verify} messages.
         * @param message FlowRate message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IFlowRate, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a FlowRate message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns FlowRate
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.FlowRate;

        /**
         * Decodes a FlowRate message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns FlowRate
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.FlowRate;

        /**
         * Verifies a FlowRate message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a FlowRate message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns FlowRate
         */
        public static fromObject(object: { [k: string]: any }): ord.FlowRate;

        /**
         * Creates a plain object from a FlowRate message. Also converts values to other types if specified.
         * @param message FlowRate
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.FlowRate, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this FlowRate to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for FlowRate
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace FlowRate {

        /** FlowRateUnit enum. */
        enum FlowRateUnit {
            UNSPECIFIED = 0,
            MICROLITER_PER_MINUTE = 1,
            MICROLITER_PER_SECOND = 2,
            MILLILITER_PER_MINUTE = 3,
            MILLILITER_PER_SECOND = 4,
            MICROLITER_PER_HOUR = 5
        }
    }

    /** Properties of a Percentage. */
    interface IPercentage {

        /** Percentage value */
        value?: (number|null);

        /** Percentage precision */
        precision?: (number|null);
    }

    /** Represents a Percentage. */
    class Percentage implements IPercentage {

        /**
         * Constructs a new Percentage.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IPercentage);

        /** Percentage value. */
        public value?: (number|null);

        /** Percentage precision. */
        public precision?: (number|null);

        /**
         * Creates a new Percentage instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Percentage instance
         */
        public static create(properties?: ord.IPercentage): ord.Percentage;

        /**
         * Encodes the specified Percentage message. Does not implicitly {@link ord.Percentage.verify|verify} messages.
         * @param message Percentage message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IPercentage, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Percentage message, length delimited. Does not implicitly {@link ord.Percentage.verify|verify} messages.
         * @param message Percentage message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IPercentage, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Percentage message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Percentage
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Percentage;

        /**
         * Decodes a Percentage message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Percentage
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Percentage;

        /**
         * Verifies a Percentage message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Percentage message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Percentage
         */
        public static fromObject(object: { [k: string]: any }): ord.Percentage;

        /**
         * Creates a plain object from a Percentage message. Also converts values to other types if specified.
         * @param message Percentage
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Percentage, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Percentage to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Percentage
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a FloatValue. */
    interface IFloatValue {

        /** FloatValue value */
        value?: (number|null);

        /** FloatValue precision */
        precision?: (number|null);
    }

    /** Represents a FloatValue. */
    class FloatValue implements IFloatValue {

        /**
         * Constructs a new FloatValue.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IFloatValue);

        /** FloatValue value. */
        public value?: (number|null);

        /** FloatValue precision. */
        public precision?: (number|null);

        /**
         * Creates a new FloatValue instance using the specified properties.
         * @param [properties] Properties to set
         * @returns FloatValue instance
         */
        public static create(properties?: ord.IFloatValue): ord.FloatValue;

        /**
         * Encodes the specified FloatValue message. Does not implicitly {@link ord.FloatValue.verify|verify} messages.
         * @param message FloatValue message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IFloatValue, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified FloatValue message, length delimited. Does not implicitly {@link ord.FloatValue.verify|verify} messages.
         * @param message FloatValue message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IFloatValue, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a FloatValue message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns FloatValue
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.FloatValue;

        /**
         * Decodes a FloatValue message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns FloatValue
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.FloatValue;

        /**
         * Verifies a FloatValue message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a FloatValue message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns FloatValue
         */
        public static fromObject(object: { [k: string]: any }): ord.FloatValue;

        /**
         * Creates a plain object from a FloatValue message. Also converts values to other types if specified.
         * @param message FloatValue
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.FloatValue, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this FloatValue to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for FloatValue
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a Data. */
    interface IData {

        /** Data floatValue */
        floatValue?: (number|null);

        /** Data integerValue */
        integerValue?: (number|null);

        /** Data bytesValue */
        bytesValue?: (Uint8Array|null);

        /** Data stringValue */
        stringValue?: (string|null);

        /** Data url */
        url?: (string|null);

        /** Data description */
        description?: (string|null);

        /** Data format */
        format?: (string|null);
    }

    /** Represents a Data. */
    class Data implements IData {

        /**
         * Constructs a new Data.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord.IData);

        /** Data floatValue. */
        public floatValue?: (number|null);

        /** Data integerValue. */
        public integerValue?: (number|null);

        /** Data bytesValue. */
        public bytesValue?: (Uint8Array|null);

        /** Data stringValue. */
        public stringValue?: (string|null);

        /** Data url. */
        public url?: (string|null);

        /** Data description. */
        public description: string;

        /** Data format. */
        public format: string;

        /** Data kind. */
        public kind?: ("floatValue"|"integerValue"|"bytesValue"|"stringValue"|"url");

        /**
         * Creates a new Data instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Data instance
         */
        public static create(properties?: ord.IData): ord.Data;

        /**
         * Encodes the specified Data message. Does not implicitly {@link ord.Data.verify|verify} messages.
         * @param message Data message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord.IData, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Data message, length delimited. Does not implicitly {@link ord.Data.verify|verify} messages.
         * @param message Data message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord.IData, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Data message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Data
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord.Data;

        /**
         * Decodes a Data message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Data
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord.Data;

        /**
         * Verifies a Data message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Data message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Data
         */
        public static fromObject(object: { [k: string]: any }): ord.Data;

        /**
         * Creates a plain object from a Data message. Also converts values to other types if specified.
         * @param message Data
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord.Data, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Data to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Data
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }
}

/** Namespace ord_test. */
export namespace ord_test {

    /** Properties of a Scalar. */
    interface IScalar {

        /** Scalar int32Value */
        int32Value?: (number|null);

        /** Scalar int64Value */
        int64Value?: (number|Long|null);

        /** Scalar floatValue */
        floatValue?: (number|null);

        /** Scalar stringValue */
        stringValue?: (string|null);

        /** Scalar bytesValue */
        bytesValue?: (Uint8Array|null);

        /** Scalar boolValue */
        boolValue?: (boolean|null);
    }

    /** Represents a Scalar. */
    class Scalar implements IScalar {

        /**
         * Constructs a new Scalar.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IScalar);

        /** Scalar int32Value. */
        public int32Value: number;

        /** Scalar int64Value. */
        public int64Value: (number|Long);

        /** Scalar floatValue. */
        public floatValue?: (number|null);

        /** Scalar stringValue. */
        public stringValue: string;

        /** Scalar bytesValue. */
        public bytesValue: Uint8Array;

        /** Scalar boolValue. */
        public boolValue?: (boolean|null);

        /**
         * Creates a new Scalar instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Scalar instance
         */
        public static create(properties?: ord_test.IScalar): ord_test.Scalar;

        /**
         * Encodes the specified Scalar message. Does not implicitly {@link ord_test.Scalar.verify|verify} messages.
         * @param message Scalar message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IScalar, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Scalar message, length delimited. Does not implicitly {@link ord_test.Scalar.verify|verify} messages.
         * @param message Scalar message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IScalar, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Scalar message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Scalar
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.Scalar;

        /**
         * Decodes a Scalar message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Scalar
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.Scalar;

        /**
         * Verifies a Scalar message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Scalar message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Scalar
         */
        public static fromObject(object: { [k: string]: any }): ord_test.Scalar;

        /**
         * Creates a plain object from a Scalar message. Also converts values to other types if specified.
         * @param message Scalar
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.Scalar, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Scalar to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Scalar
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a RepeatedScalar. */
    interface IRepeatedScalar {

        /** RepeatedScalar values */
        values?: (number[]|null);
    }

    /** Represents a RepeatedScalar. */
    class RepeatedScalar implements IRepeatedScalar {

        /**
         * Constructs a new RepeatedScalar.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IRepeatedScalar);

        /** RepeatedScalar values. */
        public values: number[];

        /**
         * Creates a new RepeatedScalar instance using the specified properties.
         * @param [properties] Properties to set
         * @returns RepeatedScalar instance
         */
        public static create(properties?: ord_test.IRepeatedScalar): ord_test.RepeatedScalar;

        /**
         * Encodes the specified RepeatedScalar message. Does not implicitly {@link ord_test.RepeatedScalar.verify|verify} messages.
         * @param message RepeatedScalar message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IRepeatedScalar, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified RepeatedScalar message, length delimited. Does not implicitly {@link ord_test.RepeatedScalar.verify|verify} messages.
         * @param message RepeatedScalar message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IRepeatedScalar, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a RepeatedScalar message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns RepeatedScalar
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.RepeatedScalar;

        /**
         * Decodes a RepeatedScalar message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns RepeatedScalar
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.RepeatedScalar;

        /**
         * Verifies a RepeatedScalar message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a RepeatedScalar message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns RepeatedScalar
         */
        public static fromObject(object: { [k: string]: any }): ord_test.RepeatedScalar;

        /**
         * Creates a plain object from a RepeatedScalar message. Also converts values to other types if specified.
         * @param message RepeatedScalar
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.RepeatedScalar, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this RepeatedScalar to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for RepeatedScalar
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of an Enum. */
    interface IEnum {

        /** Enum value */
        value?: (ord_test.Enum.EnumValues|null);
    }

    /** Represents an Enum. */
    class Enum implements IEnum {

        /**
         * Constructs a new Enum.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IEnum);

        /** Enum value. */
        public value: ord_test.Enum.EnumValues;

        /**
         * Creates a new Enum instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Enum instance
         */
        public static create(properties?: ord_test.IEnum): ord_test.Enum;

        /**
         * Encodes the specified Enum message. Does not implicitly {@link ord_test.Enum.verify|verify} messages.
         * @param message Enum message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IEnum, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Enum message, length delimited. Does not implicitly {@link ord_test.Enum.verify|verify} messages.
         * @param message Enum message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IEnum, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes an Enum message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Enum
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.Enum;

        /**
         * Decodes an Enum message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Enum
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.Enum;

        /**
         * Verifies an Enum message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates an Enum message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Enum
         */
        public static fromObject(object: { [k: string]: any }): ord_test.Enum;

        /**
         * Creates a plain object from an Enum message. Also converts values to other types if specified.
         * @param message Enum
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.Enum, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Enum to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Enum
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Enum {

        /** EnumValues enum. */
        enum EnumValues {
            UNSPECIFIED = 0,
            FIRST = 1,
            SECOND = 2
        }
    }

    /** Properties of a RepeatedEnum. */
    interface IRepeatedEnum {

        /** RepeatedEnum values */
        values?: (ord_test.RepeatedEnum.EnumValues[]|null);
    }

    /** Represents a RepeatedEnum. */
    class RepeatedEnum implements IRepeatedEnum {

        /**
         * Constructs a new RepeatedEnum.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IRepeatedEnum);

        /** RepeatedEnum values. */
        public values: ord_test.RepeatedEnum.EnumValues[];

        /**
         * Creates a new RepeatedEnum instance using the specified properties.
         * @param [properties] Properties to set
         * @returns RepeatedEnum instance
         */
        public static create(properties?: ord_test.IRepeatedEnum): ord_test.RepeatedEnum;

        /**
         * Encodes the specified RepeatedEnum message. Does not implicitly {@link ord_test.RepeatedEnum.verify|verify} messages.
         * @param message RepeatedEnum message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IRepeatedEnum, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified RepeatedEnum message, length delimited. Does not implicitly {@link ord_test.RepeatedEnum.verify|verify} messages.
         * @param message RepeatedEnum message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IRepeatedEnum, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a RepeatedEnum message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns RepeatedEnum
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.RepeatedEnum;

        /**
         * Decodes a RepeatedEnum message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns RepeatedEnum
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.RepeatedEnum;

        /**
         * Verifies a RepeatedEnum message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a RepeatedEnum message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns RepeatedEnum
         */
        public static fromObject(object: { [k: string]: any }): ord_test.RepeatedEnum;

        /**
         * Creates a plain object from a RepeatedEnum message. Also converts values to other types if specified.
         * @param message RepeatedEnum
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.RepeatedEnum, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this RepeatedEnum to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for RepeatedEnum
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace RepeatedEnum {

        /** EnumValues enum. */
        enum EnumValues {
            UNSPECIFIED = 0,
            FIRST = 1,
            SECOND = 2
        }
    }

    /** Properties of a Nested. */
    interface INested {

        /** Nested child */
        child?: (ord_test.Nested.IChild|null);
    }

    /** Represents a Nested. */
    class Nested implements INested {

        /**
         * Constructs a new Nested.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.INested);

        /** Nested child. */
        public child?: (ord_test.Nested.IChild|null);

        /**
         * Creates a new Nested instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Nested instance
         */
        public static create(properties?: ord_test.INested): ord_test.Nested;

        /**
         * Encodes the specified Nested message. Does not implicitly {@link ord_test.Nested.verify|verify} messages.
         * @param message Nested message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.INested, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Nested message, length delimited. Does not implicitly {@link ord_test.Nested.verify|verify} messages.
         * @param message Nested message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.INested, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Nested message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Nested
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.Nested;

        /**
         * Decodes a Nested message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Nested
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.Nested;

        /**
         * Verifies a Nested message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Nested message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Nested
         */
        public static fromObject(object: { [k: string]: any }): ord_test.Nested;

        /**
         * Creates a plain object from a Nested message. Also converts values to other types if specified.
         * @param message Nested
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.Nested, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Nested to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Nested
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace Nested {

        /** Properties of a Child. */
        interface IChild {

            /** Child value */
            value?: (number|null);
        }

        /** Represents a Child. */
        class Child implements IChild {

            /**
             * Constructs a new Child.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord_test.Nested.IChild);

            /** Child value. */
            public value?: (number|null);

            /**
             * Creates a new Child instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Child instance
             */
            public static create(properties?: ord_test.Nested.IChild): ord_test.Nested.Child;

            /**
             * Encodes the specified Child message. Does not implicitly {@link ord_test.Nested.Child.verify|verify} messages.
             * @param message Child message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord_test.Nested.IChild, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Child message, length delimited. Does not implicitly {@link ord_test.Nested.Child.verify|verify} messages.
             * @param message Child message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord_test.Nested.IChild, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a Child message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Child
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.Nested.Child;

            /**
             * Decodes a Child message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Child
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.Nested.Child;

            /**
             * Verifies a Child message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a Child message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Child
             */
            public static fromObject(object: { [k: string]: any }): ord_test.Nested.Child;

            /**
             * Creates a plain object from a Child message. Also converts values to other types if specified.
             * @param message Child
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord_test.Nested.Child, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Child to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Child
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }
    }

    /** Properties of a RepeatedNested. */
    interface IRepeatedNested {

        /** RepeatedNested children */
        children?: (ord_test.RepeatedNested.IChild[]|null);
    }

    /** Represents a RepeatedNested. */
    class RepeatedNested implements IRepeatedNested {

        /**
         * Constructs a new RepeatedNested.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IRepeatedNested);

        /** RepeatedNested children. */
        public children: ord_test.RepeatedNested.IChild[];

        /**
         * Creates a new RepeatedNested instance using the specified properties.
         * @param [properties] Properties to set
         * @returns RepeatedNested instance
         */
        public static create(properties?: ord_test.IRepeatedNested): ord_test.RepeatedNested;

        /**
         * Encodes the specified RepeatedNested message. Does not implicitly {@link ord_test.RepeatedNested.verify|verify} messages.
         * @param message RepeatedNested message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IRepeatedNested, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified RepeatedNested message, length delimited. Does not implicitly {@link ord_test.RepeatedNested.verify|verify} messages.
         * @param message RepeatedNested message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IRepeatedNested, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a RepeatedNested message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns RepeatedNested
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.RepeatedNested;

        /**
         * Decodes a RepeatedNested message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns RepeatedNested
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.RepeatedNested;

        /**
         * Verifies a RepeatedNested message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a RepeatedNested message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns RepeatedNested
         */
        public static fromObject(object: { [k: string]: any }): ord_test.RepeatedNested;

        /**
         * Creates a plain object from a RepeatedNested message. Also converts values to other types if specified.
         * @param message RepeatedNested
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.RepeatedNested, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this RepeatedNested to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for RepeatedNested
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace RepeatedNested {

        /** Properties of a Child. */
        interface IChild {

            /** Child value */
            value?: (number|null);
        }

        /** Represents a Child. */
        class Child implements IChild {

            /**
             * Constructs a new Child.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord_test.RepeatedNested.IChild);

            /** Child value. */
            public value?: (number|null);

            /**
             * Creates a new Child instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Child instance
             */
            public static create(properties?: ord_test.RepeatedNested.IChild): ord_test.RepeatedNested.Child;

            /**
             * Encodes the specified Child message. Does not implicitly {@link ord_test.RepeatedNested.Child.verify|verify} messages.
             * @param message Child message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord_test.RepeatedNested.IChild, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Child message, length delimited. Does not implicitly {@link ord_test.RepeatedNested.Child.verify|verify} messages.
             * @param message Child message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord_test.RepeatedNested.IChild, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a Child message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Child
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.RepeatedNested.Child;

            /**
             * Decodes a Child message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Child
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.RepeatedNested.Child;

            /**
             * Verifies a Child message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a Child message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Child
             */
            public static fromObject(object: { [k: string]: any }): ord_test.RepeatedNested.Child;

            /**
             * Creates a plain object from a Child message. Also converts values to other types if specified.
             * @param message Child
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord_test.RepeatedNested.Child, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Child to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Child
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }
    }

    /** Properties of a Map. */
    interface IMap {

        /** Map values */
        values?: ({ [k: string]: number }|null);
    }

    /** Represents a Map. */
    class Map implements IMap {

        /**
         * Constructs a new Map.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IMap);

        /** Map values. */
        public values: { [k: string]: number };

        /**
         * Creates a new Map instance using the specified properties.
         * @param [properties] Properties to set
         * @returns Map instance
         */
        public static create(properties?: ord_test.IMap): ord_test.Map;

        /**
         * Encodes the specified Map message. Does not implicitly {@link ord_test.Map.verify|verify} messages.
         * @param message Map message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IMap, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified Map message, length delimited. Does not implicitly {@link ord_test.Map.verify|verify} messages.
         * @param message Map message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IMap, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a Map message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns Map
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.Map;

        /**
         * Decodes a Map message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns Map
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.Map;

        /**
         * Verifies a Map message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a Map message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns Map
         */
        public static fromObject(object: { [k: string]: any }): ord_test.Map;

        /**
         * Creates a plain object from a Map message. Also converts values to other types if specified.
         * @param message Map
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.Map, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this Map to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for Map
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    /** Properties of a MapNested. */
    interface IMapNested {

        /** MapNested children */
        children?: ({ [k: string]: ord_test.MapNested.IChild }|null);
    }

    /** Represents a MapNested. */
    class MapNested implements IMapNested {

        /**
         * Constructs a new MapNested.
         * @param [properties] Properties to set
         */
        constructor(properties?: ord_test.IMapNested);

        /** MapNested children. */
        public children: { [k: string]: ord_test.MapNested.IChild };

        /**
         * Creates a new MapNested instance using the specified properties.
         * @param [properties] Properties to set
         * @returns MapNested instance
         */
        public static create(properties?: ord_test.IMapNested): ord_test.MapNested;

        /**
         * Encodes the specified MapNested message. Does not implicitly {@link ord_test.MapNested.verify|verify} messages.
         * @param message MapNested message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encode(message: ord_test.IMapNested, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Encodes the specified MapNested message, length delimited. Does not implicitly {@link ord_test.MapNested.verify|verify} messages.
         * @param message MapNested message or plain object to encode
         * @param [writer] Writer to encode to
         * @returns Writer
         */
        public static encodeDelimited(message: ord_test.IMapNested, writer?: $protobuf.Writer): $protobuf.Writer;

        /**
         * Decodes a MapNested message from the specified reader or buffer.
         * @param reader Reader or buffer to decode from
         * @param [length] Message length if known beforehand
         * @returns MapNested
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.MapNested;

        /**
         * Decodes a MapNested message from the specified reader or buffer, length delimited.
         * @param reader Reader or buffer to decode from
         * @returns MapNested
         * @throws {Error} If the payload is not a reader or valid buffer
         * @throws {$protobuf.util.ProtocolError} If required fields are missing
         */
        public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.MapNested;

        /**
         * Verifies a MapNested message.
         * @param message Plain object to verify
         * @returns `null` if valid, otherwise the reason why it is not
         */
        public static verify(message: { [k: string]: any }): (string|null);

        /**
         * Creates a MapNested message from a plain object. Also converts values to their respective internal types.
         * @param object Plain object
         * @returns MapNested
         */
        public static fromObject(object: { [k: string]: any }): ord_test.MapNested;

        /**
         * Creates a plain object from a MapNested message. Also converts values to other types if specified.
         * @param message MapNested
         * @param [options] Conversion options
         * @returns Plain object
         */
        public static toObject(message: ord_test.MapNested, options?: $protobuf.IConversionOptions): { [k: string]: any };

        /**
         * Converts this MapNested to JSON.
         * @returns JSON object
         */
        public toJSON(): { [k: string]: any };

        /**
         * Gets the default type url for MapNested
         * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
         * @returns The default type url
         */
        public static getTypeUrl(typeUrlPrefix?: string): string;
    }

    namespace MapNested {

        /** Properties of a Child. */
        interface IChild {

            /** Child value */
            value?: (number|null);
        }

        /** Represents a Child. */
        class Child implements IChild {

            /**
             * Constructs a new Child.
             * @param [properties] Properties to set
             */
            constructor(properties?: ord_test.MapNested.IChild);

            /** Child value. */
            public value?: (number|null);

            /**
             * Creates a new Child instance using the specified properties.
             * @param [properties] Properties to set
             * @returns Child instance
             */
            public static create(properties?: ord_test.MapNested.IChild): ord_test.MapNested.Child;

            /**
             * Encodes the specified Child message. Does not implicitly {@link ord_test.MapNested.Child.verify|verify} messages.
             * @param message Child message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encode(message: ord_test.MapNested.IChild, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Encodes the specified Child message, length delimited. Does not implicitly {@link ord_test.MapNested.Child.verify|verify} messages.
             * @param message Child message or plain object to encode
             * @param [writer] Writer to encode to
             * @returns Writer
             */
            public static encodeDelimited(message: ord_test.MapNested.IChild, writer?: $protobuf.Writer): $protobuf.Writer;

            /**
             * Decodes a Child message from the specified reader or buffer.
             * @param reader Reader or buffer to decode from
             * @param [length] Message length if known beforehand
             * @returns Child
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decode(reader: ($protobuf.Reader|Uint8Array), length?: number): ord_test.MapNested.Child;

            /**
             * Decodes a Child message from the specified reader or buffer, length delimited.
             * @param reader Reader or buffer to decode from
             * @returns Child
             * @throws {Error} If the payload is not a reader or valid buffer
             * @throws {$protobuf.util.ProtocolError} If required fields are missing
             */
            public static decodeDelimited(reader: ($protobuf.Reader|Uint8Array)): ord_test.MapNested.Child;

            /**
             * Verifies a Child message.
             * @param message Plain object to verify
             * @returns `null` if valid, otherwise the reason why it is not
             */
            public static verify(message: { [k: string]: any }): (string|null);

            /**
             * Creates a Child message from a plain object. Also converts values to their respective internal types.
             * @param object Plain object
             * @returns Child
             */
            public static fromObject(object: { [k: string]: any }): ord_test.MapNested.Child;

            /**
             * Creates a plain object from a Child message. Also converts values to other types if specified.
             * @param message Child
             * @param [options] Conversion options
             * @returns Plain object
             */
            public static toObject(message: ord_test.MapNested.Child, options?: $protobuf.IConversionOptions): { [k: string]: any };

            /**
             * Converts this Child to JSON.
             * @returns JSON object
             */
            public toJSON(): { [k: string]: any };

            /**
             * Gets the default type url for Child
             * @param [typeUrlPrefix] your custom typeUrlPrefix(default "type.googleapis.com")
             * @returns The default type url
             */
            public static getTypeUrl(typeUrlPrefix?: string): string;
        }
    }
}
