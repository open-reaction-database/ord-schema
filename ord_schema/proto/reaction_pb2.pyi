from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class Reaction(_message.Message):
    __slots__ = ["identifiers", "inputs", "setup", "conditions", "notes", "observations", "workups", "outcomes", "provenance", "reaction_id"]
    class InputsEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: ReactionInput
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[ReactionInput, _Mapping]] = ...) -> None: ...
    IDENTIFIERS_FIELD_NUMBER: _ClassVar[int]
    INPUTS_FIELD_NUMBER: _ClassVar[int]
    SETUP_FIELD_NUMBER: _ClassVar[int]
    CONDITIONS_FIELD_NUMBER: _ClassVar[int]
    NOTES_FIELD_NUMBER: _ClassVar[int]
    OBSERVATIONS_FIELD_NUMBER: _ClassVar[int]
    WORKUPS_FIELD_NUMBER: _ClassVar[int]
    OUTCOMES_FIELD_NUMBER: _ClassVar[int]
    PROVENANCE_FIELD_NUMBER: _ClassVar[int]
    REACTION_ID_FIELD_NUMBER: _ClassVar[int]
    identifiers: _containers.RepeatedCompositeFieldContainer[ReactionIdentifier]
    inputs: _containers.MessageMap[str, ReactionInput]
    setup: ReactionSetup
    conditions: ReactionConditions
    notes: ReactionNotes
    observations: _containers.RepeatedCompositeFieldContainer[ReactionObservation]
    workups: _containers.RepeatedCompositeFieldContainer[ReactionWorkup]
    outcomes: _containers.RepeatedCompositeFieldContainer[ReactionOutcome]
    provenance: ReactionProvenance
    reaction_id: str
    def __init__(self, identifiers: _Optional[_Iterable[_Union[ReactionIdentifier, _Mapping]]] = ..., inputs: _Optional[_Mapping[str, ReactionInput]] = ..., setup: _Optional[_Union[ReactionSetup, _Mapping]] = ..., conditions: _Optional[_Union[ReactionConditions, _Mapping]] = ..., notes: _Optional[_Union[ReactionNotes, _Mapping]] = ..., observations: _Optional[_Iterable[_Union[ReactionObservation, _Mapping]]] = ..., workups: _Optional[_Iterable[_Union[ReactionWorkup, _Mapping]]] = ..., outcomes: _Optional[_Iterable[_Union[ReactionOutcome, _Mapping]]] = ..., provenance: _Optional[_Union[ReactionProvenance, _Mapping]] = ..., reaction_id: _Optional[str] = ...) -> None: ...

class ReactionIdentifier(_message.Message):
    __slots__ = ["type", "details", "value", "is_mapped"]
    class ReactionIdentifierType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
        CUSTOM: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
        REACTION_SMILES: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
        REACTION_CXSMILES: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
        RDFILE: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
        RINCHI: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
        REACTION_TYPE: _ClassVar[ReactionIdentifier.ReactionIdentifierType]
    UNSPECIFIED: ReactionIdentifier.ReactionIdentifierType
    CUSTOM: ReactionIdentifier.ReactionIdentifierType
    REACTION_SMILES: ReactionIdentifier.ReactionIdentifierType
    REACTION_CXSMILES: ReactionIdentifier.ReactionIdentifierType
    RDFILE: ReactionIdentifier.ReactionIdentifierType
    RINCHI: ReactionIdentifier.ReactionIdentifierType
    REACTION_TYPE: ReactionIdentifier.ReactionIdentifierType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    IS_MAPPED_FIELD_NUMBER: _ClassVar[int]
    type: ReactionIdentifier.ReactionIdentifierType
    details: str
    value: str
    is_mapped: bool
    def __init__(self, type: _Optional[_Union[ReactionIdentifier.ReactionIdentifierType, str]] = ..., details: _Optional[str] = ..., value: _Optional[str] = ..., is_mapped: bool = ...) -> None: ...

class ReactionInput(_message.Message):
    __slots__ = ["components", "crude_components", "addition_order", "addition_time", "addition_speed", "addition_duration", "flow_rate", "addition_device", "addition_temperature", "texture"]
    class AdditionSpeed(_message.Message):
        __slots__ = ["type", "details"]
        class AdditionSpeedType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
            ALL_AT_ONCE: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
            FAST: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
            SLOW: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
            DROPWISE: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
            CONTINUOUS: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
            PORTIONWISE: _ClassVar[ReactionInput.AdditionSpeed.AdditionSpeedType]
        UNSPECIFIED: ReactionInput.AdditionSpeed.AdditionSpeedType
        ALL_AT_ONCE: ReactionInput.AdditionSpeed.AdditionSpeedType
        FAST: ReactionInput.AdditionSpeed.AdditionSpeedType
        SLOW: ReactionInput.AdditionSpeed.AdditionSpeedType
        DROPWISE: ReactionInput.AdditionSpeed.AdditionSpeedType
        CONTINUOUS: ReactionInput.AdditionSpeed.AdditionSpeedType
        PORTIONWISE: ReactionInput.AdditionSpeed.AdditionSpeedType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: ReactionInput.AdditionSpeed.AdditionSpeedType
        details: str
        def __init__(self, type: _Optional[_Union[ReactionInput.AdditionSpeed.AdditionSpeedType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    class AdditionDevice(_message.Message):
        __slots__ = ["type", "details"]
        class AdditionDeviceType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            CUSTOM: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            NONE: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            SYRINGE: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            CANNULA: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            ADDITION_FUNNEL: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            PIPETTE: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
            POSITIVE_DISPLACEMENT_PIPETTE: _ClassVar[ReactionInput.AdditionDevice.AdditionDeviceType]
        UNSPECIFIED: ReactionInput.AdditionDevice.AdditionDeviceType
        CUSTOM: ReactionInput.AdditionDevice.AdditionDeviceType
        NONE: ReactionInput.AdditionDevice.AdditionDeviceType
        SYRINGE: ReactionInput.AdditionDevice.AdditionDeviceType
        CANNULA: ReactionInput.AdditionDevice.AdditionDeviceType
        ADDITION_FUNNEL: ReactionInput.AdditionDevice.AdditionDeviceType
        PIPETTE: ReactionInput.AdditionDevice.AdditionDeviceType
        POSITIVE_DISPLACEMENT_PIPETTE: ReactionInput.AdditionDevice.AdditionDeviceType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: ReactionInput.AdditionDevice.AdditionDeviceType
        details: str
        def __init__(self, type: _Optional[_Union[ReactionInput.AdditionDevice.AdditionDeviceType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    COMPONENTS_FIELD_NUMBER: _ClassVar[int]
    CRUDE_COMPONENTS_FIELD_NUMBER: _ClassVar[int]
    ADDITION_ORDER_FIELD_NUMBER: _ClassVar[int]
    ADDITION_TIME_FIELD_NUMBER: _ClassVar[int]
    ADDITION_SPEED_FIELD_NUMBER: _ClassVar[int]
    ADDITION_DURATION_FIELD_NUMBER: _ClassVar[int]
    FLOW_RATE_FIELD_NUMBER: _ClassVar[int]
    ADDITION_DEVICE_FIELD_NUMBER: _ClassVar[int]
    ADDITION_TEMPERATURE_FIELD_NUMBER: _ClassVar[int]
    TEXTURE_FIELD_NUMBER: _ClassVar[int]
    components: _containers.RepeatedCompositeFieldContainer[Compound]
    crude_components: _containers.RepeatedCompositeFieldContainer[CrudeComponent]
    addition_order: int
    addition_time: Time
    addition_speed: ReactionInput.AdditionSpeed
    addition_duration: Time
    flow_rate: FlowRate
    addition_device: ReactionInput.AdditionDevice
    addition_temperature: Temperature
    texture: Texture
    def __init__(self, components: _Optional[_Iterable[_Union[Compound, _Mapping]]] = ..., crude_components: _Optional[_Iterable[_Union[CrudeComponent, _Mapping]]] = ..., addition_order: _Optional[int] = ..., addition_time: _Optional[_Union[Time, _Mapping]] = ..., addition_speed: _Optional[_Union[ReactionInput.AdditionSpeed, _Mapping]] = ..., addition_duration: _Optional[_Union[Time, _Mapping]] = ..., flow_rate: _Optional[_Union[FlowRate, _Mapping]] = ..., addition_device: _Optional[_Union[ReactionInput.AdditionDevice, _Mapping]] = ..., addition_temperature: _Optional[_Union[Temperature, _Mapping]] = ..., texture: _Optional[_Union[Texture, _Mapping]] = ...) -> None: ...

class Amount(_message.Message):
    __slots__ = ["mass", "moles", "volume", "unmeasured", "volume_includes_solutes"]
    MASS_FIELD_NUMBER: _ClassVar[int]
    MOLES_FIELD_NUMBER: _ClassVar[int]
    VOLUME_FIELD_NUMBER: _ClassVar[int]
    UNMEASURED_FIELD_NUMBER: _ClassVar[int]
    VOLUME_INCLUDES_SOLUTES_FIELD_NUMBER: _ClassVar[int]
    mass: Mass
    moles: Moles
    volume: Volume
    unmeasured: UnmeasuredAmount
    volume_includes_solutes: bool
    def __init__(self, mass: _Optional[_Union[Mass, _Mapping]] = ..., moles: _Optional[_Union[Moles, _Mapping]] = ..., volume: _Optional[_Union[Volume, _Mapping]] = ..., unmeasured: _Optional[_Union[UnmeasuredAmount, _Mapping]] = ..., volume_includes_solutes: bool = ...) -> None: ...

class UnmeasuredAmount(_message.Message):
    __slots__ = ["type", "details"]
    class UnmeasuredAmountType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[UnmeasuredAmount.UnmeasuredAmountType]
        CUSTOM: _ClassVar[UnmeasuredAmount.UnmeasuredAmountType]
        SATURATED: _ClassVar[UnmeasuredAmount.UnmeasuredAmountType]
        CATALYTIC: _ClassVar[UnmeasuredAmount.UnmeasuredAmountType]
        TITRATED: _ClassVar[UnmeasuredAmount.UnmeasuredAmountType]
    UNSPECIFIED: UnmeasuredAmount.UnmeasuredAmountType
    CUSTOM: UnmeasuredAmount.UnmeasuredAmountType
    SATURATED: UnmeasuredAmount.UnmeasuredAmountType
    CATALYTIC: UnmeasuredAmount.UnmeasuredAmountType
    TITRATED: UnmeasuredAmount.UnmeasuredAmountType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    type: UnmeasuredAmount.UnmeasuredAmountType
    details: str
    def __init__(self, type: _Optional[_Union[UnmeasuredAmount.UnmeasuredAmountType, str]] = ..., details: _Optional[str] = ...) -> None: ...

class Texture(_message.Message):
    __slots__ = ["type", "details"]
    class TextureType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Texture.TextureType]
        CUSTOM: _ClassVar[Texture.TextureType]
        POWDER: _ClassVar[Texture.TextureType]
        CRYSTAL: _ClassVar[Texture.TextureType]
        OIL: _ClassVar[Texture.TextureType]
        AMORPHOUS_SOLID: _ClassVar[Texture.TextureType]
        FOAM: _ClassVar[Texture.TextureType]
        WAX: _ClassVar[Texture.TextureType]
        SEMI_SOLID: _ClassVar[Texture.TextureType]
        SOLID: _ClassVar[Texture.TextureType]
        LIQUID: _ClassVar[Texture.TextureType]
        GAS: _ClassVar[Texture.TextureType]
    UNSPECIFIED: Texture.TextureType
    CUSTOM: Texture.TextureType
    POWDER: Texture.TextureType
    CRYSTAL: Texture.TextureType
    OIL: Texture.TextureType
    AMORPHOUS_SOLID: Texture.TextureType
    FOAM: Texture.TextureType
    WAX: Texture.TextureType
    SEMI_SOLID: Texture.TextureType
    SOLID: Texture.TextureType
    LIQUID: Texture.TextureType
    GAS: Texture.TextureType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    type: Texture.TextureType
    details: str
    def __init__(self, type: _Optional[_Union[Texture.TextureType, str]] = ..., details: _Optional[str] = ...) -> None: ...

class CrudeComponent(_message.Message):
    __slots__ = ["reaction_id", "includes_workup", "has_derived_amount", "amount", "texture"]
    REACTION_ID_FIELD_NUMBER: _ClassVar[int]
    INCLUDES_WORKUP_FIELD_NUMBER: _ClassVar[int]
    HAS_DERIVED_AMOUNT_FIELD_NUMBER: _ClassVar[int]
    AMOUNT_FIELD_NUMBER: _ClassVar[int]
    TEXTURE_FIELD_NUMBER: _ClassVar[int]
    reaction_id: str
    includes_workup: bool
    has_derived_amount: bool
    amount: Amount
    texture: Texture
    def __init__(self, reaction_id: _Optional[str] = ..., includes_workup: bool = ..., has_derived_amount: bool = ..., amount: _Optional[_Union[Amount, _Mapping]] = ..., texture: _Optional[_Union[Texture, _Mapping]] = ...) -> None: ...

class Compound(_message.Message):
    __slots__ = ["identifiers", "amount", "reaction_role", "is_limiting", "preparations", "source", "features", "analyses", "texture"]
    class Source(_message.Message):
        __slots__ = ["vendor", "catalog_id", "lot"]
        VENDOR_FIELD_NUMBER: _ClassVar[int]
        CATALOG_ID_FIELD_NUMBER: _ClassVar[int]
        LOT_FIELD_NUMBER: _ClassVar[int]
        vendor: str
        catalog_id: str
        lot: str
        def __init__(self, vendor: _Optional[str] = ..., catalog_id: _Optional[str] = ..., lot: _Optional[str] = ...) -> None: ...
    class FeaturesEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Data
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Data, _Mapping]] = ...) -> None: ...
    class AnalysesEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Analysis
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Analysis, _Mapping]] = ...) -> None: ...
    IDENTIFIERS_FIELD_NUMBER: _ClassVar[int]
    AMOUNT_FIELD_NUMBER: _ClassVar[int]
    REACTION_ROLE_FIELD_NUMBER: _ClassVar[int]
    IS_LIMITING_FIELD_NUMBER: _ClassVar[int]
    PREPARATIONS_FIELD_NUMBER: _ClassVar[int]
    SOURCE_FIELD_NUMBER: _ClassVar[int]
    FEATURES_FIELD_NUMBER: _ClassVar[int]
    ANALYSES_FIELD_NUMBER: _ClassVar[int]
    TEXTURE_FIELD_NUMBER: _ClassVar[int]
    identifiers: _containers.RepeatedCompositeFieldContainer[CompoundIdentifier]
    amount: Amount
    reaction_role: ReactionRole.ReactionRoleType
    is_limiting: bool
    preparations: _containers.RepeatedCompositeFieldContainer[CompoundPreparation]
    source: Compound.Source
    features: _containers.MessageMap[str, Data]
    analyses: _containers.MessageMap[str, Analysis]
    texture: Texture
    def __init__(self, identifiers: _Optional[_Iterable[_Union[CompoundIdentifier, _Mapping]]] = ..., amount: _Optional[_Union[Amount, _Mapping]] = ..., reaction_role: _Optional[_Union[ReactionRole.ReactionRoleType, str]] = ..., is_limiting: bool = ..., preparations: _Optional[_Iterable[_Union[CompoundPreparation, _Mapping]]] = ..., source: _Optional[_Union[Compound.Source, _Mapping]] = ..., features: _Optional[_Mapping[str, Data]] = ..., analyses: _Optional[_Mapping[str, Analysis]] = ..., texture: _Optional[_Union[Texture, _Mapping]] = ...) -> None: ...

class ReactionRole(_message.Message):
    __slots__ = []
    class ReactionRoleType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[ReactionRole.ReactionRoleType]
        REACTANT: _ClassVar[ReactionRole.ReactionRoleType]
        REAGENT: _ClassVar[ReactionRole.ReactionRoleType]
        SOLVENT: _ClassVar[ReactionRole.ReactionRoleType]
        CATALYST: _ClassVar[ReactionRole.ReactionRoleType]
        WORKUP: _ClassVar[ReactionRole.ReactionRoleType]
        INTERNAL_STANDARD: _ClassVar[ReactionRole.ReactionRoleType]
        AUTHENTIC_STANDARD: _ClassVar[ReactionRole.ReactionRoleType]
        PRODUCT: _ClassVar[ReactionRole.ReactionRoleType]
        BYPRODUCT: _ClassVar[ReactionRole.ReactionRoleType]
        SIDE_PRODUCT: _ClassVar[ReactionRole.ReactionRoleType]
    UNSPECIFIED: ReactionRole.ReactionRoleType
    REACTANT: ReactionRole.ReactionRoleType
    REAGENT: ReactionRole.ReactionRoleType
    SOLVENT: ReactionRole.ReactionRoleType
    CATALYST: ReactionRole.ReactionRoleType
    WORKUP: ReactionRole.ReactionRoleType
    INTERNAL_STANDARD: ReactionRole.ReactionRoleType
    AUTHENTIC_STANDARD: ReactionRole.ReactionRoleType
    PRODUCT: ReactionRole.ReactionRoleType
    BYPRODUCT: ReactionRole.ReactionRoleType
    SIDE_PRODUCT: ReactionRole.ReactionRoleType
    def __init__(self) -> None: ...

class CompoundPreparation(_message.Message):
    __slots__ = ["type", "details", "reaction_id"]
    class CompoundPreparationType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[CompoundPreparation.CompoundPreparationType]
        CUSTOM: _ClassVar[CompoundPreparation.CompoundPreparationType]
        NONE: _ClassVar[CompoundPreparation.CompoundPreparationType]
        REPURIFIED: _ClassVar[CompoundPreparation.CompoundPreparationType]
        SPARGED: _ClassVar[CompoundPreparation.CompoundPreparationType]
        DRIED: _ClassVar[CompoundPreparation.CompoundPreparationType]
        SYNTHESIZED: _ClassVar[CompoundPreparation.CompoundPreparationType]
    UNSPECIFIED: CompoundPreparation.CompoundPreparationType
    CUSTOM: CompoundPreparation.CompoundPreparationType
    NONE: CompoundPreparation.CompoundPreparationType
    REPURIFIED: CompoundPreparation.CompoundPreparationType
    SPARGED: CompoundPreparation.CompoundPreparationType
    DRIED: CompoundPreparation.CompoundPreparationType
    SYNTHESIZED: CompoundPreparation.CompoundPreparationType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    REACTION_ID_FIELD_NUMBER: _ClassVar[int]
    type: CompoundPreparation.CompoundPreparationType
    details: str
    reaction_id: str
    def __init__(self, type: _Optional[_Union[CompoundPreparation.CompoundPreparationType, str]] = ..., details: _Optional[str] = ..., reaction_id: _Optional[str] = ...) -> None: ...

class CompoundIdentifier(_message.Message):
    __slots__ = ["type", "details", "value"]
    class CompoundIdentifierType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        CUSTOM: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        SMILES: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        INCHI: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        MOLBLOCK: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        IUPAC_NAME: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        NAME: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        CAS_NUMBER: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        PUBCHEM_CID: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        CHEMSPIDER_ID: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        CXSMILES: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        INCHI_KEY: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        XYZ: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        UNIPROT_ID: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        PDB_ID: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        AMINO_ACID_SEQUENCE: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        HELM: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
        MDL: _ClassVar[CompoundIdentifier.CompoundIdentifierType]
    UNSPECIFIED: CompoundIdentifier.CompoundIdentifierType
    CUSTOM: CompoundIdentifier.CompoundIdentifierType
    SMILES: CompoundIdentifier.CompoundIdentifierType
    INCHI: CompoundIdentifier.CompoundIdentifierType
    MOLBLOCK: CompoundIdentifier.CompoundIdentifierType
    IUPAC_NAME: CompoundIdentifier.CompoundIdentifierType
    NAME: CompoundIdentifier.CompoundIdentifierType
    CAS_NUMBER: CompoundIdentifier.CompoundIdentifierType
    PUBCHEM_CID: CompoundIdentifier.CompoundIdentifierType
    CHEMSPIDER_ID: CompoundIdentifier.CompoundIdentifierType
    CXSMILES: CompoundIdentifier.CompoundIdentifierType
    INCHI_KEY: CompoundIdentifier.CompoundIdentifierType
    XYZ: CompoundIdentifier.CompoundIdentifierType
    UNIPROT_ID: CompoundIdentifier.CompoundIdentifierType
    PDB_ID: CompoundIdentifier.CompoundIdentifierType
    AMINO_ACID_SEQUENCE: CompoundIdentifier.CompoundIdentifierType
    HELM: CompoundIdentifier.CompoundIdentifierType
    MDL: CompoundIdentifier.CompoundIdentifierType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    type: CompoundIdentifier.CompoundIdentifierType
    details: str
    value: str
    def __init__(self, type: _Optional[_Union[CompoundIdentifier.CompoundIdentifierType, str]] = ..., details: _Optional[str] = ..., value: _Optional[str] = ...) -> None: ...

class Vessel(_message.Message):
    __slots__ = ["type", "details", "material", "preparations", "attachments", "volume", "vessel_id", "position", "row", "col"]
    class VesselType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Vessel.VesselType]
        CUSTOM: _ClassVar[Vessel.VesselType]
        ROUND_BOTTOM_FLASK: _ClassVar[Vessel.VesselType]
        VIAL: _ClassVar[Vessel.VesselType]
        WELL_PLATE: _ClassVar[Vessel.VesselType]
        MICROWAVE_VIAL: _ClassVar[Vessel.VesselType]
        TUBE: _ClassVar[Vessel.VesselType]
        CONTINUOUS_STIRRED_TANK_REACTOR: _ClassVar[Vessel.VesselType]
        PACKED_BED_REACTOR: _ClassVar[Vessel.VesselType]
        NMR_TUBE: _ClassVar[Vessel.VesselType]
        PRESSURE_FLASK: _ClassVar[Vessel.VesselType]
        PRESSURE_REACTOR: _ClassVar[Vessel.VesselType]
        ELECTROCHEMICAL_CELL: _ClassVar[Vessel.VesselType]
    UNSPECIFIED: Vessel.VesselType
    CUSTOM: Vessel.VesselType
    ROUND_BOTTOM_FLASK: Vessel.VesselType
    VIAL: Vessel.VesselType
    WELL_PLATE: Vessel.VesselType
    MICROWAVE_VIAL: Vessel.VesselType
    TUBE: Vessel.VesselType
    CONTINUOUS_STIRRED_TANK_REACTOR: Vessel.VesselType
    PACKED_BED_REACTOR: Vessel.VesselType
    NMR_TUBE: Vessel.VesselType
    PRESSURE_FLASK: Vessel.VesselType
    PRESSURE_REACTOR: Vessel.VesselType
    ELECTROCHEMICAL_CELL: Vessel.VesselType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    MATERIAL_FIELD_NUMBER: _ClassVar[int]
    PREPARATIONS_FIELD_NUMBER: _ClassVar[int]
    ATTACHMENTS_FIELD_NUMBER: _ClassVar[int]
    VOLUME_FIELD_NUMBER: _ClassVar[int]
    VESSEL_ID_FIELD_NUMBER: _ClassVar[int]
    POSITION_FIELD_NUMBER: _ClassVar[int]
    ROW_FIELD_NUMBER: _ClassVar[int]
    COL_FIELD_NUMBER: _ClassVar[int]
    type: Vessel.VesselType
    details: str
    material: VesselMaterial
    preparations: _containers.RepeatedCompositeFieldContainer[VesselPreparation]
    attachments: _containers.RepeatedCompositeFieldContainer[VesselAttachment]
    volume: Volume
    vessel_id: str
    position: str
    row: str
    col: str
    def __init__(self, type: _Optional[_Union[Vessel.VesselType, str]] = ..., details: _Optional[str] = ..., material: _Optional[_Union[VesselMaterial, _Mapping]] = ..., preparations: _Optional[_Iterable[_Union[VesselPreparation, _Mapping]]] = ..., attachments: _Optional[_Iterable[_Union[VesselAttachment, _Mapping]]] = ..., volume: _Optional[_Union[Volume, _Mapping]] = ..., vessel_id: _Optional[str] = ..., position: _Optional[str] = ..., row: _Optional[str] = ..., col: _Optional[str] = ...) -> None: ...

class VesselMaterial(_message.Message):
    __slots__ = ["type", "details"]
    class VesselMaterialType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[VesselMaterial.VesselMaterialType]
        CUSTOM: _ClassVar[VesselMaterial.VesselMaterialType]
        GLASS: _ClassVar[VesselMaterial.VesselMaterialType]
        POLYPROPYLENE: _ClassVar[VesselMaterial.VesselMaterialType]
        PLASTIC: _ClassVar[VesselMaterial.VesselMaterialType]
        METAL: _ClassVar[VesselMaterial.VesselMaterialType]
        QUARTZ: _ClassVar[VesselMaterial.VesselMaterialType]
    UNSPECIFIED: VesselMaterial.VesselMaterialType
    CUSTOM: VesselMaterial.VesselMaterialType
    GLASS: VesselMaterial.VesselMaterialType
    POLYPROPYLENE: VesselMaterial.VesselMaterialType
    PLASTIC: VesselMaterial.VesselMaterialType
    METAL: VesselMaterial.VesselMaterialType
    QUARTZ: VesselMaterial.VesselMaterialType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    type: VesselMaterial.VesselMaterialType
    details: str
    def __init__(self, type: _Optional[_Union[VesselMaterial.VesselMaterialType, str]] = ..., details: _Optional[str] = ...) -> None: ...

class VesselAttachment(_message.Message):
    __slots__ = ["type", "details"]
    class VesselAttachmentType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[VesselAttachment.VesselAttachmentType]
        NONE: _ClassVar[VesselAttachment.VesselAttachmentType]
        CUSTOM: _ClassVar[VesselAttachment.VesselAttachmentType]
        SEPTUM: _ClassVar[VesselAttachment.VesselAttachmentType]
        CAP: _ClassVar[VesselAttachment.VesselAttachmentType]
        MAT: _ClassVar[VesselAttachment.VesselAttachmentType]
        REFLUX_CONDENSER: _ClassVar[VesselAttachment.VesselAttachmentType]
        VENT_NEEDLE: _ClassVar[VesselAttachment.VesselAttachmentType]
        DEAN_STARK: _ClassVar[VesselAttachment.VesselAttachmentType]
        VACUUM_TUBE: _ClassVar[VesselAttachment.VesselAttachmentType]
        ADDITION_FUNNEL: _ClassVar[VesselAttachment.VesselAttachmentType]
        DRYING_TUBE: _ClassVar[VesselAttachment.VesselAttachmentType]
        ALUMINUM_FOIL: _ClassVar[VesselAttachment.VesselAttachmentType]
        THERMOCOUPLE: _ClassVar[VesselAttachment.VesselAttachmentType]
        BALLOON: _ClassVar[VesselAttachment.VesselAttachmentType]
        GAS_ADAPTER: _ClassVar[VesselAttachment.VesselAttachmentType]
        PRESSURE_REGULATOR: _ClassVar[VesselAttachment.VesselAttachmentType]
        RELEASE_VALVE: _ClassVar[VesselAttachment.VesselAttachmentType]
    UNSPECIFIED: VesselAttachment.VesselAttachmentType
    NONE: VesselAttachment.VesselAttachmentType
    CUSTOM: VesselAttachment.VesselAttachmentType
    SEPTUM: VesselAttachment.VesselAttachmentType
    CAP: VesselAttachment.VesselAttachmentType
    MAT: VesselAttachment.VesselAttachmentType
    REFLUX_CONDENSER: VesselAttachment.VesselAttachmentType
    VENT_NEEDLE: VesselAttachment.VesselAttachmentType
    DEAN_STARK: VesselAttachment.VesselAttachmentType
    VACUUM_TUBE: VesselAttachment.VesselAttachmentType
    ADDITION_FUNNEL: VesselAttachment.VesselAttachmentType
    DRYING_TUBE: VesselAttachment.VesselAttachmentType
    ALUMINUM_FOIL: VesselAttachment.VesselAttachmentType
    THERMOCOUPLE: VesselAttachment.VesselAttachmentType
    BALLOON: VesselAttachment.VesselAttachmentType
    GAS_ADAPTER: VesselAttachment.VesselAttachmentType
    PRESSURE_REGULATOR: VesselAttachment.VesselAttachmentType
    RELEASE_VALVE: VesselAttachment.VesselAttachmentType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    type: VesselAttachment.VesselAttachmentType
    details: str
    def __init__(self, type: _Optional[_Union[VesselAttachment.VesselAttachmentType, str]] = ..., details: _Optional[str] = ...) -> None: ...

class VesselPreparation(_message.Message):
    __slots__ = ["type", "details"]
    class VesselPreparationType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[VesselPreparation.VesselPreparationType]
        CUSTOM: _ClassVar[VesselPreparation.VesselPreparationType]
        NONE: _ClassVar[VesselPreparation.VesselPreparationType]
        OVEN_DRIED: _ClassVar[VesselPreparation.VesselPreparationType]
        FLAME_DRIED: _ClassVar[VesselPreparation.VesselPreparationType]
        EVACUATED_BACKFILLED: _ClassVar[VesselPreparation.VesselPreparationType]
        PURGED: _ClassVar[VesselPreparation.VesselPreparationType]
    UNSPECIFIED: VesselPreparation.VesselPreparationType
    CUSTOM: VesselPreparation.VesselPreparationType
    NONE: VesselPreparation.VesselPreparationType
    OVEN_DRIED: VesselPreparation.VesselPreparationType
    FLAME_DRIED: VesselPreparation.VesselPreparationType
    EVACUATED_BACKFILLED: VesselPreparation.VesselPreparationType
    PURGED: VesselPreparation.VesselPreparationType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    type: VesselPreparation.VesselPreparationType
    details: str
    def __init__(self, type: _Optional[_Union[VesselPreparation.VesselPreparationType, str]] = ..., details: _Optional[str] = ...) -> None: ...

class ReactionSetup(_message.Message):
    __slots__ = ["vessel", "is_automated", "automation_platform", "automation_code", "environment"]
    class AutomationCodeEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Data
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Data, _Mapping]] = ...) -> None: ...
    class ReactionEnvironment(_message.Message):
        __slots__ = ["type", "details"]
        class ReactionEnvironmentType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType]
            CUSTOM: _ClassVar[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType]
            FUME_HOOD: _ClassVar[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType]
            BENCH_TOP: _ClassVar[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType]
            GLOVE_BOX: _ClassVar[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType]
            GLOVE_BAG: _ClassVar[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType]
        UNSPECIFIED: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        CUSTOM: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        FUME_HOOD: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        BENCH_TOP: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        GLOVE_BOX: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        GLOVE_BAG: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
        details: str
        def __init__(self, type: _Optional[_Union[ReactionSetup.ReactionEnvironment.ReactionEnvironmentType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    VESSEL_FIELD_NUMBER: _ClassVar[int]
    IS_AUTOMATED_FIELD_NUMBER: _ClassVar[int]
    AUTOMATION_PLATFORM_FIELD_NUMBER: _ClassVar[int]
    AUTOMATION_CODE_FIELD_NUMBER: _ClassVar[int]
    ENVIRONMENT_FIELD_NUMBER: _ClassVar[int]
    vessel: Vessel
    is_automated: bool
    automation_platform: str
    automation_code: _containers.MessageMap[str, Data]
    environment: ReactionSetup.ReactionEnvironment
    def __init__(self, vessel: _Optional[_Union[Vessel, _Mapping]] = ..., is_automated: bool = ..., automation_platform: _Optional[str] = ..., automation_code: _Optional[_Mapping[str, Data]] = ..., environment: _Optional[_Union[ReactionSetup.ReactionEnvironment, _Mapping]] = ...) -> None: ...

class ReactionConditions(_message.Message):
    __slots__ = ["temperature", "pressure", "stirring", "illumination", "electrochemistry", "flow", "reflux", "ph", "conditions_are_dynamic", "details"]
    TEMPERATURE_FIELD_NUMBER: _ClassVar[int]
    PRESSURE_FIELD_NUMBER: _ClassVar[int]
    STIRRING_FIELD_NUMBER: _ClassVar[int]
    ILLUMINATION_FIELD_NUMBER: _ClassVar[int]
    ELECTROCHEMISTRY_FIELD_NUMBER: _ClassVar[int]
    FLOW_FIELD_NUMBER: _ClassVar[int]
    REFLUX_FIELD_NUMBER: _ClassVar[int]
    PH_FIELD_NUMBER: _ClassVar[int]
    CONDITIONS_ARE_DYNAMIC_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    temperature: TemperatureConditions
    pressure: PressureConditions
    stirring: StirringConditions
    illumination: IlluminationConditions
    electrochemistry: ElectrochemistryConditions
    flow: FlowConditions
    reflux: bool
    ph: float
    conditions_are_dynamic: bool
    details: str
    def __init__(self, temperature: _Optional[_Union[TemperatureConditions, _Mapping]] = ..., pressure: _Optional[_Union[PressureConditions, _Mapping]] = ..., stirring: _Optional[_Union[StirringConditions, _Mapping]] = ..., illumination: _Optional[_Union[IlluminationConditions, _Mapping]] = ..., electrochemistry: _Optional[_Union[ElectrochemistryConditions, _Mapping]] = ..., flow: _Optional[_Union[FlowConditions, _Mapping]] = ..., reflux: bool = ..., ph: _Optional[float] = ..., conditions_are_dynamic: bool = ..., details: _Optional[str] = ...) -> None: ...

class TemperatureConditions(_message.Message):
    __slots__ = ["control", "setpoint", "measurements"]
    class TemperatureControl(_message.Message):
        __slots__ = ["type", "details"]
        class TemperatureControlType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            CUSTOM: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            AMBIENT: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            OIL_BATH: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            WATER_BATH: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            SAND_BATH: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            ICE_BATH: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            DRY_ALUMINUM_PLATE: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            MICROWAVE: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            DRY_ICE_BATH: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            AIR_FAN: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
            LIQUID_NITROGEN: _ClassVar[TemperatureConditions.TemperatureControl.TemperatureControlType]
        UNSPECIFIED: TemperatureConditions.TemperatureControl.TemperatureControlType
        CUSTOM: TemperatureConditions.TemperatureControl.TemperatureControlType
        AMBIENT: TemperatureConditions.TemperatureControl.TemperatureControlType
        OIL_BATH: TemperatureConditions.TemperatureControl.TemperatureControlType
        WATER_BATH: TemperatureConditions.TemperatureControl.TemperatureControlType
        SAND_BATH: TemperatureConditions.TemperatureControl.TemperatureControlType
        ICE_BATH: TemperatureConditions.TemperatureControl.TemperatureControlType
        DRY_ALUMINUM_PLATE: TemperatureConditions.TemperatureControl.TemperatureControlType
        MICROWAVE: TemperatureConditions.TemperatureControl.TemperatureControlType
        DRY_ICE_BATH: TemperatureConditions.TemperatureControl.TemperatureControlType
        AIR_FAN: TemperatureConditions.TemperatureControl.TemperatureControlType
        LIQUID_NITROGEN: TemperatureConditions.TemperatureControl.TemperatureControlType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: TemperatureConditions.TemperatureControl.TemperatureControlType
        details: str
        def __init__(self, type: _Optional[_Union[TemperatureConditions.TemperatureControl.TemperatureControlType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    class TemperatureMeasurement(_message.Message):
        __slots__ = ["type", "details", "time", "temperature"]
        class TemperatureMeasurementType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType]
            CUSTOM: _ClassVar[TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType]
            THERMOCOUPLE_INTERNAL: _ClassVar[TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType]
            THERMOCOUPLE_EXTERNAL: _ClassVar[TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType]
            INFRARED: _ClassVar[TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType]
        UNSPECIFIED: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType
        CUSTOM: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType
        THERMOCOUPLE_INTERNAL: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType
        THERMOCOUPLE_EXTERNAL: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType
        INFRARED: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        TIME_FIELD_NUMBER: _ClassVar[int]
        TEMPERATURE_FIELD_NUMBER: _ClassVar[int]
        type: TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType
        details: str
        time: Time
        temperature: Temperature
        def __init__(self, type: _Optional[_Union[TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType, str]] = ..., details: _Optional[str] = ..., time: _Optional[_Union[Time, _Mapping]] = ..., temperature: _Optional[_Union[Temperature, _Mapping]] = ...) -> None: ...
    CONTROL_FIELD_NUMBER: _ClassVar[int]
    SETPOINT_FIELD_NUMBER: _ClassVar[int]
    MEASUREMENTS_FIELD_NUMBER: _ClassVar[int]
    control: TemperatureConditions.TemperatureControl
    setpoint: Temperature
    measurements: _containers.RepeatedCompositeFieldContainer[TemperatureConditions.TemperatureMeasurement]
    def __init__(self, control: _Optional[_Union[TemperatureConditions.TemperatureControl, _Mapping]] = ..., setpoint: _Optional[_Union[Temperature, _Mapping]] = ..., measurements: _Optional[_Iterable[_Union[TemperatureConditions.TemperatureMeasurement, _Mapping]]] = ...) -> None: ...

class PressureConditions(_message.Message):
    __slots__ = ["control", "setpoint", "atmosphere", "measurements"]
    class PressureControl(_message.Message):
        __slots__ = ["type", "details"]
        class PressureControlType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[PressureConditions.PressureControl.PressureControlType]
            CUSTOM: _ClassVar[PressureConditions.PressureControl.PressureControlType]
            AMBIENT: _ClassVar[PressureConditions.PressureControl.PressureControlType]
            SLIGHT_POSITIVE: _ClassVar[PressureConditions.PressureControl.PressureControlType]
            SEALED: _ClassVar[PressureConditions.PressureControl.PressureControlType]
            PRESSURIZED: _ClassVar[PressureConditions.PressureControl.PressureControlType]
        UNSPECIFIED: PressureConditions.PressureControl.PressureControlType
        CUSTOM: PressureConditions.PressureControl.PressureControlType
        AMBIENT: PressureConditions.PressureControl.PressureControlType
        SLIGHT_POSITIVE: PressureConditions.PressureControl.PressureControlType
        SEALED: PressureConditions.PressureControl.PressureControlType
        PRESSURIZED: PressureConditions.PressureControl.PressureControlType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: PressureConditions.PressureControl.PressureControlType
        details: str
        def __init__(self, type: _Optional[_Union[PressureConditions.PressureControl.PressureControlType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    class Atmosphere(_message.Message):
        __slots__ = ["type", "details"]
        class AtmosphereType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            CUSTOM: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            AIR: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            NITROGEN: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            ARGON: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            OXYGEN: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            HYDROGEN: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            CARBON_MONOXIDE: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            CARBON_DIOXIDE: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            METHANE: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            AMMONIA: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            OZONE: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            ETHYLENE: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
            ACETYLENE: _ClassVar[PressureConditions.Atmosphere.AtmosphereType]
        UNSPECIFIED: PressureConditions.Atmosphere.AtmosphereType
        CUSTOM: PressureConditions.Atmosphere.AtmosphereType
        AIR: PressureConditions.Atmosphere.AtmosphereType
        NITROGEN: PressureConditions.Atmosphere.AtmosphereType
        ARGON: PressureConditions.Atmosphere.AtmosphereType
        OXYGEN: PressureConditions.Atmosphere.AtmosphereType
        HYDROGEN: PressureConditions.Atmosphere.AtmosphereType
        CARBON_MONOXIDE: PressureConditions.Atmosphere.AtmosphereType
        CARBON_DIOXIDE: PressureConditions.Atmosphere.AtmosphereType
        METHANE: PressureConditions.Atmosphere.AtmosphereType
        AMMONIA: PressureConditions.Atmosphere.AtmosphereType
        OZONE: PressureConditions.Atmosphere.AtmosphereType
        ETHYLENE: PressureConditions.Atmosphere.AtmosphereType
        ACETYLENE: PressureConditions.Atmosphere.AtmosphereType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: PressureConditions.Atmosphere.AtmosphereType
        details: str
        def __init__(self, type: _Optional[_Union[PressureConditions.Atmosphere.AtmosphereType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    class PressureMeasurement(_message.Message):
        __slots__ = ["type", "details", "time", "pressure"]
        class PressureMeasurementType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[PressureConditions.PressureMeasurement.PressureMeasurementType]
            CUSTOM: _ClassVar[PressureConditions.PressureMeasurement.PressureMeasurementType]
            PRESSURE_TRANSDUCER: _ClassVar[PressureConditions.PressureMeasurement.PressureMeasurementType]
        UNSPECIFIED: PressureConditions.PressureMeasurement.PressureMeasurementType
        CUSTOM: PressureConditions.PressureMeasurement.PressureMeasurementType
        PRESSURE_TRANSDUCER: PressureConditions.PressureMeasurement.PressureMeasurementType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        TIME_FIELD_NUMBER: _ClassVar[int]
        PRESSURE_FIELD_NUMBER: _ClassVar[int]
        type: PressureConditions.PressureMeasurement.PressureMeasurementType
        details: str
        time: Time
        pressure: Pressure
        def __init__(self, type: _Optional[_Union[PressureConditions.PressureMeasurement.PressureMeasurementType, str]] = ..., details: _Optional[str] = ..., time: _Optional[_Union[Time, _Mapping]] = ..., pressure: _Optional[_Union[Pressure, _Mapping]] = ...) -> None: ...
    CONTROL_FIELD_NUMBER: _ClassVar[int]
    SETPOINT_FIELD_NUMBER: _ClassVar[int]
    ATMOSPHERE_FIELD_NUMBER: _ClassVar[int]
    MEASUREMENTS_FIELD_NUMBER: _ClassVar[int]
    control: PressureConditions.PressureControl
    setpoint: Pressure
    atmosphere: PressureConditions.Atmosphere
    measurements: _containers.RepeatedCompositeFieldContainer[PressureConditions.PressureMeasurement]
    def __init__(self, control: _Optional[_Union[PressureConditions.PressureControl, _Mapping]] = ..., setpoint: _Optional[_Union[Pressure, _Mapping]] = ..., atmosphere: _Optional[_Union[PressureConditions.Atmosphere, _Mapping]] = ..., measurements: _Optional[_Iterable[_Union[PressureConditions.PressureMeasurement, _Mapping]]] = ...) -> None: ...

class StirringConditions(_message.Message):
    __slots__ = ["type", "details", "rate"]
    class StirringMethodType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[StirringConditions.StirringMethodType]
        CUSTOM: _ClassVar[StirringConditions.StirringMethodType]
        NONE: _ClassVar[StirringConditions.StirringMethodType]
        STIR_BAR: _ClassVar[StirringConditions.StirringMethodType]
        OVERHEAD_MIXER: _ClassVar[StirringConditions.StirringMethodType]
        AGITATION: _ClassVar[StirringConditions.StirringMethodType]
        BALL_MILLING: _ClassVar[StirringConditions.StirringMethodType]
        SONICATION: _ClassVar[StirringConditions.StirringMethodType]
    UNSPECIFIED: StirringConditions.StirringMethodType
    CUSTOM: StirringConditions.StirringMethodType
    NONE: StirringConditions.StirringMethodType
    STIR_BAR: StirringConditions.StirringMethodType
    OVERHEAD_MIXER: StirringConditions.StirringMethodType
    AGITATION: StirringConditions.StirringMethodType
    BALL_MILLING: StirringConditions.StirringMethodType
    SONICATION: StirringConditions.StirringMethodType
    class StirringRate(_message.Message):
        __slots__ = ["type", "details", "rpm"]
        class StirringRateType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[StirringConditions.StirringRate.StirringRateType]
            HIGH: _ClassVar[StirringConditions.StirringRate.StirringRateType]
            MEDIUM: _ClassVar[StirringConditions.StirringRate.StirringRateType]
            LOW: _ClassVar[StirringConditions.StirringRate.StirringRateType]
        UNSPECIFIED: StirringConditions.StirringRate.StirringRateType
        HIGH: StirringConditions.StirringRate.StirringRateType
        MEDIUM: StirringConditions.StirringRate.StirringRateType
        LOW: StirringConditions.StirringRate.StirringRateType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        RPM_FIELD_NUMBER: _ClassVar[int]
        type: StirringConditions.StirringRate.StirringRateType
        details: str
        rpm: int
        def __init__(self, type: _Optional[_Union[StirringConditions.StirringRate.StirringRateType, str]] = ..., details: _Optional[str] = ..., rpm: _Optional[int] = ...) -> None: ...
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    RATE_FIELD_NUMBER: _ClassVar[int]
    type: StirringConditions.StirringMethodType
    details: str
    rate: StirringConditions.StirringRate
    def __init__(self, type: _Optional[_Union[StirringConditions.StirringMethodType, str]] = ..., details: _Optional[str] = ..., rate: _Optional[_Union[StirringConditions.StirringRate, _Mapping]] = ...) -> None: ...

class IlluminationConditions(_message.Message):
    __slots__ = ["type", "details", "peak_wavelength", "color", "distance_to_vessel"]
    class IlluminationType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[IlluminationConditions.IlluminationType]
        CUSTOM: _ClassVar[IlluminationConditions.IlluminationType]
        AMBIENT: _ClassVar[IlluminationConditions.IlluminationType]
        DARK: _ClassVar[IlluminationConditions.IlluminationType]
        LED: _ClassVar[IlluminationConditions.IlluminationType]
        HALOGEN_LAMP: _ClassVar[IlluminationConditions.IlluminationType]
        DEUTERIUM_LAMP: _ClassVar[IlluminationConditions.IlluminationType]
        SOLAR_SIMULATOR: _ClassVar[IlluminationConditions.IlluminationType]
        BROAD_SPECTRUM: _ClassVar[IlluminationConditions.IlluminationType]
    UNSPECIFIED: IlluminationConditions.IlluminationType
    CUSTOM: IlluminationConditions.IlluminationType
    AMBIENT: IlluminationConditions.IlluminationType
    DARK: IlluminationConditions.IlluminationType
    LED: IlluminationConditions.IlluminationType
    HALOGEN_LAMP: IlluminationConditions.IlluminationType
    DEUTERIUM_LAMP: IlluminationConditions.IlluminationType
    SOLAR_SIMULATOR: IlluminationConditions.IlluminationType
    BROAD_SPECTRUM: IlluminationConditions.IlluminationType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    PEAK_WAVELENGTH_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    DISTANCE_TO_VESSEL_FIELD_NUMBER: _ClassVar[int]
    type: IlluminationConditions.IlluminationType
    details: str
    peak_wavelength: Wavelength
    color: str
    distance_to_vessel: Length
    def __init__(self, type: _Optional[_Union[IlluminationConditions.IlluminationType, str]] = ..., details: _Optional[str] = ..., peak_wavelength: _Optional[_Union[Wavelength, _Mapping]] = ..., color: _Optional[str] = ..., distance_to_vessel: _Optional[_Union[Length, _Mapping]] = ...) -> None: ...

class ElectrochemistryConditions(_message.Message):
    __slots__ = ["type", "details", "current", "voltage", "anode_material", "cathode_material", "electrode_separation", "measurements", "cell"]
    class ElectrochemistryType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[ElectrochemistryConditions.ElectrochemistryType]
        CUSTOM: _ClassVar[ElectrochemistryConditions.ElectrochemistryType]
        CONSTANT_CURRENT: _ClassVar[ElectrochemistryConditions.ElectrochemistryType]
        CONSTANT_VOLTAGE: _ClassVar[ElectrochemistryConditions.ElectrochemistryType]
    UNSPECIFIED: ElectrochemistryConditions.ElectrochemistryType
    CUSTOM: ElectrochemistryConditions.ElectrochemistryType
    CONSTANT_CURRENT: ElectrochemistryConditions.ElectrochemistryType
    CONSTANT_VOLTAGE: ElectrochemistryConditions.ElectrochemistryType
    class ElectrochemistryMeasurement(_message.Message):
        __slots__ = ["time", "current", "voltage"]
        TIME_FIELD_NUMBER: _ClassVar[int]
        CURRENT_FIELD_NUMBER: _ClassVar[int]
        VOLTAGE_FIELD_NUMBER: _ClassVar[int]
        time: Time
        current: Current
        voltage: Voltage
        def __init__(self, time: _Optional[_Union[Time, _Mapping]] = ..., current: _Optional[_Union[Current, _Mapping]] = ..., voltage: _Optional[_Union[Voltage, _Mapping]] = ...) -> None: ...
    class ElectrochemistryCell(_message.Message):
        __slots__ = ["type", "details"]
        class ElectrochemistryCellType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType]
            CUSTOM: _ClassVar[ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType]
            DIVIDED_CELL: _ClassVar[ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType]
            UNDIVIDED_CELL: _ClassVar[ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType]
        UNSPECIFIED: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType
        CUSTOM: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType
        DIVIDED_CELL: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType
        UNDIVIDED_CELL: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType
        details: str
        def __init__(self, type: _Optional[_Union[ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    CURRENT_FIELD_NUMBER: _ClassVar[int]
    VOLTAGE_FIELD_NUMBER: _ClassVar[int]
    ANODE_MATERIAL_FIELD_NUMBER: _ClassVar[int]
    CATHODE_MATERIAL_FIELD_NUMBER: _ClassVar[int]
    ELECTRODE_SEPARATION_FIELD_NUMBER: _ClassVar[int]
    MEASUREMENTS_FIELD_NUMBER: _ClassVar[int]
    CELL_FIELD_NUMBER: _ClassVar[int]
    type: ElectrochemistryConditions.ElectrochemistryType
    details: str
    current: Current
    voltage: Voltage
    anode_material: str
    cathode_material: str
    electrode_separation: Length
    measurements: _containers.RepeatedCompositeFieldContainer[ElectrochemistryConditions.ElectrochemistryMeasurement]
    cell: ElectrochemistryConditions.ElectrochemistryCell
    def __init__(self, type: _Optional[_Union[ElectrochemistryConditions.ElectrochemistryType, str]] = ..., details: _Optional[str] = ..., current: _Optional[_Union[Current, _Mapping]] = ..., voltage: _Optional[_Union[Voltage, _Mapping]] = ..., anode_material: _Optional[str] = ..., cathode_material: _Optional[str] = ..., electrode_separation: _Optional[_Union[Length, _Mapping]] = ..., measurements: _Optional[_Iterable[_Union[ElectrochemistryConditions.ElectrochemistryMeasurement, _Mapping]]] = ..., cell: _Optional[_Union[ElectrochemistryConditions.ElectrochemistryCell, _Mapping]] = ...) -> None: ...

class FlowConditions(_message.Message):
    __slots__ = ["type", "details", "pump_type", "tubing"]
    class FlowType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[FlowConditions.FlowType]
        CUSTOM: _ClassVar[FlowConditions.FlowType]
        PLUG_FLOW_REACTOR: _ClassVar[FlowConditions.FlowType]
        CONTINUOUS_STIRRED_TANK_REACTOR: _ClassVar[FlowConditions.FlowType]
        PACKED_BED_REACTOR: _ClassVar[FlowConditions.FlowType]
    UNSPECIFIED: FlowConditions.FlowType
    CUSTOM: FlowConditions.FlowType
    PLUG_FLOW_REACTOR: FlowConditions.FlowType
    CONTINUOUS_STIRRED_TANK_REACTOR: FlowConditions.FlowType
    PACKED_BED_REACTOR: FlowConditions.FlowType
    class Tubing(_message.Message):
        __slots__ = ["type", "details", "diameter"]
        class TubingType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[FlowConditions.Tubing.TubingType]
            CUSTOM: _ClassVar[FlowConditions.Tubing.TubingType]
            STEEL: _ClassVar[FlowConditions.Tubing.TubingType]
            COPPER: _ClassVar[FlowConditions.Tubing.TubingType]
            PFA: _ClassVar[FlowConditions.Tubing.TubingType]
            FEP: _ClassVar[FlowConditions.Tubing.TubingType]
            TEFLONAF: _ClassVar[FlowConditions.Tubing.TubingType]
            PTFE: _ClassVar[FlowConditions.Tubing.TubingType]
            GLASS: _ClassVar[FlowConditions.Tubing.TubingType]
            QUARTZ: _ClassVar[FlowConditions.Tubing.TubingType]
            SILICON: _ClassVar[FlowConditions.Tubing.TubingType]
            PDMS: _ClassVar[FlowConditions.Tubing.TubingType]
        UNSPECIFIED: FlowConditions.Tubing.TubingType
        CUSTOM: FlowConditions.Tubing.TubingType
        STEEL: FlowConditions.Tubing.TubingType
        COPPER: FlowConditions.Tubing.TubingType
        PFA: FlowConditions.Tubing.TubingType
        FEP: FlowConditions.Tubing.TubingType
        TEFLONAF: FlowConditions.Tubing.TubingType
        PTFE: FlowConditions.Tubing.TubingType
        GLASS: FlowConditions.Tubing.TubingType
        QUARTZ: FlowConditions.Tubing.TubingType
        SILICON: FlowConditions.Tubing.TubingType
        PDMS: FlowConditions.Tubing.TubingType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        DIAMETER_FIELD_NUMBER: _ClassVar[int]
        type: FlowConditions.Tubing.TubingType
        details: str
        diameter: Length
        def __init__(self, type: _Optional[_Union[FlowConditions.Tubing.TubingType, str]] = ..., details: _Optional[str] = ..., diameter: _Optional[_Union[Length, _Mapping]] = ...) -> None: ...
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    PUMP_TYPE_FIELD_NUMBER: _ClassVar[int]
    TUBING_FIELD_NUMBER: _ClassVar[int]
    type: FlowConditions.FlowType
    details: str
    pump_type: str
    tubing: FlowConditions.Tubing
    def __init__(self, type: _Optional[_Union[FlowConditions.FlowType, str]] = ..., details: _Optional[str] = ..., pump_type: _Optional[str] = ..., tubing: _Optional[_Union[FlowConditions.Tubing, _Mapping]] = ...) -> None: ...

class ReactionNotes(_message.Message):
    __slots__ = ["is_heterogeneous", "forms_precipitate", "is_exothermic", "offgasses", "is_sensitive_to_moisture", "is_sensitive_to_oxygen", "is_sensitive_to_light", "safety_notes", "procedure_details"]
    IS_HETEROGENEOUS_FIELD_NUMBER: _ClassVar[int]
    FORMS_PRECIPITATE_FIELD_NUMBER: _ClassVar[int]
    IS_EXOTHERMIC_FIELD_NUMBER: _ClassVar[int]
    OFFGASSES_FIELD_NUMBER: _ClassVar[int]
    IS_SENSITIVE_TO_MOISTURE_FIELD_NUMBER: _ClassVar[int]
    IS_SENSITIVE_TO_OXYGEN_FIELD_NUMBER: _ClassVar[int]
    IS_SENSITIVE_TO_LIGHT_FIELD_NUMBER: _ClassVar[int]
    SAFETY_NOTES_FIELD_NUMBER: _ClassVar[int]
    PROCEDURE_DETAILS_FIELD_NUMBER: _ClassVar[int]
    is_heterogeneous: bool
    forms_precipitate: bool
    is_exothermic: bool
    offgasses: bool
    is_sensitive_to_moisture: bool
    is_sensitive_to_oxygen: bool
    is_sensitive_to_light: bool
    safety_notes: str
    procedure_details: str
    def __init__(self, is_heterogeneous: bool = ..., forms_precipitate: bool = ..., is_exothermic: bool = ..., offgasses: bool = ..., is_sensitive_to_moisture: bool = ..., is_sensitive_to_oxygen: bool = ..., is_sensitive_to_light: bool = ..., safety_notes: _Optional[str] = ..., procedure_details: _Optional[str] = ...) -> None: ...

class ReactionObservation(_message.Message):
    __slots__ = ["time", "comment", "image"]
    TIME_FIELD_NUMBER: _ClassVar[int]
    COMMENT_FIELD_NUMBER: _ClassVar[int]
    IMAGE_FIELD_NUMBER: _ClassVar[int]
    time: Time
    comment: str
    image: Data
    def __init__(self, time: _Optional[_Union[Time, _Mapping]] = ..., comment: _Optional[str] = ..., image: _Optional[_Union[Data, _Mapping]] = ...) -> None: ...

class ReactionWorkup(_message.Message):
    __slots__ = ["type", "details", "duration", "input", "amount", "temperature", "keep_phase", "stirring", "target_ph", "is_automated"]
    class ReactionWorkupType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[ReactionWorkup.ReactionWorkupType]
        CUSTOM: _ClassVar[ReactionWorkup.ReactionWorkupType]
        ADDITION: _ClassVar[ReactionWorkup.ReactionWorkupType]
        ALIQUOT: _ClassVar[ReactionWorkup.ReactionWorkupType]
        TEMPERATURE: _ClassVar[ReactionWorkup.ReactionWorkupType]
        CONCENTRATION: _ClassVar[ReactionWorkup.ReactionWorkupType]
        EXTRACTION: _ClassVar[ReactionWorkup.ReactionWorkupType]
        FILTRATION: _ClassVar[ReactionWorkup.ReactionWorkupType]
        WASH: _ClassVar[ReactionWorkup.ReactionWorkupType]
        DRY_IN_VACUUM: _ClassVar[ReactionWorkup.ReactionWorkupType]
        DRY_WITH_MATERIAL: _ClassVar[ReactionWorkup.ReactionWorkupType]
        FLASH_CHROMATOGRAPHY: _ClassVar[ReactionWorkup.ReactionWorkupType]
        OTHER_CHROMATOGRAPHY: _ClassVar[ReactionWorkup.ReactionWorkupType]
        SCAVENGING: _ClassVar[ReactionWorkup.ReactionWorkupType]
        WAIT: _ClassVar[ReactionWorkup.ReactionWorkupType]
        STIRRING: _ClassVar[ReactionWorkup.ReactionWorkupType]
        PH_ADJUST: _ClassVar[ReactionWorkup.ReactionWorkupType]
        DISSOLUTION: _ClassVar[ReactionWorkup.ReactionWorkupType]
        DISTILLATION: _ClassVar[ReactionWorkup.ReactionWorkupType]
    UNSPECIFIED: ReactionWorkup.ReactionWorkupType
    CUSTOM: ReactionWorkup.ReactionWorkupType
    ADDITION: ReactionWorkup.ReactionWorkupType
    ALIQUOT: ReactionWorkup.ReactionWorkupType
    TEMPERATURE: ReactionWorkup.ReactionWorkupType
    CONCENTRATION: ReactionWorkup.ReactionWorkupType
    EXTRACTION: ReactionWorkup.ReactionWorkupType
    FILTRATION: ReactionWorkup.ReactionWorkupType
    WASH: ReactionWorkup.ReactionWorkupType
    DRY_IN_VACUUM: ReactionWorkup.ReactionWorkupType
    DRY_WITH_MATERIAL: ReactionWorkup.ReactionWorkupType
    FLASH_CHROMATOGRAPHY: ReactionWorkup.ReactionWorkupType
    OTHER_CHROMATOGRAPHY: ReactionWorkup.ReactionWorkupType
    SCAVENGING: ReactionWorkup.ReactionWorkupType
    WAIT: ReactionWorkup.ReactionWorkupType
    STIRRING: ReactionWorkup.ReactionWorkupType
    PH_ADJUST: ReactionWorkup.ReactionWorkupType
    DISSOLUTION: ReactionWorkup.ReactionWorkupType
    DISTILLATION: ReactionWorkup.ReactionWorkupType
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    DURATION_FIELD_NUMBER: _ClassVar[int]
    INPUT_FIELD_NUMBER: _ClassVar[int]
    AMOUNT_FIELD_NUMBER: _ClassVar[int]
    TEMPERATURE_FIELD_NUMBER: _ClassVar[int]
    KEEP_PHASE_FIELD_NUMBER: _ClassVar[int]
    STIRRING_FIELD_NUMBER: _ClassVar[int]
    TARGET_PH_FIELD_NUMBER: _ClassVar[int]
    IS_AUTOMATED_FIELD_NUMBER: _ClassVar[int]
    type: ReactionWorkup.ReactionWorkupType
    details: str
    duration: Time
    input: ReactionInput
    amount: Amount
    temperature: TemperatureConditions
    keep_phase: str
    stirring: StirringConditions
    target_ph: float
    is_automated: bool
    def __init__(self, type: _Optional[_Union[ReactionWorkup.ReactionWorkupType, str]] = ..., details: _Optional[str] = ..., duration: _Optional[_Union[Time, _Mapping]] = ..., input: _Optional[_Union[ReactionInput, _Mapping]] = ..., amount: _Optional[_Union[Amount, _Mapping]] = ..., temperature: _Optional[_Union[TemperatureConditions, _Mapping]] = ..., keep_phase: _Optional[str] = ..., stirring: _Optional[_Union[StirringConditions, _Mapping]] = ..., target_ph: _Optional[float] = ..., is_automated: bool = ...) -> None: ...

class ReactionOutcome(_message.Message):
    __slots__ = ["reaction_time", "conversion", "products", "analyses"]
    class AnalysesEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Analysis
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Analysis, _Mapping]] = ...) -> None: ...
    REACTION_TIME_FIELD_NUMBER: _ClassVar[int]
    CONVERSION_FIELD_NUMBER: _ClassVar[int]
    PRODUCTS_FIELD_NUMBER: _ClassVar[int]
    ANALYSES_FIELD_NUMBER: _ClassVar[int]
    reaction_time: Time
    conversion: Percentage
    products: _containers.RepeatedCompositeFieldContainer[ProductCompound]
    analyses: _containers.MessageMap[str, Analysis]
    def __init__(self, reaction_time: _Optional[_Union[Time, _Mapping]] = ..., conversion: _Optional[_Union[Percentage, _Mapping]] = ..., products: _Optional[_Iterable[_Union[ProductCompound, _Mapping]]] = ..., analyses: _Optional[_Mapping[str, Analysis]] = ...) -> None: ...

class ProductCompound(_message.Message):
    __slots__ = ["identifiers", "is_desired_product", "measurements", "isolated_color", "texture", "features", "reaction_role"]
    class FeaturesEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Data
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Data, _Mapping]] = ...) -> None: ...
    IDENTIFIERS_FIELD_NUMBER: _ClassVar[int]
    IS_DESIRED_PRODUCT_FIELD_NUMBER: _ClassVar[int]
    MEASUREMENTS_FIELD_NUMBER: _ClassVar[int]
    ISOLATED_COLOR_FIELD_NUMBER: _ClassVar[int]
    TEXTURE_FIELD_NUMBER: _ClassVar[int]
    FEATURES_FIELD_NUMBER: _ClassVar[int]
    REACTION_ROLE_FIELD_NUMBER: _ClassVar[int]
    identifiers: _containers.RepeatedCompositeFieldContainer[CompoundIdentifier]
    is_desired_product: bool
    measurements: _containers.RepeatedCompositeFieldContainer[ProductMeasurement]
    isolated_color: str
    texture: Texture
    features: _containers.MessageMap[str, Data]
    reaction_role: ReactionRole.ReactionRoleType
    def __init__(self, identifiers: _Optional[_Iterable[_Union[CompoundIdentifier, _Mapping]]] = ..., is_desired_product: bool = ..., measurements: _Optional[_Iterable[_Union[ProductMeasurement, _Mapping]]] = ..., isolated_color: _Optional[str] = ..., texture: _Optional[_Union[Texture, _Mapping]] = ..., features: _Optional[_Mapping[str, Data]] = ..., reaction_role: _Optional[_Union[ReactionRole.ReactionRoleType, str]] = ...) -> None: ...

class ProductMeasurement(_message.Message):
    __slots__ = ["analysis_key", "type", "details", "uses_internal_standard", "is_normalized", "uses_authentic_standard", "authentic_standard", "percentage", "float_value", "string_value", "amount", "retention_time", "mass_spec_details", "selectivity", "wavelength"]
    class ProductMeasurementType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[ProductMeasurement.ProductMeasurementType]
        CUSTOM: _ClassVar[ProductMeasurement.ProductMeasurementType]
        IDENTITY: _ClassVar[ProductMeasurement.ProductMeasurementType]
        YIELD: _ClassVar[ProductMeasurement.ProductMeasurementType]
        SELECTIVITY: _ClassVar[ProductMeasurement.ProductMeasurementType]
        PURITY: _ClassVar[ProductMeasurement.ProductMeasurementType]
        AREA: _ClassVar[ProductMeasurement.ProductMeasurementType]
        COUNTS: _ClassVar[ProductMeasurement.ProductMeasurementType]
        INTENSITY: _ClassVar[ProductMeasurement.ProductMeasurementType]
        AMOUNT: _ClassVar[ProductMeasurement.ProductMeasurementType]
    UNSPECIFIED: ProductMeasurement.ProductMeasurementType
    CUSTOM: ProductMeasurement.ProductMeasurementType
    IDENTITY: ProductMeasurement.ProductMeasurementType
    YIELD: ProductMeasurement.ProductMeasurementType
    SELECTIVITY: ProductMeasurement.ProductMeasurementType
    PURITY: ProductMeasurement.ProductMeasurementType
    AREA: ProductMeasurement.ProductMeasurementType
    COUNTS: ProductMeasurement.ProductMeasurementType
    INTENSITY: ProductMeasurement.ProductMeasurementType
    AMOUNT: ProductMeasurement.ProductMeasurementType
    class MassSpecMeasurementDetails(_message.Message):
        __slots__ = ["type", "details", "tic_minimum_mz", "tic_maximum_mz", "eic_masses"]
        class MassSpecMeasurementType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType]
            CUSTOM: _ClassVar[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType]
            TIC: _ClassVar[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType]
            TIC_POSITIVE: _ClassVar[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType]
            TIC_NEGATIVE: _ClassVar[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType]
            EIC: _ClassVar[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType]
        UNSPECIFIED: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        CUSTOM: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        TIC: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        TIC_POSITIVE: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        TIC_NEGATIVE: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        EIC: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        TIC_MINIMUM_MZ_FIELD_NUMBER: _ClassVar[int]
        TIC_MAXIMUM_MZ_FIELD_NUMBER: _ClassVar[int]
        EIC_MASSES_FIELD_NUMBER: _ClassVar[int]
        type: ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType
        details: str
        tic_minimum_mz: float
        tic_maximum_mz: float
        eic_masses: _containers.RepeatedScalarFieldContainer[float]
        def __init__(self, type: _Optional[_Union[ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType, str]] = ..., details: _Optional[str] = ..., tic_minimum_mz: _Optional[float] = ..., tic_maximum_mz: _Optional[float] = ..., eic_masses: _Optional[_Iterable[float]] = ...) -> None: ...
    class Selectivity(_message.Message):
        __slots__ = ["type", "details"]
        class SelectivityType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
            __slots__ = []
            UNSPECIFIED: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
            CUSTOM: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
            EE: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
            ER: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
            DR: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
            EZ: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
            ZE: _ClassVar[ProductMeasurement.Selectivity.SelectivityType]
        UNSPECIFIED: ProductMeasurement.Selectivity.SelectivityType
        CUSTOM: ProductMeasurement.Selectivity.SelectivityType
        EE: ProductMeasurement.Selectivity.SelectivityType
        ER: ProductMeasurement.Selectivity.SelectivityType
        DR: ProductMeasurement.Selectivity.SelectivityType
        EZ: ProductMeasurement.Selectivity.SelectivityType
        ZE: ProductMeasurement.Selectivity.SelectivityType
        TYPE_FIELD_NUMBER: _ClassVar[int]
        DETAILS_FIELD_NUMBER: _ClassVar[int]
        type: ProductMeasurement.Selectivity.SelectivityType
        details: str
        def __init__(self, type: _Optional[_Union[ProductMeasurement.Selectivity.SelectivityType, str]] = ..., details: _Optional[str] = ...) -> None: ...
    ANALYSIS_KEY_FIELD_NUMBER: _ClassVar[int]
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    USES_INTERNAL_STANDARD_FIELD_NUMBER: _ClassVar[int]
    IS_NORMALIZED_FIELD_NUMBER: _ClassVar[int]
    USES_AUTHENTIC_STANDARD_FIELD_NUMBER: _ClassVar[int]
    AUTHENTIC_STANDARD_FIELD_NUMBER: _ClassVar[int]
    PERCENTAGE_FIELD_NUMBER: _ClassVar[int]
    FLOAT_VALUE_FIELD_NUMBER: _ClassVar[int]
    STRING_VALUE_FIELD_NUMBER: _ClassVar[int]
    AMOUNT_FIELD_NUMBER: _ClassVar[int]
    RETENTION_TIME_FIELD_NUMBER: _ClassVar[int]
    MASS_SPEC_DETAILS_FIELD_NUMBER: _ClassVar[int]
    SELECTIVITY_FIELD_NUMBER: _ClassVar[int]
    WAVELENGTH_FIELD_NUMBER: _ClassVar[int]
    analysis_key: str
    type: ProductMeasurement.ProductMeasurementType
    details: str
    uses_internal_standard: bool
    is_normalized: bool
    uses_authentic_standard: bool
    authentic_standard: Compound
    percentage: Percentage
    float_value: FloatValue
    string_value: str
    amount: Amount
    retention_time: Time
    mass_spec_details: ProductMeasurement.MassSpecMeasurementDetails
    selectivity: ProductMeasurement.Selectivity
    wavelength: Wavelength
    def __init__(self, analysis_key: _Optional[str] = ..., type: _Optional[_Union[ProductMeasurement.ProductMeasurementType, str]] = ..., details: _Optional[str] = ..., uses_internal_standard: bool = ..., is_normalized: bool = ..., uses_authentic_standard: bool = ..., authentic_standard: _Optional[_Union[Compound, _Mapping]] = ..., percentage: _Optional[_Union[Percentage, _Mapping]] = ..., float_value: _Optional[_Union[FloatValue, _Mapping]] = ..., string_value: _Optional[str] = ..., amount: _Optional[_Union[Amount, _Mapping]] = ..., retention_time: _Optional[_Union[Time, _Mapping]] = ..., mass_spec_details: _Optional[_Union[ProductMeasurement.MassSpecMeasurementDetails, _Mapping]] = ..., selectivity: _Optional[_Union[ProductMeasurement.Selectivity, _Mapping]] = ..., wavelength: _Optional[_Union[Wavelength, _Mapping]] = ...) -> None: ...

class DateTime(_message.Message):
    __slots__ = ["value"]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    value: str
    def __init__(self, value: _Optional[str] = ...) -> None: ...

class Analysis(_message.Message):
    __slots__ = ["type", "details", "chmo_id", "is_of_isolated_species", "data", "instrument_manufacturer", "instrument_last_calibrated"]
    class AnalysisType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Analysis.AnalysisType]
        CUSTOM: _ClassVar[Analysis.AnalysisType]
        LC: _ClassVar[Analysis.AnalysisType]
        GC: _ClassVar[Analysis.AnalysisType]
        IR: _ClassVar[Analysis.AnalysisType]
        NMR_1H: _ClassVar[Analysis.AnalysisType]
        NMR_13C: _ClassVar[Analysis.AnalysisType]
        NMR_OTHER: _ClassVar[Analysis.AnalysisType]
        MP: _ClassVar[Analysis.AnalysisType]
        UV: _ClassVar[Analysis.AnalysisType]
        TLC: _ClassVar[Analysis.AnalysisType]
        MS: _ClassVar[Analysis.AnalysisType]
        HRMS: _ClassVar[Analysis.AnalysisType]
        MSMS: _ClassVar[Analysis.AnalysisType]
        WEIGHT: _ClassVar[Analysis.AnalysisType]
        LCMS: _ClassVar[Analysis.AnalysisType]
        GCMS: _ClassVar[Analysis.AnalysisType]
        ELSD: _ClassVar[Analysis.AnalysisType]
        CD: _ClassVar[Analysis.AnalysisType]
        SFC: _ClassVar[Analysis.AnalysisType]
        EPR: _ClassVar[Analysis.AnalysisType]
        XRD: _ClassVar[Analysis.AnalysisType]
        RAMAN: _ClassVar[Analysis.AnalysisType]
        ED: _ClassVar[Analysis.AnalysisType]
        OPTICAL_ROTATION: _ClassVar[Analysis.AnalysisType]
        CAD: _ClassVar[Analysis.AnalysisType]
    UNSPECIFIED: Analysis.AnalysisType
    CUSTOM: Analysis.AnalysisType
    LC: Analysis.AnalysisType
    GC: Analysis.AnalysisType
    IR: Analysis.AnalysisType
    NMR_1H: Analysis.AnalysisType
    NMR_13C: Analysis.AnalysisType
    NMR_OTHER: Analysis.AnalysisType
    MP: Analysis.AnalysisType
    UV: Analysis.AnalysisType
    TLC: Analysis.AnalysisType
    MS: Analysis.AnalysisType
    HRMS: Analysis.AnalysisType
    MSMS: Analysis.AnalysisType
    WEIGHT: Analysis.AnalysisType
    LCMS: Analysis.AnalysisType
    GCMS: Analysis.AnalysisType
    ELSD: Analysis.AnalysisType
    CD: Analysis.AnalysisType
    SFC: Analysis.AnalysisType
    EPR: Analysis.AnalysisType
    XRD: Analysis.AnalysisType
    RAMAN: Analysis.AnalysisType
    ED: Analysis.AnalysisType
    OPTICAL_ROTATION: Analysis.AnalysisType
    CAD: Analysis.AnalysisType
    class DataEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Data
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Data, _Mapping]] = ...) -> None: ...
    TYPE_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    CHMO_ID_FIELD_NUMBER: _ClassVar[int]
    IS_OF_ISOLATED_SPECIES_FIELD_NUMBER: _ClassVar[int]
    DATA_FIELD_NUMBER: _ClassVar[int]
    INSTRUMENT_MANUFACTURER_FIELD_NUMBER: _ClassVar[int]
    INSTRUMENT_LAST_CALIBRATED_FIELD_NUMBER: _ClassVar[int]
    type: Analysis.AnalysisType
    details: str
    chmo_id: int
    is_of_isolated_species: bool
    data: _containers.MessageMap[str, Data]
    instrument_manufacturer: str
    instrument_last_calibrated: DateTime
    def __init__(self, type: _Optional[_Union[Analysis.AnalysisType, str]] = ..., details: _Optional[str] = ..., chmo_id: _Optional[int] = ..., is_of_isolated_species: bool = ..., data: _Optional[_Mapping[str, Data]] = ..., instrument_manufacturer: _Optional[str] = ..., instrument_last_calibrated: _Optional[_Union[DateTime, _Mapping]] = ...) -> None: ...

class ReactionProvenance(_message.Message):
    __slots__ = ["experimenter", "city", "experiment_start", "doi", "patent", "publication_url", "record_created", "record_modified", "reaction_metadata", "is_mined"]
    class ReactionMetadataEntry(_message.Message):
        __slots__ = ["key", "value"]
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: Data
        def __init__(self, key: _Optional[str] = ..., value: _Optional[_Union[Data, _Mapping]] = ...) -> None: ...
    EXPERIMENTER_FIELD_NUMBER: _ClassVar[int]
    CITY_FIELD_NUMBER: _ClassVar[int]
    EXPERIMENT_START_FIELD_NUMBER: _ClassVar[int]
    DOI_FIELD_NUMBER: _ClassVar[int]
    PATENT_FIELD_NUMBER: _ClassVar[int]
    PUBLICATION_URL_FIELD_NUMBER: _ClassVar[int]
    RECORD_CREATED_FIELD_NUMBER: _ClassVar[int]
    RECORD_MODIFIED_FIELD_NUMBER: _ClassVar[int]
    REACTION_METADATA_FIELD_NUMBER: _ClassVar[int]
    IS_MINED_FIELD_NUMBER: _ClassVar[int]
    experimenter: Person
    city: str
    experiment_start: DateTime
    doi: str
    patent: str
    publication_url: str
    record_created: RecordEvent
    record_modified: _containers.RepeatedCompositeFieldContainer[RecordEvent]
    reaction_metadata: _containers.MessageMap[str, Data]
    is_mined: bool
    def __init__(self, experimenter: _Optional[_Union[Person, _Mapping]] = ..., city: _Optional[str] = ..., experiment_start: _Optional[_Union[DateTime, _Mapping]] = ..., doi: _Optional[str] = ..., patent: _Optional[str] = ..., publication_url: _Optional[str] = ..., record_created: _Optional[_Union[RecordEvent, _Mapping]] = ..., record_modified: _Optional[_Iterable[_Union[RecordEvent, _Mapping]]] = ..., reaction_metadata: _Optional[_Mapping[str, Data]] = ..., is_mined: bool = ...) -> None: ...

class Person(_message.Message):
    __slots__ = ["username", "name", "orcid", "organization", "email"]
    USERNAME_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    ORCID_FIELD_NUMBER: _ClassVar[int]
    ORGANIZATION_FIELD_NUMBER: _ClassVar[int]
    EMAIL_FIELD_NUMBER: _ClassVar[int]
    username: str
    name: str
    orcid: str
    organization: str
    email: str
    def __init__(self, username: _Optional[str] = ..., name: _Optional[str] = ..., orcid: _Optional[str] = ..., organization: _Optional[str] = ..., email: _Optional[str] = ...) -> None: ...

class RecordEvent(_message.Message):
    __slots__ = ["time", "person", "details"]
    TIME_FIELD_NUMBER: _ClassVar[int]
    PERSON_FIELD_NUMBER: _ClassVar[int]
    DETAILS_FIELD_NUMBER: _ClassVar[int]
    time: DateTime
    person: Person
    details: str
    def __init__(self, time: _Optional[_Union[DateTime, _Mapping]] = ..., person: _Optional[_Union[Person, _Mapping]] = ..., details: _Optional[str] = ...) -> None: ...

class Time(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class TimeUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Time.TimeUnit]
        DAY: _ClassVar[Time.TimeUnit]
        HOUR: _ClassVar[Time.TimeUnit]
        MINUTE: _ClassVar[Time.TimeUnit]
        SECOND: _ClassVar[Time.TimeUnit]
    UNSPECIFIED: Time.TimeUnit
    DAY: Time.TimeUnit
    HOUR: Time.TimeUnit
    MINUTE: Time.TimeUnit
    SECOND: Time.TimeUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Time.TimeUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Time.TimeUnit, str]] = ...) -> None: ...

class Mass(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class MassUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Mass.MassUnit]
        KILOGRAM: _ClassVar[Mass.MassUnit]
        GRAM: _ClassVar[Mass.MassUnit]
        MILLIGRAM: _ClassVar[Mass.MassUnit]
        MICROGRAM: _ClassVar[Mass.MassUnit]
    UNSPECIFIED: Mass.MassUnit
    KILOGRAM: Mass.MassUnit
    GRAM: Mass.MassUnit
    MILLIGRAM: Mass.MassUnit
    MICROGRAM: Mass.MassUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Mass.MassUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Mass.MassUnit, str]] = ...) -> None: ...

class Moles(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class MolesUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Moles.MolesUnit]
        MOLE: _ClassVar[Moles.MolesUnit]
        MILLIMOLE: _ClassVar[Moles.MolesUnit]
        MICROMOLE: _ClassVar[Moles.MolesUnit]
        NANOMOLE: _ClassVar[Moles.MolesUnit]
    UNSPECIFIED: Moles.MolesUnit
    MOLE: Moles.MolesUnit
    MILLIMOLE: Moles.MolesUnit
    MICROMOLE: Moles.MolesUnit
    NANOMOLE: Moles.MolesUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Moles.MolesUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Moles.MolesUnit, str]] = ...) -> None: ...

class Volume(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class VolumeUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Volume.VolumeUnit]
        LITER: _ClassVar[Volume.VolumeUnit]
        MILLILITER: _ClassVar[Volume.VolumeUnit]
        MICROLITER: _ClassVar[Volume.VolumeUnit]
        NANOLITER: _ClassVar[Volume.VolumeUnit]
    UNSPECIFIED: Volume.VolumeUnit
    LITER: Volume.VolumeUnit
    MILLILITER: Volume.VolumeUnit
    MICROLITER: Volume.VolumeUnit
    NANOLITER: Volume.VolumeUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Volume.VolumeUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Volume.VolumeUnit, str]] = ...) -> None: ...

class Concentration(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class ConcentrationUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Concentration.ConcentrationUnit]
        MOLAR: _ClassVar[Concentration.ConcentrationUnit]
        MILLIMOLAR: _ClassVar[Concentration.ConcentrationUnit]
        MICROMOLAR: _ClassVar[Concentration.ConcentrationUnit]
    UNSPECIFIED: Concentration.ConcentrationUnit
    MOLAR: Concentration.ConcentrationUnit
    MILLIMOLAR: Concentration.ConcentrationUnit
    MICROMOLAR: Concentration.ConcentrationUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Concentration.ConcentrationUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Concentration.ConcentrationUnit, str]] = ...) -> None: ...

class Pressure(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class PressureUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Pressure.PressureUnit]
        BAR: _ClassVar[Pressure.PressureUnit]
        ATMOSPHERE: _ClassVar[Pressure.PressureUnit]
        PSI: _ClassVar[Pressure.PressureUnit]
        KPSI: _ClassVar[Pressure.PressureUnit]
        PASCAL: _ClassVar[Pressure.PressureUnit]
        KILOPASCAL: _ClassVar[Pressure.PressureUnit]
        TORR: _ClassVar[Pressure.PressureUnit]
        MM_HG: _ClassVar[Pressure.PressureUnit]
    UNSPECIFIED: Pressure.PressureUnit
    BAR: Pressure.PressureUnit
    ATMOSPHERE: Pressure.PressureUnit
    PSI: Pressure.PressureUnit
    KPSI: Pressure.PressureUnit
    PASCAL: Pressure.PressureUnit
    KILOPASCAL: Pressure.PressureUnit
    TORR: Pressure.PressureUnit
    MM_HG: Pressure.PressureUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Pressure.PressureUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Pressure.PressureUnit, str]] = ...) -> None: ...

class Temperature(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class TemperatureUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Temperature.TemperatureUnit]
        CELSIUS: _ClassVar[Temperature.TemperatureUnit]
        FAHRENHEIT: _ClassVar[Temperature.TemperatureUnit]
        KELVIN: _ClassVar[Temperature.TemperatureUnit]
    UNSPECIFIED: Temperature.TemperatureUnit
    CELSIUS: Temperature.TemperatureUnit
    FAHRENHEIT: Temperature.TemperatureUnit
    KELVIN: Temperature.TemperatureUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Temperature.TemperatureUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Temperature.TemperatureUnit, str]] = ...) -> None: ...

class Current(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class CurrentUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Current.CurrentUnit]
        AMPERE: _ClassVar[Current.CurrentUnit]
        MILLIAMPERE: _ClassVar[Current.CurrentUnit]
    UNSPECIFIED: Current.CurrentUnit
    AMPERE: Current.CurrentUnit
    MILLIAMPERE: Current.CurrentUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Current.CurrentUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Current.CurrentUnit, str]] = ...) -> None: ...

class Voltage(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class VoltageUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Voltage.VoltageUnit]
        VOLT: _ClassVar[Voltage.VoltageUnit]
        MILLIVOLT: _ClassVar[Voltage.VoltageUnit]
    UNSPECIFIED: Voltage.VoltageUnit
    VOLT: Voltage.VoltageUnit
    MILLIVOLT: Voltage.VoltageUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Voltage.VoltageUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Voltage.VoltageUnit, str]] = ...) -> None: ...

class Length(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class LengthUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Length.LengthUnit]
        CENTIMETER: _ClassVar[Length.LengthUnit]
        MILLIMETER: _ClassVar[Length.LengthUnit]
        METER: _ClassVar[Length.LengthUnit]
        INCH: _ClassVar[Length.LengthUnit]
        FOOT: _ClassVar[Length.LengthUnit]
    UNSPECIFIED: Length.LengthUnit
    CENTIMETER: Length.LengthUnit
    MILLIMETER: Length.LengthUnit
    METER: Length.LengthUnit
    INCH: Length.LengthUnit
    FOOT: Length.LengthUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Length.LengthUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Length.LengthUnit, str]] = ...) -> None: ...

class Wavelength(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class WavelengthUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[Wavelength.WavelengthUnit]
        NANOMETER: _ClassVar[Wavelength.WavelengthUnit]
        WAVENUMBER: _ClassVar[Wavelength.WavelengthUnit]
    UNSPECIFIED: Wavelength.WavelengthUnit
    NANOMETER: Wavelength.WavelengthUnit
    WAVENUMBER: Wavelength.WavelengthUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: Wavelength.WavelengthUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[Wavelength.WavelengthUnit, str]] = ...) -> None: ...

class FlowRate(_message.Message):
    __slots__ = ["value", "precision", "units"]
    class FlowRateUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = []
        UNSPECIFIED: _ClassVar[FlowRate.FlowRateUnit]
        MICROLITER_PER_MINUTE: _ClassVar[FlowRate.FlowRateUnit]
        MICROLITER_PER_SECOND: _ClassVar[FlowRate.FlowRateUnit]
        MILLILITER_PER_MINUTE: _ClassVar[FlowRate.FlowRateUnit]
        MILLILITER_PER_SECOND: _ClassVar[FlowRate.FlowRateUnit]
        MICROLITER_PER_HOUR: _ClassVar[FlowRate.FlowRateUnit]
    UNSPECIFIED: FlowRate.FlowRateUnit
    MICROLITER_PER_MINUTE: FlowRate.FlowRateUnit
    MICROLITER_PER_SECOND: FlowRate.FlowRateUnit
    MILLILITER_PER_MINUTE: FlowRate.FlowRateUnit
    MILLILITER_PER_SECOND: FlowRate.FlowRateUnit
    MICROLITER_PER_HOUR: FlowRate.FlowRateUnit
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    UNITS_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    units: FlowRate.FlowRateUnit
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ..., units: _Optional[_Union[FlowRate.FlowRateUnit, str]] = ...) -> None: ...

class Percentage(_message.Message):
    __slots__ = ["value", "precision"]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ...) -> None: ...

class FloatValue(_message.Message):
    __slots__ = ["value", "precision"]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    PRECISION_FIELD_NUMBER: _ClassVar[int]
    value: float
    precision: float
    def __init__(self, value: _Optional[float] = ..., precision: _Optional[float] = ...) -> None: ...

class Data(_message.Message):
    __slots__ = ["float_value", "integer_value", "bytes_value", "string_value", "url", "description", "format"]
    FLOAT_VALUE_FIELD_NUMBER: _ClassVar[int]
    INTEGER_VALUE_FIELD_NUMBER: _ClassVar[int]
    BYTES_VALUE_FIELD_NUMBER: _ClassVar[int]
    STRING_VALUE_FIELD_NUMBER: _ClassVar[int]
    URL_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    FORMAT_FIELD_NUMBER: _ClassVar[int]
    float_value: float
    integer_value: int
    bytes_value: bytes
    string_value: str
    url: str
    description: str
    format: str
    def __init__(self, float_value: _Optional[float] = ..., integer_value: _Optional[int] = ..., bytes_value: _Optional[bytes] = ..., string_value: _Optional[str] = ..., url: _Optional[str] = ..., description: _Optional[str] = ..., format: _Optional[str] = ...) -> None: ...
