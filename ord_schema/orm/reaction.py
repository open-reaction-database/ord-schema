"""Table mappings for Reaction protos.

Notes:
    * Foreign keys to the `reaction` table are done using the `id` column, not the ORD reaction ID (`name`).
    * Every table has its own `id` column, even if it has a one-to-one mapping with another table.
    * When a message type is used in multiple places, use Joined Table Inheritance; see
      https://docs.sqlalchemy.org/en/14/orm/inheritance.html#joined-table-inheritance.
    * The naming might be a bit confusing: classes for single-use and parent protos have the same name
      as the corresponding message, while classes for multi-use protos (child classes) are named as
      <ContainingClass><Attribute>. For example, Reaction.identifiers uses the ReactionIdentifier class,
      while Reaction.inputs uses the ReactionInputs class (since ReactionInput is used in multiple places).
      The effect is that all tables storing proto information are named after their corresponding proto
      messages, with extra child tables for storing relationships.
"""
from __future__ import annotations

from inspect import getmro

from inflection import underscore
from google.protobuf.message import Message
from sqlalchemy import (
    Boolean,
    Column,
    Enum,
    Float,
    Integer,
    ForeignKey,
    LargeBinary,
    String,
    Text,
)
from sqlalchemy.orm import declarative_base, declarative_mixin, declared_attr, relationship

from ord_schema.proto import reaction_pb2


def enum_to_str(message: Message, field: str) -> str:
    assert len(message.DESCRIPTOR.enum_types) == 1, message.DESCRIPTOR.name
    return message.DESCRIPTOR.enum_types[0].values_by_number[getattr(message, field)].name


class Base:
    """See https://docs.sqlalchemy.org/en/14/orm/declarative_mixins.html#augmenting-the-base."""

    @declared_attr
    def __tablename__(cls):
        return underscore(cls.__name__)

    id = Column(Integer, primary_key=True)


Base = declarative_base(cls=Base)


@declarative_mixin
class Parent:
    """Mixin class for parent classes."""

    _type = Column(String(255), nullable=False)

    @declared_attr
    def __mapper_args__(cls):
        return {
            "polymorphic_identity": underscore(cls.__name__),
            "polymorphic_on": cls._type,
        }


@declarative_mixin
class Child:
    """Mixin class for child classes."""

    child_id = Column(Integer, primary_key=True)  # Avoid conflict with parent `id` attr.

    @declared_attr
    def parent_id(cls):
        parent_table = None
        for parent_class in getmro(cls)[1:]:
            if issubclass(parent_class, Base):
                parent_table = parent_class.__tablename__
                break
        assert parent_table is not None, (cls, getmro(cls))
        return Column(f"{parent_table}_id", Integer, ForeignKey(f"{parent_table}.id"), unique=True)

    @declared_attr
    def __mapper_args__(cls):
        return {"polymorphic_identity": underscore(cls.__name__)}


class Reaction(Base):
    PROTO = reaction_pb2.Reaction

    name = Column(String(32), nullable=False)
    proto = Column(LargeBinary, nullable=False)

    identifiers = relationship("ReactionIdentifier")
    inputs = relationship("ReactionInputs")
    setup = relationship("ReactionSetup", uselist=False)
    conditions = relationship("ReactionConditions", uselist=False)
    notes = relationship("ReactionNotes", uselist=False)
    observations = relationship("ReactionObservation")
    workups = relationship("ReactionWorkup")
    outcomes = relationship("ReactionOutcome")
    provenance = relationship("ReactionProvenance", uselist=False)

    @classmethod
    def from_proto(cls, message: reaction_pb2.Reaction) -> Reaction:
        pass

    def to_proto(self) -> reaction_pb2.Reaction:
        pass


class ReactionIdentifier(Base):
    PROTO = reaction_pb2.ReactionIdentifier
    reaction_id = Column(Integer, ForeignKey("reaction.id"))

    type = Column(
        Enum(*reaction_pb2.ReactionIdentifier.IdentifierType.keys(), name="ReactionIdentifier.IdentifierType")
    )
    details = Column(Text)
    value = Column(Text)
    is_mapped = Column(Boolean)


class ReactionInput(Parent, Base):
    PROTO = reaction_pb2.ReactionInput

    components = relationship("ReactionInputComponents")
    crude_components = relationship("CrudeComponent")
    addition_order = Column(Integer)
    addition_time = relationship("ReactionInputAdditionTime", uselist=False)
    addition_speed = relationship("AdditionSpeed", uselist=False)
    addition_duration = relationship("ReactionInputAdditionDuration", uselist=False)
    flow_rate = relationship("FlowRate", uselist=False)
    addition_device = relationship("AdditionDevice", uselist=False)
    addition_temperature = relationship("ReactionInputAdditionTemperature", uselist=False)


class ReactionInputs(Child, ReactionInput):
    reaction_id = Column(Integer, ForeignKey("reaction.id"))
    name = Column(String(255))  # Map key.


class ReactionWorkupInput(Child, ReactionInput):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), unique=True)


class AdditionSpeed(Base):
    PROTO = reaction_pb2.ReactionInput.AdditionSpeed
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ReactionInput.AdditionSpeed.AdditionSpeedType.keys(),
            name="ReactionInput.AdditionSpeed.AdditionSpeedType",
        )
    )
    details = Column(Text)


class AdditionDevice(Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ReactionInput.AdditionDevice.AdditionDeviceType.keys(),
            name="ReactionInput.AdditionDevice.AdditionDeviceType",
        )
    )
    details = Column(Text)


class Amount(Parent, Base):
    PROTO = reaction_pb2.Amount

    mass = relationship("Mass", uselist=False)
    moles = relationship("Moles", uselist=False)
    volume = relationship("AmountVolume", uselist=False)
    unmeasured = relationship("UnmeasuredAmount", uselist=False)
    volume_includes_solutes = Column(Boolean)


class CrudeComponentAmount(Child, Amount):
    crude_component_id = Column(Integer, ForeignKey("crude_component.id"))


class CompoundAmount(Child, Amount):
    compound_id = Column(Integer, ForeignKey("compound.id"))


class ReactionWorkupAmount(Child, Amount):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"))


class ProductMeasurementAmount(Child, Amount):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"))


class UnmeasuredAmount(Base):
    PROTO = reaction_pb2.UnmeasuredAmount
    amount_id = Column(Integer, ForeignKey("amount.id"), unique=True)

    type = Column(
        Enum(*reaction_pb2.UnmeasuredAmount.UnmeasuredAmountType.keys(), name="UnmeasuredAmount.UnmeasuredAmountType")
    )
    details = Column(Text)


class CrudeComponent(Base):
    PROTO = reaction_pb2.CrudeComponent
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"))

    includes_workup = Column(Boolean)
    has_derived_amount = Column(Boolean)
    amount = relationship("CrudeComponentAmount", uselist=False)


# Shared enum; see https://docs.sqlalchemy.org/en/14/dialects/postgresql.html#sqlalchemy.dialects.postgresql.ENUM.
ReactionRoleType = Enum(
    *reaction_pb2.ReactionRole.ReactionRoleType.keys(), name="ReactionRole.ReactionRoleType", metadata=Base.metadata
)


class Compound(Parent, Base):
    PROTO = reaction_pb2.Compound

    identifiers = relationship("CompoundIdentifiers")
    amount = relationship("CompoundAmount")
    reaction_role = Column(ReactionRoleType)
    is_limiting = Column(Boolean)
    preparations = relationship("CompoundPreparation")
    source = relationship("Source", uselist=False)
    features = relationship("CompoundFeatures")
    analyses = relationship("CompoundAnalyses")


class Source(Base):
    PROTO = reaction_pb2.Compound.Source
    compound_id = Column(Integer, ForeignKey("compound.id"), unique=True)

    vendor = Column(String(255))
    vendor_id = Column(String(255))
    lot = Column(String(255))


class ReactionInputComponents(Child, Compound):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"))


class ProductMeasurementAuthenticStandard(Child, Compound):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)


class CompoundPreparation(Base):
    PROTO = reaction_pb2.CompoundPreparation
    compound_id = Column(Integer, ForeignKey("compound.id"))

    type = Column(
        Enum(*reaction_pb2.CompoundPreparation.PreparationType.keys(), name="CompoundPreparation.PreparationType")
    )
    details = Column(Text)
    reaction_id = Column(Integer, ForeignKey("reaction.id"))


class CompoundIdentifier(Parent, Base):
    PROTO = reaction_pb2.CompoundIdentifier

    type = Column(
        Enum(*reaction_pb2.CompoundIdentifier.IdentifierType.keys(), name="CompoundIdentifier.IdentifierType")
    )
    details = Column(Text)
    value = Column(Text)


class CompoundIdentifiers(Child, CompoundIdentifier):
    compound_id = Column(Integer, ForeignKey("compound.id"))


class ProductCompoundIdentifiers(Child, CompoundIdentifier):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"))


class Vessel(Base):
    PROTO = reaction_pb2.Vessel
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), unique=True)

    type = Column(Enum(*reaction_pb2.Vessel.VesselType.keys(), name="Vessel.VesselType"))
    details = Column(Text)
    material = relationship("VesselMaterial", uselist=False)
    preparations = relationship("VesselPreparation")
    attachments = relationship("VesselAttachment")
    volume = relationship("VesselVolume", uselist=False)
    plate_id = Column(String(255))
    plate_position = Column(String(32))


class VesselMaterial(Base):
    PROTO = reaction_pb2.VesselMaterial
    vessel_id = Column(Integer, ForeignKey("vessel.id"), unique=True)

    type = Column(
        Enum(*reaction_pb2.VesselMaterial.VesselMaterialType.keys(), name="VesselMaterial.VesselMaterialType")
    )
    details = Column(Text)


class VesselAttachment(Base):
    PROTO = reaction_pb2.VesselAttachment
    vessel_id = Column(Integer, ForeignKey("vessel.id"))

    type = Column(
        Enum(*reaction_pb2.VesselAttachment.VesselAttachmentType.keys(), name="VesselAttachment.VesselAttachmentType")
    )
    details = Column(Text)


class VesselPreparation(Base):
    PROTO = reaction_pb2.VesselPreparation
    vessel_id = Column(Integer, ForeignKey("vessel.id"))

    type = Column(
        Enum(
            *reaction_pb2.VesselPreparation.VesselPreparationType.keys(), name="VesselPreparation.VesselPreparationType"
        )
    )
    details = Column(Text)


class ReactionSetup(Base):
    PROTO = reaction_pb2.ReactionSetup
    reaction_id = Column(Integer, ForeignKey("reaction.id"), unique=True)

    vessel = relationship("Vessel", uselist=False)
    is_automated = Column(Boolean)
    automation_platform = Column(String(255))
    automation_code = relationship("ReactionSetupAutomationCode")
    environment = relationship("ReactionEnvironment", uselist=False)


class ReactionEnvironment(Base):
    PROTO = reaction_pb2.ReactionSetup.ReactionEnvironment
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType.keys(),
            name="ReactionSetup.ReactionEnvironment.ReactionEnvironmentType",
        )
    )
    details = Column(Text)


class ReactionConditions(Base):
    PROTO = reaction_pb2.ReactionConditions
    reaction_id = Column(Integer, ForeignKey("reaction.id"), unique=True)

    temperature = relationship("ReactionConditionsTemperature", uselist=False)
    pressure = relationship("PressureConditions", uselist=False)
    stirring = relationship("ReactionConditionsStirring", uselist=False)
    illumination = relationship("IlluminationConditions", uselist=False)
    electrochemistry = relationship("ElectrochemistryConditions", uselist=False)
    flow = relationship("FlowConditions", uselist=False)
    reflux = Column(Boolean)
    ph = Column(Float)
    conditions_are_dynamic = Column(Boolean)
    details = Column(Text)


class TemperatureConditions(Parent, Base):
    PROTO = reaction_pb2.TemperatureConditions

    control = relationship("TemperatureControl", uselist=False)
    setpoint = relationship("TemperatureConditionsSetpoint", uselist=False)
    measurements = relationship("TemperatureConditionsMeasurement")


class ReactionConditionsTemperature(Child, TemperatureConditions):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), unique=True)


class ReactionWorkupTemperature(Child, TemperatureConditions):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), unique=True)


class TemperatureControl(Base):
    PROTO = reaction_pb2.TemperatureConditions.TemperatureControl
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.TemperatureConditions.TemperatureControl.TemperatureControlType.keys(),
            name="TemperatureConditions.TemperatureControl.TemperatureControlType",
        )
    )
    details = Column(Text)


class TemperatureConditionsMeasurement(Base):
    PROTO = reaction_pb2.TemperatureConditions.Measurement
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"))

    type = Column(
        Enum(
            *reaction_pb2.TemperatureConditions.Measurement.MeasurementType.keys(),
            name="TemperatureConditions.Measurement.MeasurementType",
        )
    )
    details = Column(Text)
    time = relationship("TemperatureConditionsMeasurementTime", uselist=False)
    temperature = relationship("TemperatureConditionsMeasurementTemperature", uselist=False)


class PressureConditions(Base):
    PROTO = reaction_pb2.PressureConditions
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), unique=True)

    control = relationship("PressureControl", uselist=False)
    setpoint = relationship("PressureConditionsSetpoint", uselist=False)
    atmosphere = relationship("Atmosphere", uselist=False)
    measurements = relationship("PressureConditionsMeasurement")


class PressureControl(Base):
    PROTO = reaction_pb2.PressureConditions.PressureControl
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.PressureConditions.PressureControl.PressureControlType.keys(),
            name="PressureConditions.PressureControl.PressureControlType",
        )
    )
    details = Column(Text)


class Atmosphere(Base):
    PROTO = reaction_pb2.PressureConditions.Atmosphere
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.PressureConditions.Atmosphere.AtmosphereType.keys(),
            name="PressureConditions.Atmosphere.AtmosphereType",
        )
    )
    details = Column(Text)


class PressureConditionsMeasurement(Base):
    PROTO = reaction_pb2.PressureConditions.Measurement
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"))

    type = Column(
        Enum(
            *reaction_pb2.PressureConditions.Measurement.MeasurementType.keys(),
            name="PressureConditions.Measurement.MeasurementType",
        )
    )
    details = Column(Text)
    time = relationship("PressureConditionsMeasurementTime", uselist=False)
    pressure = relationship("PressureConditionsMeasurementPressure", uselist=False)


class StirringConditions(Parent, Base):
    PROTO = reaction_pb2.StirringConditions

    type = Column(
        Enum(*reaction_pb2.StirringConditions.StirringMethodType.keys(), name="StirringConditions.StirringMethodType")
    )
    details = Column(Text)
    rate = relationship("StirringRate", uselist=False)


class ReactionConditionsStirring(Child, StirringConditions):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), unique=True)


class ReactionWorkupStirring(Child, StirringConditions):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), unique=True)


class StirringRate(Base):
    PROTO = reaction_pb2.StirringConditions.StirringRate
    stirring_conditions_id = Column(Integer, ForeignKey("stirring_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.StirringConditions.StirringRate.StirringRateType.keys(),
            name="StirringConditions.StirringRate.StirringRateType",
        )
    )
    details = Column(Text)
    rpm = Column(Integer)


class IlluminationConditions(Base):
    PROTO = reaction_pb2.IlluminationConditions
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.IlluminationConditions.IlluminationType.keys(), name="IlluminationConditions.IlluminationType"
        )
    )
    details = Column(Text)
    peak_wavelength = relationship("IlluminationConditionsPeakWavelength", uselist=False)
    color = Column(String(255))
    distance_to_vessel = relationship("IlluminationConditionsDistanceToVessel", uselist=False)


class ElectrochemistryConditions(Base):
    PROTO = reaction_pb2.ElectrochemistryConditions
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ElectrochemistryConditions.ElectrochemistryType.keys(),
            name="ElectrochemistryConditions.ElectrochemistryType",
        )
    )
    details = Column(Text)
    current = relationship("ElectrochemistryConditionsCurrent", uselist=False)
    voltage = relationship("ElectrochemistryConditionsVoltage", uselist=False)
    anode_material = Column(String(255))
    cathode_material = Column(String(255))
    electrode_separation = relationship("ElectrochemistryConditionsElectrodeSeparation", uselist=False)
    measurements = relationship("ElectrochemistryConditionsMeasurement")
    cell = relationship("ElectrochemistryCell", uselist=False)


class ElectrochemistryConditionsMeasurement(Base):
    PROTO = reaction_pb2.ElectrochemistryConditions.Measurement
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"))

    time = relationship("ElectrochemistryConditionsMeasurementTime", uselist=False)
    current = relationship("ElectrochemistryConditionsMeasurementCurrent", uselist=False)
    voltage = relationship("ElectrochemistryConditionsMeasurementVoltage", uselist=False)


class ElectrochemistryCell(Base):
    PROTO = reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType.keys(),
            name="ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType",
        )
    )
    details = Column(Text)


class FlowConditions(Base):
    PROTO = reaction_pb2.FlowConditions
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), unique=True)

    type = Column(Enum(*reaction_pb2.FlowConditions.FlowType.keys(), name="FlowConditions.FlowType"))
    details = Column(Text)
    pump_type = Column(Text)
    tubing = relationship("Tubing", uselist=False)


class Tubing(Base):
    PROTO = reaction_pb2.FlowConditions.Tubing
    flow_conditions_id = Column(Integer, ForeignKey("flow_conditions.id"), unique=True)

    type = Column(Enum(*reaction_pb2.FlowConditions.Tubing.TubingType.keys(), name="FlowConditions.Tubing.TubingType"))
    details = Column(Text)
    diameter = relationship("TubingDiameter", uselist=False)


class ReactionNotes(Base):
    PROTO = reaction_pb2.ReactionNotes
    reaction_id = Column(Integer, ForeignKey("reaction.id"), unique=True)

    is_heterogeneous = Column(Boolean)
    forms_precipitate = Column(Boolean)
    is_exothermic = Column(Boolean)
    offgasses = Column(Boolean)
    is_sensitive_to_moisture = Column(Boolean)
    is_sensitive_to_oxygen = Column(Boolean)
    is_sensitive_to_light = Column(Boolean)
    safety_notes = Column(Text)
    procedure_details = Column(Text)


class ReactionObservation(Base):
    PROTO = reaction_pb2.ReactionObservation
    reaction_id = Column(Integer, ForeignKey("reaction.id"))

    time = relationship("ReactionObservationTime", uselist=False)
    comment = Column(Text)
    image = relationship("ReactionObservationImage", uselist=False)


class ReactionWorkup(Base):
    PROTO = reaction_pb2.ReactionWorkup
    reaction_id = Column(Integer, ForeignKey("reaction.id"))

    type = Column(Enum(*reaction_pb2.ReactionWorkup.WorkupType.keys(), name="ReactionWorkup.WorkupType"))
    details = Column(Text)
    duration = relationship("ReactionWorkupDuration", uselist=False)
    input = relationship("ReactionWorkupInput", uselist=False)
    amount = relationship("ReactionWorkupAmount", uselist=False)
    temperature = relationship("ReactionWorkupTemperature", uselist=False)
    keep_phase = Column(String(255))
    stirring = relationship("ReactionWorkupStirring", uselist=False)
    target_ph = Column(Float)
    is_automated = Column(Boolean)


class ReactionOutcome(Base):
    PROTO = reaction_pb2.ReactionOutcome
    reaction_id = Column(Integer, ForeignKey("reaction.id"))

    reaction_time = relationship("ReactionOutcomeReactionTime", uselist=False)
    conversion = relationship("ReactionOutcomeConversion", uselist=False)
    products = relationship("ProductCompound")
    analyses = relationship("ReactionOutcomeAnalyses")


class ProductCompound(Base):
    PROTO = reaction_pb2.ProductCompound
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"))

    identifiers = relationship("ProductCompoundIdentifiers")
    is_desired_product = Column(Boolean)
    measurements = relationship("ProductMeasurement")
    isolated_color = Column(String(255))
    texture = relationship("Texture", uselist=False)
    features = relationship("ProductCompoundFeatures")
    reaction_role = Column(ReactionRoleType)


class Texture(Base):
    PROTO = reaction_pb2.ProductCompound.Texture
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), unique=True)

    type = Column(
        Enum(*reaction_pb2.ProductCompound.Texture.TextureType.keys(), name="ProductCompound.Texture.TextureType")
    )
    details = Column(Text)


class ProductMeasurement(Base):
    PROTO = reaction_pb2.ProductMeasurement
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"))

    analysis_key = Column(String(255))
    type = Column(
        Enum(*reaction_pb2.ProductMeasurement.MeasurementType.keys(), name="ProductMeasurement.MeasurementType")
    )
    details = Column(Text)
    uses_internal_standard = Column(Boolean)
    is_normalized = Column(Boolean)
    uses_authentic_standard = Column(Boolean)
    authentic_standard = relationship("ProductMeasurementAuthenticStandard", uselist=False)
    percentage = relationship("ProductMeasurementPercentage", uselist=False)
    float_value = relationship("ProductMeasurementFloatValue", uselist=False)
    string_value = Column(Text)
    amount = relationship("ProductMeasurementAmount", uselist=False)
    retention_time = relationship("ProductMeasurementRetentionTime", uselist=False)
    mass_spec_details = relationship("MassSpecMeasurementDetails", uselist=False)
    selectivity = relationship("Selectivity", uselist=False)
    wavelength = relationship("ProductMeasurementWavelength", uselist=False)


class MassSpecMeasurementDetails(Base):
    PROTO = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType.keys(),
            name="ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType",
        )
    )
    details = Column(Text)
    tic_minimum_mz = Column(Float)
    tic_maximum_mz = Column(Float)
    eic_masses = relationship("MassSpecMeasurementDetailsEicMasses")


class Selectivity(Base):
    PROTO = reaction_pb2.ProductMeasurement.Selectivity
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ProductMeasurement.Selectivity.SelectivityType.keys(),
            name="ProductMeasurement.Selectivity.SelectivityType",
        )
    )
    details = Column(Text)


class DateTime(Parent, Base):
    PROTO = reaction_pb2.DateTime

    value = Column(String(255))


class AnalysisInstrumentLastCalibrated(Child, DateTime):
    analysis_id = Column(Integer, ForeignKey("analysis.id"), unique=True)


class ReactionProvenanceExperimentStart(Child, DateTime):
    analysis_id = Column(Integer, ForeignKey("reaction_provenance.id"), unique=True)


class RecordEventTime(Child, DateTime):
    analysis_id = Column(Integer, ForeignKey("record_event.id"), unique=True)


class Analysis(Parent, Base):
    PROTO = reaction_pb2.Analysis

    type = Column(Enum(*reaction_pb2.Analysis.AnalysisType.keys(), name="Analysis.AnalysisType"))
    details = Column(Text)
    chmo_id = Column(Integer)
    is_of_isolated_species = Column(Boolean)
    data = relationship("AnalysisData")
    instrument_manufacturer = Column(String(255))
    instrument_last_calibrated = relationship("AnalysisInstrumentLastCalibrated", uselist=False)


class CompoundAnalyses(Child, Analysis):
    compound_id = Column(Integer, ForeignKey("compound.id"))
    name = Column(String(255))  # Map key.


class ReactionOutcomeAnalyses(Child, Analysis):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"))
    name = Column(String(255))  # Map key.


class ReactionProvenance(Base):
    PROTO = reaction_pb2.ReactionProvenance
    reaction_id = Column(Integer, ForeignKey("reaction.id"), unique=True)

    experimenter = relationship("ReactionProvenanceExperimenter", uselist=False)
    city = Column(String(255))
    experiment_start = relationship("ReactionProvenanceExperimentStart", uselist=False)
    doi = Column(String(255))
    patent = Column(String(255))
    publication_url = Column(Text)
    record_created = relationship("ReactionProvenanceRecordCreated", uselist=False)
    record_modified = relationship("ReactionProvenanceRecordModified")

    @classmethod
    def from_proto(cls, message: reaction_pb2.ReactionProvenance) -> ReactionProvenance:
        return ReactionProvenance(
            experimenter=Person.from_proto(message.experimenter),
            city=message.city,
            experiment_start=DateTime.from_proto(message.experiment_start),
            doi=message.doi,
            patent=message.patent,
            publication_url=message.publication_url,
            record_created=RecordEvent.from_proto(message.record_created),
            record_modified=[RecordEvent.from_proto(event) for event in message.record_modified],
        )

    def to_proto(self) -> reaction_pb2.ReactionProvenance:
        return reaction_pb2.ReactionProvenance(username=self.username,
                                   name=self.name,
                                   orcid=self.orcid,
                                   organization=self.organization,
                                   email=self.email)


class Person(Parent, Base):
    PROTO = reaction_pb2.Person

    username = Column(String(255))
    name = Column(String(255))
    orcid = Column(String(19))
    organization = Column(String(255))
    email = Column(String(255))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Person) -> Person:
        return Person(
            username=message.username, name=message.name, orcid=message.orcid, organization=message.organization,
            email=message.email
        )

    def to_proto(self) -> reaction_pb2.Person:
        return reaction_pb2.Person(username=self.username,
                                   name=self.name,
                                   orcid=self.orcid,
                                   organization=self.organization,
                                   email=self.email)


class ReactionProvenanceExperimenter(Child, Person):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), unique=True)


class RecordEventPerson(Child, Person):
    record_event_id = Column(Integer, ForeignKey("record_event.id"), unique=True)


class RecordEvent(Parent, Base):
    PROTO = reaction_pb2.RecordEvent

    time = relationship("RecordEventTime", uselist=False)
    person = relationship("RecordEventPerson", uselist=False)
    details = Column(Text)

    @classmethod
    def from_proto(cls, message: reaction_pb2.RecordEvent) -> RecordEvent:
        return RecordEvent(
            time=Time.from_proto(message.time), person=Person.from_proto(message.person), details=message.details
        )

    def to_proto(self) -> reaction_pb2.RecordEvent:
        return reaction_pb2.RecordEvent(time=self.time.to_proto(), person=self.person.to_proto(), details=self.details)


class ReactionProvenanceRecordCreated(Child, RecordEvent):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), unique=True)


class ReactionProvenanceRecordModified(Child, RecordEvent):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"))


class Time(Parent, Base):
    PROTO = reaction_pb2.Time

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Time.TimeUnit.keys(), name="Time.TimeUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Time) -> Time:
        return Time(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Time:
        return reaction_pb2.Time(value=self.value, precision=self.precision, units=self.units)


class ReactionInputAdditionTime(Child, Time):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), unique=True)


class ReactionInputAdditionDuration(Child, Time):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), unique=True)


class TemperatureConditionsMeasurementTime(Child, Time):
    temperature_conditions_measurement_id = Column(
        Integer, ForeignKey("temperature_conditions_measurement.id"), unique=True
    )


class PressureConditionsMeasurementTime(Child, Time):
    pressure_conditions_measurement_id = Column(Integer, ForeignKey("pressure_conditions_measurement.id"), unique=True)


class ElectrochemistryConditionsMeasurementTime(Child, Time):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), unique=True
    )


class ReactionObservationTime(Child, Time):
    reaction_observation_id = Column(Integer, ForeignKey("reaction_observation.id"), unique=True)


class ReactionWorkupDuration(Child, Time):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), unique=True)


class ReactionOutcomeReactionTime(Child, Time):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), unique=True)


class ProductMeasurementRetentionTime(Child, Time):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)


class Mass(Base):
    PROTO = reaction_pb2.Mass
    amount_id = Column(Integer, ForeignKey("amount.id"), unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Mass.MassUnit.keys(), name="Mass.MassUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Mass) -> Mass:
        return Mass(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Mass:
        return reaction_pb2.Mass(value=self.value, precision=self.precision, units=self.units)


class Moles(Base):
    PROTO = reaction_pb2.Moles
    amount_id = Column(Integer, ForeignKey("amount.id"), unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Moles.MolesUnit.keys(), name="Moles.MolesUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Moles) -> Moles:
        return Moles(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Moles:
        return reaction_pb2.Moles(value=self.value, precision=self.precision, units=self.units)


class Volume(Parent, Base):
    PROTO = reaction_pb2.Volume

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Volume.VolumeUnit.keys(), name="Volume.VolumeUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Volume) -> Volume:
        return Volume(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Volume:
        return reaction_pb2.Volume(value=self.value, precision=self.precision, units=self.units)


class AmountVolume(Child, Volume):
    amount_id = Column(Integer, ForeignKey("amount.id"), unique=True)


class VesselVolume(Child, Volume):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), unique=True)


class Concentration(Base):
    PROTO = reaction_pb2.Concentration

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Concentration.ConcentrationUnit.keys(), name="Concentration.ConcentrationUnit"))


class Pressure(Parent, Base):
    PROTO = reaction_pb2.Pressure

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Pressure.PressureUnit.keys(), name="Pressure.PressureUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Pressure) -> Pressure:
        return Pressure(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Pressure:
        return reaction_pb2.Pressure(value=self.value, precision=self.precision, units=self.units)


class PressureConditionsSetpoint(Child, Pressure):
    pressure_conditions_id = Column(Integer, ForeignKey("pressure_conditions.id"), unique=True)


class PressureConditionsMeasurementPressure(Child, Pressure):
    pressure_conditions_measurement_id = Column(Integer, ForeignKey("pressure_conditions_measurement.id"), unique=True)


class Temperature(Parent, Base):
    PROTO = reaction_pb2.Temperature

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Temperature.TemperatureUnit.keys(), name="Temperature.TemperatureUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Temperature) -> Temperature:
        return Temperature(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Temperature:
        return reaction_pb2.Temperature(value=self.value, precision=self.precision, units=self.units)


class ReactionInputAdditionTemperature(Child, Temperature):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), unique=True)


class TemperatureConditionsSetpoint(Child, Temperature):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), unique=True)


class TemperatureConditionsMeasurementTemperature(Child, Temperature):
    temperature_conditions_measurement_id = Column(
        Integer, ForeignKey("temperature_conditions_measurement.id"), unique=True
    )


class Current(Parent, Base):
    PROTO = reaction_pb2.Current

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Current.CurrentUnit.keys(), name="Current.CurrentUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Current) -> Current:
        return Current(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Current:
        return reaction_pb2.Current(value=self.value, precision=self.precision, units=self.units)


class ElectrochemistryConditionsCurrent(Child, Current):
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), unique=True)


class ElectrochemistryConditionsMeasurementCurrent(Child, Current):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), unique=True
    )


class Voltage(Parent, Base):
    PROTO = reaction_pb2.Voltage

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Voltage.VoltageUnit.keys(), name="Voltage.VoltageUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Voltage) -> Voltage:
        return Voltage(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Voltage:
        return reaction_pb2.Voltage(value=self.value, precision=self.precision, units=self.units)


class ElectrochemistryConditionsVoltage(Child, Voltage):
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), unique=True)


class ElectrochemistryConditionsMeasurementVoltage(Child, Voltage):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), unique=True
    )


class Length(Parent, Base):
    PROTO = reaction_pb2.Length

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Length.LengthUnit.keys(), name="Length.LengthUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Length) -> Length:
        return Length(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Length:
        return reaction_pb2.Length(value=self.value, precision=self.precision, units=self.units)


class IlluminationConditionsDistanceToVessel(Child, Length):
    illumination_conditions_id = Column(Integer, ForeignKey("illumination_conditions.id"), unique=True)


class ElectrochemistryConditionsElectrodeSeparation(Child, Length):
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), unique=True)


class TubingDiameter(Child, Length):
    tubing_id = Column(Integer, ForeignKey("tubing.id"), unique=True)


class Wavelength(Parent, Base):
    PROTO = reaction_pb2.Wavelength

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Wavelength.WavelengthUnit.keys(), name="Wavelength.WavelengthUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Wavelength) -> Wavelength:
        return Wavelength(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.Wavelength:
        return reaction_pb2.Wavelength(value=self.value, precision=self.precision, units=self.units)


class IlluminationConditionsPeakWavelength(Child, Wavelength):
    illumination_conditions_id = Column(Integer, ForeignKey("illumination_conditions.id"), unique=True)


class ProductMeasurementWavelength(Child, Wavelength):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)


class FlowRate(Base):
    PROTO = reaction_pb2.FlowRate
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.FlowRate.FlowRateUnit.keys(), name="FlowRate.FlowRateUnit"))

    @classmethod
    def from_proto(cls, message: reaction_pb2.FlowRate) -> FlowRate:
        return FlowRate(value=message.value, precision=message.precision, units=enum_to_str(message, "units"))

    def to_proto(self) -> reaction_pb2.FlowRate:
        return reaction_pb2.FlowRate(value=self.value, precision=self.precision, units=self.units)


class Percentage(Parent, Base):
    PROTO = reaction_pb2.Percentage

    value = Column(Float)
    precision = Column(Float)

    @classmethod
    def from_proto(cls, message: reaction_pb2.Percentage) -> Percentage:
        return Percentage(value=message.value, precision=message.precision)

    def to_proto(self) -> reaction_pb2.Percentage:
        return reaction_pb2.Percentage(value=self.value, precision=self.precision)


class ReactionOutcomeConversion(Child, Percentage):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), unique=True)


class ProductMeasurementPercentage(Child, Percentage):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)


class FloatValue(Parent, Base):
    PROTO = reaction_pb2.FloatValue

    value = Column(Float)
    precision = Column(Float)

    @classmethod
    def from_proto(cls, message: reaction_pb2.FloatValue) -> FloatValue:
        return FloatValue(value=message.value, precision=message.precision)

    def to_proto(self) -> reaction_pb2.FloatValue:
        return reaction_pb2.FloatValue(value=self.value, precision=self.precision)


class ProductMeasurementFloatValue(Child, FloatValue):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), unique=True)


class MassSpecMeasurementDetailsEicMasses(Child, FloatValue):
    mass_spec_measurement_details_id = Column(Integer, ForeignKey("mass_spec_measurement_details.id"))


class Data(Parent, Base):
    PROTO = reaction_pb2.Data

    float_value = Column(Float)
    integer_value = Column(Integer)
    bytes_value = Column(LargeBinary)
    string_value = Column(Text)
    url = Column(Text)
    description = Column(Text)
    format = Column(String(255))

    @classmethod
    def from_proto(cls, message: reaction_pb2.Data) -> Data:
        kind = message.WhichOneof("kind")
        kwargs = {kind: getattr(message, kind), "description": message.description, "format": message.format}
        return Data(**kwargs)

    def to_proto(self) -> reaction_pb2.Data:
        return reaction_pb2.Data(
            float_value=self.float_value,
            integer_value=self.integer_value,
            bytes_value=self.bytes_value,
            string_value=self.string_value,
            url=self.url,
            description=self.description,
            format=self.format,
        )


class CompoundFeatures(Child, Data):
    compound_id = Column(Integer, ForeignKey("compound.id"))
    name = Column(String(255))  # Map key.


class ReactionSetupAutomationCode(Child, Data):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"))
    name = Column(String(255))  # Map key.


class ReactionObservationImage(Child, Data):
    reaction_observation_id = Column(Integer, ForeignKey("reaction_observation.id"), unique=True)


class ProductCompoundFeatures(Child, Data):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"))
    name = Column(String(255))  # Map key.


class AnalysisData(Child, Data):
    analysis_id = Column(Integer, ForeignKey("analysis.id"))
    name = Column(String(255))  # Map key.


if __name__ == "__main__":
    from sqlalchemy import create_engine

    engine = create_engine("postgresql://postgres:postgres@localhost:5433/test", echo=True, future=True)
    Base.metadata.create_all(engine)
