# Copyright 2022 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Table mappings for Reaction protos.

Notes:
    * Foreign keys to the `reaction` table are done using the `id` column, not the ORD reaction ID (`reaction_id`).
      However, the `reaction_id` column is used when the reaction ID is specifically called for, as with crude inputs.
    * When a message type is used in multiple places, use Joined Table Inheritance; see
      https://docs.sqlalchemy.org/en/14/orm/inheritance.html#joined-table-inheritance.
    * The naming might be a bit confusing: classes for single-use and parent protos have the same name
      as the corresponding message, while classes for multi-use protos (child classes) are named as
      <ContainingClass><Attribute>. For example, Reaction.identifiers uses the ReactionIdentifier class,
      while Reaction.inputs uses the ReactionInputs class (since ReactionInput is used in multiple places).
      The effect is that all tables storing proto information are named after their corresponding proto
      messages, with extra child tables for storing relationships.
    * Only message types are allowed for repeated/mapped values in the ORM (not scalar types). Specifically:
        * MassSpecMeasurementDetails.eic_masses is converted from repeated float to repeated FloatValue.
        * Dataset.reaction_ids is converted from repeated string to repeated DatasetId.
      # TODO(skearnes): Update the proto definition?
    * Source.id was renamed in the ORM to Source.vendor_id to avoid conflicting with the `id` column.
      # TODO(skearnes): Update the proto definition?
"""
from __future__ import annotations

from inspect import getmro
from typing import Optional, Type

from google.protobuf.descriptor import FieldDescriptor
from google.protobuf.message import Message
from google.protobuf.pyext._message import MessageMapContainer
from inflection import underscore
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

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


class Base:
    """See https://docs.sqlalchemy.org/en/14/orm/declarative_mixins.html#augmenting-the-base.

    * The table name is a snake_case rendering of the class name.
    * Every table has an `id` column used as the primary key. Note that Child tables do not have this column.
    """

    id = Column(Integer, primary_key=True)

    @declared_attr
    def __tablename__(cls):  # pylint: disable=no-self-argument
        return underscore(cls.__name__)  # pylint: disable=no-member


Base = declarative_base(cls=Base)


@declarative_mixin
class Parent:
    """Mixin class for parent classes.

    * Polymorphic types are identified in the `_type` column.
    """

    _type = Column(String(255), nullable=False)

    @declared_attr
    def __mapper_args__(cls):  # pylint: disable=no-self-argument
        return {
            "polymorphic_identity": underscore(cls.__name__),  # pylint: disable=no-member
            "polymorphic_on": cls._type,
        }


@declarative_mixin
class Child:
    """Mixin class for child classes.

    * The parent ID is stored in the `parent_id` column (also used as the primary key).
    """

    @declared_attr
    def parent_id(cls):  # pylint: disable=no-self-argument
        """Creates a `parent_id` column referring to the parent table."""
        parent_table = None
        for parent_class in getmro(cls)[1:]:
            if issubclass(parent_class, Base):
                parent_table = parent_class.__tablename__
                break
        assert parent_table is not None, (cls, getmro(cls))
        return Column(f"{parent_table}_id", Integer, ForeignKey(f"{parent_table}.id"), primary_key=True)

    @declared_attr
    def __mapper_args__(cls):  # pylint: disable=no-self-argument
        return {"polymorphic_identity": underscore(cls.__name__)}  # pylint: disable=no-member


# pylint:disable=missing-class-docstring


class Dataset(Base):
    name = Column(Text)
    description = Column(Text)
    reactions = relationship("Reaction")
    reaction_ids = relationship("ReactionId")
    dataset_id = Column(String(255), nullable=False, unique=True)


class ReactionId(Base):
    dataset_id = Column(Integer, ForeignKey("dataset.id"), nullable=False)

    reaction_id = Column(String(255), ForeignKey("reaction.reaction_id"), nullable=False)


class DatasetExample(Base):
    dataset_id = Column(String(255), ForeignKey("dataset.dataset_id"), nullable=False)
    description = Column(Text)
    url = Column(Text)
    created = relationship("DatasetExampleCreated", uselist=False)


class Reaction(Base):
    dataset_id = Column(String(255), ForeignKey("dataset.dataset_id"), nullable=False)
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
    reaction_id = Column(String(255), nullable=False, unique=True)


class ReactionIdentifier(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

    type = Column(
        Enum(*reaction_pb2.ReactionIdentifier.IdentifierType.keys(), name="ReactionIdentifier.IdentifierType")
    )
    details = Column(Text)
    value = Column(Text)
    is_mapped = Column(Boolean)


class ReactionInput(Parent, Base):
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
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)
    key = Column(String(255))  # Map key.


class ReactionWorkupInput(Child, ReactionInput):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False, unique=True)


class AdditionSpeed(Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ReactionInput.AdditionSpeed.AdditionSpeedType.keys(),
            name="ReactionInput.AdditionSpeed.AdditionSpeedType",
        )
    )
    details = Column(Text)


class AdditionDevice(Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ReactionInput.AdditionDevice.AdditionDeviceType.keys(),
            name="ReactionInput.AdditionDevice.AdditionDeviceType",
        )
    )
    details = Column(Text)


class Amount(Parent, Base):
    mass = relationship("Mass", uselist=False)
    moles = relationship("Moles", uselist=False)
    volume = relationship("AmountVolume", uselist=False)
    unmeasured = relationship("UnmeasuredAmount", uselist=False)
    volume_includes_solutes = Column(Boolean)


class CrudeComponentAmount(Child, Amount):
    crude_component_id = Column(Integer, ForeignKey("crude_component.id"), nullable=False)


class CompoundAmount(Child, Amount):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)


class ReactionWorkupAmount(Child, Amount):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False)


class ProductMeasurementAmount(Child, Amount):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False)


class UnmeasuredAmount(Base):
    amount_id = Column(Integer, ForeignKey("amount.id"), nullable=False, unique=True)

    type = Column(
        Enum(*reaction_pb2.UnmeasuredAmount.UnmeasuredAmountType.keys(), name="UnmeasuredAmount.UnmeasuredAmountType")
    )
    details = Column(Text)


class CrudeComponent(Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False)

    reaction_id = Column(String(255), ForeignKey("reaction.reaction_id"), nullable=False)
    includes_workup = Column(Boolean)
    has_derived_amount = Column(Boolean)
    amount = relationship("CrudeComponentAmount", uselist=False)


# Shared enum; see https://docs.sqlalchemy.org/en/14/dialects/postgresql.html#sqlalchemy.dialects.postgresql.ENUM.
ReactionRoleType = Enum(
    *reaction_pb2.ReactionRole.ReactionRoleType.keys(), name="ReactionRole.ReactionRoleType", metadata=Base.metadata
)


class Compound(Parent, Base):
    identifiers = relationship("CompoundIdentifiers")
    amount = relationship("CompoundAmount", uselist=False)
    reaction_role = Column(ReactionRoleType)
    is_limiting = Column(Boolean)
    preparations = relationship("CompoundPreparation")
    source = relationship("Source", uselist=False)
    features = relationship("CompoundFeatures")
    analyses = relationship("CompoundAnalyses")


class Source(Base):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False, unique=True)

    vendor = Column(String(255))
    vendor_id = Column(String(255))
    lot = Column(String(255))


class ReactionInputComponents(Child, Compound):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False)


class ProductMeasurementAuthenticStandard(Child, Compound):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class CompoundPreparation(Base):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)

    type = Column(
        Enum(*reaction_pb2.CompoundPreparation.PreparationType.keys(), name="CompoundPreparation.PreparationType")
    )
    details = Column(Text)
    reaction_id = Column(String(255), ForeignKey("reaction.reaction_id"))


class CompoundIdentifier(Parent, Base):
    type = Column(
        Enum(*reaction_pb2.CompoundIdentifier.IdentifierType.keys(), name="CompoundIdentifier.IdentifierType")
    )
    details = Column(Text)
    value = Column(Text)


class CompoundIdentifiers(Child, CompoundIdentifier):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)


class ProductCompoundIdentifiers(Child, CompoundIdentifier):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False)


class Vessel(Base):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), nullable=False, unique=True)

    type = Column(Enum(*reaction_pb2.Vessel.VesselType.keys(), name="Vessel.VesselType"))
    details = Column(Text)
    material = relationship("VesselMaterial", uselist=False)
    preparations = relationship("VesselPreparation")
    attachments = relationship("VesselAttachment")
    volume = relationship("VesselVolume", uselist=False)
    plate_id = Column(String(255))
    plate_position = Column(String(32))


class VesselMaterial(Base):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), nullable=False, unique=True)

    type = Column(
        Enum(*reaction_pb2.VesselMaterial.VesselMaterialType.keys(), name="VesselMaterial.VesselMaterialType")
    )
    details = Column(Text)


class VesselAttachment(Base):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), nullable=False)

    type = Column(
        Enum(*reaction_pb2.VesselAttachment.VesselAttachmentType.keys(), name="VesselAttachment.VesselAttachmentType")
    )
    details = Column(Text)


class VesselPreparation(Base):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), nullable=False)

    type = Column(
        Enum(
            *reaction_pb2.VesselPreparation.VesselPreparationType.keys(), name="VesselPreparation.VesselPreparationType"
        )
    )
    details = Column(Text)


class ReactionSetup(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False, unique=True)

    vessel = relationship("Vessel", uselist=False)
    is_automated = Column(Boolean)
    automation_platform = Column(String(255))
    automation_code = relationship("ReactionSetupAutomationCode")
    environment = relationship("ReactionEnvironment", uselist=False)


class ReactionEnvironment(Base):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType.keys(),
            name="ReactionSetup.ReactionEnvironment.ReactionEnvironmentType",
        )
    )
    details = Column(Text)


class ReactionConditions(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False, unique=True)

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
    control = relationship("TemperatureControl", uselist=False)
    setpoint = relationship("TemperatureConditionsSetpoint", uselist=False)
    measurements = relationship("TemperatureConditionsMeasurement")


class ReactionConditionsTemperature(Child, TemperatureConditions):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)


class ReactionWorkupTemperature(Child, TemperatureConditions):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False, unique=True)


class TemperatureControl(Base):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.TemperatureConditions.TemperatureControl.TemperatureControlType.keys(),
            name="TemperatureConditions.TemperatureControl.TemperatureControlType",
        )
    )
    details = Column(Text)


class TemperatureConditionsMeasurement(Base):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), nullable=False)

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
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)

    control = relationship("PressureControl", uselist=False)
    setpoint = relationship("PressureConditionsSetpoint", uselist=False)
    atmosphere = relationship("Atmosphere", uselist=False)
    measurements = relationship("PressureConditionsMeasurement")


class PressureControl(Base):
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.PressureConditions.PressureControl.PressureControlType.keys(),
            name="PressureConditions.PressureControl.PressureControlType",
        )
    )
    details = Column(Text)


class Atmosphere(Base):
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.PressureConditions.Atmosphere.AtmosphereType.keys(),
            name="PressureConditions.Atmosphere.AtmosphereType",
        )
    )
    details = Column(Text)


class PressureConditionsMeasurement(Base):
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), nullable=False)

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
    type = Column(
        Enum(*reaction_pb2.StirringConditions.StirringMethodType.keys(), name="StirringConditions.StirringMethodType")
    )
    details = Column(Text)
    rate = relationship("StirringRate", uselist=False)


class ReactionConditionsStirring(Child, StirringConditions):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)


class ReactionWorkupStirring(Child, StirringConditions):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False, unique=True)


class StirringRate(Base):
    stirring_conditions_id = Column(Integer, ForeignKey("stirring_conditions.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.StirringConditions.StirringRate.StirringRateType.keys(),
            name="StirringConditions.StirringRate.StirringRateType",
        )
    )
    details = Column(Text)
    rpm = Column(Integer)


class IlluminationConditions(Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)

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
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)

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
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False)

    time = relationship("ElectrochemistryConditionsMeasurementTime", uselist=False)
    current = relationship("ElectrochemistryConditionsMeasurementCurrent", uselist=False)
    voltage = relationship("ElectrochemistryConditionsMeasurementVoltage", uselist=False)


class ElectrochemistryCell(Base):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )

    type = Column(
        Enum(
            *reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType.keys(),
            name="ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType",
        )
    )
    details = Column(Text)


class FlowConditions(Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)

    type = Column(Enum(*reaction_pb2.FlowConditions.FlowType.keys(), name="FlowConditions.FlowType"))
    details = Column(Text)
    pump_type = Column(Text)
    tubing = relationship("Tubing", uselist=False)


class Tubing(Base):
    flow_conditions_id = Column(Integer, ForeignKey("flow_conditions.id"), nullable=False, unique=True)

    type = Column(Enum(*reaction_pb2.FlowConditions.Tubing.TubingType.keys(), name="FlowConditions.Tubing.TubingType"))
    details = Column(Text)
    diameter = relationship("TubingDiameter", uselist=False)


class ReactionNotes(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False, unique=True)

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
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

    time = relationship("ReactionObservationTime", uselist=False)
    comment = Column(Text)
    image = relationship("ReactionObservationImage", uselist=False)


class ReactionWorkup(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

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
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

    reaction_time = relationship("ReactionOutcomeReactionTime", uselist=False)
    conversion = relationship("ReactionOutcomeConversion", uselist=False)
    products = relationship("ProductCompound")
    analyses = relationship("ReactionOutcomeAnalyses")


class ProductCompound(Base):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False)

    identifiers = relationship("ProductCompoundIdentifiers")
    is_desired_product = Column(Boolean)
    measurements = relationship("ProductMeasurement")
    isolated_color = Column(String(255))
    texture = relationship("Texture", uselist=False)
    features = relationship("ProductCompoundFeatures")
    reaction_role = Column(ReactionRoleType)


class Texture(Base):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False, unique=True)

    type = Column(
        Enum(*reaction_pb2.ProductCompound.Texture.TextureType.keys(), name="ProductCompound.Texture.TextureType")
    )
    details = Column(Text)


class ProductMeasurement(Base):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False)

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
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)

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
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ProductMeasurement.Selectivity.SelectivityType.keys(),
            name="ProductMeasurement.Selectivity.SelectivityType",
        )
    )
    details = Column(Text)


class DateTime(Parent, Base):
    value = Column(String(255))


class AnalysisInstrumentLastCalibrated(Child, DateTime):
    analysis_id = Column(Integer, ForeignKey("analysis.id"), nullable=False, unique=True)


class ReactionProvenanceExperimentStart(Child, DateTime):
    analysis_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False, unique=True)


class RecordEventTime(Child, DateTime):
    analysis_id = Column(Integer, ForeignKey("record_event.id"), nullable=False, unique=True)


class Analysis(Parent, Base):
    type = Column(Enum(*reaction_pb2.Analysis.AnalysisType.keys(), name="Analysis.AnalysisType"))
    details = Column(Text)
    chmo_id = Column(Integer)
    is_of_isolated_species = Column(Boolean)
    data = relationship("AnalysisData")
    instrument_manufacturer = Column(String(255))
    instrument_last_calibrated = relationship("AnalysisInstrumentLastCalibrated", uselist=False)


class CompoundAnalyses(Child, Analysis):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)
    key = Column(String(255))  # Map key.


class ReactionOutcomeAnalyses(Child, Analysis):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False)
    key = Column(String(255))  # Map key.


class ReactionProvenance(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False, unique=True)

    experimenter = relationship("ReactionProvenanceExperimenter", uselist=False)
    city = Column(String(255))
    experiment_start = relationship("ReactionProvenanceExperimentStart", uselist=False)
    doi = Column(String(255))
    patent = Column(String(255))
    publication_url = Column(Text)
    record_created = relationship("ReactionProvenanceRecordCreated", uselist=False)
    record_modified = relationship("ReactionProvenanceRecordModified")


class Person(Parent, Base):
    username = Column(String(255))
    name = Column(String(255))
    orcid = Column(String(19))
    organization = Column(String(255))
    email = Column(String(255))


class ReactionProvenanceExperimenter(Child, Person):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False, unique=True)


class RecordEventPerson(Child, Person):
    record_event_id = Column(Integer, ForeignKey("record_event.id"), nullable=False, unique=True)


class RecordEvent(Parent, Base):
    time = relationship("RecordEventTime", uselist=False)
    person = relationship("RecordEventPerson", uselist=False)
    details = Column(Text)


class ReactionProvenanceRecordCreated(Child, RecordEvent):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False, unique=True)


class ReactionProvenanceRecordModified(Child, RecordEvent):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False)


class DatasetExampleCreated(Child, RecordEvent):
    dataset_example_id = Column(Integer, ForeignKey("dataset_example.id"), nullable=False, unique=True)


class Time(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Time.TimeUnit.keys(), name="Time.TimeUnit"))


class ReactionInputAdditionTime(Child, Time):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)


class ReactionInputAdditionDuration(Child, Time):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)


class TemperatureConditionsMeasurementTime(Child, Time):
    temperature_conditions_measurement_id = Column(
        Integer, ForeignKey("temperature_conditions_measurement.id"), nullable=False, unique=True
    )


class PressureConditionsMeasurementTime(Child, Time):
    pressure_conditions_measurement_id = Column(
        Integer, ForeignKey("pressure_conditions_measurement.id"), nullable=False, unique=True
    )


class ElectrochemistryConditionsMeasurementTime(Child, Time):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), nullable=False, unique=True
    )


class ReactionObservationTime(Child, Time):
    reaction_observation_id = Column(Integer, ForeignKey("reaction_observation.id"), nullable=False, unique=True)


class ReactionWorkupDuration(Child, Time):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False, unique=True)


class ReactionOutcomeReactionTime(Child, Time):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False, unique=True)


class ProductMeasurementRetentionTime(Child, Time):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class Mass(Base):
    amount_id = Column(Integer, ForeignKey("amount.id"), nullable=False, unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Mass.MassUnit.keys(), name="Mass.MassUnit"))


class Moles(Base):
    amount_id = Column(Integer, ForeignKey("amount.id"), nullable=False, unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Moles.MolesUnit.keys(), name="Moles.MolesUnit"))


class Volume(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Volume.VolumeUnit.keys(), name="Volume.VolumeUnit"))


class AmountVolume(Child, Volume):
    amount_id = Column(Integer, ForeignKey("amount.id"), nullable=False, unique=True)


class VesselVolume(Child, Volume):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), nullable=False, unique=True)


class Concentration(Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Concentration.ConcentrationUnit.keys(), name="Concentration.ConcentrationUnit"))


class Pressure(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Pressure.PressureUnit.keys(), name="Pressure.PressureUnit"))


class PressureConditionsSetpoint(Child, Pressure):
    pressure_conditions_id = Column(Integer, ForeignKey("pressure_conditions.id"), nullable=False, unique=True)


class PressureConditionsMeasurementPressure(Child, Pressure):
    pressure_conditions_measurement_id = Column(
        Integer, ForeignKey("pressure_conditions_measurement.id"), nullable=False, unique=True
    )


class Temperature(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Temperature.TemperatureUnit.keys(), name="Temperature.TemperatureUnit"))


class ReactionInputAdditionTemperature(Child, Temperature):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)


class TemperatureConditionsSetpoint(Child, Temperature):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), nullable=False, unique=True)


class TemperatureConditionsMeasurementTemperature(Child, Temperature):
    temperature_conditions_measurement_id = Column(
        Integer, ForeignKey("temperature_conditions_measurement.id"), nullable=False, unique=True
    )


class Current(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Current.CurrentUnit.keys(), name="Current.CurrentUnit"))


class ElectrochemistryConditionsCurrent(Child, Current):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )


class ElectrochemistryConditionsMeasurementCurrent(Child, Current):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), nullable=False, unique=True
    )


class Voltage(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Voltage.VoltageUnit.keys(), name="Voltage.VoltageUnit"))


class ElectrochemistryConditionsVoltage(Child, Voltage):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )


class ElectrochemistryConditionsMeasurementVoltage(Child, Voltage):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), nullable=False, unique=True
    )


class Length(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Length.LengthUnit.keys(), name="Length.LengthUnit"))


class IlluminationConditionsDistanceToVessel(Child, Length):
    illumination_conditions_id = Column(Integer, ForeignKey("illumination_conditions.id"), nullable=False, unique=True)


class ElectrochemistryConditionsElectrodeSeparation(Child, Length):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )


class TubingDiameter(Child, Length):
    tubing_id = Column(Integer, ForeignKey("tubing.id"), nullable=False, unique=True)


class Wavelength(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Wavelength.WavelengthUnit.keys(), name="Wavelength.WavelengthUnit"))


class IlluminationConditionsPeakWavelength(Child, Wavelength):
    illumination_conditions_id = Column(Integer, ForeignKey("illumination_conditions.id"), nullable=False, unique=True)


class ProductMeasurementWavelength(Child, Wavelength):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class FlowRate(Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.FlowRate.FlowRateUnit.keys(), name="FlowRate.FlowRateUnit"))


class Percentage(Parent, Base):
    value = Column(Float)
    precision = Column(Float)


class ReactionOutcomeConversion(Child, Percentage):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False, unique=True)


class ProductMeasurementPercentage(Child, Percentage):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class FloatValue(Parent, Base):
    value = Column(Float)
    precision = Column(Float)


class ProductMeasurementFloatValue(Child, FloatValue):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class MassSpecMeasurementDetailsEicMasses(Child, FloatValue):
    mass_spec_measurement_details_id = Column(Integer, ForeignKey("mass_spec_measurement_details.id"), nullable=False)


class Data(Parent, Base):
    float_value = Column(Float)
    integer_value = Column(Integer)
    bytes_value = Column(LargeBinary)
    string_value = Column(Text)
    url = Column(Text)
    description = Column(Text)
    format = Column(String(255))


class CompoundFeatures(Child, Data):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)
    key = Column(String(255))  # Map key.


class ReactionSetupAutomationCode(Child, Data):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), nullable=False)
    key = Column(String(255))  # Map key.


class ReactionObservationImage(Child, Data):
    reaction_observation_id = Column(Integer, ForeignKey("reaction_observation.id"), nullable=False, unique=True)


class ProductCompoundFeatures(Child, Data):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False)
    key = Column(String(255))  # Map key.


class AnalysisData(Child, Data):
    analysis_id = Column(Integer, ForeignKey("analysis.id"), nullable=False)
    key = Column(String(255))  # Map key.


MAPPERS: dict[Type[Message], Type[Base]] = {
    dataset_pb2.Dataset: Dataset,
    dataset_pb2.DatasetExample: DatasetExample,
    reaction_pb2.Reaction: Reaction,
    reaction_pb2.ReactionIdentifier: ReactionIdentifier,
    reaction_pb2.ReactionInput: ReactionInput,
    reaction_pb2.ReactionInput.AdditionSpeed: AdditionSpeed,
    reaction_pb2.ReactionInput.AdditionDevice: AdditionDevice,
    reaction_pb2.Amount: Amount,
    reaction_pb2.UnmeasuredAmount: UnmeasuredAmount,
    reaction_pb2.CrudeComponent: CrudeComponent,
    reaction_pb2.Compound: Compound,
    reaction_pb2.Compound.Source: Source,
    reaction_pb2.CompoundPreparation: CompoundPreparation,
    reaction_pb2.CompoundIdentifier: CompoundIdentifier,
    reaction_pb2.Vessel: Vessel,
    reaction_pb2.VesselMaterial: VesselMaterial,
    reaction_pb2.VesselAttachment: VesselAttachment,
    reaction_pb2.VesselPreparation: VesselPreparation,
    reaction_pb2.ReactionSetup: ReactionSetup,
    reaction_pb2.ReactionSetup.ReactionEnvironment: ReactionEnvironment,
    reaction_pb2.ReactionConditions: ReactionConditions,
    reaction_pb2.TemperatureConditions: TemperatureConditions,
    reaction_pb2.TemperatureConditions.TemperatureControl: TemperatureControl,
    reaction_pb2.TemperatureConditions.Measurement: TemperatureConditionsMeasurement,
    reaction_pb2.PressureConditions: PressureConditions,
    reaction_pb2.PressureConditions.PressureControl: PressureControl,
    reaction_pb2.PressureConditions.Atmosphere: Atmosphere,
    reaction_pb2.PressureConditions.Measurement: PressureConditionsMeasurement,
    reaction_pb2.StirringConditions: StirringConditions,
    reaction_pb2.StirringConditions.StirringRate: StirringRate,
    reaction_pb2.IlluminationConditions: IlluminationConditions,
    reaction_pb2.ElectrochemistryConditions: ElectrochemistryConditions,
    reaction_pb2.ElectrochemistryConditions.Measurement: ElectrochemistryConditionsMeasurement,
    reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell: ElectrochemistryCell,
    reaction_pb2.FlowConditions: FlowConditions,
    reaction_pb2.FlowConditions.Tubing: Tubing,
    reaction_pb2.ReactionNotes: ReactionNotes,
    reaction_pb2.ReactionObservation: ReactionObservation,
    reaction_pb2.ReactionWorkup: ReactionWorkup,
    reaction_pb2.ReactionOutcome: ReactionOutcome,
    reaction_pb2.ProductCompound: ProductCompound,
    reaction_pb2.ProductCompound.Texture: Texture,
    reaction_pb2.ProductMeasurement: ProductMeasurement,
    reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails: MassSpecMeasurementDetails,
    reaction_pb2.ProductMeasurement.Selectivity: Selectivity,
    reaction_pb2.DateTime: DateTime,
    reaction_pb2.Analysis: Analysis,
    reaction_pb2.ReactionProvenance: ReactionProvenance,
    reaction_pb2.Person: Person,
    reaction_pb2.RecordEvent: RecordEvent,
    reaction_pb2.Time: Time,
    reaction_pb2.Mass: Mass,
    reaction_pb2.Moles: Moles,
    reaction_pb2.Volume: Volume,
    reaction_pb2.Concentration: Concentration,
    reaction_pb2.Pressure: Pressure,
    reaction_pb2.Temperature: Temperature,
    reaction_pb2.Current: Current,
    reaction_pb2.Voltage: Voltage,
    reaction_pb2.Length: Length,
    reaction_pb2.Wavelength: Wavelength,
    reaction_pb2.FlowRate: FlowRate,
    reaction_pb2.Percentage: Percentage,
    reaction_pb2.FloatValue: FloatValue,
    reaction_pb2.Data: Data,
}
PROTOS: dict[Type[Base], Type[Message]] = {value: key for key, value in MAPPERS.items()}

MAPPER_RENAMES: dict[tuple[Type[Message], str], str] = {
    (reaction_pb2.Compound.Source, "id"): "vendor_id",
}
PROTO_RENAMES: dict[tuple[Type[Base], str], str] = {
    (MAPPERS[key[0]], value): key[1] for key, value in MAPPER_RENAMES.items()
}


def from_proto(message: Message, mapper: Optional[Type[Base]] = None, key: Optional[str] = None) -> Base:
    """Converts a protobuf message into an ORM object.

    Args:
        message: Protobuf message.
        mapper: ORM mapper class. For top-level protos like Dataset and Reaction this can be left as None; it must
            be provided for Child subclasses to properly handle polymorphism.
        key: Map key (we store maps as rows of (key, value) tuples).

    Returns:
        ORM object.
    """
    if mapper is None:
        mapper = MAPPERS[type(message)]
    kwargs = {}
    if key is not None:
        kwargs["key"] = key
    if mapper == Reaction:
        kwargs["proto"] = message.SerializeToString()
    for field, value in message.ListFields():
        field_name = MAPPER_RENAMES.get((type(message), field.name), field.name)
        if field_name == "eic_masses":
            # Convert repeated float to repeated FloatValue.
            kwargs[field_name] = [MassSpecMeasurementDetailsEicMasses(value=v) for v in value]
        elif field_name == "reaction_ids":
            # Convert repeated string to repeated ReactionId.
            kwargs[field_name] = [ReactionId(reaction_id=v) for v in value]
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


def to_proto(base: Base) -> Message:
    """Converts an ORM object into a protobuf message.

    Args:
        base: ORM object.

    Returns:
        Protobuf message.
    """
    kwargs = {}
    proto = None
    for mapper in getmro(type(base)):
        if not issubclass(mapper, Child):
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
            if len(value) == 0:
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
