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
      _<ContainingClass><Attribute>. For example, Reaction.identifiers uses the ReactionIdentifier class,
      while Reaction.inputs uses the _ReactionInputs class (since _ReactionInput is used in multiple places).
      The effect is that all tables storing proto information are named after their corresponding proto
      messages, with extra child tables for storing relationships.
    * Parent and child class names are prepended with underscores; e.g. _ReactionInput. Every parent class has an
      associated with_polymorphic wrapper for query construction, named without the leading underscore;
      e.g. ReactionInput.
    * Only message types are allowed for repeated/mapped values in the ORM (not scalar types). Specifically:
        * MassSpecMeasurementDetails.eic_masses is converted from repeated float to repeated FloatValue.
        * Dataset.reaction_ids is converted from repeated string to repeated DatasetId.
      # TODO(skearnes): Update the proto definition?
    * Source.id was renamed in the ORM to Source.vendor_id to avoid conflicting with the `id` column.
      # TODO(skearnes): Update the proto definition?
"""
from __future__ import annotations

from inspect import getmro
from typing import Mapping, Optional, Type

from google.protobuf.descriptor import FieldDescriptor
from google.protobuf.message import Message
from sqlalchemy import Boolean, Column, Enum, Float, Integer, ForeignKey, LargeBinary, String, Text
from sqlalchemy.orm import relationship, with_polymorphic

import ord_schema.orm.structure  # pylint: disable=unused-import
from ord_schema import message_helpers
from ord_schema.orm import Base, Child, Parent
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

# pylint:disable=missing-class-docstring


class Dataset(Base):
    name = Column(Text)
    description = Column(Text)
    reactions = relationship("Reaction", cascade="all, delete-orphan")
    reaction_ids = relationship("ReactionId", cascade="all, delete-orphan")
    dataset_id = Column(String(255), nullable=False, unique=True)


class ReactionId(Base):
    dataset_id = Column(Integer, ForeignKey("dataset.id"), nullable=False)

    reaction_id = Column(String(255), ForeignKey("reaction.reaction_id"), nullable=False)


class DatasetExample(Base):
    dataset_id = Column(String(255), ForeignKey("dataset.dataset_id"), nullable=False)
    description = Column(Text)
    url = Column(Text)
    created = relationship("_DatasetExampleCreated", uselist=False, cascade="all, delete-orphan")


class Reaction(Base):
    dataset_id = Column(String(255), ForeignKey("dataset.dataset_id"), nullable=False)
    proto = Column(LargeBinary, nullable=False)

    identifiers = relationship("ReactionIdentifier", cascade="all, delete-orphan")
    inputs = relationship("_ReactionInputs", cascade="all, delete-orphan")
    setup = relationship("ReactionSetup", uselist=False, cascade="all, delete-orphan")
    conditions = relationship("ReactionConditions", uselist=False, cascade="all, delete-orphan")
    notes = relationship("ReactionNotes", uselist=False, cascade="all, delete-orphan")
    observations = relationship("ReactionObservation", cascade="all, delete-orphan")
    workups = relationship("ReactionWorkup", cascade="all, delete-orphan")
    outcomes = relationship("ReactionOutcome", cascade="all, delete-orphan")
    provenance = relationship("ReactionProvenance", uselist=False, cascade="all, delete-orphan")
    reaction_id = Column(String(255), nullable=False, unique=True)


class ReactionIdentifier(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

    type = Column(
        Enum(*reaction_pb2.ReactionIdentifier.IdentifierType.keys(), name="ReactionIdentifier.IdentifierType")
    )
    details = Column(Text)
    value = Column(Text)
    is_mapped = Column(Boolean)


class _ReactionInput(Parent, Base):
    components = relationship("_ReactionInputComponents", cascade="all, delete-orphan")
    crude_components = relationship("CrudeComponent", cascade="all, delete-orphan")
    addition_order = Column(Integer)
    addition_time = relationship("_ReactionInputAdditionTime", uselist=False, cascade="all, delete-orphan")
    addition_speed = relationship("AdditionSpeed", uselist=False, cascade="all, delete-orphan")
    addition_duration = relationship("_ReactionInputAdditionDuration", uselist=False, cascade="all, delete-orphan")
    flow_rate = relationship("FlowRate", uselist=False, cascade="all, delete-orphan")
    addition_device = relationship("AdditionDevice", uselist=False, cascade="all, delete-orphan")
    addition_temperature = relationship(
        "_ReactionInputAdditionTemperature", uselist=False, cascade="all, delete-orphan"
    )


class _ReactionInputs(Child, _ReactionInput):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)
    key = Column(String(255))  # Map key.


class _ReactionWorkupInput(Child, _ReactionInput):
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


class _Amount(Parent, Base):
    mass = relationship("Mass", uselist=False, cascade="all, delete-orphan")
    moles = relationship("Moles", uselist=False, cascade="all, delete-orphan")
    volume = relationship("_AmountVolume", uselist=False, cascade="all, delete-orphan")
    unmeasured = relationship("UnmeasuredAmount", uselist=False, cascade="all, delete-orphan")
    volume_includes_solutes = Column(Boolean)


class _CrudeComponentAmount(Child, _Amount):
    crude_component_id = Column(Integer, ForeignKey("crude_component.id"), nullable=False)


class _CompoundAmount(Child, _Amount):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)


class _ReactionWorkupAmount(Child, _Amount):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False)


class _ProductMeasurementAmount(Child, _Amount):
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
    amount = relationship("_CrudeComponentAmount", uselist=False, cascade="all, delete-orphan")


# Shared enum; see https://docs.sqlalchemy.org/en/14/dialects/postgresql.html#sqlalchemy.dialects.postgresql.ENUM.
ReactionRoleType = Enum(
    *reaction_pb2.ReactionRole.ReactionRoleType.keys(), name="ReactionRole.ReactionRoleType", metadata=Base.metadata
)


class _Compound(Parent, Base):
    identifiers = relationship("_CompoundIdentifiers", cascade="all, delete-orphan")
    amount = relationship("_CompoundAmount", uselist=False, cascade="all, delete-orphan")
    reaction_role = Column(ReactionRoleType)
    is_limiting = Column(Boolean)
    preparations = relationship("CompoundPreparation", cascade="all, delete-orphan")
    source = relationship("Source", uselist=False, cascade="all, delete-orphan")
    features = relationship("_CompoundFeatures", cascade="all, delete-orphan")
    analyses = relationship("_CompoundAnalyses", cascade="all, delete-orphan")
    structure = relationship("_CompoundStructure", uselist=False, cascade="all, delete-orphan")


class _ReactionInputComponents(Child, _Compound):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False)


class _ProductMeasurementAuthenticStandard(Child, _Compound):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class Source(Base):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False, unique=True)

    vendor = Column(String(255))
    vendor_id = Column(String(255))
    lot = Column(String(255))


class CompoundPreparation(Base):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)

    type = Column(
        Enum(*reaction_pb2.CompoundPreparation.PreparationType.keys(), name="CompoundPreparation.PreparationType")
    )
    details = Column(Text)
    reaction_id = Column(String(255), ForeignKey("reaction.reaction_id"))


class _CompoundIdentifier(Parent, Base):
    type = Column(
        Enum(*reaction_pb2.CompoundIdentifier.IdentifierType.keys(), name="CompoundIdentifier.IdentifierType")
    )
    details = Column(Text)
    value = Column(Text)


class _CompoundIdentifiers(Child, _CompoundIdentifier):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)


class _ProductCompoundIdentifiers(Child, _CompoundIdentifier):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False)


class Vessel(Base):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), nullable=False, unique=True)

    type = Column(Enum(*reaction_pb2.Vessel.VesselType.keys(), name="Vessel.VesselType"))
    details = Column(Text)
    material = relationship("VesselMaterial", uselist=False, cascade="all, delete-orphan")
    preparations = relationship("VesselPreparation", cascade="all, delete-orphan")
    attachments = relationship("VesselAttachment", cascade="all, delete-orphan")
    volume = relationship("_VesselVolume", uselist=False, cascade="all, delete-orphan")
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

    vessel = relationship("Vessel", uselist=False, cascade="all, delete-orphan")
    is_automated = Column(Boolean)
    automation_platform = Column(String(255))
    automation_code = relationship("_ReactionSetupAutomationCode", cascade="all, delete-orphan")
    environment = relationship("ReactionEnvironment", uselist=False, cascade="all, delete-orphan")


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

    temperature = relationship("_ReactionConditionsTemperature", uselist=False, cascade="all, delete-orphan")
    pressure = relationship("PressureConditions", uselist=False, cascade="all, delete-orphan")
    stirring = relationship("_ReactionConditionsStirring", uselist=False, cascade="all, delete-orphan")
    illumination = relationship("IlluminationConditions", uselist=False, cascade="all, delete-orphan")
    electrochemistry = relationship("ElectrochemistryConditions", uselist=False, cascade="all, delete-orphan")
    flow = relationship("FlowConditions", uselist=False, cascade="all, delete-orphan")
    reflux = Column(Boolean)
    ph = Column(Float)
    conditions_are_dynamic = Column(Boolean)
    details = Column(Text)


class _TemperatureConditions(Parent, Base):
    control = relationship("TemperatureControl", uselist=False, cascade="all, delete-orphan")
    setpoint = relationship("_TemperatureConditionsSetpoint", uselist=False, cascade="all, delete-orphan")
    measurements = relationship("TemperatureConditionsMeasurement", cascade="all, delete-orphan")


class _ReactionConditionsTemperature(Child, _TemperatureConditions):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)


class _ReactionWorkupTemperature(Child, _TemperatureConditions):
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
    time = relationship("_TemperatureConditionsMeasurementTime", uselist=False, cascade="all, delete-orphan")
    temperature = relationship(
        "_TemperatureConditionsMeasurementTemperature", uselist=False, cascade="all, delete-orphan"
    )


class PressureConditions(Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)

    control = relationship("PressureControl", uselist=False, cascade="all, delete-orphan")
    setpoint = relationship("_PressureConditionsSetpoint", uselist=False, cascade="all, delete-orphan")
    atmosphere = relationship("Atmosphere", uselist=False, cascade="all, delete-orphan")
    measurements = relationship("PressureConditionsMeasurement", cascade="all, delete-orphan")


class PressureControl(Base):
    pressure_conditions_id = Column(Integer, ForeignKey("pressure_conditions.id"), nullable=False, unique=True)

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
    time = relationship("_PressureConditionsMeasurementTime", uselist=False, cascade="all, delete-orphan")
    pressure = relationship("_PressureConditionsMeasurementPressure", uselist=False, cascade="all, delete-orphan")


class _StirringConditions(Parent, Base):
    type = Column(
        Enum(*reaction_pb2.StirringConditions.StirringMethodType.keys(), name="StirringConditions.StirringMethodType")
    )
    details = Column(Text)
    rate = relationship("StirringRate", uselist=False, cascade="all, delete-orphan")


class _ReactionConditionsStirring(Child, _StirringConditions):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)


class _ReactionWorkupStirring(Child, _StirringConditions):
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
    peak_wavelength = relationship("_IlluminationConditionsPeakWavelength", uselist=False, cascade="all, delete-orphan")
    color = Column(String(255))
    distance_to_vessel = relationship(
        "_IlluminationConditionsDistanceToVessel", uselist=False, cascade="all, delete-orphan"
    )


class ElectrochemistryConditions(Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ElectrochemistryConditions.ElectrochemistryType.keys(),
            name="ElectrochemistryConditions.ElectrochemistryType",
        )
    )
    details = Column(Text)
    current = relationship("_ElectrochemistryConditionsCurrent", uselist=False, cascade="all, delete-orphan")
    voltage = relationship("_ElectrochemistryConditionsVoltage", uselist=False, cascade="all, delete-orphan")
    anode_material = Column(String(255))
    cathode_material = Column(String(255))
    electrode_separation = relationship(
        "_ElectrochemistryConditionsElectrodeSeparation", uselist=False, cascade="all, delete-orphan"
    )
    measurements = relationship("ElectrochemistryConditionsMeasurement", cascade="all, delete-orphan")
    cell = relationship("ElectrochemistryCell", uselist=False, cascade="all, delete-orphan")


class ElectrochemistryConditionsMeasurement(Base):
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False)

    time = relationship("_ElectrochemistryConditionsMeasurementTime", uselist=False, cascade="all, delete-orphan")
    current = relationship("_ElectrochemistryConditionsMeasurementCurrent", uselist=False, cascade="all, delete-orphan")
    voltage = relationship("_ElectrochemistryConditionsMeasurementVoltage", uselist=False, cascade="all, delete-orphan")


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
    tubing = relationship("Tubing", uselist=False, cascade="all, delete-orphan")


class Tubing(Base):
    flow_conditions_id = Column(Integer, ForeignKey("flow_conditions.id"), nullable=False, unique=True)

    type = Column(Enum(*reaction_pb2.FlowConditions.Tubing.TubingType.keys(), name="FlowConditions.Tubing.TubingType"))
    details = Column(Text)
    diameter = relationship("_TubingDiameter", uselist=False, cascade="all, delete-orphan")


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

    time = relationship("_ReactionObservationTime", uselist=False, cascade="all, delete-orphan")
    comment = Column(Text)
    image = relationship("_ReactionObservationImage", uselist=False, cascade="all, delete-orphan")


class ReactionWorkup(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

    type = Column(Enum(*reaction_pb2.ReactionWorkup.WorkupType.keys(), name="ReactionWorkup.WorkupType"))
    details = Column(Text)
    duration = relationship("_ReactionWorkupDuration", uselist=False, cascade="all, delete-orphan")
    input = relationship("_ReactionWorkupInput", uselist=False, cascade="all, delete-orphan")
    amount = relationship("_ReactionWorkupAmount", uselist=False, cascade="all, delete-orphan")
    temperature = relationship("_ReactionWorkupTemperature", uselist=False, cascade="all, delete-orphan")
    keep_phase = Column(String(255))
    stirring = relationship("_ReactionWorkupStirring", uselist=False, cascade="all, delete-orphan")
    target_ph = Column(Float)
    is_automated = Column(Boolean)


class ReactionOutcome(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False)

    reaction_time = relationship("_ReactionOutcomeReactionTime", uselist=False, cascade="all, delete-orphan")
    conversion = relationship("_ReactionOutcomeConversion", uselist=False, cascade="all, delete-orphan")
    products = relationship("ProductCompound", cascade="all, delete-orphan")
    analyses = relationship("_ReactionOutcomeAnalyses", cascade="all, delete-orphan")


class ProductCompound(Base):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False)

    identifiers = relationship("_ProductCompoundIdentifiers", cascade="all, delete-orphan")
    is_desired_product = Column(Boolean)
    measurements = relationship("ProductMeasurement", cascade="all, delete-orphan")
    isolated_color = Column(String(255))
    texture = relationship("Texture", uselist=False, cascade="all, delete-orphan")
    features = relationship("_ProductCompoundFeatures", cascade="all, delete-orphan")
    reaction_role = Column(ReactionRoleType)
    structure = relationship("_ProductCompoundStructure", uselist=False, cascade="all, delete-orphan")


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
    authentic_standard = relationship(
        "_ProductMeasurementAuthenticStandard", uselist=False, cascade="all, delete-orphan"
    )
    percentage = relationship("_ProductMeasurementPercentage", uselist=False, cascade="all, delete-orphan")
    float_value = relationship("_ProductMeasurementFloatValue", uselist=False, cascade="all, delete-orphan")
    string_value = Column(Text)
    amount = relationship("_ProductMeasurementAmount", uselist=False, cascade="all, delete-orphan")
    retention_time = relationship("_ProductMeasurementRetentionTime", uselist=False, cascade="all, delete-orphan")
    mass_spec_details = relationship("MassSpecMeasurementDetails", uselist=False, cascade="all, delete-orphan")
    selectivity = relationship("Selectivity", uselist=False, cascade="all, delete-orphan")
    wavelength = relationship("_ProductMeasurementWavelength", uselist=False, cascade="all, delete-orphan")


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
    eic_masses = relationship("_MassSpecMeasurementDetailsEicMasses", cascade="all, delete-orphan")


class Selectivity(Base):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)

    type = Column(
        Enum(
            *reaction_pb2.ProductMeasurement.Selectivity.SelectivityType.keys(),
            name="ProductMeasurement.Selectivity.SelectivityType",
        )
    )
    details = Column(Text)


class _DateTime(Parent, Base):
    value = Column(String(255))


class _AnalysisInstrumentLastCalibrated(Child, _DateTime):
    analysis_id = Column(Integer, ForeignKey("analysis.id"), nullable=False, unique=True)


class _ReactionProvenanceExperimentStart(Child, _DateTime):
    analysis_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False, unique=True)


class _RecordEventTime(Child, _DateTime):
    analysis_id = Column(Integer, ForeignKey("record_event.id"), nullable=False, unique=True)


class _Analysis(Parent, Base):
    type = Column(Enum(*reaction_pb2.Analysis.AnalysisType.keys(), name="Analysis.AnalysisType"))
    details = Column(Text)
    chmo_id = Column(Integer)
    is_of_isolated_species = Column(Boolean)
    data = relationship("_AnalysisData", cascade="all, delete-orphan")
    instrument_manufacturer = Column(String(255))
    instrument_last_calibrated = relationship(
        "_AnalysisInstrumentLastCalibrated", uselist=False, cascade="all, delete-orphan"
    )


class _CompoundAnalyses(Child, _Analysis):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)
    key = Column(String(255))  # Map key.


class _ReactionOutcomeAnalyses(Child, _Analysis):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False)
    key = Column(String(255))  # Map key.


class ReactionProvenance(Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), nullable=False, unique=True)

    experimenter = relationship("_ReactionProvenanceExperimenter", uselist=False, cascade="all, delete-orphan")
    city = Column(String(255))
    experiment_start = relationship("_ReactionProvenanceExperimentStart", uselist=False, cascade="all, delete-orphan")
    doi = Column(String(255))
    patent = Column(String(255))
    publication_url = Column(Text)
    record_created = relationship("_ReactionProvenanceRecordCreated", uselist=False, cascade="all, delete-orphan")
    record_modified = relationship("_ReactionProvenanceRecordModified", cascade="all, delete-orphan")


class _Person(Parent, Base):
    username = Column(String(255))
    name = Column(String(255))
    orcid = Column(String(19))
    organization = Column(String(255))
    email = Column(String(255))


class _ReactionProvenanceExperimenter(Child, _Person):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False, unique=True)


class _RecordEventPerson(Child, _Person):
    record_event_id = Column(Integer, ForeignKey("record_event.id"), nullable=False, unique=True)


class _RecordEvent(Parent, Base):
    time = relationship("_RecordEventTime", uselist=False, cascade="all, delete-orphan")
    person = relationship("_RecordEventPerson", uselist=False, cascade="all, delete-orphan")
    details = Column(Text)


class _ReactionProvenanceRecordCreated(Child, _RecordEvent):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False, unique=True)


class _ReactionProvenanceRecordModified(Child, _RecordEvent):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), nullable=False)


class _DatasetExampleCreated(Child, _RecordEvent):
    dataset_example_id = Column(Integer, ForeignKey("dataset_example.id"), nullable=False, unique=True)


class _Time(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Time.TimeUnit.keys(), name="Time.TimeUnit"))


class _ReactionInputAdditionTime(Child, _Time):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)


class _ReactionInputAdditionDuration(Child, _Time):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)


class _TemperatureConditionsMeasurementTime(Child, _Time):
    temperature_conditions_measurement_id = Column(
        Integer, ForeignKey("temperature_conditions_measurement.id"), nullable=False, unique=True
    )


class _PressureConditionsMeasurementTime(Child, _Time):
    pressure_conditions_measurement_id = Column(
        Integer, ForeignKey("pressure_conditions_measurement.id"), nullable=False, unique=True
    )


class _ElectrochemistryConditionsMeasurementTime(Child, _Time):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), nullable=False, unique=True
    )


class _ReactionObservationTime(Child, _Time):
    reaction_observation_id = Column(Integer, ForeignKey("reaction_observation.id"), nullable=False, unique=True)


class _ReactionWorkupDuration(Child, _Time):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"), nullable=False, unique=True)


class _ReactionOutcomeReactionTime(Child, _Time):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False, unique=True)


class _ProductMeasurementRetentionTime(Child, _Time):
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


class _Volume(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Volume.VolumeUnit.keys(), name="Volume.VolumeUnit"))


class _AmountVolume(Child, _Volume):
    amount_id = Column(Integer, ForeignKey("amount.id"), nullable=False, unique=True)


class _VesselVolume(Child, _Volume):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), nullable=False, unique=True)


class Concentration(Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Concentration.ConcentrationUnit.keys(), name="Concentration.ConcentrationUnit"))


class _Pressure(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Pressure.PressureUnit.keys(), name="Pressure.PressureUnit"))


class _PressureConditionsSetpoint(Child, _Pressure):
    pressure_conditions_id = Column(Integer, ForeignKey("pressure_conditions.id"), nullable=False, unique=True)


class _PressureConditionsMeasurementPressure(Child, _Pressure):
    pressure_conditions_measurement_id = Column(
        Integer, ForeignKey("pressure_conditions_measurement.id"), nullable=False, unique=True
    )


class _Temperature(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Temperature.TemperatureUnit.keys(), name="Temperature.TemperatureUnit"))


class _ReactionInputAdditionTemperature(Child, _Temperature):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)


class _TemperatureConditionsSetpoint(Child, _Temperature):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), nullable=False, unique=True)


class _TemperatureConditionsMeasurementTemperature(Child, _Temperature):
    temperature_conditions_measurement_id = Column(
        Integer, ForeignKey("temperature_conditions_measurement.id"), nullable=False, unique=True
    )


class _Current(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Current.CurrentUnit.keys(), name="Current.CurrentUnit"))


class _ElectrochemistryConditionsCurrent(Child, _Current):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )


class _ElectrochemistryConditionsMeasurementCurrent(Child, _Current):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), nullable=False, unique=True
    )


class _Voltage(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Voltage.VoltageUnit.keys(), name="Voltage.VoltageUnit"))


class _ElectrochemistryConditionsVoltage(Child, _Voltage):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )


class _ElectrochemistryConditionsMeasurementVoltage(Child, _Voltage):
    electrochemistry_conditions_measurement_id = Column(
        Integer, ForeignKey("electrochemistry_conditions_measurement.id"), nullable=False, unique=True
    )


class _Length(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Length.LengthUnit.keys(), name="Length.LengthUnit"))


class _IlluminationConditionsDistanceToVessel(Child, _Length):
    illumination_conditions_id = Column(Integer, ForeignKey("illumination_conditions.id"), nullable=False, unique=True)


class _ElectrochemistryConditionsElectrodeSeparation(Child, _Length):
    electrochemistry_conditions_id = Column(
        Integer, ForeignKey("electrochemistry_conditions.id"), nullable=False, unique=True
    )


class _TubingDiameter(Child, _Length):
    tubing_id = Column(Integer, ForeignKey("tubing.id"), nullable=False, unique=True)


class _Wavelength(Parent, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Wavelength.WavelengthUnit.keys(), name="Wavelength.WavelengthUnit"))


class _IlluminationConditionsPeakWavelength(Child, _Wavelength):
    illumination_conditions_id = Column(Integer, ForeignKey("illumination_conditions.id"), nullable=False, unique=True)


class _ProductMeasurementWavelength(Child, _Wavelength):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class FlowRate(Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), nullable=False, unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.FlowRate.FlowRateUnit.keys(), name="FlowRate.FlowRateUnit"))


class _Percentage(Parent, Base):
    value = Column(Float)
    precision = Column(Float)


class _ReactionOutcomeConversion(Child, _Percentage):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), nullable=False, unique=True)


class _ProductMeasurementPercentage(Child, _Percentage):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class _FloatValue(Parent, Base):
    value = Column(Float)
    precision = Column(Float)


class _ProductMeasurementFloatValue(Child, _FloatValue):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), nullable=False, unique=True)


class _MassSpecMeasurementDetailsEicMasses(Child, _FloatValue):
    mass_spec_measurement_details_id = Column(Integer, ForeignKey("mass_spec_measurement_details.id"), nullable=False)


class _Data(Parent, Base):
    float_value = Column(Float)
    integer_value = Column(Integer)
    bytes_value = Column(LargeBinary)
    string_value = Column(Text)
    url = Column(Text)
    description = Column(Text)
    format = Column(String(255))


class _CompoundFeatures(Child, _Data):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)
    key = Column(String(255))  # Map key.


class _ReactionSetupAutomationCode(Child, _Data):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), nullable=False)
    key = Column(String(255))  # Map key.


class _ReactionObservationImage(Child, _Data):
    reaction_observation_id = Column(Integer, ForeignKey("reaction_observation.id"), nullable=False, unique=True)


class _ProductCompoundFeatures(Child, _Data):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False)
    key = Column(String(255))  # Map key.


class _AnalysisData(Child, _Data):
    analysis_id = Column(Integer, ForeignKey("analysis.id"), nullable=False)
    key = Column(String(255))  # Map key.


_MAPPERS: dict[Type[Message], Type[Base]] = {
    dataset_pb2.Dataset: Dataset,
    dataset_pb2.DatasetExample: DatasetExample,
    reaction_pb2.Reaction: Reaction,
    reaction_pb2.ReactionIdentifier: ReactionIdentifier,
    reaction_pb2.ReactionInput: _ReactionInput,
    reaction_pb2.ReactionInput.AdditionSpeed: AdditionSpeed,
    reaction_pb2.ReactionInput.AdditionDevice: AdditionDevice,
    reaction_pb2.Amount: _Amount,
    reaction_pb2.UnmeasuredAmount: UnmeasuredAmount,
    reaction_pb2.CrudeComponent: CrudeComponent,
    reaction_pb2.Compound: _Compound,
    reaction_pb2.Compound.Source: Source,
    reaction_pb2.CompoundPreparation: CompoundPreparation,
    reaction_pb2.CompoundIdentifier: _CompoundIdentifier,
    reaction_pb2.Vessel: Vessel,
    reaction_pb2.VesselMaterial: VesselMaterial,
    reaction_pb2.VesselAttachment: VesselAttachment,
    reaction_pb2.VesselPreparation: VesselPreparation,
    reaction_pb2.ReactionSetup: ReactionSetup,
    reaction_pb2.ReactionSetup.ReactionEnvironment: ReactionEnvironment,
    reaction_pb2.ReactionConditions: ReactionConditions,
    reaction_pb2.TemperatureConditions: _TemperatureConditions,
    reaction_pb2.TemperatureConditions.TemperatureControl: TemperatureControl,
    reaction_pb2.TemperatureConditions.Measurement: TemperatureConditionsMeasurement,
    reaction_pb2.PressureConditions: PressureConditions,
    reaction_pb2.PressureConditions.PressureControl: PressureControl,
    reaction_pb2.PressureConditions.Atmosphere: Atmosphere,
    reaction_pb2.PressureConditions.Measurement: PressureConditionsMeasurement,
    reaction_pb2.StirringConditions: _StirringConditions,
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
    reaction_pb2.DateTime: _DateTime,
    reaction_pb2.Analysis: _Analysis,
    reaction_pb2.ReactionProvenance: ReactionProvenance,
    reaction_pb2.Person: _Person,
    reaction_pb2.RecordEvent: _RecordEvent,
    reaction_pb2.Time: _Time,
    reaction_pb2.Mass: Mass,
    reaction_pb2.Moles: Moles,
    reaction_pb2.Volume: _Volume,
    reaction_pb2.Concentration: Concentration,
    reaction_pb2.Pressure: _Pressure,
    reaction_pb2.Temperature: _Temperature,
    reaction_pb2.Current: _Current,
    reaction_pb2.Voltage: _Voltage,
    reaction_pb2.Length: _Length,
    reaction_pb2.Wavelength: _Wavelength,
    reaction_pb2.FlowRate: FlowRate,
    reaction_pb2.Percentage: _Percentage,
    reaction_pb2.FloatValue: _FloatValue,
    reaction_pb2.Data: _Data,
}
_PROTOS: dict[Type[Base], Type[Message]] = {value: key for key, value in _MAPPERS.items()}

_MAPPER_RENAMES: dict[tuple[Type[Message], str], str] = {
    (reaction_pb2.Compound.Source, "id"): "vendor_id",
}
_PROTO_RENAMES: dict[tuple[Type[Base], str], str] = {
    (_MAPPERS[key[0]], value): key[1] for key, value in _MAPPER_RENAMES.items()
}


# Polymorphic matchers for constructing queries.
Amount = with_polymorphic(_Amount, "*")
Analysis = with_polymorphic(_Analysis, "*")
Compound = with_polymorphic(_Compound, "*")
CompoundIdentifier = with_polymorphic(_CompoundIdentifier, "*")
Current = with_polymorphic(_Current, "*")
Data = with_polymorphic(_Data, "*")
DateTime = with_polymorphic(_DateTime, "*")
FloatValue = with_polymorphic(_FloatValue, "*")
Length = with_polymorphic(_Length, "*")
Percentage = with_polymorphic(_Percentage, "*")
Person = with_polymorphic(_Person, "*")
Pressure = with_polymorphic(_Pressure, "*")
ReactionInput = with_polymorphic(_ReactionInput, "*")
RecordEvent = with_polymorphic(_RecordEvent, "*")
StirringConditions = with_polymorphic(_StirringConditions, "*")
Temperature = with_polymorphic(_Temperature, "*")
TemperatureConditions = with_polymorphic(_TemperatureConditions, "*")
Time = with_polymorphic(_Time, "*")
Voltage = with_polymorphic(_Voltage, "*")
Volume = with_polymorphic(_Volume, "*")
Wavelength = with_polymorphic(_Wavelength, "*")


def from_proto(  # pylint: disable=too-many-branches
    message: Message, mapper: Optional[Type[Base]] = None, key: Optional[str] = None
) -> Base:
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
        mapper = _MAPPERS[type(message)]
    kwargs = {}
    if key is not None:
        kwargs["key"] = key
    if mapper == Reaction:
        kwargs["proto"] = message.SerializeToString()
    for field, value in message.ListFields():
        field_name = _MAPPER_RENAMES.get((type(message), field.name), field.name)
        if field_name == "eic_masses":
            # Convert repeated float to repeated FloatValue.
            kwargs[field_name] = [_MassSpecMeasurementDetailsEicMasses(value=v) for v in value]
        elif field_name == "reaction_ids":
            # Convert repeated string to repeated ReactionId.
            kwargs[field_name] = [ReactionId(reaction_id=v) for v in value]
        elif field.type == FieldDescriptor.TYPE_MESSAGE:
            field_mapper = getattr(mapper, field_name).mapper.class_
            if isinstance(value, Mapping):
                kwargs[field_name] = [from_proto(v, mapper=field_mapper, key=k) for k, v in value.items()]
            elif field.label == FieldDescriptor.LABEL_REPEATED:
                kwargs[field_name] = [from_proto(v, mapper=field_mapper) for v in value]
            else:
                kwargs[field_name] = from_proto(value, mapper=field_mapper)
        elif field.type == FieldDescriptor.TYPE_ENUM:
            kwargs[field_name] = field.enum_type.values_by_number[value].name
        else:
            kwargs[field_name] = value
    if isinstance(message, (reaction_pb2.Compound, reaction_pb2.ProductCompound)):
        # Add RDKit cartridge functionality.
        field_mapper = getattr(mapper, "structure").mapper.class_
        try:
            kwargs["structure"] = field_mapper(smiles=message_helpers.smiles_from_compound(message))
        except ValueError:
            pass
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
            proto = _PROTOS[mapper]
            break
    assert issubclass(proto, Message)
    for field in proto.DESCRIPTOR.fields:
        mapper_field_name = _MAPPER_RENAMES.get((proto, field.name), field.name)
        value = getattr(base, mapper_field_name)
        field_name = _PROTO_RENAMES.get((type(base), mapper_field_name), mapper_field_name)
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
