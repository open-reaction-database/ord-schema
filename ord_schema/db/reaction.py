"""Table mappings for Reaction protos.

Notes:
    * Foreign keys to the `reaction` table are done using the `id` column, not the ORD reaction ID (`name`).
    * Every table has its own `id` column, even if it has a one-to-one mapping with another table.
    * When a message type is used in multiple places, use Single Table Inheritance; see
      https://docs.sqlalchemy.org/en/14/orm/inheritance.html#single-table-inheritance.
"""
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

from ord_schema.proto import reaction_pb2

Base = declarative_base()

@declarative_mixin
class OrdMixin:
    """See https://docs.sqlalchemy.org/en/14/orm/declarative_mixins.html."""

    @classmethod
    @declared_attr
    def __tablename__(cls):
        return underscore(cls.__name__)

    id = Column(Integer, primary_key=True)

class Reaction(OrdMixin, Base):
    name = Column(String(32), nullable=False)
    proto = Column(LargeBinary, nullable=False)

    identifiers = relationship("ReactionIdentifier")
    inputs = relationship("ReactionInput")
    setup = relationship("ReactionSetup", uselist=False)
    conditions = relationship("ReactionConditions", uselist=False)
    notes = relationship("ReactionNotes", uselist=False)
    observations = relationship("ReactionObservation")
    workups = relationship("ReactionWorkup")
    outcomes = relationship("ReactionOutcome")
    provenance = relationship("ReactionProvenance", uselist=False)

class ReactionIdentifier(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ReactionIdentifier.IdentifierType.keys()))
    details = Column(Text)
    value = Column(Text)
    is_mapped = Column(Boolean)

class ReactionInput(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    name = Column(Text)  # Map key.
    components = relationship("Compound")
    crude_components = relationship("CrudeComponent")
    addition_order = Column(Integer)
    addition_time = relationship("Time", uselist=False)
    addition_speed = relationship("AdditionSpeed", uselist=False)
    addition_duration = relationship("Time", uselist=False)
    flow_rate = relationship("FlowRate", uselist=False)
    addition_device = relationship("AdditionDevice", uselist=False)
    addition_temperature = relationship("Temperature", uselist=False)

class AdditionSpeed(OrdMixin, Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ReactionInput.AdditionSpeed.AdditionSpeedType.keys()))
    details = Column(Text)


class AdditionDevice(OrdMixin, Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ReactionInput.AdditionDevice.AdditionDeviceType.keys()))
    details = Column(Text)


class Amount(OrdMixin, Base):
    type = Column(String(255), nullable=False)

    mass = relationship("Mass", uselist=False)
    moles = relationship("Moles", uselist=False)
    volume = relationship("Volume", uselist=False)
    unmeasured = relationship("UnmeasuredAmount", uselist=False)
    volume_includes_solutes = Column(Boolean)

    __mapper_args__ = {
        "polymorphic_identity": "amount",
        "polymorphic_on": type,
    }


class CrudeComponentAmount(Amount):
    crude_component_id = Column(Integer, ForeignKey("crude_component.id"))

    __mapper_args__ = {
        "polymorphic_identity": "crude_component_amount",
    }

class CompoundAmount(Amount):
    compound_id = Column(Integer, ForeignKey("compound.id"))

    __mapper_args__ = {
        "polymorphic_identity": "compound_amount",
    }

class ReactionWorkupAmount(Amount):
    reaction_workup_id = Column(Integer, ForeignKey("reaction_workup.id"))

    __mapper_args__ = {
        "polymorphic_identity": "reaction_workup_amount",
    }


class ProductMeasurementAmount(Amount):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"))

    __mapper_args__ = {
        "polymorphic_identity": "crude_component_amount",
    }


class UnmeasuredAmount(OrdMixin, Base):
    amount_id = Column(Integer, ForeignKey("amount.id"), primary_key=True)
    type = Column(Enum(reaction_pb2.UnmeasuredAmount.UnmeasuredAmountType.keys()))
    details = Column(Text)

class CrudeComponent(OrdMixin, Base):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"), primary_key=True)
    includes_workup = Column(Boolean)
    has_derived_amount = Column(Boolean)
    amount = relationship("Amount", uselist=False)

class Compound(OrdMixin, Base):
    type = Column(String(255), nullable=False)

    identifiers = relationship("CompoundIdentifier")
    amount = relationship("Amount")
    reaction_role = Column(Enum(reaction_pb2.ReactionRole.ReactionRoleType.keys()))
    is_limiting = Column(Boolean)
    preparations = relationship("CompoundPreparation")
    source = relationship("Source", uselist=False)
    features = relationship("Data")
    analyses = relationship("Analysis")

    __mapper_args__ = {
        "polymorphic_identity": "compound",
        "polymorphic_on": type,
    }


class ReactionInputCompound(Compound):
    reaction_input_id = Column(Integer, ForeignKey("reaction_input.id"))

    __mapper_args__ = {
        "polymorphic_identity": "reaction_input_compound",
    }


class ProductMeasurementCompound(Compound):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"))

    __mapper_args__ = {
        "polymorphic_identity": "product_measurement_compound",
    }


class CompoundPreparation(OrdMixin, Base):
    compound_id = Column(Integer, ForeignKey("compound.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.CompoundPreparation.PreparationType.keys()))
    details = Column(Text)
    reaction_id = Column(Integer, ForeignKey("reaction.id"))


class CompoundIdentifier(OrdMixin, Base):
    compound_id = Column(Integer, ForeignKey("compound.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.CompoundIdentifier.IdentifierType.keys()))
    details = Column(Text)
    value = Column(Text)


class Vessel(OrdMixin, Base):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.Vessel.VesselType.keys()))
    details = Column(Text)
    material = relationship("VesselMaterial", uselist=False)
    preparations = relationship("VesselPreparation")
    attachments = relationship("VesselAttachment")
    volume = relationship("Volume", uselist=False)
    plate_id = Column(String(255))
    plate_position = Column(String(32))


class VesselMaterial(OrdMixin, Base):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.VesselMaterial.VesselMaterialType.keys()))
    details = Column(Text)

class VesselAttachment(OrdMixin, Base):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.VesselAttachment.VesselAttachmentType.keys()))
    details = Column(Text)


class VesselPreparation(OrdMixin, Base):
    vessel_id = Column(Integer, ForeignKey("vessel.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.VesselPreparation.VesselPreparationType.keys()))
    details = Column(Text)

class ReactionSetup(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    vessel = relationship("Vessel", uselist=False)
    is_automated = Column(Boolean)
    automation_platform = Column(String(255))
    automation_code = relationship("Data")
    environment = relationship("ReactionEnvironment", uselist=False)


class ReactionEnvironment(OrdMixin, Base):
    reaction_setup_id = Column(Integer, ForeignKey("reaction_setup.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType.keys()))
    details = Column(Text)


class ReactionConditions(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    temperature = relationship("TemperatureConditions", uselist=False)
    pressure = relationship("PressureConditions", uselist=False)
    stirring = relationship("StirringConditions", uselist=False)
    illumination = relationship("IlluminationConditions", uselist=False)
    electrochemistry = relationship("ElectrochemistryConditions", uselist=False)
    flow = relationship("FlowConditions", uselist=False)
    reflux = Column(Boolean)
    ph = Column(Float)
    conditions_are_dynamic = Column(Boolean)
    details = Column(Text)


class TemperatureConditions(OrdMixin, Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), primary_key=True)

    control = relationship("TemperatureControl", uselist=False)
    setpoint = relationship("Temperature", uselist=False)
    measurements = relationship("TemperatureConditionsMeasurement")


class TemperatureControl(OrdMixin, Base):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.TemperatureConditions.TemperatureControl.TemperatureControlType.keys()))
    details = Column(Text)

class TemperatureControlMeasurement(OrdMixin, Base):
    temperature_conditions_id = Column(Integer, ForeignKey("temperature_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.TemperatureConditions.Measurement.MeasurementType.keys()))
    details = Column(Text)
    time = relationship("Time", uselist=False)
    temperature = relationship("Temperature", uselist=False)


class PressureConditions(OrdMixin, Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), primary_key=True)

    control = relationship("PressureControl", uselist=False)
    setpoint = relationship("Pressure", uselist=False)
    atmosphere = relationship("Atmosphere", uselist=False)
    measurements = relationship("PressureConditionsMeasurement")


class PressureControl(OrdMixin, Base):
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.PressureConditions.PressureControl.PressureControlType.keys()))
    details = Column(Text)

class Atmosphere(OrdMixin, Base):
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.PressureConditions.Atmosphere.AtmosphereType.keys()))
    details = Column(Text)


class PressureControlMeasurement(OrdMixin, Base):
    pressure_conditionss_id = Column(Integer, ForeignKey("pressure_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.PressureConditions.Measurement.MeasurementType.keys()))
    details = Column(Text)
    time = relationship("Time", uselist=False)
    pressure = relationship("Pressure", uselist=False)


class StirringConditions(OrdMixin, Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.StirringConditions.StirringConditionsType.keys()))
    details = Column(Text)
    rate = relationship("StirringRate", uselist=False)


class StirringRate(OrdMixin, Base):
    stirring_conditions_id = Column(Integer, ForeignKey("stirring_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.StirringConditions.StirringRate.StirringRateType.keys()))
    details = Column(Text)
    rpm = Column(Integer)


class IlluminationConditions(OrdMixin, Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.IlluminationConditions.IlluminationType.keys()))
    details = Column(Text)
    peak_wavelength = relationship("Wavelength", uselist=False)
    color = Column(String(255))
    distance_to_vessel = relationship("Length", uselist=False)

class ElectrochemistryConditions(OrdMixin, Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ElectrochemistryConditions.ElectrochemistryType.keys()))
    details = Column(Text)
    current = relationship("Current", uselist=False)
    voltage = relationship("Voltage", uselist=False)
    anode_material = Column(String(255))
    cathode_material = Column(String(255))
    electrode_separation = relationship("Length", uselist=False)
    measurements = relationship("ElectrochemistryConditionsMeasurement")
    cell = relationship("ElectrochemistryCell", uselist=False)


class ElectrochemistryConditionsMeasurement(OrdMixin, Base):
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), primary_key=True)

    time = relationship("Time", uselist=False)
    current = relationship("Current", uselist=False)
    voltage = relationship("Voltage", uselist=False)


class ElectrochemistryCell(OrdMixin, Base):
    electrochemistry_conditions_id = Column(Integer, ForeignKey("electrochemistry_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType.keys()))
    details = Column(Text)


class FlowConditions(OrdMixin, Base):
    reaction_conditions_id = Column(Integer, ForeignKey("reaction_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.FlowConditions.FlowType.keys()))
    details = Column(Text)
    pump_type = Column(Text)
    tubing = relationship("Tubing", uselist=False)


class Tubing(OrdMixin, Base):
    flow_conditions_id = Column(Integer, ForeignKey("flow_conditions.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.FlowConditions.Tubing.TubingType.keys()))
    details = Column(Text)
    diameter = relationship("Length", uselist=False)

class ReactionNotes(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    is_heterogeneous = Column(Boolean)
    forms_precipitate = Column(Boolean)
    is_exothermic = Column(Boolean)
    offgasses = Column(Boolean)
    is_sensitive_to_moisture = Column(Boolean)
    is_sensitive_to_oxygen = Column(Boolean)
    is_sensitive_to_light = Column(Boolean)
    safety_notes = Column(Text)
    procedure_details = Column(Text)

class ReactionObservation(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    time = relationship("Time", uselist=False)
    comment = Column(Text)
    image = relationship("Data", uselist=False)

class ReactionWorkup(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ReactionWorkup.WorkupType.keys()))
    details = Column(Text)
    duration = relationship("Time", uselist=False)
    input = relationship("ReactionInput", uselist=False)
    amount = relationship("Amount", uselist=False)
    temperature = relationship("TemperatureConditions", uselist=False)
    keep_phase = Column(String(255))
    stirring = relationship("StirringConditions", uselist=False)
    target_ph = Column(Float)
    is_automated = Column(Boolean)


class ReactionOutcome(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    reaction_time = relationship("Time", uselist=False)
    conversion = relationship("Percentage", uselist=False)
    products = relationship("ProductCompound")
    analyses = relationship("Analysis")


class ProductCompound(OrdMixin, Base):
    reaction_outcome_id = Column(Integer, ForeignKey("reaction_outcome.id"), primary_key=True)

    identifiers = relationship("CompoundIdentifier")
    is_desired_product = Column(Boolean)
    measurements = relationship("ProductMeasurement")
    isolated_color = Column(String(255))
    texture = relationship("Texture", uselist=False)
    features = relationship("Data")
    reaction_role = Column(Enum(reaction_pb2.ReactionRole.ReactionRoleType.keys()))


class Texture(OrdMixin, Base):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ProductCompound.Texture.TextureType.keys()))
    details = Column(Text)

class ProductMeasurement(OrdMixin, Base):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), primary_key=True)

    analysis_key = Column(String(255))
    type = Column(Enum(reaction_pb2.ProductMeasurement.MeasurementType.keys()))
    details = Column(Text)
    uses_internal_standard = Column(Boolean)
    is_normalized = Column(Boolean)
    uses_authentic_standard = Column(Boolean)
    authentic_standard = relationship("Compound", uselist=False)
    percentage = relationship("Percentage", uselist=False)
    float_value = relationship("FloatValue", uselist=False)
    string_value = Column(Text)
    amount = relationship("Amount", uselist=False)
    retention_time = relationship("Time", uselist=False)
    mass_spec_details = relationship("MassSpecMeasurementDetails", uselist=False)
    selectivity = relationship("Selectivity", uselist=False)
    wavelength = relationship("Wavelength", uselist=False)

class MassSpecMeasurementDetails(OrdMixin, Base):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType.keys()))
    details = Column(Text)
    tic_minimum_mz = Column(Float)
    tic_maximum_mz = Column(Float)
    eic_masses = relationship("FloatValue")


class Selectivity(OrdMixin, Base):
    product_measurement_id = Column(Integer, ForeignKey("product_measurement.id"), primary_key=True)

    type = Column(Enum(reaction_pb2.ProductMeasurement.Selectivity.SelectivityType.keys()))
    details = Column(Text)


class DateTime(OrdMixin, Base):
    value = Column(String(255))

class Analysis(OrdMixin, Base):
    type = Column(Enum(reaction_pb2.Analysis.AnalysisType.keys()))
    details = Column(Text)
    chmo_id = Column(Integer)
    is_of_isolated_species = Column(Boolean)
    data = relationship("Data")
    instrument_manufacturer = Column(String(255))
    instrument_last_calibrated = relationship("DateTime", uselist=False)


class ReactionProvenance(OrdMixin, Base):
    reaction_id = Column(Integer, ForeignKey("reaction.id"), primary_key=True)

    experimenter = relationship("Person", uselist=False)
    city = Column(String(255))
    experiment_start = relationship("DateTime", uselist=False)
    doi = Column(String(255))
    patent = Column(String(255))
    publication_url = Column(Text)
    record_created = relationship("RecordEvent", uselist=False)
    record_modified = relationship("RecordEvent")

class Person(OrdMixin, Base):
    reaction_provenance_id = Column(Integer, ForeignKey("reaction_provenance.id"), primary_key=True)

    username = Column(String(255))
    name = Column(String(255))
    orcid = Column(String(19))
    organization = Column(String(255))
    email = Column(String(255))


class RecordEvent(OrdMixin, Base):
    time = relationship("DateTime", uselist=False)
    person = relationship("Person", uselist=False)
    details = Column(Text)

class Time(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Time.TimeUnit.keys()))


class Mass(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Mass.MassUnit.keys()))


class Moles(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Moles.MolesUnit.keys()))



class Volume(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Volume.VolumeUnit.keys()))


class Concentration(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Concentration.ConcentrationUnit.keys()))


class Pressure(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Pressure.PressureUnit.keys()))


class Temperature(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Temperature.TemperatureUnit.keys()))


class Current(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Current.CurrentUnit.keys()))


class Voltage(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Voltage.VoltageUnit.keys()))


class Length(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Length.LengthUnit.keys()))


class Wavelength(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.Wavelength.WavelengthUnit.keys()))


class FlowRate(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(reaction_pb2.FlowRate.FlowRateUnit.keys()))


class Percentage(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)


class FloatValue(OrdMixin, Base):
    value = Column(Float)
    precision = Column(Float)


class Data(OrdMixin, Base):
    float_value = Column(Float)
    integer_value = Column(Integer)
    bytes_value = Column(LargeBinary)
    string_value = Column(Text)
    url = Column(Text)
    description = Column(Text)
    format = Column(String(255))