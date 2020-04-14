# Protocol Documentation
<a name="top"></a>

## Table of Contents

- [reaction.proto](#reaction.proto)
    - [Compound](#ord.Compound)
    - [Compound.Feature](#ord.Compound.Feature)
    - [Compound.ReactionRole](#ord.Compound.ReactionRole)
    - [CompoundIdentifier](#ord.CompoundIdentifier)
    - [CompoundPreparation](#ord.CompoundPreparation)
    - [Concentration](#ord.Concentration)
    - [Current](#ord.Current)
    - [Data](#ord.Data)
    - [DateTime](#ord.DateTime)
    - [ElectrochemistryConditions](#ord.ElectrochemistryConditions)
    - [ElectrochemistryConditions.Measurement](#ord.ElectrochemistryConditions.Measurement)
    - [FlowConditions](#ord.FlowConditions)
    - [FlowConditions.Tubing](#ord.FlowConditions.Tubing)
    - [FlowRate](#ord.FlowRate)
    - [IlluminationConditions](#ord.IlluminationConditions)
    - [Length](#ord.Length)
    - [Mass](#ord.Mass)
    - [Moles](#ord.Moles)
    - [Percentage](#ord.Percentage)
    - [Person](#ord.Person)
    - [Pressure](#ord.Pressure)
    - [PressureConditions](#ord.PressureConditions)
    - [PressureConditions.Atmosphere](#ord.PressureConditions.Atmosphere)
    - [PressureConditions.Measurement](#ord.PressureConditions.Measurement)
    - [PressureConditions.PressureControl](#ord.PressureConditions.PressureControl)
    - [Reaction](#ord.Reaction)
    - [Reaction.InputsEntry](#ord.Reaction.InputsEntry)
    - [ReactionAnalysis](#ord.ReactionAnalysis)
    - [ReactionAnalysis.ProcessedDataEntry](#ord.ReactionAnalysis.ProcessedDataEntry)
    - [ReactionAnalysis.RawDataEntry](#ord.ReactionAnalysis.RawDataEntry)
    - [ReactionConditions](#ord.ReactionConditions)
    - [ReactionIdentifier](#ord.ReactionIdentifier)
    - [ReactionInput](#ord.ReactionInput)
    - [ReactionInput.AdditionSpeed](#ord.ReactionInput.AdditionSpeed)
    - [ReactionNotes](#ord.ReactionNotes)
    - [ReactionObservation](#ord.ReactionObservation)
    - [ReactionOutcome](#ord.ReactionOutcome)
    - [ReactionOutcome.AnalysesEntry](#ord.ReactionOutcome.AnalysesEntry)
    - [ReactionProduct](#ord.ReactionProduct)
    - [ReactionProduct.Texture](#ord.ReactionProduct.Texture)
    - [ReactionProvenance](#ord.ReactionProvenance)
    - [ReactionProvenance.RecordEvent](#ord.ReactionProvenance.RecordEvent)
    - [ReactionSetup](#ord.ReactionSetup)
    - [ReactionSetup.AutomationCodeEntry](#ord.ReactionSetup.AutomationCodeEntry)
    - [ReactionWorkup](#ord.ReactionWorkup)
    - [Selectivity](#ord.Selectivity)
    - [StirringConditions](#ord.StirringConditions)
    - [StirringConditions.StirringMethod](#ord.StirringConditions.StirringMethod)
    - [StirringConditions.StirringRate](#ord.StirringConditions.StirringRate)
    - [Temperature](#ord.Temperature)
    - [TemperatureConditions](#ord.TemperatureConditions)
    - [TemperatureConditions.Measurement](#ord.TemperatureConditions.Measurement)
    - [TemperatureConditions.TemperatureControl](#ord.TemperatureConditions.TemperatureControl)
    - [Time](#ord.Time)
    - [Vessel](#ord.Vessel)
    - [Vessel.VesselMaterial](#ord.Vessel.VesselMaterial)
    - [Vessel.VesselPreparation](#ord.Vessel.VesselPreparation)
    - [Vessel.VesselType](#ord.Vessel.VesselType)
    - [Voltage](#ord.Voltage)
    - [Volume](#ord.Volume)
    - [Wavelength](#ord.Wavelength)
  
    - [Compound.ReactionRole.ReactionRoleType](#ord.Compound.ReactionRole.ReactionRoleType)
    - [CompoundIdentifier.IdentifierType](#ord.CompoundIdentifier.IdentifierType)
    - [CompoundPreparation.PreparationType](#ord.CompoundPreparation.PreparationType)
    - [Concentration.ConcentrationUnit](#ord.Concentration.ConcentrationUnit)
    - [Current.CurrentUnit](#ord.Current.CurrentUnit)
    - [ElectrochemistryConditions.ElectrochemistryType](#ord.ElectrochemistryConditions.ElectrochemistryType)
    - [FlowConditions.FlowType](#ord.FlowConditions.FlowType)
    - [FlowConditions.Tubing.TubingMaterialType](#ord.FlowConditions.Tubing.TubingMaterialType)
    - [FlowRate.FlowRateUnit](#ord.FlowRate.FlowRateUnit)
    - [IlluminationConditions.IlluminationType](#ord.IlluminationConditions.IlluminationType)
    - [Length.LengthUnit](#ord.Length.LengthUnit)
    - [Mass.MassUnit](#ord.Mass.MassUnit)
    - [Moles.MolesUnit](#ord.Moles.MolesUnit)
    - [Pressure.PressureUnit](#ord.Pressure.PressureUnit)
    - [PressureConditions.Atmosphere.AtmosphereType](#ord.PressureConditions.Atmosphere.AtmosphereType)
    - [PressureConditions.Measurement.MeasurementType](#ord.PressureConditions.Measurement.MeasurementType)
    - [PressureConditions.PressureControl.PressureControlType](#ord.PressureConditions.PressureControl.PressureControlType)
    - [ReactionAnalysis.AnalysisType](#ord.ReactionAnalysis.AnalysisType)
    - [ReactionIdentifier.IdentifierType](#ord.ReactionIdentifier.IdentifierType)
    - [ReactionInput.AdditionSpeed.AdditionSpeedType](#ord.ReactionInput.AdditionSpeed.AdditionSpeedType)
    - [ReactionProduct.Texture.TextureType](#ord.ReactionProduct.Texture.TextureType)
    - [ReactionWorkup.WorkupType](#ord.ReactionWorkup.WorkupType)
    - [Selectivity.SelectivityType](#ord.Selectivity.SelectivityType)
    - [StirringConditions.StirringMethod.StirringMethodType](#ord.StirringConditions.StirringMethod.StirringMethodType)
    - [StirringConditions.StirringRate.StirringRateType](#ord.StirringConditions.StirringRate.StirringRateType)
    - [Temperature.TemperatureUnit](#ord.Temperature.TemperatureUnit)
    - [TemperatureConditions.Measurement.MeasurementType](#ord.TemperatureConditions.Measurement.MeasurementType)
    - [TemperatureConditions.TemperatureControl.TemperatureControlType](#ord.TemperatureConditions.TemperatureControl.TemperatureControlType)
    - [Time.TimeUnit](#ord.Time.TimeUnit)
    - [Vessel.VesselMaterial.VesselMaterialType](#ord.Vessel.VesselMaterial.VesselMaterialType)
    - [Vessel.VesselPreparation.VesselPreparationType](#ord.Vessel.VesselPreparation.VesselPreparationType)
    - [Vessel.VesselType.VesselTypeEnum](#ord.Vessel.VesselType.VesselTypeEnum)
    - [Voltage.VoltageUnit](#ord.Voltage.VoltageUnit)
    - [Volume.VolumeUnit](#ord.Volume.VolumeUnit)
    - [Wavelength.WavelengthUnit](#ord.Wavelength.WavelengthUnit)
  
  
  

- [Scalar Value Types](#scalar-value-types)



<a name="reaction.proto"></a>
<p align="right"><a href="#top">Top</a></p>

## reaction.proto
Schema for the Open Reaction Database.


<a name="ord.Compound"></a>

### Compound



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| identifiers | [CompoundIdentifier](#ord.CompoundIdentifier) | repeated | Set of identifiers used to uniquely define this compound. Solutions or mixed compounds should use the NAME identifier and list all constituent compounds in the &#34;components&#34; field. |
| mass | [Mass](#ord.Mass) |  |  |
| moles | [Moles](#ord.Moles) |  |  |
| volume | [Volume](#ord.Volume) |  |  |
| reaction_role | [Compound.ReactionRole.ReactionRoleType](#ord.Compound.ReactionRole.ReactionRoleType) |  |  |
| is_limiting | [bool](#bool) |  | Whether this species was intended to be a limiting reactant. |
| preparation | [CompoundPreparation](#ord.CompoundPreparation) |  |  |
| vendor_source | [string](#string) |  | Name of the vendor or supplier the compound was purchased from. |
| vendor_id | [string](#string) |  | Compound ID in the vendor database or catalog. |
| vendor_lot | [string](#string) |  | Batch/lot identification. |
| features | [Compound.Feature](#ord.Compound.Feature) | repeated |  |






<a name="ord.Compound.Feature"></a>

### Compound.Feature
Compounds can accommodate any number of features. These may include simple
properties of the compound (e.g., molecular weight), heuristic estimates
of physical properties (e.g., ClogP), optimized geometries (e.g., through
DFT), and calculated stereoselectronic descriptors.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [string](#string) |  |  |
| string_value | [string](#string) |  |  |
| float_value | [float](#float) |  |  |
| how_computed | [string](#string) |  |  |






<a name="ord.Compound.ReactionRole"></a>

### Compound.ReactionRole







<a name="ord.CompoundIdentifier"></a>

### CompoundIdentifier
Compound identifiers uniquely define a single (pure) chemical species.
While we encourage the use of SMILES strings, these do not work well in
all cases (e.g., handling tautomerism, axial chirality). Multiple
identifiers may be specified for a single compound to avoid ambiguity.
We discourage chemicals from being defined only by a name. For compounds
that are prepared or isolated as salts, the identifier should include
specification of which salt.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [CompoundIdentifier.IdentifierType](#ord.CompoundIdentifier.IdentifierType) |  |  |
| details | [string](#string) |  |  |
| value | [string](#string) |  |  |
| bytes_value | [bytes](#bytes) |  |  |






<a name="ord.CompoundPreparation"></a>

### CompoundPreparation
Compounds may undergo additional preparation before being used in a
reaction after being received from a supplier or vendor. We encourage
the use of the &#39;preparation&#39; enum when possible, even if the description
is an oversimplification of the full procedure, which can be described
in the &#39;details&#39; field.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [CompoundPreparation.PreparationType](#ord.CompoundPreparation.PreparationType) |  |  |
| details | [string](#string) |  | Full description of how the received compound was prepared. |






<a name="ord.Concentration"></a>

### Concentration



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Concentration.ConcentrationUnit](#ord.Concentration.ConcentrationUnit) |  |  |






<a name="ord.Current"></a>

### Current



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Current.CurrentUnit](#ord.Current.CurrentUnit) |  |  |






<a name="ord.Data"></a>

### Data
Data is a container for arbitrary string or bytes data.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [string](#string) |  |  |
| bytes_value | [bytes](#bytes) |  |  |
| url | [string](#string) |  | URL for data stored elsewhere. |
| description | [string](#string) |  |  |
| format | [string](#string) |  | Description of the file format (if applicable); usually the file extension. For example, &#39;png&#39; or &#39;tiff&#39; for images. If empty, we assume string data. |






<a name="ord.DateTime"></a>

### DateTime
TODO(ccoley): If we want the DateTime to be a string that we parse as
needed, should it simply be &#34;string datetime&#34; when used? Or is there any 
benefit to having a separate message type that could be changed in the 
future if needed?


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [string](#string) |  |  |






<a name="ord.ElectrochemistryConditions"></a>

### ElectrochemistryConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [ElectrochemistryConditions.ElectrochemistryType](#ord.ElectrochemistryConditions.ElectrochemistryType) |  |  |
| details | [string](#string) |  |  |
| current | [Current](#ord.Current) |  |  |
| voltage | [Voltage](#ord.Voltage) |  |  |
| anode_material | [string](#string) |  |  |
| cathode_material | [string](#string) |  |  |
| electrode_separation | [Length](#ord.Length) |  |  |
| measurements | [ElectrochemistryConditions.Measurement](#ord.ElectrochemistryConditions.Measurement) | repeated |  |






<a name="ord.ElectrochemistryConditions.Measurement"></a>

### ElectrochemistryConditions.Measurement



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| time | [Time](#ord.Time) |  |  |
| current | [Current](#ord.Current) |  |  |
| voltage | [Voltage](#ord.Voltage) |  |  |






<a name="ord.FlowConditions"></a>

### FlowConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [FlowConditions.FlowType](#ord.FlowConditions.FlowType) |  |  |
| details | [string](#string) |  |  |
| pump_type | [string](#string) |  |  |
| tubing | [FlowConditions.Tubing](#ord.FlowConditions.Tubing) |  |  |






<a name="ord.FlowConditions.Tubing"></a>

### FlowConditions.Tubing



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [FlowConditions.Tubing.TubingMaterialType](#ord.FlowConditions.Tubing.TubingMaterialType) |  |  |
| details | [string](#string) |  |  |
| diameter | [Length](#ord.Length) |  |  |






<a name="ord.FlowRate"></a>

### FlowRate



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [FlowRate.FlowRateUnit](#ord.FlowRate.FlowRateUnit) |  |  |






<a name="ord.IlluminationConditions"></a>

### IlluminationConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [IlluminationConditions.IlluminationType](#ord.IlluminationConditions.IlluminationType) |  |  |
| details | [string](#string) |  |  |
| peak_wavelength | [Wavelength](#ord.Wavelength) |  |  |
| color | [string](#string) |  |  |
| distance_to_vessel | [Length](#ord.Length) |  |  |






<a name="ord.Length"></a>

### Length



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Length.LengthUnit](#ord.Length.LengthUnit) |  |  |






<a name="ord.Mass"></a>

### Mass



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Mass.MassUnit](#ord.Mass.MassUnit) |  |  |






<a name="ord.Moles"></a>

### Moles



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Moles.MolesUnit](#ord.Moles.MolesUnit) |  |  |






<a name="ord.Percentage"></a>

### Percentage
Used for things like conversion and yield.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |






<a name="ord.Person"></a>

### Person



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| username | [string](#string) |  |  |
| name | [string](#string) |  |  |
| orcid | [string](#string) |  |  |
| organization | [string](#string) |  |  |






<a name="ord.Pressure"></a>

### Pressure



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Pressure.PressureUnit](#ord.Pressure.PressureUnit) |  |  |






<a name="ord.PressureConditions"></a>

### PressureConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [PressureConditions.PressureControl.PressureControlType](#ord.PressureConditions.PressureControl.PressureControlType) |  |  |
| details | [string](#string) |  |  |
| setpoint | [Pressure](#ord.Pressure) |  |  |
| atmosphere | [PressureConditions.Atmosphere.AtmosphereType](#ord.PressureConditions.Atmosphere.AtmosphereType) |  |  |
| atmosphere_details | [string](#string) |  |  |
| measurements | [PressureConditions.Measurement](#ord.PressureConditions.Measurement) | repeated |  |






<a name="ord.PressureConditions.Atmosphere"></a>

### PressureConditions.Atmosphere







<a name="ord.PressureConditions.Measurement"></a>

### PressureConditions.Measurement



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [PressureConditions.Measurement.MeasurementType](#ord.PressureConditions.Measurement.MeasurementType) |  |  |
| details | [string](#string) |  |  |
| time | [Time](#ord.Time) |  |  |
| pressure | [Pressure](#ord.Pressure) |  |  |






<a name="ord.PressureConditions.PressureControl"></a>

### PressureConditions.PressureControl







<a name="ord.Reaction"></a>

### Reaction
Throughout this schema, we introduce enums to encourage consistency in
nomenclature and to avoid unnecessary downstream data processing that would
otherwise be required to consolidate equivalent entries. However, we do
not wish to restrict what users are able to specify if their synthesis
does not fit cleanly into a pre-existing enum field. For that reason, many
enums contain a CUSTOM field, which must be accompanied by setting the
&#39;details&#39; field (or &#39;&lt;field_name&gt;_details&#39;, where appropriate).

NOTE(kearnes): In many places, we deliberately violate the style guide for
enums by nesting instead of prefixing; this is not done lightly. The primary
consideration is API consistency and the ability to use unqualified strings
as enum values. For instance, we want &#39;CUSTOM&#39; to be a valid value for all
enums that support custom types.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| identifiers | [ReactionIdentifier](#ord.ReactionIdentifier) | repeated |  |
| inputs | [Reaction.InputsEntry](#ord.Reaction.InputsEntry) | repeated | List of pure substances or mixtures that were added to the reaction vessel. This is a map, not a repeated, to simplify reaction templating through the use of keys. String keys are simple descriptions and are present only for convenience. |
| setup | [ReactionSetup](#ord.ReactionSetup) |  |  |
| conditions | [ReactionConditions](#ord.ReactionConditions) |  |  |
| notes | [ReactionNotes](#ord.ReactionNotes) |  | Reaction notes largely pertain to safety considerations. |
| observations | [ReactionObservation](#ord.ReactionObservation) | repeated |  |
| workup | [ReactionWorkup](#ord.ReactionWorkup) | repeated | Workup steps are listed in the order they are performed. |
| outcomes | [ReactionOutcome](#ord.ReactionOutcome) | repeated |  |
| provenance | [ReactionProvenance](#ord.ReactionProvenance) |  |  |






<a name="ord.Reaction.InputsEntry"></a>

### Reaction.InputsEntry



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| key | [string](#string) |  |  |
| value | [ReactionInput](#ord.ReactionInput) |  |  |






<a name="ord.ReactionAnalysis"></a>

### ReactionAnalysis



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [ReactionAnalysis.AnalysisType](#ord.ReactionAnalysis.AnalysisType) |  |  |
| details | [string](#string) |  | Any details about analysis (e.g., NMR type, columns, gradients, conditions) |
| processed_data | [ReactionAnalysis.ProcessedDataEntry](#ord.ReactionAnalysis.ProcessedDataEntry) | repeated | Data files (processed or annotated). |
| raw_data | [ReactionAnalysis.RawDataEntry](#ord.ReactionAnalysis.RawDataEntry) | repeated | Data files (raw) obtained directly from the instrument |
| instrument_manufacturer | [string](#string) |  |  |
| instrument_last_calibrated | [DateTime](#ord.DateTime) |  |  |






<a name="ord.ReactionAnalysis.ProcessedDataEntry"></a>

### ReactionAnalysis.ProcessedDataEntry



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| key | [string](#string) |  |  |
| value | [Data](#ord.Data) |  |  |






<a name="ord.ReactionAnalysis.RawDataEntry"></a>

### ReactionAnalysis.RawDataEntry



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| key | [string](#string) |  |  |
| value | [Data](#ord.Data) |  |  |






<a name="ord.ReactionConditions"></a>

### ReactionConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| temperature | [TemperatureConditions](#ord.TemperatureConditions) |  |  |
| pressure | [PressureConditions](#ord.PressureConditions) |  |  |
| stirring | [StirringConditions](#ord.StirringConditions) |  |  |
| illumination | [IlluminationConditions](#ord.IlluminationConditions) |  |  |
| electrochemistry | [ElectrochemistryConditions](#ord.ElectrochemistryConditions) |  |  |
| flow | [FlowConditions](#ord.FlowConditions) |  |  |
| reflux | [bool](#bool) |  |  |
| pH | [float](#float) |  |  |
| conditions_are_dynamic | [bool](#bool) |  | Boolean to describe whether the conditions cannot be represented by the static, single-step schema. |
| details | [string](#string) |  | A catch-all string field for providing more information about the conditions (e.g., multiple stages) |






<a name="ord.ReactionIdentifier"></a>

### ReactionIdentifier
Reaction identifiers define descriptions of the overall reaction.
While we encourage the use of SMILES strings, these do not work well in
all cases. The &lt;reaction_smiles&gt; field should be able to be derived
from the information present in the ReactionInput and ReactionOutcome
fields of any Reaction message.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [ReactionIdentifier.IdentifierType](#ord.ReactionIdentifier.IdentifierType) |  |  |
| details | [string](#string) |  |  |
| value | [string](#string) |  |  |
| bytes_value | [bytes](#bytes) |  |  |






<a name="ord.ReactionInput"></a>

### ReactionInput



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| components | [Compound](#ord.Compound) | repeated | We use the components field for pure substances and mixtures.

For example, suppose we are adding 3 mL of a 4 M solution of NaOH in water.

input { description: &#34;3 mL of 4M NaOH solution in water&#34; components: [ { identifiers: [ {type: IDENTIFIER_SMILES, value: &#34;O&#34;}, {type: IDENTIFIER_NAME, value: &#34;water&#34;} ] amount: { volume: {value: 3, units: MILLILITER} } } components: [ { identifiers: [ {type: IDENTIFIER_SMILES, value: &#34;[Na&#43;].[OH-]&#34;}, {type: IDENTIFIER_NAME, value: &#34;sodium hydroxide&#34;} ] amount { moles: {value: 12, units: MILLIMOLES} } } ] } |
| addition_order | [int32](#int32) |  | Used to define order of addition. ReactionInputs with the same addition_order were added simultaneously. One ReactionInput with a lower addition_order than another was added earlier in the procedure. This field is 1-indexed. |
| addition_time | [Time](#ord.Time) |  | When the addition event took place in terms of the reaction time (or, in the case of flow chemistry, the residence time). |
| addition_speed | [ReactionInput.AdditionSpeed.AdditionSpeedType](#ord.ReactionInput.AdditionSpeed.AdditionSpeedType) |  | The qualitative rate of addition. |
| addition_duration | [Time](#ord.Time) |  | Quantitatively, how long addition took |
| flow_rate | [FlowRate](#ord.FlowRate) |  | For continuous synthesis, we instead specify a flow rate. |






<a name="ord.ReactionInput.AdditionSpeed"></a>

### ReactionInput.AdditionSpeed







<a name="ord.ReactionNotes"></a>

### ReactionNotes



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| is_heterogeneous | [bool](#bool) |  | Equivalent to &#34;not single phase&#34;. |
| is_exothermic | [bool](#bool) |  | Qualitative exothermicity (primarily for safety). |
| is_offgasses | [bool](#bool) |  | Qualitative offgassing (primarily for safety). |
| is_sensitive_to_moisture | [bool](#bool) |  |  |
| is_sensitive_to_oxygen | [bool](#bool) |  |  |
| is_sensitive_to_light | [bool](#bool) |  |  |
| safety_notes | [string](#string) |  |  |
| procedure_details | [string](#string) |  | Overflow field for full procedure details |






<a name="ord.ReactionObservation"></a>

### ReactionObservation



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| time | [Time](#ord.Time) |  |  |
| comment | [string](#string) |  | e.g. what color is the reaction? |
| image | [Data](#ord.Data) |  |  |






<a name="ord.ReactionOutcome"></a>

### ReactionOutcome
The outcomes of a reaction describe the conversion, yield, and/or other
analyses of the resulting product mixture after workup step(s). Each
outcome is associated with a reaction/residence time. To allow for
one Reaction message to contain the results of a full kinetic profiling
experiment, this is a repeated field of the Reaction message.

It is the parent message for product characterization and any analytical
data.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| reaction_time | [Time](#ord.Time) |  | Reaction time (for flow, equivalent to residence time or spacetime). |
| conversion | [Percentage](#ord.Percentage) |  | Conversion with respect to the limiting reactant. |
| products | [ReactionProduct](#ord.ReactionProduct) | repeated |  |
| analyses | [ReactionOutcome.AnalysesEntry](#ord.ReactionOutcome.AnalysesEntry) | repeated | Analyses are stored in a map to associate each with a unique key. The key is cross-referenced in ReactionProduct messages to indicate which analyses were used to derive which performance values/metrics. The string used for the key carries no meaning outside of this cross-referencing. |






<a name="ord.ReactionOutcome.AnalysesEntry"></a>

### ReactionOutcome.AnalysesEntry



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| key | [string](#string) |  |  |
| value | [ReactionAnalysis](#ord.ReactionAnalysis) |  |  |






<a name="ord.ReactionProduct"></a>

### ReactionProduct



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| compound | [Compound](#ord.Compound) |  |  |
| is_desired_product | [bool](#bool) |  |  |
| compound_yield | [Percentage](#ord.Percentage) |  |  |
| purity | [Percentage](#ord.Percentage) |  |  |
| selectivity | [Selectivity](#ord.Selectivity) |  |  |
| analysis_identity | [string](#string) | repeated | Key(s) of the analysis used to confirm identity. |
| analysis_yield | [string](#string) | repeated | Key(s) of the analysis used to assess yield. |
| analysis_purity | [string](#string) | repeated | Key(s) of the analysis used to assess purity. |
| analysis_selectivity | [string](#string) | repeated | Key(s) of the analysis used to assess selectivity |
| isolated_color | [string](#string) |  | TODO(ccoley): How to allow specification of the state of matter of the purified compound? For example, &#34;___ was recovered as a white powder in x% yield (y.z mg)&#34;. Or oils, crystal texture, etc. This is only relevant for compounds that are isolated. TODO(kearnes): Should this be an Observation message? |
| texture | [ReactionProduct.Texture.TextureType](#ord.ReactionProduct.Texture.TextureType) |  |  |
| texture_details | [string](#string) |  |  |






<a name="ord.ReactionProduct.Texture"></a>

### ReactionProduct.Texture







<a name="ord.ReactionProvenance"></a>

### ReactionProvenance



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| experimenter | [Person](#ord.Person) |  |  |
| city | [string](#string) |  |  |
| experiment_start | [DateTime](#ord.DateTime) |  |  |
| doi | [string](#string) |  |  |
| patent | [string](#string) |  |  |
| publication_url | [string](#string) |  |  |
| record_created | [ReactionProvenance.RecordEvent](#ord.ReactionProvenance.RecordEvent) |  |  |
| record_modified | [ReactionProvenance.RecordEvent](#ord.ReactionProvenance.RecordEvent) | repeated |  |
| record_id | [string](#string) |  | This is a unique ID field that the centralized database will write to. |






<a name="ord.ReactionProvenance.RecordEvent"></a>

### ReactionProvenance.RecordEvent
Metadata for the public database.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| time | [DateTime](#ord.DateTime) |  |  |
| person | [Person](#ord.Person) |  |  |
| details | [string](#string) |  |  |






<a name="ord.ReactionSetup"></a>

### ReactionSetup



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| vessel | [Vessel](#ord.Vessel) |  |  |
| is_automated | [bool](#bool) |  | Specification of automated protocols. |
| automation_platform | [string](#string) |  | Automated platform name, brand, or model number. |
| automation_code | [ReactionSetup.AutomationCodeEntry](#ord.ReactionSetup.AutomationCodeEntry) | repeated | Raw automation code or synthetic recipe definition. |






<a name="ord.ReactionSetup.AutomationCodeEntry"></a>

### ReactionSetup.AutomationCodeEntry



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| key | [string](#string) |  |  |
| value | [Data](#ord.Data) |  |  |






<a name="ord.ReactionWorkup"></a>

### ReactionWorkup



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [ReactionWorkup.WorkupType](#ord.ReactionWorkup.WorkupType) |  |  |
| details | [string](#string) |  |  |
| duration | [Time](#ord.Time) |  |  |
| components | [Compound](#ord.Compound) | repeated |  |
| temperature | [TemperatureConditions](#ord.TemperatureConditions) |  |  |
| keep_phase | [string](#string) |  |  |
| stirring | [StirringConditions](#ord.StirringConditions) |  |  |
| target_ph | [float](#float) |  |  |






<a name="ord.Selectivity"></a>

### Selectivity



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [Selectivity.SelectivityType](#ord.Selectivity.SelectivityType) |  |  |
| details | [string](#string) |  |  |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | TODO(kearnes): What does precision mean in this context? |






<a name="ord.StirringConditions"></a>

### StirringConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [StirringConditions.StirringMethod.StirringMethodType](#ord.StirringConditions.StirringMethod.StirringMethodType) |  |  |
| details | [string](#string) |  |  |
| rate | [StirringConditions.StirringRate.StirringRateType](#ord.StirringConditions.StirringRate.StirringRateType) |  |  |
| rpm | [int32](#int32) |  |  |






<a name="ord.StirringConditions.StirringMethod"></a>

### StirringConditions.StirringMethod







<a name="ord.StirringConditions.StirringRate"></a>

### StirringConditions.StirringRate







<a name="ord.Temperature"></a>

### Temperature



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Temperature.TemperatureUnit](#ord.Temperature.TemperatureUnit) |  |  |






<a name="ord.TemperatureConditions"></a>

### TemperatureConditions



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [TemperatureConditions.TemperatureControl.TemperatureControlType](#ord.TemperatureConditions.TemperatureControl.TemperatureControlType) |  |  |
| details | [string](#string) |  |  |
| setpoint | [Temperature](#ord.Temperature) |  |  |
| measurements | [TemperatureConditions.Measurement](#ord.TemperatureConditions.Measurement) | repeated |  |






<a name="ord.TemperatureConditions.Measurement"></a>

### TemperatureConditions.Measurement



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [TemperatureConditions.Measurement.MeasurementType](#ord.TemperatureConditions.Measurement.MeasurementType) |  |  |
| details | [string](#string) |  |  |
| time | [Time](#ord.Time) |  |  |
| temperature | [Temperature](#ord.Temperature) |  |  |






<a name="ord.TemperatureConditions.TemperatureControl"></a>

### TemperatureConditions.TemperatureControl







<a name="ord.Time"></a>

### Time
To allow users to describe synthetic processes in whatever units they find
most natural, we define a fixed list of allowable units for each measurement
type. Upon submission to a centralized database, or using a validation and
canonicalization script, we will convert all values to the default units
(the first nonzero item in each enum).

Each message also contains a `precision` field, which specifies the precision
of the measurement in the same units as the measurement itself. Often the
precision will be the standard deviation from an instrument calibration.


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Time.TimeUnit](#ord.Time.TimeUnit) |  |  |






<a name="ord.Vessel"></a>

### Vessel



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| type | [Vessel.VesselType.VesselTypeEnum](#ord.Vessel.VesselType.VesselTypeEnum) |  |  |
| details | [string](#string) |  |  |
| material | [Vessel.VesselMaterial.VesselMaterialType](#ord.Vessel.VesselMaterial.VesselMaterialType) |  |  |
| material_details | [string](#string) |  |  |
| preparation | [Vessel.VesselPreparation.VesselPreparationType](#ord.Vessel.VesselPreparation.VesselPreparationType) |  |  |
| preparation_details | [string](#string) |  |  |
| volume | [Volume](#ord.Volume) |  | Size (volume) of the vessel. |






<a name="ord.Vessel.VesselMaterial"></a>

### Vessel.VesselMaterial







<a name="ord.Vessel.VesselPreparation"></a>

### Vessel.VesselPreparation







<a name="ord.Vessel.VesselType"></a>

### Vessel.VesselType







<a name="ord.Voltage"></a>

### Voltage



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Voltage.VoltageUnit](#ord.Voltage.VoltageUnit) |  |  |






<a name="ord.Volume"></a>

### Volume



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Volume.VolumeUnit](#ord.Volume.VolumeUnit) |  |  |






<a name="ord.Wavelength"></a>

### Wavelength



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| value | [float](#float) |  |  |
| precision | [float](#float) |  | Precision of the measurement (with the same units as `value`). |
| units | [Wavelength.WavelengthUnit](#ord.Wavelength.WavelengthUnit) |  |  |





 


<a name="ord.Compound.ReactionRole.ReactionRoleType"></a>

### Compound.ReactionRole.ReactionRoleType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| REACTANT | 1 | TODO(ccoley): Do we want to use the definition of a reactant aligned with Reaxys, or say that any species that contributes heavy atoms to a desired product is a reactant? This field might be kind of a throwaway anyway... |
| REAGENT | 2 |  |
| SOLVENT | 3 |  |
| CATALYST | 4 |  |
| WORKUP | 5 |  |



<a name="ord.CompoundIdentifier.IdentifierType"></a>

### CompoundIdentifier.IdentifierType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| SMILES | 2 | Simplified molecular-input line-entry system. |
| INCHI | 3 | IUPAC International Chemical Identifier. |
| MOLBLOCK | 4 | Molblock from a MDL Molfile V3000. |
| IUPAC_NAME | 5 | Chemical name following IUPAC nomenclature recommendations. |
| NAME | 6 | Any accepted common name, trade name, etc. |
| CAS_NUMBER | 7 | Chemical Abstracts Service Registry Number (with hyphens). |
| PUBCHEM_CID | 8 | PubChem Compound ID number. |
| CHEMSPIDER_ID | 9 | ChemSpider ID number. |
| CXSMILES | 10 | ChemAxon extended SMILES |
| INCHI_KEY | 11 | IUPAC International Chemical Identifier key |
| XYZ | 12 | XYZ molecule file |
| UNIPROT_ID | 13 | UniProt ID (for enzymes) |
| PDB_ID | 14 | Protein data bank ID (for enzymes) |
| RDKIT_BINARY | 15 | RDKit binary format (for fast loading) |



<a name="ord.CompoundPreparation.PreparationType"></a>

### CompoundPreparation.PreparationType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| NONE | 2 | Compound used as received. |
| REPURIFIED | 3 | Compound repurified (e.g., recrystallized). |
| SPARGED | 4 | Compound sparged, most likely to be the case with solvents. |
| DRIED | 5 | Moisture removed, e.g., using molecular sieves. |
| SYNTHESIZED | 6 | Compound synthesized in-house |



<a name="ord.Concentration.ConcentrationUnit"></a>

### Concentration.ConcentrationUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| MOLAR | 1 |  |
| MILLIMOLAR | 2 |  |
| MICROMOLAR | 3 |  |



<a name="ord.Current.CurrentUnit"></a>

### Current.CurrentUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| AMPERE | 1 |  |
| MILLIAMPERE | 2 |  |



<a name="ord.ElectrochemistryConditions.ElectrochemistryType"></a>

### ElectrochemistryConditions.ElectrochemistryType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| CONSTANT_CURRENT | 2 |  |
| CONSTANT_VOLTAGE | 3 |  |



<a name="ord.FlowConditions.FlowType"></a>

### FlowConditions.FlowType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| PLUG_FLOW_REACTOR | 2 |  |
| CONTINUOUS_STIRRED_TANK_REACTOR | 3 |  |
| PACKED_BED_REACTOR | 4 |  |



<a name="ord.FlowConditions.Tubing.TubingMaterialType"></a>

### FlowConditions.Tubing.TubingMaterialType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| STEEL | 2 |  |
| COPPER | 3 |  |
| PFA | 4 |  |
| FEP | 5 |  |
| TEFLONAF | 6 |  |
| PTFE | 7 |  |
| GLASS | 8 |  |
| QUARTZ | 9 |  |
| SILICON | 10 | e.g., a chip-based microreactor |
| PDMS | 11 |  |



<a name="ord.FlowRate.FlowRateUnit"></a>

### FlowRate.FlowRateUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| MICROLITER_PER_MINUTE | 1 |  |
| MICROLITER_PER_SECOND | 2 |  |
| MILLILITER_PER_MINUTE | 3 |  |
| MILLILITER_PER_SECOND | 4 |  |
| MICROLITER_PER_HOUR | 5 |  |



<a name="ord.IlluminationConditions.IlluminationType"></a>

### IlluminationConditions.IlluminationType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| AMBIENT | 2 |  |
| DARK | 3 |  |
| LED | 4 |  |
| HALOGEN_LAMP | 5 |  |
| DEUTERIUM_LAMP | 6 |  |
| SOLAR_SIMULATOR | 7 |  |
| BROAD_SPECTRUM | 8 |  |



<a name="ord.Length.LengthUnit"></a>

### Length.LengthUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CENTIMETER | 1 |  |
| MILLIMETER | 2 |  |
| METER | 3 |  |
| INCH | 4 |  |
| FOOT | 5 |  |



<a name="ord.Mass.MassUnit"></a>

### Mass.MassUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| GRAM | 1 |  |
| MILLIGRAM | 2 |  |
| MICROGRAM | 3 |  |
| KILOGRAM | 4 |  |



<a name="ord.Moles.MolesUnit"></a>

### Moles.MolesUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| MOLES | 1 |  |
| MILLIMOLES | 2 |  |
| MICROMOLES | 3 |  |
| NANOMOLES | 4 |  |



<a name="ord.Pressure.PressureUnit"></a>

### Pressure.PressureUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| BAR | 1 |  |
| ATMOSPHERE | 2 |  |
| PSI | 3 |  |
| KPSI | 4 |  |
| PASCAL | 5 |  |
| KILOPASCAL | 6 |  |



<a name="ord.PressureConditions.Atmosphere.AtmosphereType"></a>

### PressureConditions.Atmosphere.AtmosphereType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| AIR | 2 |  |
| NITROGEN | 3 |  |
| ARGON | 4 |  |
| OXYGEN | 5 |  |
| HYDROGEN | 6 |  |



<a name="ord.PressureConditions.Measurement.MeasurementType"></a>

### PressureConditions.Measurement.MeasurementType
TODO(ccoley) get input on how to expand this enum, among others

| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| PRESSURE_TRANSDUCER | 2 |  |



<a name="ord.PressureConditions.PressureControl.PressureControlType"></a>

### PressureConditions.PressureControl.PressureControlType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| AMBIENT | 2 |  |
| BALLOON | 3 |  |
| SEALED | 4 | Fully sealed vessel (e.g., microwave vial). |
| SEPTUM_WITH_NEEDLE | 5 | Slight positive pressure maintained |
| RELEASEVALVE | 6 |  |
| BPR | 7 | Back pressure regulator, as used in flow synthesis. |



<a name="ord.ReactionAnalysis.AnalysisType"></a>

### ReactionAnalysis.AnalysisType
TODO(ccoley): Solicit more feedback from experimentalists

| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| LC | 2 | Liquid chromatography. |
| GC | 3 | Gas chromatography. |
| IR | 4 | Infrared spectroscopy. |
| NMR | 5 | NMR spectroscopy. |
| MP | 6 | Melting point characterization. |
| UV | 7 | Ultraviolet spectroscopy. |
| TLC | 8 | Thin-layer chromatography. |
| MS | 9 | Mass spectrometry. |
| HRMS | 10 | High resolution mass spectrometry. |
| MSMS | 11 | Two-dimensional mass spectrometry. |
| WEIGHT | 12 | Weight of an isolated compound. |
| LCMS | 13 | Combined LC/MS. |
| GCMS | 14 | Combined GC/MS. |
| ELSD | 15 | Evaporative light scattering detector. |
| CD | 16 | Circular Dichroism. |
| SFC | 17 | Supercritical fluid chromatography. |



<a name="ord.ReactionIdentifier.IdentifierType"></a>

### ReactionIdentifier.IdentifierType
Possible identifier types are listed in an enum for extensibility

| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| REACTION_SMILES | 2 |  |
| ATOM_MAPPED_SMILES | 3 |  |
| RINCHI | 4 | Reaction InChI. |
| NAME | 5 | Named reaction or reaction category. |
| RDKIT_BINARY | 6 | RDKit binary format (for fast loading). |



<a name="ord.ReactionInput.AdditionSpeed.AdditionSpeedType"></a>

### ReactionInput.AdditionSpeed.AdditionSpeedType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 | Unspecified. |
| ALL_AT_ONCE | 1 |  |
| FAST | 2 |  |
| SLOW | 3 |  |
| DROPWISE | 4 |  |
| CONTINUOUS | 5 |  |



<a name="ord.ReactionProduct.Texture.TextureType"></a>

### ReactionProduct.Texture.TextureType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| POWDER | 2 |  |
| CRYSTAL | 3 |  |
| OIL | 4 |  |



<a name="ord.ReactionWorkup.WorkupType"></a>

### ReactionWorkup.WorkupType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| ADDITION | 2 | Addition (quench, dilution, extraction solvent, etc.) Specify composition/amount in &#34;components&#34;. |
| TEMPERATURE | 3 | Change of temperature. Specify conditions in &#34;temperature&#34;. |
| CONCENTRATION | 4 | Concentration step, often using a rotovap. |
| EXTRACTION | 5 | Liquid extractions are often preceded by Additions. If there are multiple distinct additions prior to an extraction, it is assumed that the kept phases are pooled. Specify which phase to keep in &#34;keep_phase&#34;. |
| FILTRATION | 6 | Filtration (can keep solid or filtrate). Specify which phase to keep in &#34;keep phase&#34;. |
| WASH | 7 | Washing a solid or liquid, keeping the original phase. Specify &#34;components&#34; of rinse. Rinses performed in multiple stages should be given multiple workup steps |
| DRY_IN_VACUUM | 8 | Dried under vacuum. |
| DRY_WITH_MATERIAL | 9 | Dried with chemical additive. Specify chemical additive in &#34;components&#34;. |
| FLASH_CHROMATOGRAPHY | 10 | Purification by flash chromatography. |
| OTHER_CHROMATOGRAPHY | 11 | Purification by other prep chromatography. |
| SCAVENGING | 12 | Scavenging step (e.g., pass through alumina pad) Specify any material additives in &#34;components&#34;. |
| WAIT | 13 | Waiting step. Specify &#34;duration&#34;. |
| STIRRING | 14 | Mixing step. Specify &#34;stirring&#34; |
| CRYSTALLIZATION | 15 |  |
| PH_ADJUST | 16 | pH adjustments should specify &#34;components&#34; to define species used as well as &#34;ph&#34; for target ph |
| DISSOLUTION | 17 | Redissolution considered to be a special form of addition. Specify &#34;components&#34; |



<a name="ord.Selectivity.SelectivityType"></a>

### Selectivity.SelectivityType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| EE | 2 | Enantiomeric excess as a percentage. |
| ER | 3 | Enantiomeric ratio. (x:1) |
| DE | 4 | Diasteromeric ratio (x:1) |



<a name="ord.StirringConditions.StirringMethod.StirringMethodType"></a>

### StirringConditions.StirringMethod.StirringMethodType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| NONE | 2 |  |
| STIR_BAR | 3 |  |
| OVERHEAD_MIXER | 4 |  |
| AGITATION | 5 |  |



<a name="ord.StirringConditions.StirringRate.StirringRateType"></a>

### StirringConditions.StirringRate.StirringRateType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| HIGH | 1 |  |
| MEDIUM | 2 |  |
| LOW | 3 |  |



<a name="ord.Temperature.TemperatureUnit"></a>

### Temperature.TemperatureUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CELSIUS | 1 |  |
| FAHRENHEIT | 2 |  |
| KELVIN | 3 |  |



<a name="ord.TemperatureConditions.Measurement.MeasurementType"></a>

### TemperatureConditions.Measurement.MeasurementType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| THERMOCOUPLE_INTERNAL | 2 | Physically in reaction solution. |
| THERMOCOUPLE_EXTERNAL | 3 | On outside of vessel or, e.g., in oil bath. |
| INFRARED | 4 | Contactless infrared probe. |



<a name="ord.TemperatureConditions.TemperatureControl.TemperatureControlType"></a>

### TemperatureConditions.TemperatureControl.TemperatureControlType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| AMBIENT | 2 |  |
| OIL_BATH | 3 |  |
| WATER_BATH | 4 |  |
| SAND_BATH | 5 |  |
| ICE_BATH | 6 |  |
| DRY_ALUMINUM_PLATE | 7 |  |
| MICROWAVE | 8 |  |
| DRY_ICE_BATH | 9 |  |
| AIR_FAN | 10 |  |
| LIQUID_NITROGEN | 11 |  |



<a name="ord.Time.TimeUnit"></a>

### Time.TimeUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| HOUR | 1 |  |
| MINUTE | 2 |  |
| SECOND | 3 |  |



<a name="ord.Vessel.VesselMaterial.VesselMaterialType"></a>

### Vessel.VesselMaterial.VesselMaterialType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| GLASS | 2 |  |
| POLYPROPYLENE | 3 |  |
| PLASTIC | 4 |  |



<a name="ord.Vessel.VesselPreparation.VesselPreparationType"></a>

### Vessel.VesselPreparation.VesselPreparationType


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| NONE | 2 |  |
| OVEN_DRIED | 3 |  |



<a name="ord.Vessel.VesselType.VesselTypeEnum"></a>

### Vessel.VesselType.VesselTypeEnum


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| CUSTOM | 1 |  |
| ROUND_BOTTOM_FLASK | 2 |  |
| VIAL | 3 |  |
| WELL_PLATE | 4 |  |
| MICROWAVE_VIAL | 5 |  |
| TUBE | 6 |  |
| CONTINUOUS_STIRRED_TANK_REACTOR | 7 |  |
| PACKED_BED_REACTOR | 8 |  |



<a name="ord.Voltage.VoltageUnit"></a>

### Voltage.VoltageUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| VOLT | 1 |  |
| MILLIVOLT | 2 |  |



<a name="ord.Volume.VolumeUnit"></a>

### Volume.VolumeUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| MILLILITER | 1 |  |
| MICROLITER | 2 |  |
| LITER | 3 |  |



<a name="ord.Wavelength.WavelengthUnit"></a>

### Wavelength.WavelengthUnit


| Name | Number | Description |
| ---- | ------ | ----------- |
| UNSPECIFIED | 0 |  |
| NANOMETER | 1 |  |
| WAVENUMBER | 2 | cm^{-1} |


 

 

 



## Scalar Value Types

| .proto Type | Notes | C++ | Java | Python | Go | C# | PHP | Ruby |
| ----------- | ----- | --- | ---- | ------ | -- | -- | --- | ---- |
| <a name="double" /> double |  | double | double | float | float64 | double | float | Float |
| <a name="float" /> float |  | float | float | float | float32 | float | float | Float |
| <a name="int32" /> int32 | Uses variable-length encoding. Inefficient for encoding negative numbers – if your field is likely to have negative values, use sint32 instead. | int32 | int | int | int32 | int | integer | Bignum or Fixnum (as required) |
| <a name="int64" /> int64 | Uses variable-length encoding. Inefficient for encoding negative numbers – if your field is likely to have negative values, use sint64 instead. | int64 | long | int/long | int64 | long | integer/string | Bignum |
| <a name="uint32" /> uint32 | Uses variable-length encoding. | uint32 | int | int/long | uint32 | uint | integer | Bignum or Fixnum (as required) |
| <a name="uint64" /> uint64 | Uses variable-length encoding. | uint64 | long | int/long | uint64 | ulong | integer/string | Bignum or Fixnum (as required) |
| <a name="sint32" /> sint32 | Uses variable-length encoding. Signed int value. These more efficiently encode negative numbers than regular int32s. | int32 | int | int | int32 | int | integer | Bignum or Fixnum (as required) |
| <a name="sint64" /> sint64 | Uses variable-length encoding. Signed int value. These more efficiently encode negative numbers than regular int64s. | int64 | long | int/long | int64 | long | integer/string | Bignum |
| <a name="fixed32" /> fixed32 | Always four bytes. More efficient than uint32 if values are often greater than 2^28. | uint32 | int | int | uint32 | uint | integer | Bignum or Fixnum (as required) |
| <a name="fixed64" /> fixed64 | Always eight bytes. More efficient than uint64 if values are often greater than 2^56. | uint64 | long | int/long | uint64 | ulong | integer/string | Bignum |
| <a name="sfixed32" /> sfixed32 | Always four bytes. | int32 | int | int | int32 | int | integer | Bignum or Fixnum (as required) |
| <a name="sfixed64" /> sfixed64 | Always eight bytes. | int64 | long | int/long | int64 | long | integer/string | Bignum |
| <a name="bool" /> bool |  | bool | boolean | boolean | bool | bool | boolean | TrueClass/FalseClass |
| <a name="string" /> string | A string must always contain UTF-8 encoded or 7-bit ASCII text. | string | String | str/unicode | string | string | string | String (UTF-8) |
| <a name="bytes" /> bytes | May contain any arbitrary sequence of bytes. | string | ByteString | str | []byte | ByteString | string | String (ASCII-8BIT) |

