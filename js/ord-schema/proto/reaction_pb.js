/**
 * Copyright 2023 Open Reaction Database Project Authors
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

// source: ord-schema/proto/reaction.proto
/**
 * @fileoverview
 * @enhanceable
 * @suppress {missingRequire} reports error on implicit type usages.
 * @suppress {messageConventions} JS Compiler reports an error if a variable or
 *     field starts with 'MSG_' and isn't a translatable message.
 * @public
 */
// GENERATED CODE -- DO NOT EDIT!
/* eslint-disable */
// @ts-nocheck

var jspb = require('google-protobuf');
var goog = jspb;
var global =
    (typeof globalThis !== 'undefined' && globalThis) ||
    (typeof window !== 'undefined' && window) ||
    (typeof global !== 'undefined' && global) ||
    (typeof self !== 'undefined' && self) ||
    (function () { return this; }).call(null) ||
    Function('return this')();

goog.exportSymbol('proto.ord.Amount', null, global);
goog.exportSymbol('proto.ord.Amount.KindCase', null, global);
goog.exportSymbol('proto.ord.Analysis', null, global);
goog.exportSymbol('proto.ord.Analysis.AnalysisType', null, global);
goog.exportSymbol('proto.ord.Compound', null, global);
goog.exportSymbol('proto.ord.Compound.Source', null, global);
goog.exportSymbol('proto.ord.CompoundIdentifier', null, global);
goog.exportSymbol('proto.ord.CompoundIdentifier.CompoundIdentifierType', null, global);
goog.exportSymbol('proto.ord.CompoundPreparation', null, global);
goog.exportSymbol('proto.ord.CompoundPreparation.CompoundPreparationType', null, global);
goog.exportSymbol('proto.ord.Concentration', null, global);
goog.exportSymbol('proto.ord.Concentration.ConcentrationUnit', null, global);
goog.exportSymbol('proto.ord.CrudeComponent', null, global);
goog.exportSymbol('proto.ord.Current', null, global);
goog.exportSymbol('proto.ord.Current.CurrentUnit', null, global);
goog.exportSymbol('proto.ord.Data', null, global);
goog.exportSymbol('proto.ord.Data.KindCase', null, global);
goog.exportSymbol('proto.ord.DateTime', null, global);
goog.exportSymbol('proto.ord.ElectrochemistryConditions', null, global);
goog.exportSymbol('proto.ord.ElectrochemistryConditions.ElectrochemistryCell', null, global);
goog.exportSymbol('proto.ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType', null, global);
goog.exportSymbol('proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement', null, global);
goog.exportSymbol('proto.ord.ElectrochemistryConditions.ElectrochemistryType', null, global);
goog.exportSymbol('proto.ord.FloatValue', null, global);
goog.exportSymbol('proto.ord.FlowConditions', null, global);
goog.exportSymbol('proto.ord.FlowConditions.FlowType', null, global);
goog.exportSymbol('proto.ord.FlowConditions.Tubing', null, global);
goog.exportSymbol('proto.ord.FlowConditions.Tubing.TubingType', null, global);
goog.exportSymbol('proto.ord.FlowRate', null, global);
goog.exportSymbol('proto.ord.FlowRate.FlowRateUnit', null, global);
goog.exportSymbol('proto.ord.IlluminationConditions', null, global);
goog.exportSymbol('proto.ord.IlluminationConditions.IlluminationType', null, global);
goog.exportSymbol('proto.ord.Length', null, global);
goog.exportSymbol('proto.ord.Length.LengthUnit', null, global);
goog.exportSymbol('proto.ord.Mass', null, global);
goog.exportSymbol('proto.ord.Mass.MassUnit', null, global);
goog.exportSymbol('proto.ord.Moles', null, global);
goog.exportSymbol('proto.ord.Moles.MolesUnit', null, global);
goog.exportSymbol('proto.ord.Percentage', null, global);
goog.exportSymbol('proto.ord.Person', null, global);
goog.exportSymbol('proto.ord.Pressure', null, global);
goog.exportSymbol('proto.ord.Pressure.PressureUnit', null, global);
goog.exportSymbol('proto.ord.PressureConditions', null, global);
goog.exportSymbol('proto.ord.PressureConditions.Atmosphere', null, global);
goog.exportSymbol('proto.ord.PressureConditions.Atmosphere.AtmosphereType', null, global);
goog.exportSymbol('proto.ord.PressureConditions.PressureControl', null, global);
goog.exportSymbol('proto.ord.PressureConditions.PressureControl.PressureControlType', null, global);
goog.exportSymbol('proto.ord.PressureConditions.PressureMeasurement', null, global);
goog.exportSymbol('proto.ord.PressureConditions.PressureMeasurement.PressureMeasurementType', null, global);
goog.exportSymbol('proto.ord.ProductCompound', null, global);
goog.exportSymbol('proto.ord.ProductCompound.Texture', null, global);
goog.exportSymbol('proto.ord.ProductCompound.Texture.TextureType', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement.MassSpecMeasurementDetails', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement.ProductMeasurementType', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement.Selectivity', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement.Selectivity.SelectivityType', null, global);
goog.exportSymbol('proto.ord.ProductMeasurement.ValueCase', null, global);
goog.exportSymbol('proto.ord.Reaction', null, global);
goog.exportSymbol('proto.ord.ReactionConditions', null, global);
goog.exportSymbol('proto.ord.ReactionIdentifier', null, global);
goog.exportSymbol('proto.ord.ReactionIdentifier.ReactionIdentifierType', null, global);
goog.exportSymbol('proto.ord.ReactionInput', null, global);
goog.exportSymbol('proto.ord.ReactionInput.AdditionDevice', null, global);
goog.exportSymbol('proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType', null, global);
goog.exportSymbol('proto.ord.ReactionInput.AdditionSpeed', null, global);
goog.exportSymbol('proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType', null, global);
goog.exportSymbol('proto.ord.ReactionNotes', null, global);
goog.exportSymbol('proto.ord.ReactionObservation', null, global);
goog.exportSymbol('proto.ord.ReactionOutcome', null, global);
goog.exportSymbol('proto.ord.ReactionProvenance', null, global);
goog.exportSymbol('proto.ord.ReactionRole', null, global);
goog.exportSymbol('proto.ord.ReactionRole.ReactionRoleType', null, global);
goog.exportSymbol('proto.ord.ReactionSetup', null, global);
goog.exportSymbol('proto.ord.ReactionSetup.ReactionEnvironment', null, global);
goog.exportSymbol('proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType', null, global);
goog.exportSymbol('proto.ord.ReactionWorkup', null, global);
goog.exportSymbol('proto.ord.ReactionWorkup.ReactionWorkupType', null, global);
goog.exportSymbol('proto.ord.RecordEvent', null, global);
goog.exportSymbol('proto.ord.StirringConditions', null, global);
goog.exportSymbol('proto.ord.StirringConditions.StirringMethodType', null, global);
goog.exportSymbol('proto.ord.StirringConditions.StirringRate', null, global);
goog.exportSymbol('proto.ord.StirringConditions.StirringRate.StirringRateType', null, global);
goog.exportSymbol('proto.ord.Temperature', null, global);
goog.exportSymbol('proto.ord.Temperature.TemperatureUnit', null, global);
goog.exportSymbol('proto.ord.TemperatureConditions', null, global);
goog.exportSymbol('proto.ord.TemperatureConditions.TemperatureControl', null, global);
goog.exportSymbol('proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType', null, global);
goog.exportSymbol('proto.ord.TemperatureConditions.TemperatureMeasurement', null, global);
goog.exportSymbol('proto.ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType', null, global);
goog.exportSymbol('proto.ord.Time', null, global);
goog.exportSymbol('proto.ord.Time.TimeUnit', null, global);
goog.exportSymbol('proto.ord.UnmeasuredAmount', null, global);
goog.exportSymbol('proto.ord.UnmeasuredAmount.UnmeasuredAmountType', null, global);
goog.exportSymbol('proto.ord.Vessel', null, global);
goog.exportSymbol('proto.ord.Vessel.VesselType', null, global);
goog.exportSymbol('proto.ord.VesselAttachment', null, global);
goog.exportSymbol('proto.ord.VesselAttachment.VesselAttachmentType', null, global);
goog.exportSymbol('proto.ord.VesselMaterial', null, global);
goog.exportSymbol('proto.ord.VesselMaterial.VesselMaterialType', null, global);
goog.exportSymbol('proto.ord.VesselPreparation', null, global);
goog.exportSymbol('proto.ord.VesselPreparation.VesselPreparationType', null, global);
goog.exportSymbol('proto.ord.Voltage', null, global);
goog.exportSymbol('proto.ord.Voltage.VoltageUnit', null, global);
goog.exportSymbol('proto.ord.Volume', null, global);
goog.exportSymbol('proto.ord.Volume.VolumeUnit', null, global);
goog.exportSymbol('proto.ord.Wavelength', null, global);
goog.exportSymbol('proto.ord.Wavelength.WavelengthUnit', null, global);
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Reaction = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.Reaction.repeatedFields_, null);
};
goog.inherits(proto.ord.Reaction, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Reaction.displayName = 'proto.ord.Reaction';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionIdentifier = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionIdentifier, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionIdentifier.displayName = 'proto.ord.ReactionIdentifier';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionInput = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.ReactionInput.repeatedFields_, null);
};
goog.inherits(proto.ord.ReactionInput, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionInput.displayName = 'proto.ord.ReactionInput';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionInput.AdditionSpeed = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionInput.AdditionSpeed, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionInput.AdditionSpeed.displayName = 'proto.ord.ReactionInput.AdditionSpeed';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionInput.AdditionDevice = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionInput.AdditionDevice, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionInput.AdditionDevice.displayName = 'proto.ord.ReactionInput.AdditionDevice';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Amount = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, proto.ord.Amount.oneofGroups_);
};
goog.inherits(proto.ord.Amount, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Amount.displayName = 'proto.ord.Amount';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.UnmeasuredAmount = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.UnmeasuredAmount, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.UnmeasuredAmount.displayName = 'proto.ord.UnmeasuredAmount';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.CrudeComponent = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.CrudeComponent, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.CrudeComponent.displayName = 'proto.ord.CrudeComponent';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Compound = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.Compound.repeatedFields_, null);
};
goog.inherits(proto.ord.Compound, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Compound.displayName = 'proto.ord.Compound';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Compound.Source = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Compound.Source, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Compound.Source.displayName = 'proto.ord.Compound.Source';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionRole = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionRole, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionRole.displayName = 'proto.ord.ReactionRole';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.CompoundPreparation = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.CompoundPreparation, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.CompoundPreparation.displayName = 'proto.ord.CompoundPreparation';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.CompoundIdentifier = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.CompoundIdentifier, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.CompoundIdentifier.displayName = 'proto.ord.CompoundIdentifier';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Vessel = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.Vessel.repeatedFields_, null);
};
goog.inherits(proto.ord.Vessel, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Vessel.displayName = 'proto.ord.Vessel';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.VesselMaterial = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.VesselMaterial, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.VesselMaterial.displayName = 'proto.ord.VesselMaterial';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.VesselAttachment = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.VesselAttachment, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.VesselAttachment.displayName = 'proto.ord.VesselAttachment';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.VesselPreparation = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.VesselPreparation, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.VesselPreparation.displayName = 'proto.ord.VesselPreparation';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionSetup = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionSetup, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionSetup.displayName = 'proto.ord.ReactionSetup';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionSetup.ReactionEnvironment = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionSetup.ReactionEnvironment, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionSetup.ReactionEnvironment.displayName = 'proto.ord.ReactionSetup.ReactionEnvironment';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionConditions.displayName = 'proto.ord.ReactionConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.TemperatureConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.TemperatureConditions.repeatedFields_, null);
};
goog.inherits(proto.ord.TemperatureConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.TemperatureConditions.displayName = 'proto.ord.TemperatureConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.TemperatureConditions.TemperatureControl = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.TemperatureConditions.TemperatureControl, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.TemperatureConditions.TemperatureControl.displayName = 'proto.ord.TemperatureConditions.TemperatureControl';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.TemperatureConditions.TemperatureMeasurement = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.TemperatureConditions.TemperatureMeasurement, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.TemperatureConditions.TemperatureMeasurement.displayName = 'proto.ord.TemperatureConditions.TemperatureMeasurement';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.PressureConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.PressureConditions.repeatedFields_, null);
};
goog.inherits(proto.ord.PressureConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.PressureConditions.displayName = 'proto.ord.PressureConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.PressureConditions.PressureControl = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.PressureConditions.PressureControl, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.PressureConditions.PressureControl.displayName = 'proto.ord.PressureConditions.PressureControl';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.PressureConditions.Atmosphere = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.PressureConditions.Atmosphere, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.PressureConditions.Atmosphere.displayName = 'proto.ord.PressureConditions.Atmosphere';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.PressureConditions.PressureMeasurement = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.PressureConditions.PressureMeasurement, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.PressureConditions.PressureMeasurement.displayName = 'proto.ord.PressureConditions.PressureMeasurement';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.StirringConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.StirringConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.StirringConditions.displayName = 'proto.ord.StirringConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.StirringConditions.StirringRate = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.StirringConditions.StirringRate, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.StirringConditions.StirringRate.displayName = 'proto.ord.StirringConditions.StirringRate';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.IlluminationConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.IlluminationConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.IlluminationConditions.displayName = 'proto.ord.IlluminationConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ElectrochemistryConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.ElectrochemistryConditions.repeatedFields_, null);
};
goog.inherits(proto.ord.ElectrochemistryConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ElectrochemistryConditions.displayName = 'proto.ord.ElectrochemistryConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.displayName = 'proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ElectrochemistryConditions.ElectrochemistryCell, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ElectrochemistryConditions.ElectrochemistryCell.displayName = 'proto.ord.ElectrochemistryConditions.ElectrochemistryCell';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.FlowConditions = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.FlowConditions, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.FlowConditions.displayName = 'proto.ord.FlowConditions';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.FlowConditions.Tubing = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.FlowConditions.Tubing, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.FlowConditions.Tubing.displayName = 'proto.ord.FlowConditions.Tubing';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionNotes = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionNotes, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionNotes.displayName = 'proto.ord.ReactionNotes';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionObservation = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionObservation, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionObservation.displayName = 'proto.ord.ReactionObservation';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionWorkup = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ReactionWorkup, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionWorkup.displayName = 'proto.ord.ReactionWorkup';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionOutcome = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.ReactionOutcome.repeatedFields_, null);
};
goog.inherits(proto.ord.ReactionOutcome, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionOutcome.displayName = 'proto.ord.ReactionOutcome';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ProductCompound = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.ProductCompound.repeatedFields_, null);
};
goog.inherits(proto.ord.ProductCompound, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ProductCompound.displayName = 'proto.ord.ProductCompound';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ProductCompound.Texture = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ProductCompound.Texture, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ProductCompound.Texture.displayName = 'proto.ord.ProductCompound.Texture';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ProductMeasurement = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, proto.ord.ProductMeasurement.oneofGroups_);
};
goog.inherits(proto.ord.ProductMeasurement, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ProductMeasurement.displayName = 'proto.ord.ProductMeasurement';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.ProductMeasurement.MassSpecMeasurementDetails.repeatedFields_, null);
};
goog.inherits(proto.ord.ProductMeasurement.MassSpecMeasurementDetails, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ProductMeasurement.MassSpecMeasurementDetails.displayName = 'proto.ord.ProductMeasurement.MassSpecMeasurementDetails';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ProductMeasurement.Selectivity = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.ProductMeasurement.Selectivity, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ProductMeasurement.Selectivity.displayName = 'proto.ord.ProductMeasurement.Selectivity';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.DateTime = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.DateTime, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.DateTime.displayName = 'proto.ord.DateTime';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Analysis = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Analysis, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Analysis.displayName = 'proto.ord.Analysis';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.ReactionProvenance = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, proto.ord.ReactionProvenance.repeatedFields_, null);
};
goog.inherits(proto.ord.ReactionProvenance, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.ReactionProvenance.displayName = 'proto.ord.ReactionProvenance';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Person = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Person, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Person.displayName = 'proto.ord.Person';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.RecordEvent = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.RecordEvent, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.RecordEvent.displayName = 'proto.ord.RecordEvent';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Time = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Time, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Time.displayName = 'proto.ord.Time';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Mass = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Mass, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Mass.displayName = 'proto.ord.Mass';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Moles = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Moles, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Moles.displayName = 'proto.ord.Moles';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Volume = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Volume, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Volume.displayName = 'proto.ord.Volume';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Concentration = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Concentration, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Concentration.displayName = 'proto.ord.Concentration';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Pressure = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Pressure, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Pressure.displayName = 'proto.ord.Pressure';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Temperature = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Temperature, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Temperature.displayName = 'proto.ord.Temperature';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Current = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Current, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Current.displayName = 'proto.ord.Current';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Voltage = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Voltage, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Voltage.displayName = 'proto.ord.Voltage';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Length = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Length, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Length.displayName = 'proto.ord.Length';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Wavelength = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Wavelength, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Wavelength.displayName = 'proto.ord.Wavelength';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.FlowRate = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.FlowRate, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.FlowRate.displayName = 'proto.ord.FlowRate';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Percentage = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.Percentage, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Percentage.displayName = 'proto.ord.Percentage';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.FloatValue = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, null);
};
goog.inherits(proto.ord.FloatValue, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.FloatValue.displayName = 'proto.ord.FloatValue';
}
/**
 * Generated by JsPbCodeGenerator.
 * @param {Array=} opt_data Optional initial data array, typically from a
 * server response, or constructed directly in Javascript. The array is used
 * in place and becomes part of the constructed object. It is not cloned.
 * If no data is provided, the constructed object will be empty, but still
 * valid.
 * @extends {jspb.Message}
 * @constructor
 */
proto.ord.Data = function(opt_data) {
  jspb.Message.initialize(this, opt_data, 0, -1, null, proto.ord.Data.oneofGroups_);
};
goog.inherits(proto.ord.Data, jspb.Message);
if (goog.DEBUG && !COMPILED) {
  /**
   * @public
   * @override
   */
  proto.ord.Data.displayName = 'proto.ord.Data';
}

/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.Reaction.repeatedFields_ = [1,6,7,8];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Reaction.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Reaction.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Reaction} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Reaction.toObject = function(includeInstance, msg) {
  var f, obj = {
    identifiersList: jspb.Message.toObjectList(msg.getIdentifiersList(),
    proto.ord.ReactionIdentifier.toObject, includeInstance),
    inputsMap: (f = msg.getInputsMap()) ? f.toObject(includeInstance, proto.ord.ReactionInput.toObject) : [],
    setup: (f = msg.getSetup()) && proto.ord.ReactionSetup.toObject(includeInstance, f),
    conditions: (f = msg.getConditions()) && proto.ord.ReactionConditions.toObject(includeInstance, f),
    notes: (f = msg.getNotes()) && proto.ord.ReactionNotes.toObject(includeInstance, f),
    observationsList: jspb.Message.toObjectList(msg.getObservationsList(),
    proto.ord.ReactionObservation.toObject, includeInstance),
    workupsList: jspb.Message.toObjectList(msg.getWorkupsList(),
    proto.ord.ReactionWorkup.toObject, includeInstance),
    outcomesList: jspb.Message.toObjectList(msg.getOutcomesList(),
    proto.ord.ReactionOutcome.toObject, includeInstance),
    provenance: (f = msg.getProvenance()) && proto.ord.ReactionProvenance.toObject(includeInstance, f),
    reactionId: jspb.Message.getFieldWithDefault(msg, 10, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Reaction}
 */
proto.ord.Reaction.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Reaction;
  return proto.ord.Reaction.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Reaction} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Reaction}
 */
proto.ord.Reaction.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.ReactionIdentifier;
      reader.readMessage(value,proto.ord.ReactionIdentifier.deserializeBinaryFromReader);
      msg.addIdentifiers(value);
      break;
    case 2:
      var value = msg.getInputsMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.ReactionInput.deserializeBinaryFromReader, "", new proto.ord.ReactionInput());
         });
      break;
    case 3:
      var value = new proto.ord.ReactionSetup;
      reader.readMessage(value,proto.ord.ReactionSetup.deserializeBinaryFromReader);
      msg.setSetup(value);
      break;
    case 4:
      var value = new proto.ord.ReactionConditions;
      reader.readMessage(value,proto.ord.ReactionConditions.deserializeBinaryFromReader);
      msg.setConditions(value);
      break;
    case 5:
      var value = new proto.ord.ReactionNotes;
      reader.readMessage(value,proto.ord.ReactionNotes.deserializeBinaryFromReader);
      msg.setNotes(value);
      break;
    case 6:
      var value = new proto.ord.ReactionObservation;
      reader.readMessage(value,proto.ord.ReactionObservation.deserializeBinaryFromReader);
      msg.addObservations(value);
      break;
    case 7:
      var value = new proto.ord.ReactionWorkup;
      reader.readMessage(value,proto.ord.ReactionWorkup.deserializeBinaryFromReader);
      msg.addWorkups(value);
      break;
    case 8:
      var value = new proto.ord.ReactionOutcome;
      reader.readMessage(value,proto.ord.ReactionOutcome.deserializeBinaryFromReader);
      msg.addOutcomes(value);
      break;
    case 9:
      var value = new proto.ord.ReactionProvenance;
      reader.readMessage(value,proto.ord.ReactionProvenance.deserializeBinaryFromReader);
      msg.setProvenance(value);
      break;
    case 10:
      var value = /** @type {string} */ (reader.readString());
      msg.setReactionId(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Reaction.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Reaction.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Reaction} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Reaction.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getIdentifiersList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      1,
      f,
      proto.ord.ReactionIdentifier.serializeBinaryToWriter
    );
  }
  f = message.getInputsMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(2, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.ReactionInput.serializeBinaryToWriter);
  }
  f = message.getSetup();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.ReactionSetup.serializeBinaryToWriter
    );
  }
  f = message.getConditions();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.ReactionConditions.serializeBinaryToWriter
    );
  }
  f = message.getNotes();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.ReactionNotes.serializeBinaryToWriter
    );
  }
  f = message.getObservationsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      6,
      f,
      proto.ord.ReactionObservation.serializeBinaryToWriter
    );
  }
  f = message.getWorkupsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      7,
      f,
      proto.ord.ReactionWorkup.serializeBinaryToWriter
    );
  }
  f = message.getOutcomesList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      8,
      f,
      proto.ord.ReactionOutcome.serializeBinaryToWriter
    );
  }
  f = message.getProvenance();
  if (f != null) {
    writer.writeMessage(
      9,
      f,
      proto.ord.ReactionProvenance.serializeBinaryToWriter
    );
  }
  f = message.getReactionId();
  if (f.length > 0) {
    writer.writeString(
      10,
      f
    );
  }
};


/**
 * repeated ReactionIdentifier identifiers = 1;
 * @return {!Array<!proto.ord.ReactionIdentifier>}
 */
proto.ord.Reaction.prototype.getIdentifiersList = function() {
  return /** @type{!Array<!proto.ord.ReactionIdentifier>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ReactionIdentifier, 1));
};


/**
 * @param {!Array<!proto.ord.ReactionIdentifier>} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setIdentifiersList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 1, value);
};


/**
 * @param {!proto.ord.ReactionIdentifier=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ReactionIdentifier}
 */
proto.ord.Reaction.prototype.addIdentifiers = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 1, opt_value, proto.ord.ReactionIdentifier, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearIdentifiersList = function() {
  return this.setIdentifiersList([]);
};


/**
 * map<string, ReactionInput> inputs = 2;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.ReactionInput>}
 */
proto.ord.Reaction.prototype.getInputsMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.ReactionInput>} */ (
      jspb.Message.getMapField(this, 2, opt_noLazyCreate,
      proto.ord.ReactionInput));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearInputsMap = function() {
  this.getInputsMap().clear();
  return this;
};


/**
 * optional ReactionSetup setup = 3;
 * @return {?proto.ord.ReactionSetup}
 */
proto.ord.Reaction.prototype.getSetup = function() {
  return /** @type{?proto.ord.ReactionSetup} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionSetup, 3));
};


/**
 * @param {?proto.ord.ReactionSetup|undefined} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setSetup = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearSetup = function() {
  return this.setSetup(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Reaction.prototype.hasSetup = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional ReactionConditions conditions = 4;
 * @return {?proto.ord.ReactionConditions}
 */
proto.ord.Reaction.prototype.getConditions = function() {
  return /** @type{?proto.ord.ReactionConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionConditions, 4));
};


/**
 * @param {?proto.ord.ReactionConditions|undefined} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setConditions = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearConditions = function() {
  return this.setConditions(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Reaction.prototype.hasConditions = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional ReactionNotes notes = 5;
 * @return {?proto.ord.ReactionNotes}
 */
proto.ord.Reaction.prototype.getNotes = function() {
  return /** @type{?proto.ord.ReactionNotes} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionNotes, 5));
};


/**
 * @param {?proto.ord.ReactionNotes|undefined} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setNotes = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearNotes = function() {
  return this.setNotes(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Reaction.prototype.hasNotes = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * repeated ReactionObservation observations = 6;
 * @return {!Array<!proto.ord.ReactionObservation>}
 */
proto.ord.Reaction.prototype.getObservationsList = function() {
  return /** @type{!Array<!proto.ord.ReactionObservation>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ReactionObservation, 6));
};


/**
 * @param {!Array<!proto.ord.ReactionObservation>} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setObservationsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 6, value);
};


/**
 * @param {!proto.ord.ReactionObservation=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ReactionObservation}
 */
proto.ord.Reaction.prototype.addObservations = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 6, opt_value, proto.ord.ReactionObservation, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearObservationsList = function() {
  return this.setObservationsList([]);
};


/**
 * repeated ReactionWorkup workups = 7;
 * @return {!Array<!proto.ord.ReactionWorkup>}
 */
proto.ord.Reaction.prototype.getWorkupsList = function() {
  return /** @type{!Array<!proto.ord.ReactionWorkup>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ReactionWorkup, 7));
};


/**
 * @param {!Array<!proto.ord.ReactionWorkup>} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setWorkupsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 7, value);
};


/**
 * @param {!proto.ord.ReactionWorkup=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ReactionWorkup}
 */
proto.ord.Reaction.prototype.addWorkups = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 7, opt_value, proto.ord.ReactionWorkup, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearWorkupsList = function() {
  return this.setWorkupsList([]);
};


/**
 * repeated ReactionOutcome outcomes = 8;
 * @return {!Array<!proto.ord.ReactionOutcome>}
 */
proto.ord.Reaction.prototype.getOutcomesList = function() {
  return /** @type{!Array<!proto.ord.ReactionOutcome>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ReactionOutcome, 8));
};


/**
 * @param {!Array<!proto.ord.ReactionOutcome>} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setOutcomesList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 8, value);
};


/**
 * @param {!proto.ord.ReactionOutcome=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ReactionOutcome}
 */
proto.ord.Reaction.prototype.addOutcomes = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 8, opt_value, proto.ord.ReactionOutcome, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearOutcomesList = function() {
  return this.setOutcomesList([]);
};


/**
 * optional ReactionProvenance provenance = 9;
 * @return {?proto.ord.ReactionProvenance}
 */
proto.ord.Reaction.prototype.getProvenance = function() {
  return /** @type{?proto.ord.ReactionProvenance} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionProvenance, 9));
};


/**
 * @param {?proto.ord.ReactionProvenance|undefined} value
 * @return {!proto.ord.Reaction} returns this
*/
proto.ord.Reaction.prototype.setProvenance = function(value) {
  return jspb.Message.setWrapperField(this, 9, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.clearProvenance = function() {
  return this.setProvenance(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Reaction.prototype.hasProvenance = function() {
  return jspb.Message.getField(this, 9) != null;
};


/**
 * optional string reaction_id = 10;
 * @return {string}
 */
proto.ord.Reaction.prototype.getReactionId = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 10, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Reaction} returns this
 */
proto.ord.Reaction.prototype.setReactionId = function(value) {
  return jspb.Message.setProto3StringField(this, 10, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionIdentifier.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionIdentifier.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionIdentifier} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionIdentifier.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    value: jspb.Message.getFieldWithDefault(msg, 3, ""),
    isMapped: jspb.Message.getBooleanFieldWithDefault(msg, 4, false)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionIdentifier}
 */
proto.ord.ReactionIdentifier.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionIdentifier;
  return proto.ord.ReactionIdentifier.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionIdentifier} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionIdentifier}
 */
proto.ord.ReactionIdentifier.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ReactionIdentifier.ReactionIdentifierType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setValue(value);
      break;
    case 4:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsMapped(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionIdentifier.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionIdentifier.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionIdentifier} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionIdentifier.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getValue();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeBool(
      4,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ReactionIdentifier.ReactionIdentifierType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  REACTION_SMILES: 2,
  REACTION_CXSMILES: 6,
  RDFILE: 3,
  RINCHI: 4,
  REACTION_TYPE: 5
};

/**
 * optional ReactionIdentifierType type = 1;
 * @return {!proto.ord.ReactionIdentifier.ReactionIdentifierType}
 */
proto.ord.ReactionIdentifier.prototype.getType = function() {
  return /** @type {!proto.ord.ReactionIdentifier.ReactionIdentifierType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ReactionIdentifier.ReactionIdentifierType} value
 * @return {!proto.ord.ReactionIdentifier} returns this
 */
proto.ord.ReactionIdentifier.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ReactionIdentifier.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionIdentifier} returns this
 */
proto.ord.ReactionIdentifier.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string value = 3;
 * @return {string}
 */
proto.ord.ReactionIdentifier.prototype.getValue = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionIdentifier} returns this
 */
proto.ord.ReactionIdentifier.prototype.setValue = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};


/**
 * optional bool is_mapped = 4;
 * @return {boolean}
 */
proto.ord.ReactionIdentifier.prototype.getIsMapped = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 4, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionIdentifier} returns this
 */
proto.ord.ReactionIdentifier.prototype.setIsMapped = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionIdentifier} returns this
 */
proto.ord.ReactionIdentifier.prototype.clearIsMapped = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionIdentifier.prototype.hasIsMapped = function() {
  return jspb.Message.getField(this, 4) != null;
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.ReactionInput.repeatedFields_ = [1,2];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionInput.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionInput.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionInput} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionInput.toObject = function(includeInstance, msg) {
  var f, obj = {
    componentsList: jspb.Message.toObjectList(msg.getComponentsList(),
    proto.ord.Compound.toObject, includeInstance),
    crudeComponentsList: jspb.Message.toObjectList(msg.getCrudeComponentsList(),
    proto.ord.CrudeComponent.toObject, includeInstance),
    additionOrder: jspb.Message.getFieldWithDefault(msg, 3, 0),
    additionTime: (f = msg.getAdditionTime()) && proto.ord.Time.toObject(includeInstance, f),
    additionSpeed: (f = msg.getAdditionSpeed()) && proto.ord.ReactionInput.AdditionSpeed.toObject(includeInstance, f),
    additionDuration: (f = msg.getAdditionDuration()) && proto.ord.Time.toObject(includeInstance, f),
    flowRate: (f = msg.getFlowRate()) && proto.ord.FlowRate.toObject(includeInstance, f),
    additionDevice: (f = msg.getAdditionDevice()) && proto.ord.ReactionInput.AdditionDevice.toObject(includeInstance, f),
    additionTemperature: (f = msg.getAdditionTemperature()) && proto.ord.Temperature.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionInput}
 */
proto.ord.ReactionInput.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionInput;
  return proto.ord.ReactionInput.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionInput} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionInput}
 */
proto.ord.ReactionInput.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Compound;
      reader.readMessage(value,proto.ord.Compound.deserializeBinaryFromReader);
      msg.addComponents(value);
      break;
    case 2:
      var value = new proto.ord.CrudeComponent;
      reader.readMessage(value,proto.ord.CrudeComponent.deserializeBinaryFromReader);
      msg.addCrudeComponents(value);
      break;
    case 3:
      var value = /** @type {number} */ (reader.readInt32());
      msg.setAdditionOrder(value);
      break;
    case 4:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setAdditionTime(value);
      break;
    case 5:
      var value = new proto.ord.ReactionInput.AdditionSpeed;
      reader.readMessage(value,proto.ord.ReactionInput.AdditionSpeed.deserializeBinaryFromReader);
      msg.setAdditionSpeed(value);
      break;
    case 6:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setAdditionDuration(value);
      break;
    case 7:
      var value = new proto.ord.FlowRate;
      reader.readMessage(value,proto.ord.FlowRate.deserializeBinaryFromReader);
      msg.setFlowRate(value);
      break;
    case 8:
      var value = new proto.ord.ReactionInput.AdditionDevice;
      reader.readMessage(value,proto.ord.ReactionInput.AdditionDevice.deserializeBinaryFromReader);
      msg.setAdditionDevice(value);
      break;
    case 9:
      var value = new proto.ord.Temperature;
      reader.readMessage(value,proto.ord.Temperature.deserializeBinaryFromReader);
      msg.setAdditionTemperature(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionInput.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionInput.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionInput} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionInput.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getComponentsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      1,
      f,
      proto.ord.Compound.serializeBinaryToWriter
    );
  }
  f = message.getCrudeComponentsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      2,
      f,
      proto.ord.CrudeComponent.serializeBinaryToWriter
    );
  }
  f = message.getAdditionOrder();
  if (f !== 0) {
    writer.writeInt32(
      3,
      f
    );
  }
  f = message.getAdditionTime();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getAdditionSpeed();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.ReactionInput.AdditionSpeed.serializeBinaryToWriter
    );
  }
  f = message.getAdditionDuration();
  if (f != null) {
    writer.writeMessage(
      6,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getFlowRate();
  if (f != null) {
    writer.writeMessage(
      7,
      f,
      proto.ord.FlowRate.serializeBinaryToWriter
    );
  }
  f = message.getAdditionDevice();
  if (f != null) {
    writer.writeMessage(
      8,
      f,
      proto.ord.ReactionInput.AdditionDevice.serializeBinaryToWriter
    );
  }
  f = message.getAdditionTemperature();
  if (f != null) {
    writer.writeMessage(
      9,
      f,
      proto.ord.Temperature.serializeBinaryToWriter
    );
  }
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionInput.AdditionSpeed.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionInput.AdditionSpeed.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionInput.AdditionSpeed} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionInput.AdditionSpeed.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionInput.AdditionSpeed}
 */
proto.ord.ReactionInput.AdditionSpeed.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionInput.AdditionSpeed;
  return proto.ord.ReactionInput.AdditionSpeed.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionInput.AdditionSpeed} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionInput.AdditionSpeed}
 */
proto.ord.ReactionInput.AdditionSpeed.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionInput.AdditionSpeed.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionInput.AdditionSpeed.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionInput.AdditionSpeed} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionInput.AdditionSpeed.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType = {
  UNSPECIFIED: 0,
  ALL_AT_ONCE: 1,
  FAST: 2,
  SLOW: 3,
  DROPWISE: 4,
  CONTINUOUS: 5,
  PORTIONWISE: 6
};

/**
 * optional AdditionSpeedType type = 1;
 * @return {!proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType}
 */
proto.ord.ReactionInput.AdditionSpeed.prototype.getType = function() {
  return /** @type {!proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType} value
 * @return {!proto.ord.ReactionInput.AdditionSpeed} returns this
 */
proto.ord.ReactionInput.AdditionSpeed.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ReactionInput.AdditionSpeed.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionInput.AdditionSpeed} returns this
 */
proto.ord.ReactionInput.AdditionSpeed.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionInput.AdditionDevice.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionInput.AdditionDevice.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionInput.AdditionDevice} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionInput.AdditionDevice.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionInput.AdditionDevice}
 */
proto.ord.ReactionInput.AdditionDevice.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionInput.AdditionDevice;
  return proto.ord.ReactionInput.AdditionDevice.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionInput.AdditionDevice} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionInput.AdditionDevice}
 */
proto.ord.ReactionInput.AdditionDevice.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionInput.AdditionDevice.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionInput.AdditionDevice.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionInput.AdditionDevice} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionInput.AdditionDevice.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  NONE: 2,
  SYRINGE: 3,
  CANNULA: 4,
  ADDITION_FUNNEL: 5
};

/**
 * optional AdditionDeviceType type = 1;
 * @return {!proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType}
 */
proto.ord.ReactionInput.AdditionDevice.prototype.getType = function() {
  return /** @type {!proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType} value
 * @return {!proto.ord.ReactionInput.AdditionDevice} returns this
 */
proto.ord.ReactionInput.AdditionDevice.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ReactionInput.AdditionDevice.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionInput.AdditionDevice} returns this
 */
proto.ord.ReactionInput.AdditionDevice.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * repeated Compound components = 1;
 * @return {!Array<!proto.ord.Compound>}
 */
proto.ord.ReactionInput.prototype.getComponentsList = function() {
  return /** @type{!Array<!proto.ord.Compound>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.Compound, 1));
};


/**
 * @param {!Array<!proto.ord.Compound>} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setComponentsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 1, value);
};


/**
 * @param {!proto.ord.Compound=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.Compound}
 */
proto.ord.ReactionInput.prototype.addComponents = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 1, opt_value, proto.ord.Compound, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearComponentsList = function() {
  return this.setComponentsList([]);
};


/**
 * repeated CrudeComponent crude_components = 2;
 * @return {!Array<!proto.ord.CrudeComponent>}
 */
proto.ord.ReactionInput.prototype.getCrudeComponentsList = function() {
  return /** @type{!Array<!proto.ord.CrudeComponent>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.CrudeComponent, 2));
};


/**
 * @param {!Array<!proto.ord.CrudeComponent>} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setCrudeComponentsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 2, value);
};


/**
 * @param {!proto.ord.CrudeComponent=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.CrudeComponent}
 */
proto.ord.ReactionInput.prototype.addCrudeComponents = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 2, opt_value, proto.ord.CrudeComponent, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearCrudeComponentsList = function() {
  return this.setCrudeComponentsList([]);
};


/**
 * optional int32 addition_order = 3;
 * @return {number}
 */
proto.ord.ReactionInput.prototype.getAdditionOrder = function() {
  return /** @type {number} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {number} value
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.setAdditionOrder = function(value) {
  return jspb.Message.setProto3IntField(this, 3, value);
};


/**
 * optional Time addition_time = 4;
 * @return {?proto.ord.Time}
 */
proto.ord.ReactionInput.prototype.getAdditionTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 4));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setAdditionTime = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearAdditionTime = function() {
  return this.setAdditionTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionInput.prototype.hasAdditionTime = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional AdditionSpeed addition_speed = 5;
 * @return {?proto.ord.ReactionInput.AdditionSpeed}
 */
proto.ord.ReactionInput.prototype.getAdditionSpeed = function() {
  return /** @type{?proto.ord.ReactionInput.AdditionSpeed} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionInput.AdditionSpeed, 5));
};


/**
 * @param {?proto.ord.ReactionInput.AdditionSpeed|undefined} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setAdditionSpeed = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearAdditionSpeed = function() {
  return this.setAdditionSpeed(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionInput.prototype.hasAdditionSpeed = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional Time addition_duration = 6;
 * @return {?proto.ord.Time}
 */
proto.ord.ReactionInput.prototype.getAdditionDuration = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 6));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setAdditionDuration = function(value) {
  return jspb.Message.setWrapperField(this, 6, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearAdditionDuration = function() {
  return this.setAdditionDuration(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionInput.prototype.hasAdditionDuration = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * optional FlowRate flow_rate = 7;
 * @return {?proto.ord.FlowRate}
 */
proto.ord.ReactionInput.prototype.getFlowRate = function() {
  return /** @type{?proto.ord.FlowRate} */ (
    jspb.Message.getWrapperField(this, proto.ord.FlowRate, 7));
};


/**
 * @param {?proto.ord.FlowRate|undefined} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setFlowRate = function(value) {
  return jspb.Message.setWrapperField(this, 7, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearFlowRate = function() {
  return this.setFlowRate(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionInput.prototype.hasFlowRate = function() {
  return jspb.Message.getField(this, 7) != null;
};


/**
 * optional AdditionDevice addition_device = 8;
 * @return {?proto.ord.ReactionInput.AdditionDevice}
 */
proto.ord.ReactionInput.prototype.getAdditionDevice = function() {
  return /** @type{?proto.ord.ReactionInput.AdditionDevice} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionInput.AdditionDevice, 8));
};


/**
 * @param {?proto.ord.ReactionInput.AdditionDevice|undefined} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setAdditionDevice = function(value) {
  return jspb.Message.setWrapperField(this, 8, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearAdditionDevice = function() {
  return this.setAdditionDevice(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionInput.prototype.hasAdditionDevice = function() {
  return jspb.Message.getField(this, 8) != null;
};


/**
 * optional Temperature addition_temperature = 9;
 * @return {?proto.ord.Temperature}
 */
proto.ord.ReactionInput.prototype.getAdditionTemperature = function() {
  return /** @type{?proto.ord.Temperature} */ (
    jspb.Message.getWrapperField(this, proto.ord.Temperature, 9));
};


/**
 * @param {?proto.ord.Temperature|undefined} value
 * @return {!proto.ord.ReactionInput} returns this
*/
proto.ord.ReactionInput.prototype.setAdditionTemperature = function(value) {
  return jspb.Message.setWrapperField(this, 9, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionInput} returns this
 */
proto.ord.ReactionInput.prototype.clearAdditionTemperature = function() {
  return this.setAdditionTemperature(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionInput.prototype.hasAdditionTemperature = function() {
  return jspb.Message.getField(this, 9) != null;
};



/**
 * Oneof group definitions for this message. Each group defines the field
 * numbers belonging to that group. When of these fields' value is set, all
 * other fields in the group are cleared. During deserialization, if multiple
 * fields are encountered for a group, only the last value seen will be kept.
 * @private {!Array<!Array<number>>}
 * @const
 */
proto.ord.Amount.oneofGroups_ = [[1,2,3,5]];

/**
 * @enum {number}
 */
proto.ord.Amount.KindCase = {
  KIND_NOT_SET: 0,
  MASS: 1,
  MOLES: 2,
  VOLUME: 3,
  UNMEASURED: 5
};

/**
 * @return {proto.ord.Amount.KindCase}
 */
proto.ord.Amount.prototype.getKindCase = function() {
  return /** @type {proto.ord.Amount.KindCase} */(jspb.Message.computeOneofCase(this, proto.ord.Amount.oneofGroups_[0]));
};



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Amount.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Amount.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Amount} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Amount.toObject = function(includeInstance, msg) {
  var f, obj = {
    mass: (f = msg.getMass()) && proto.ord.Mass.toObject(includeInstance, f),
    moles: (f = msg.getMoles()) && proto.ord.Moles.toObject(includeInstance, f),
    volume: (f = msg.getVolume()) && proto.ord.Volume.toObject(includeInstance, f),
    unmeasured: (f = msg.getUnmeasured()) && proto.ord.UnmeasuredAmount.toObject(includeInstance, f),
    volumeIncludesSolutes: jspb.Message.getBooleanFieldWithDefault(msg, 4, false)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Amount}
 */
proto.ord.Amount.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Amount;
  return proto.ord.Amount.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Amount} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Amount}
 */
proto.ord.Amount.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Mass;
      reader.readMessage(value,proto.ord.Mass.deserializeBinaryFromReader);
      msg.setMass(value);
      break;
    case 2:
      var value = new proto.ord.Moles;
      reader.readMessage(value,proto.ord.Moles.deserializeBinaryFromReader);
      msg.setMoles(value);
      break;
    case 3:
      var value = new proto.ord.Volume;
      reader.readMessage(value,proto.ord.Volume.deserializeBinaryFromReader);
      msg.setVolume(value);
      break;
    case 5:
      var value = new proto.ord.UnmeasuredAmount;
      reader.readMessage(value,proto.ord.UnmeasuredAmount.deserializeBinaryFromReader);
      msg.setUnmeasured(value);
      break;
    case 4:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setVolumeIncludesSolutes(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Amount.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Amount.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Amount} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Amount.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getMass();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.Mass.serializeBinaryToWriter
    );
  }
  f = message.getMoles();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Moles.serializeBinaryToWriter
    );
  }
  f = message.getVolume();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Volume.serializeBinaryToWriter
    );
  }
  f = message.getUnmeasured();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.UnmeasuredAmount.serializeBinaryToWriter
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeBool(
      4,
      f
    );
  }
};


/**
 * optional Mass mass = 1;
 * @return {?proto.ord.Mass}
 */
proto.ord.Amount.prototype.getMass = function() {
  return /** @type{?proto.ord.Mass} */ (
    jspb.Message.getWrapperField(this, proto.ord.Mass, 1));
};


/**
 * @param {?proto.ord.Mass|undefined} value
 * @return {!proto.ord.Amount} returns this
*/
proto.ord.Amount.prototype.setMass = function(value) {
  return jspb.Message.setOneofWrapperField(this, 1, proto.ord.Amount.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Amount} returns this
 */
proto.ord.Amount.prototype.clearMass = function() {
  return this.setMass(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Amount.prototype.hasMass = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional Moles moles = 2;
 * @return {?proto.ord.Moles}
 */
proto.ord.Amount.prototype.getMoles = function() {
  return /** @type{?proto.ord.Moles} */ (
    jspb.Message.getWrapperField(this, proto.ord.Moles, 2));
};


/**
 * @param {?proto.ord.Moles|undefined} value
 * @return {!proto.ord.Amount} returns this
*/
proto.ord.Amount.prototype.setMoles = function(value) {
  return jspb.Message.setOneofWrapperField(this, 2, proto.ord.Amount.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Amount} returns this
 */
proto.ord.Amount.prototype.clearMoles = function() {
  return this.setMoles(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Amount.prototype.hasMoles = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional Volume volume = 3;
 * @return {?proto.ord.Volume}
 */
proto.ord.Amount.prototype.getVolume = function() {
  return /** @type{?proto.ord.Volume} */ (
    jspb.Message.getWrapperField(this, proto.ord.Volume, 3));
};


/**
 * @param {?proto.ord.Volume|undefined} value
 * @return {!proto.ord.Amount} returns this
*/
proto.ord.Amount.prototype.setVolume = function(value) {
  return jspb.Message.setOneofWrapperField(this, 3, proto.ord.Amount.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Amount} returns this
 */
proto.ord.Amount.prototype.clearVolume = function() {
  return this.setVolume(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Amount.prototype.hasVolume = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional UnmeasuredAmount unmeasured = 5;
 * @return {?proto.ord.UnmeasuredAmount}
 */
proto.ord.Amount.prototype.getUnmeasured = function() {
  return /** @type{?proto.ord.UnmeasuredAmount} */ (
    jspb.Message.getWrapperField(this, proto.ord.UnmeasuredAmount, 5));
};


/**
 * @param {?proto.ord.UnmeasuredAmount|undefined} value
 * @return {!proto.ord.Amount} returns this
*/
proto.ord.Amount.prototype.setUnmeasured = function(value) {
  return jspb.Message.setOneofWrapperField(this, 5, proto.ord.Amount.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Amount} returns this
 */
proto.ord.Amount.prototype.clearUnmeasured = function() {
  return this.setUnmeasured(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Amount.prototype.hasUnmeasured = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional bool volume_includes_solutes = 4;
 * @return {boolean}
 */
proto.ord.Amount.prototype.getVolumeIncludesSolutes = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 4, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.Amount} returns this
 */
proto.ord.Amount.prototype.setVolumeIncludesSolutes = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Amount} returns this
 */
proto.ord.Amount.prototype.clearVolumeIncludesSolutes = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Amount.prototype.hasVolumeIncludesSolutes = function() {
  return jspb.Message.getField(this, 4) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.UnmeasuredAmount.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.UnmeasuredAmount.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.UnmeasuredAmount} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.UnmeasuredAmount.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.UnmeasuredAmount}
 */
proto.ord.UnmeasuredAmount.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.UnmeasuredAmount;
  return proto.ord.UnmeasuredAmount.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.UnmeasuredAmount} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.UnmeasuredAmount}
 */
proto.ord.UnmeasuredAmount.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.UnmeasuredAmount.UnmeasuredAmountType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.UnmeasuredAmount.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.UnmeasuredAmount.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.UnmeasuredAmount} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.UnmeasuredAmount.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.UnmeasuredAmount.UnmeasuredAmountType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  SATURATED: 2,
  CATALYTIC: 3,
  TITRATED: 4
};

/**
 * optional UnmeasuredAmountType type = 1;
 * @return {!proto.ord.UnmeasuredAmount.UnmeasuredAmountType}
 */
proto.ord.UnmeasuredAmount.prototype.getType = function() {
  return /** @type {!proto.ord.UnmeasuredAmount.UnmeasuredAmountType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.UnmeasuredAmount.UnmeasuredAmountType} value
 * @return {!proto.ord.UnmeasuredAmount} returns this
 */
proto.ord.UnmeasuredAmount.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.UnmeasuredAmount.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.UnmeasuredAmount} returns this
 */
proto.ord.UnmeasuredAmount.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.CrudeComponent.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.CrudeComponent.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.CrudeComponent} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.CrudeComponent.toObject = function(includeInstance, msg) {
  var f, obj = {
    reactionId: jspb.Message.getFieldWithDefault(msg, 1, ""),
    includesWorkup: jspb.Message.getBooleanFieldWithDefault(msg, 2, false),
    hasDerivedAmount: jspb.Message.getBooleanFieldWithDefault(msg, 3, false),
    amount: (f = msg.getAmount()) && proto.ord.Amount.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.CrudeComponent}
 */
proto.ord.CrudeComponent.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.CrudeComponent;
  return proto.ord.CrudeComponent.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.CrudeComponent} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.CrudeComponent}
 */
proto.ord.CrudeComponent.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {string} */ (reader.readString());
      msg.setReactionId(value);
      break;
    case 2:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIncludesWorkup(value);
      break;
    case 3:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setHasDerivedAmount(value);
      break;
    case 4:
      var value = new proto.ord.Amount;
      reader.readMessage(value,proto.ord.Amount.deserializeBinaryFromReader);
      msg.setAmount(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.CrudeComponent.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.CrudeComponent.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.CrudeComponent} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.CrudeComponent.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getReactionId();
  if (f.length > 0) {
    writer.writeString(
      1,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeBool(
      2,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 3));
  if (f != null) {
    writer.writeBool(
      3,
      f
    );
  }
  f = message.getAmount();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.Amount.serializeBinaryToWriter
    );
  }
};


/**
 * optional string reaction_id = 1;
 * @return {string}
 */
proto.ord.CrudeComponent.prototype.getReactionId = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 1, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.CrudeComponent} returns this
 */
proto.ord.CrudeComponent.prototype.setReactionId = function(value) {
  return jspb.Message.setProto3StringField(this, 1, value);
};


/**
 * optional bool includes_workup = 2;
 * @return {boolean}
 */
proto.ord.CrudeComponent.prototype.getIncludesWorkup = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 2, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.CrudeComponent} returns this
 */
proto.ord.CrudeComponent.prototype.setIncludesWorkup = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.CrudeComponent} returns this
 */
proto.ord.CrudeComponent.prototype.clearIncludesWorkup = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.CrudeComponent.prototype.hasIncludesWorkup = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional bool has_derived_amount = 3;
 * @return {boolean}
 */
proto.ord.CrudeComponent.prototype.getHasDerivedAmount = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 3, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.CrudeComponent} returns this
 */
proto.ord.CrudeComponent.prototype.setHasDerivedAmount = function(value) {
  return jspb.Message.setField(this, 3, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.CrudeComponent} returns this
 */
proto.ord.CrudeComponent.prototype.clearHasDerivedAmount = function() {
  return jspb.Message.setField(this, 3, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.CrudeComponent.prototype.hasHasDerivedAmount = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional Amount amount = 4;
 * @return {?proto.ord.Amount}
 */
proto.ord.CrudeComponent.prototype.getAmount = function() {
  return /** @type{?proto.ord.Amount} */ (
    jspb.Message.getWrapperField(this, proto.ord.Amount, 4));
};


/**
 * @param {?proto.ord.Amount|undefined} value
 * @return {!proto.ord.CrudeComponent} returns this
*/
proto.ord.CrudeComponent.prototype.setAmount = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.CrudeComponent} returns this
 */
proto.ord.CrudeComponent.prototype.clearAmount = function() {
  return this.setAmount(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.CrudeComponent.prototype.hasAmount = function() {
  return jspb.Message.getField(this, 4) != null;
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.Compound.repeatedFields_ = [1,5];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Compound.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Compound.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Compound} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Compound.toObject = function(includeInstance, msg) {
  var f, obj = {
    identifiersList: jspb.Message.toObjectList(msg.getIdentifiersList(),
    proto.ord.CompoundIdentifier.toObject, includeInstance),
    amount: (f = msg.getAmount()) && proto.ord.Amount.toObject(includeInstance, f),
    reactionRole: jspb.Message.getFieldWithDefault(msg, 3, 0),
    isLimiting: jspb.Message.getBooleanFieldWithDefault(msg, 4, false),
    preparationsList: jspb.Message.toObjectList(msg.getPreparationsList(),
    proto.ord.CompoundPreparation.toObject, includeInstance),
    source: (f = msg.getSource()) && proto.ord.Compound.Source.toObject(includeInstance, f),
    featuresMap: (f = msg.getFeaturesMap()) ? f.toObject(includeInstance, proto.ord.Data.toObject) : [],
    analysesMap: (f = msg.getAnalysesMap()) ? f.toObject(includeInstance, proto.ord.Analysis.toObject) : []
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Compound}
 */
proto.ord.Compound.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Compound;
  return proto.ord.Compound.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Compound} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Compound}
 */
proto.ord.Compound.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.CompoundIdentifier;
      reader.readMessage(value,proto.ord.CompoundIdentifier.deserializeBinaryFromReader);
      msg.addIdentifiers(value);
      break;
    case 2:
      var value = new proto.ord.Amount;
      reader.readMessage(value,proto.ord.Amount.deserializeBinaryFromReader);
      msg.setAmount(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.ReactionRole.ReactionRoleType} */ (reader.readEnum());
      msg.setReactionRole(value);
      break;
    case 4:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsLimiting(value);
      break;
    case 5:
      var value = new proto.ord.CompoundPreparation;
      reader.readMessage(value,proto.ord.CompoundPreparation.deserializeBinaryFromReader);
      msg.addPreparations(value);
      break;
    case 6:
      var value = new proto.ord.Compound.Source;
      reader.readMessage(value,proto.ord.Compound.Source.deserializeBinaryFromReader);
      msg.setSource(value);
      break;
    case 7:
      var value = msg.getFeaturesMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Data.deserializeBinaryFromReader, "", new proto.ord.Data());
         });
      break;
    case 8:
      var value = msg.getAnalysesMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Analysis.deserializeBinaryFromReader, "", new proto.ord.Analysis());
         });
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Compound.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Compound.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Compound} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Compound.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getIdentifiersList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      1,
      f,
      proto.ord.CompoundIdentifier.serializeBinaryToWriter
    );
  }
  f = message.getAmount();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Amount.serializeBinaryToWriter
    );
  }
  f = message.getReactionRole();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeBool(
      4,
      f
    );
  }
  f = message.getPreparationsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      5,
      f,
      proto.ord.CompoundPreparation.serializeBinaryToWriter
    );
  }
  f = message.getSource();
  if (f != null) {
    writer.writeMessage(
      6,
      f,
      proto.ord.Compound.Source.serializeBinaryToWriter
    );
  }
  f = message.getFeaturesMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(7, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Data.serializeBinaryToWriter);
  }
  f = message.getAnalysesMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(8, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Analysis.serializeBinaryToWriter);
  }
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Compound.Source.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Compound.Source.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Compound.Source} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Compound.Source.toObject = function(includeInstance, msg) {
  var f, obj = {
    vendor: jspb.Message.getFieldWithDefault(msg, 1, ""),
    catalogId: jspb.Message.getFieldWithDefault(msg, 2, ""),
    lot: jspb.Message.getFieldWithDefault(msg, 3, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Compound.Source}
 */
proto.ord.Compound.Source.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Compound.Source;
  return proto.ord.Compound.Source.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Compound.Source} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Compound.Source}
 */
proto.ord.Compound.Source.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {string} */ (reader.readString());
      msg.setVendor(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setCatalogId(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setLot(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Compound.Source.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Compound.Source.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Compound.Source} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Compound.Source.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getVendor();
  if (f.length > 0) {
    writer.writeString(
      1,
      f
    );
  }
  f = message.getCatalogId();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getLot();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
};


/**
 * optional string vendor = 1;
 * @return {string}
 */
proto.ord.Compound.Source.prototype.getVendor = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 1, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Compound.Source} returns this
 */
proto.ord.Compound.Source.prototype.setVendor = function(value) {
  return jspb.Message.setProto3StringField(this, 1, value);
};


/**
 * optional string catalog_id = 2;
 * @return {string}
 */
proto.ord.Compound.Source.prototype.getCatalogId = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Compound.Source} returns this
 */
proto.ord.Compound.Source.prototype.setCatalogId = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string lot = 3;
 * @return {string}
 */
proto.ord.Compound.Source.prototype.getLot = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Compound.Source} returns this
 */
proto.ord.Compound.Source.prototype.setLot = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};


/**
 * repeated CompoundIdentifier identifiers = 1;
 * @return {!Array<!proto.ord.CompoundIdentifier>}
 */
proto.ord.Compound.prototype.getIdentifiersList = function() {
  return /** @type{!Array<!proto.ord.CompoundIdentifier>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.CompoundIdentifier, 1));
};


/**
 * @param {!Array<!proto.ord.CompoundIdentifier>} value
 * @return {!proto.ord.Compound} returns this
*/
proto.ord.Compound.prototype.setIdentifiersList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 1, value);
};


/**
 * @param {!proto.ord.CompoundIdentifier=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.CompoundIdentifier}
 */
proto.ord.Compound.prototype.addIdentifiers = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 1, opt_value, proto.ord.CompoundIdentifier, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearIdentifiersList = function() {
  return this.setIdentifiersList([]);
};


/**
 * optional Amount amount = 2;
 * @return {?proto.ord.Amount}
 */
proto.ord.Compound.prototype.getAmount = function() {
  return /** @type{?proto.ord.Amount} */ (
    jspb.Message.getWrapperField(this, proto.ord.Amount, 2));
};


/**
 * @param {?proto.ord.Amount|undefined} value
 * @return {!proto.ord.Compound} returns this
*/
proto.ord.Compound.prototype.setAmount = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearAmount = function() {
  return this.setAmount(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Compound.prototype.hasAmount = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional ReactionRole.ReactionRoleType reaction_role = 3;
 * @return {!proto.ord.ReactionRole.ReactionRoleType}
 */
proto.ord.Compound.prototype.getReactionRole = function() {
  return /** @type {!proto.ord.ReactionRole.ReactionRoleType} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.ReactionRole.ReactionRoleType} value
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.setReactionRole = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};


/**
 * optional bool is_limiting = 4;
 * @return {boolean}
 */
proto.ord.Compound.prototype.getIsLimiting = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 4, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.setIsLimiting = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearIsLimiting = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Compound.prototype.hasIsLimiting = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * repeated CompoundPreparation preparations = 5;
 * @return {!Array<!proto.ord.CompoundPreparation>}
 */
proto.ord.Compound.prototype.getPreparationsList = function() {
  return /** @type{!Array<!proto.ord.CompoundPreparation>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.CompoundPreparation, 5));
};


/**
 * @param {!Array<!proto.ord.CompoundPreparation>} value
 * @return {!proto.ord.Compound} returns this
*/
proto.ord.Compound.prototype.setPreparationsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 5, value);
};


/**
 * @param {!proto.ord.CompoundPreparation=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.CompoundPreparation}
 */
proto.ord.Compound.prototype.addPreparations = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 5, opt_value, proto.ord.CompoundPreparation, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearPreparationsList = function() {
  return this.setPreparationsList([]);
};


/**
 * optional Source source = 6;
 * @return {?proto.ord.Compound.Source}
 */
proto.ord.Compound.prototype.getSource = function() {
  return /** @type{?proto.ord.Compound.Source} */ (
    jspb.Message.getWrapperField(this, proto.ord.Compound.Source, 6));
};


/**
 * @param {?proto.ord.Compound.Source|undefined} value
 * @return {!proto.ord.Compound} returns this
*/
proto.ord.Compound.prototype.setSource = function(value) {
  return jspb.Message.setWrapperField(this, 6, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearSource = function() {
  return this.setSource(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Compound.prototype.hasSource = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * map<string, Data> features = 7;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Data>}
 */
proto.ord.Compound.prototype.getFeaturesMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Data>} */ (
      jspb.Message.getMapField(this, 7, opt_noLazyCreate,
      proto.ord.Data));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearFeaturesMap = function() {
  this.getFeaturesMap().clear();
  return this;
};


/**
 * map<string, Analysis> analyses = 8;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Analysis>}
 */
proto.ord.Compound.prototype.getAnalysesMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Analysis>} */ (
      jspb.Message.getMapField(this, 8, opt_noLazyCreate,
      proto.ord.Analysis));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.Compound} returns this
 */
proto.ord.Compound.prototype.clearAnalysesMap = function() {
  this.getAnalysesMap().clear();
  return this;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionRole.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionRole.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionRole} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionRole.toObject = function(includeInstance, msg) {
  var f, obj = {

  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionRole}
 */
proto.ord.ReactionRole.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionRole;
  return proto.ord.ReactionRole.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionRole} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionRole}
 */
proto.ord.ReactionRole.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionRole.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionRole.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionRole} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionRole.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
};


/**
 * @enum {number}
 */
proto.ord.ReactionRole.ReactionRoleType = {
  UNSPECIFIED: 0,
  REACTANT: 1,
  REAGENT: 2,
  SOLVENT: 3,
  CATALYST: 4,
  WORKUP: 5,
  INTERNAL_STANDARD: 6,
  AUTHENTIC_STANDARD: 7,
  PRODUCT: 8,
  BYPRODUCT: 9,
  SIDE_PRODUCT: 10
};




if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.CompoundPreparation.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.CompoundPreparation.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.CompoundPreparation} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.CompoundPreparation.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    reactionId: jspb.Message.getFieldWithDefault(msg, 3, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.CompoundPreparation}
 */
proto.ord.CompoundPreparation.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.CompoundPreparation;
  return proto.ord.CompoundPreparation.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.CompoundPreparation} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.CompoundPreparation}
 */
proto.ord.CompoundPreparation.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.CompoundPreparation.CompoundPreparationType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setReactionId(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.CompoundPreparation.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.CompoundPreparation.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.CompoundPreparation} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.CompoundPreparation.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getReactionId();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.CompoundPreparation.CompoundPreparationType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  NONE: 2,
  REPURIFIED: 3,
  SPARGED: 4,
  DRIED: 5,
  SYNTHESIZED: 6
};

/**
 * optional CompoundPreparationType type = 1;
 * @return {!proto.ord.CompoundPreparation.CompoundPreparationType}
 */
proto.ord.CompoundPreparation.prototype.getType = function() {
  return /** @type {!proto.ord.CompoundPreparation.CompoundPreparationType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.CompoundPreparation.CompoundPreparationType} value
 * @return {!proto.ord.CompoundPreparation} returns this
 */
proto.ord.CompoundPreparation.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.CompoundPreparation.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.CompoundPreparation} returns this
 */
proto.ord.CompoundPreparation.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string reaction_id = 3;
 * @return {string}
 */
proto.ord.CompoundPreparation.prototype.getReactionId = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.CompoundPreparation} returns this
 */
proto.ord.CompoundPreparation.prototype.setReactionId = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.CompoundIdentifier.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.CompoundIdentifier.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.CompoundIdentifier} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.CompoundIdentifier.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    value: jspb.Message.getFieldWithDefault(msg, 3, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.CompoundIdentifier}
 */
proto.ord.CompoundIdentifier.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.CompoundIdentifier;
  return proto.ord.CompoundIdentifier.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.CompoundIdentifier} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.CompoundIdentifier}
 */
proto.ord.CompoundIdentifier.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.CompoundIdentifier.CompoundIdentifierType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setValue(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.CompoundIdentifier.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.CompoundIdentifier.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.CompoundIdentifier} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.CompoundIdentifier.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getValue();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.CompoundIdentifier.CompoundIdentifierType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  SMILES: 2,
  INCHI: 3,
  MOLBLOCK: 4,
  IUPAC_NAME: 5,
  NAME: 6,
  CAS_NUMBER: 7,
  PUBCHEM_CID: 8,
  CHEMSPIDER_ID: 9,
  CXSMILES: 10,
  INCHI_KEY: 11,
  XYZ: 12,
  UNIPROT_ID: 13,
  PDB_ID: 14,
  AMINO_ACID_SEQUENCE: 15,
  HELM: 16
};

/**
 * optional CompoundIdentifierType type = 1;
 * @return {!proto.ord.CompoundIdentifier.CompoundIdentifierType}
 */
proto.ord.CompoundIdentifier.prototype.getType = function() {
  return /** @type {!proto.ord.CompoundIdentifier.CompoundIdentifierType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.CompoundIdentifier.CompoundIdentifierType} value
 * @return {!proto.ord.CompoundIdentifier} returns this
 */
proto.ord.CompoundIdentifier.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.CompoundIdentifier.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.CompoundIdentifier} returns this
 */
proto.ord.CompoundIdentifier.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string value = 3;
 * @return {string}
 */
proto.ord.CompoundIdentifier.prototype.getValue = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.CompoundIdentifier} returns this
 */
proto.ord.CompoundIdentifier.prototype.setValue = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.Vessel.repeatedFields_ = [4,5];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Vessel.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Vessel.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Vessel} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Vessel.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    material: (f = msg.getMaterial()) && proto.ord.VesselMaterial.toObject(includeInstance, f),
    preparationsList: jspb.Message.toObjectList(msg.getPreparationsList(),
    proto.ord.VesselPreparation.toObject, includeInstance),
    attachmentsList: jspb.Message.toObjectList(msg.getAttachmentsList(),
    proto.ord.VesselAttachment.toObject, includeInstance),
    volume: (f = msg.getVolume()) && proto.ord.Volume.toObject(includeInstance, f),
    vesselId: jspb.Message.getFieldWithDefault(msg, 7, ""),
    position: jspb.Message.getFieldWithDefault(msg, 8, ""),
    row: jspb.Message.getFieldWithDefault(msg, 9, ""),
    col: jspb.Message.getFieldWithDefault(msg, 10, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Vessel}
 */
proto.ord.Vessel.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Vessel;
  return proto.ord.Vessel.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Vessel} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Vessel}
 */
proto.ord.Vessel.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.Vessel.VesselType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.VesselMaterial;
      reader.readMessage(value,proto.ord.VesselMaterial.deserializeBinaryFromReader);
      msg.setMaterial(value);
      break;
    case 4:
      var value = new proto.ord.VesselPreparation;
      reader.readMessage(value,proto.ord.VesselPreparation.deserializeBinaryFromReader);
      msg.addPreparations(value);
      break;
    case 5:
      var value = new proto.ord.VesselAttachment;
      reader.readMessage(value,proto.ord.VesselAttachment.deserializeBinaryFromReader);
      msg.addAttachments(value);
      break;
    case 6:
      var value = new proto.ord.Volume;
      reader.readMessage(value,proto.ord.Volume.deserializeBinaryFromReader);
      msg.setVolume(value);
      break;
    case 7:
      var value = /** @type {string} */ (reader.readString());
      msg.setVesselId(value);
      break;
    case 8:
      var value = /** @type {string} */ (reader.readString());
      msg.setPosition(value);
      break;
    case 9:
      var value = /** @type {string} */ (reader.readString());
      msg.setRow(value);
      break;
    case 10:
      var value = /** @type {string} */ (reader.readString());
      msg.setCol(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Vessel.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Vessel.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Vessel} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Vessel.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getMaterial();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.VesselMaterial.serializeBinaryToWriter
    );
  }
  f = message.getPreparationsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      4,
      f,
      proto.ord.VesselPreparation.serializeBinaryToWriter
    );
  }
  f = message.getAttachmentsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      5,
      f,
      proto.ord.VesselAttachment.serializeBinaryToWriter
    );
  }
  f = message.getVolume();
  if (f != null) {
    writer.writeMessage(
      6,
      f,
      proto.ord.Volume.serializeBinaryToWriter
    );
  }
  f = message.getVesselId();
  if (f.length > 0) {
    writer.writeString(
      7,
      f
    );
  }
  f = message.getPosition();
  if (f.length > 0) {
    writer.writeString(
      8,
      f
    );
  }
  f = message.getRow();
  if (f.length > 0) {
    writer.writeString(
      9,
      f
    );
  }
  f = message.getCol();
  if (f.length > 0) {
    writer.writeString(
      10,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Vessel.VesselType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  ROUND_BOTTOM_FLASK: 2,
  VIAL: 3,
  WELL_PLATE: 4,
  MICROWAVE_VIAL: 5,
  TUBE: 6,
  CONTINUOUS_STIRRED_TANK_REACTOR: 7,
  PACKED_BED_REACTOR: 8,
  NMR_TUBE: 9,
  PRESSURE_FLASK: 10,
  PRESSURE_REACTOR: 11,
  ELECTROCHEMICAL_CELL: 12
};

/**
 * optional VesselType type = 1;
 * @return {!proto.ord.Vessel.VesselType}
 */
proto.ord.Vessel.prototype.getType = function() {
  return /** @type {!proto.ord.Vessel.VesselType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.Vessel.VesselType} value
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.Vessel.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional VesselMaterial material = 3;
 * @return {?proto.ord.VesselMaterial}
 */
proto.ord.Vessel.prototype.getMaterial = function() {
  return /** @type{?proto.ord.VesselMaterial} */ (
    jspb.Message.getWrapperField(this, proto.ord.VesselMaterial, 3));
};


/**
 * @param {?proto.ord.VesselMaterial|undefined} value
 * @return {!proto.ord.Vessel} returns this
*/
proto.ord.Vessel.prototype.setMaterial = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.clearMaterial = function() {
  return this.setMaterial(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Vessel.prototype.hasMaterial = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * repeated VesselPreparation preparations = 4;
 * @return {!Array<!proto.ord.VesselPreparation>}
 */
proto.ord.Vessel.prototype.getPreparationsList = function() {
  return /** @type{!Array<!proto.ord.VesselPreparation>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.VesselPreparation, 4));
};


/**
 * @param {!Array<!proto.ord.VesselPreparation>} value
 * @return {!proto.ord.Vessel} returns this
*/
proto.ord.Vessel.prototype.setPreparationsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 4, value);
};


/**
 * @param {!proto.ord.VesselPreparation=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.VesselPreparation}
 */
proto.ord.Vessel.prototype.addPreparations = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 4, opt_value, proto.ord.VesselPreparation, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.clearPreparationsList = function() {
  return this.setPreparationsList([]);
};


/**
 * repeated VesselAttachment attachments = 5;
 * @return {!Array<!proto.ord.VesselAttachment>}
 */
proto.ord.Vessel.prototype.getAttachmentsList = function() {
  return /** @type{!Array<!proto.ord.VesselAttachment>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.VesselAttachment, 5));
};


/**
 * @param {!Array<!proto.ord.VesselAttachment>} value
 * @return {!proto.ord.Vessel} returns this
*/
proto.ord.Vessel.prototype.setAttachmentsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 5, value);
};


/**
 * @param {!proto.ord.VesselAttachment=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.VesselAttachment}
 */
proto.ord.Vessel.prototype.addAttachments = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 5, opt_value, proto.ord.VesselAttachment, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.clearAttachmentsList = function() {
  return this.setAttachmentsList([]);
};


/**
 * optional Volume volume = 6;
 * @return {?proto.ord.Volume}
 */
proto.ord.Vessel.prototype.getVolume = function() {
  return /** @type{?proto.ord.Volume} */ (
    jspb.Message.getWrapperField(this, proto.ord.Volume, 6));
};


/**
 * @param {?proto.ord.Volume|undefined} value
 * @return {!proto.ord.Vessel} returns this
*/
proto.ord.Vessel.prototype.setVolume = function(value) {
  return jspb.Message.setWrapperField(this, 6, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.clearVolume = function() {
  return this.setVolume(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Vessel.prototype.hasVolume = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * optional string vessel_id = 7;
 * @return {string}
 */
proto.ord.Vessel.prototype.getVesselId = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 7, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.setVesselId = function(value) {
  return jspb.Message.setProto3StringField(this, 7, value);
};


/**
 * optional string position = 8;
 * @return {string}
 */
proto.ord.Vessel.prototype.getPosition = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 8, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.setPosition = function(value) {
  return jspb.Message.setProto3StringField(this, 8, value);
};


/**
 * optional string row = 9;
 * @return {string}
 */
proto.ord.Vessel.prototype.getRow = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 9, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.setRow = function(value) {
  return jspb.Message.setProto3StringField(this, 9, value);
};


/**
 * optional string col = 10;
 * @return {string}
 */
proto.ord.Vessel.prototype.getCol = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 10, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Vessel} returns this
 */
proto.ord.Vessel.prototype.setCol = function(value) {
  return jspb.Message.setProto3StringField(this, 10, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.VesselMaterial.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.VesselMaterial.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.VesselMaterial} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.VesselMaterial.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.VesselMaterial}
 */
proto.ord.VesselMaterial.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.VesselMaterial;
  return proto.ord.VesselMaterial.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.VesselMaterial} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.VesselMaterial}
 */
proto.ord.VesselMaterial.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.VesselMaterial.VesselMaterialType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.VesselMaterial.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.VesselMaterial.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.VesselMaterial} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.VesselMaterial.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.VesselMaterial.VesselMaterialType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  GLASS: 2,
  POLYPROPYLENE: 3,
  PLASTIC: 4,
  METAL: 5,
  QUARTZ: 6
};

/**
 * optional VesselMaterialType type = 1;
 * @return {!proto.ord.VesselMaterial.VesselMaterialType}
 */
proto.ord.VesselMaterial.prototype.getType = function() {
  return /** @type {!proto.ord.VesselMaterial.VesselMaterialType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.VesselMaterial.VesselMaterialType} value
 * @return {!proto.ord.VesselMaterial} returns this
 */
proto.ord.VesselMaterial.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.VesselMaterial.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.VesselMaterial} returns this
 */
proto.ord.VesselMaterial.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.VesselAttachment.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.VesselAttachment.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.VesselAttachment} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.VesselAttachment.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.VesselAttachment}
 */
proto.ord.VesselAttachment.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.VesselAttachment;
  return proto.ord.VesselAttachment.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.VesselAttachment} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.VesselAttachment}
 */
proto.ord.VesselAttachment.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.VesselAttachment.VesselAttachmentType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.VesselAttachment.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.VesselAttachment.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.VesselAttachment} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.VesselAttachment.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.VesselAttachment.VesselAttachmentType = {
  UNSPECIFIED: 0,
  NONE: 1,
  CUSTOM: 2,
  SEPTUM: 3,
  CAP: 4,
  MAT: 5,
  REFLUX_CONDENSER: 6,
  VENT_NEEDLE: 7,
  DEAN_STARK: 8,
  VACUUM_TUBE: 9,
  ADDITION_FUNNEL: 10,
  DRYING_TUBE: 11,
  ALUMINUM_FOIL: 12,
  THERMOCOUPLE: 13,
  BALLOON: 14,
  GAS_ADAPTER: 15,
  PRESSURE_REGULATOR: 16,
  RELEASE_VALVE: 17
};

/**
 * optional VesselAttachmentType type = 1;
 * @return {!proto.ord.VesselAttachment.VesselAttachmentType}
 */
proto.ord.VesselAttachment.prototype.getType = function() {
  return /** @type {!proto.ord.VesselAttachment.VesselAttachmentType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.VesselAttachment.VesselAttachmentType} value
 * @return {!proto.ord.VesselAttachment} returns this
 */
proto.ord.VesselAttachment.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.VesselAttachment.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.VesselAttachment} returns this
 */
proto.ord.VesselAttachment.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.VesselPreparation.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.VesselPreparation.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.VesselPreparation} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.VesselPreparation.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.VesselPreparation}
 */
proto.ord.VesselPreparation.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.VesselPreparation;
  return proto.ord.VesselPreparation.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.VesselPreparation} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.VesselPreparation}
 */
proto.ord.VesselPreparation.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.VesselPreparation.VesselPreparationType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.VesselPreparation.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.VesselPreparation.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.VesselPreparation} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.VesselPreparation.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.VesselPreparation.VesselPreparationType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  NONE: 2,
  OVEN_DRIED: 3,
  FLAME_DRIED: 4,
  EVACUATED_BACKFILLED: 5,
  PURGED: 6
};

/**
 * optional VesselPreparationType type = 1;
 * @return {!proto.ord.VesselPreparation.VesselPreparationType}
 */
proto.ord.VesselPreparation.prototype.getType = function() {
  return /** @type {!proto.ord.VesselPreparation.VesselPreparationType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.VesselPreparation.VesselPreparationType} value
 * @return {!proto.ord.VesselPreparation} returns this
 */
proto.ord.VesselPreparation.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.VesselPreparation.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.VesselPreparation} returns this
 */
proto.ord.VesselPreparation.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionSetup.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionSetup.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionSetup} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionSetup.toObject = function(includeInstance, msg) {
  var f, obj = {
    vessel: (f = msg.getVessel()) && proto.ord.Vessel.toObject(includeInstance, f),
    isAutomated: jspb.Message.getBooleanFieldWithDefault(msg, 2, false),
    automationPlatform: jspb.Message.getFieldWithDefault(msg, 3, ""),
    automationCodeMap: (f = msg.getAutomationCodeMap()) ? f.toObject(includeInstance, proto.ord.Data.toObject) : [],
    environment: (f = msg.getEnvironment()) && proto.ord.ReactionSetup.ReactionEnvironment.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionSetup}
 */
proto.ord.ReactionSetup.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionSetup;
  return proto.ord.ReactionSetup.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionSetup} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionSetup}
 */
proto.ord.ReactionSetup.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Vessel;
      reader.readMessage(value,proto.ord.Vessel.deserializeBinaryFromReader);
      msg.setVessel(value);
      break;
    case 2:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsAutomated(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setAutomationPlatform(value);
      break;
    case 4:
      var value = msg.getAutomationCodeMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Data.deserializeBinaryFromReader, "", new proto.ord.Data());
         });
      break;
    case 5:
      var value = new proto.ord.ReactionSetup.ReactionEnvironment;
      reader.readMessage(value,proto.ord.ReactionSetup.ReactionEnvironment.deserializeBinaryFromReader);
      msg.setEnvironment(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionSetup.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionSetup.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionSetup} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionSetup.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getVessel();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.Vessel.serializeBinaryToWriter
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeBool(
      2,
      f
    );
  }
  f = message.getAutomationPlatform();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
  f = message.getAutomationCodeMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(4, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Data.serializeBinaryToWriter);
  }
  f = message.getEnvironment();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.ReactionSetup.ReactionEnvironment.serializeBinaryToWriter
    );
  }
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionSetup.ReactionEnvironment.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionSetup.ReactionEnvironment.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionSetup.ReactionEnvironment} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionSetup.ReactionEnvironment.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionSetup.ReactionEnvironment}
 */
proto.ord.ReactionSetup.ReactionEnvironment.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionSetup.ReactionEnvironment;
  return proto.ord.ReactionSetup.ReactionEnvironment.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionSetup.ReactionEnvironment} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionSetup.ReactionEnvironment}
 */
proto.ord.ReactionSetup.ReactionEnvironment.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionSetup.ReactionEnvironment.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionSetup.ReactionEnvironment.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionSetup.ReactionEnvironment} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionSetup.ReactionEnvironment.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  FUME_HOOD: 2,
  BENCH_TOP: 3,
  GLOVE_BOX: 4,
  GLOVE_BAG: 5
};

/**
 * optional ReactionEnvironmentType type = 1;
 * @return {!proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType}
 */
proto.ord.ReactionSetup.ReactionEnvironment.prototype.getType = function() {
  return /** @type {!proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType} value
 * @return {!proto.ord.ReactionSetup.ReactionEnvironment} returns this
 */
proto.ord.ReactionSetup.ReactionEnvironment.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ReactionSetup.ReactionEnvironment.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionSetup.ReactionEnvironment} returns this
 */
proto.ord.ReactionSetup.ReactionEnvironment.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Vessel vessel = 1;
 * @return {?proto.ord.Vessel}
 */
proto.ord.ReactionSetup.prototype.getVessel = function() {
  return /** @type{?proto.ord.Vessel} */ (
    jspb.Message.getWrapperField(this, proto.ord.Vessel, 1));
};


/**
 * @param {?proto.ord.Vessel|undefined} value
 * @return {!proto.ord.ReactionSetup} returns this
*/
proto.ord.ReactionSetup.prototype.setVessel = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionSetup} returns this
 */
proto.ord.ReactionSetup.prototype.clearVessel = function() {
  return this.setVessel(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionSetup.prototype.hasVessel = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional bool is_automated = 2;
 * @return {boolean}
 */
proto.ord.ReactionSetup.prototype.getIsAutomated = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 2, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionSetup} returns this
 */
proto.ord.ReactionSetup.prototype.setIsAutomated = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionSetup} returns this
 */
proto.ord.ReactionSetup.prototype.clearIsAutomated = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionSetup.prototype.hasIsAutomated = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional string automation_platform = 3;
 * @return {string}
 */
proto.ord.ReactionSetup.prototype.getAutomationPlatform = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionSetup} returns this
 */
proto.ord.ReactionSetup.prototype.setAutomationPlatform = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};


/**
 * map<string, Data> automation_code = 4;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Data>}
 */
proto.ord.ReactionSetup.prototype.getAutomationCodeMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Data>} */ (
      jspb.Message.getMapField(this, 4, opt_noLazyCreate,
      proto.ord.Data));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.ReactionSetup} returns this
 */
proto.ord.ReactionSetup.prototype.clearAutomationCodeMap = function() {
  this.getAutomationCodeMap().clear();
  return this;
};


/**
 * optional ReactionEnvironment environment = 5;
 * @return {?proto.ord.ReactionSetup.ReactionEnvironment}
 */
proto.ord.ReactionSetup.prototype.getEnvironment = function() {
  return /** @type{?proto.ord.ReactionSetup.ReactionEnvironment} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionSetup.ReactionEnvironment, 5));
};


/**
 * @param {?proto.ord.ReactionSetup.ReactionEnvironment|undefined} value
 * @return {!proto.ord.ReactionSetup} returns this
*/
proto.ord.ReactionSetup.prototype.setEnvironment = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionSetup} returns this
 */
proto.ord.ReactionSetup.prototype.clearEnvironment = function() {
  return this.setEnvironment(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionSetup.prototype.hasEnvironment = function() {
  return jspb.Message.getField(this, 5) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    temperature: (f = msg.getTemperature()) && proto.ord.TemperatureConditions.toObject(includeInstance, f),
    pressure: (f = msg.getPressure()) && proto.ord.PressureConditions.toObject(includeInstance, f),
    stirring: (f = msg.getStirring()) && proto.ord.StirringConditions.toObject(includeInstance, f),
    illumination: (f = msg.getIllumination()) && proto.ord.IlluminationConditions.toObject(includeInstance, f),
    electrochemistry: (f = msg.getElectrochemistry()) && proto.ord.ElectrochemistryConditions.toObject(includeInstance, f),
    flow: (f = msg.getFlow()) && proto.ord.FlowConditions.toObject(includeInstance, f),
    reflux: jspb.Message.getBooleanFieldWithDefault(msg, 7, false),
    ph: jspb.Message.getFloatingPointFieldWithDefault(msg, 8, 0.0),
    conditionsAreDynamic: jspb.Message.getBooleanFieldWithDefault(msg, 9, false),
    details: jspb.Message.getFieldWithDefault(msg, 10, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionConditions}
 */
proto.ord.ReactionConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionConditions;
  return proto.ord.ReactionConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionConditions}
 */
proto.ord.ReactionConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.TemperatureConditions;
      reader.readMessage(value,proto.ord.TemperatureConditions.deserializeBinaryFromReader);
      msg.setTemperature(value);
      break;
    case 2:
      var value = new proto.ord.PressureConditions;
      reader.readMessage(value,proto.ord.PressureConditions.deserializeBinaryFromReader);
      msg.setPressure(value);
      break;
    case 3:
      var value = new proto.ord.StirringConditions;
      reader.readMessage(value,proto.ord.StirringConditions.deserializeBinaryFromReader);
      msg.setStirring(value);
      break;
    case 4:
      var value = new proto.ord.IlluminationConditions;
      reader.readMessage(value,proto.ord.IlluminationConditions.deserializeBinaryFromReader);
      msg.setIllumination(value);
      break;
    case 5:
      var value = new proto.ord.ElectrochemistryConditions;
      reader.readMessage(value,proto.ord.ElectrochemistryConditions.deserializeBinaryFromReader);
      msg.setElectrochemistry(value);
      break;
    case 6:
      var value = new proto.ord.FlowConditions;
      reader.readMessage(value,proto.ord.FlowConditions.deserializeBinaryFromReader);
      msg.setFlow(value);
      break;
    case 7:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setReflux(value);
      break;
    case 8:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPh(value);
      break;
    case 9:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setConditionsAreDynamic(value);
      break;
    case 10:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getTemperature();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.TemperatureConditions.serializeBinaryToWriter
    );
  }
  f = message.getPressure();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.PressureConditions.serializeBinaryToWriter
    );
  }
  f = message.getStirring();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.StirringConditions.serializeBinaryToWriter
    );
  }
  f = message.getIllumination();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.IlluminationConditions.serializeBinaryToWriter
    );
  }
  f = message.getElectrochemistry();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.ElectrochemistryConditions.serializeBinaryToWriter
    );
  }
  f = message.getFlow();
  if (f != null) {
    writer.writeMessage(
      6,
      f,
      proto.ord.FlowConditions.serializeBinaryToWriter
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 7));
  if (f != null) {
    writer.writeBool(
      7,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 8));
  if (f != null) {
    writer.writeFloat(
      8,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 9));
  if (f != null) {
    writer.writeBool(
      9,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      10,
      f
    );
  }
};


/**
 * optional TemperatureConditions temperature = 1;
 * @return {?proto.ord.TemperatureConditions}
 */
proto.ord.ReactionConditions.prototype.getTemperature = function() {
  return /** @type{?proto.ord.TemperatureConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.TemperatureConditions, 1));
};


/**
 * @param {?proto.ord.TemperatureConditions|undefined} value
 * @return {!proto.ord.ReactionConditions} returns this
*/
proto.ord.ReactionConditions.prototype.setTemperature = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearTemperature = function() {
  return this.setTemperature(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasTemperature = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional PressureConditions pressure = 2;
 * @return {?proto.ord.PressureConditions}
 */
proto.ord.ReactionConditions.prototype.getPressure = function() {
  return /** @type{?proto.ord.PressureConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.PressureConditions, 2));
};


/**
 * @param {?proto.ord.PressureConditions|undefined} value
 * @return {!proto.ord.ReactionConditions} returns this
*/
proto.ord.ReactionConditions.prototype.setPressure = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearPressure = function() {
  return this.setPressure(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasPressure = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional StirringConditions stirring = 3;
 * @return {?proto.ord.StirringConditions}
 */
proto.ord.ReactionConditions.prototype.getStirring = function() {
  return /** @type{?proto.ord.StirringConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.StirringConditions, 3));
};


/**
 * @param {?proto.ord.StirringConditions|undefined} value
 * @return {!proto.ord.ReactionConditions} returns this
*/
proto.ord.ReactionConditions.prototype.setStirring = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearStirring = function() {
  return this.setStirring(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasStirring = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional IlluminationConditions illumination = 4;
 * @return {?proto.ord.IlluminationConditions}
 */
proto.ord.ReactionConditions.prototype.getIllumination = function() {
  return /** @type{?proto.ord.IlluminationConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.IlluminationConditions, 4));
};


/**
 * @param {?proto.ord.IlluminationConditions|undefined} value
 * @return {!proto.ord.ReactionConditions} returns this
*/
proto.ord.ReactionConditions.prototype.setIllumination = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearIllumination = function() {
  return this.setIllumination(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasIllumination = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional ElectrochemistryConditions electrochemistry = 5;
 * @return {?proto.ord.ElectrochemistryConditions}
 */
proto.ord.ReactionConditions.prototype.getElectrochemistry = function() {
  return /** @type{?proto.ord.ElectrochemistryConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.ElectrochemistryConditions, 5));
};


/**
 * @param {?proto.ord.ElectrochemistryConditions|undefined} value
 * @return {!proto.ord.ReactionConditions} returns this
*/
proto.ord.ReactionConditions.prototype.setElectrochemistry = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearElectrochemistry = function() {
  return this.setElectrochemistry(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasElectrochemistry = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional FlowConditions flow = 6;
 * @return {?proto.ord.FlowConditions}
 */
proto.ord.ReactionConditions.prototype.getFlow = function() {
  return /** @type{?proto.ord.FlowConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.FlowConditions, 6));
};


/**
 * @param {?proto.ord.FlowConditions|undefined} value
 * @return {!proto.ord.ReactionConditions} returns this
*/
proto.ord.ReactionConditions.prototype.setFlow = function(value) {
  return jspb.Message.setWrapperField(this, 6, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearFlow = function() {
  return this.setFlow(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasFlow = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * optional bool reflux = 7;
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.getReflux = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 7, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.setReflux = function(value) {
  return jspb.Message.setField(this, 7, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearReflux = function() {
  return jspb.Message.setField(this, 7, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasReflux = function() {
  return jspb.Message.getField(this, 7) != null;
};


/**
 * optional float ph = 8;
 * @return {number}
 */
proto.ord.ReactionConditions.prototype.getPh = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 8, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.setPh = function(value) {
  return jspb.Message.setField(this, 8, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearPh = function() {
  return jspb.Message.setField(this, 8, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasPh = function() {
  return jspb.Message.getField(this, 8) != null;
};


/**
 * optional bool conditions_are_dynamic = 9;
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.getConditionsAreDynamic = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 9, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.setConditionsAreDynamic = function(value) {
  return jspb.Message.setField(this, 9, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.clearConditionsAreDynamic = function() {
  return jspb.Message.setField(this, 9, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionConditions.prototype.hasConditionsAreDynamic = function() {
  return jspb.Message.getField(this, 9) != null;
};


/**
 * optional string details = 10;
 * @return {string}
 */
proto.ord.ReactionConditions.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 10, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionConditions} returns this
 */
proto.ord.ReactionConditions.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 10, value);
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.TemperatureConditions.repeatedFields_ = [3];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.TemperatureConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.TemperatureConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.TemperatureConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.TemperatureConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    control: (f = msg.getControl()) && proto.ord.TemperatureConditions.TemperatureControl.toObject(includeInstance, f),
    setpoint: (f = msg.getSetpoint()) && proto.ord.Temperature.toObject(includeInstance, f),
    measurementsList: jspb.Message.toObjectList(msg.getMeasurementsList(),
    proto.ord.TemperatureConditions.TemperatureMeasurement.toObject, includeInstance)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.TemperatureConditions}
 */
proto.ord.TemperatureConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.TemperatureConditions;
  return proto.ord.TemperatureConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.TemperatureConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.TemperatureConditions}
 */
proto.ord.TemperatureConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.TemperatureConditions.TemperatureControl;
      reader.readMessage(value,proto.ord.TemperatureConditions.TemperatureControl.deserializeBinaryFromReader);
      msg.setControl(value);
      break;
    case 2:
      var value = new proto.ord.Temperature;
      reader.readMessage(value,proto.ord.Temperature.deserializeBinaryFromReader);
      msg.setSetpoint(value);
      break;
    case 3:
      var value = new proto.ord.TemperatureConditions.TemperatureMeasurement;
      reader.readMessage(value,proto.ord.TemperatureConditions.TemperatureMeasurement.deserializeBinaryFromReader);
      msg.addMeasurements(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.TemperatureConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.TemperatureConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.TemperatureConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.TemperatureConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getControl();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.TemperatureConditions.TemperatureControl.serializeBinaryToWriter
    );
  }
  f = message.getSetpoint();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Temperature.serializeBinaryToWriter
    );
  }
  f = message.getMeasurementsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      3,
      f,
      proto.ord.TemperatureConditions.TemperatureMeasurement.serializeBinaryToWriter
    );
  }
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.TemperatureConditions.TemperatureControl.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.TemperatureConditions.TemperatureControl.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.TemperatureConditions.TemperatureControl} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.TemperatureConditions.TemperatureControl.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.TemperatureConditions.TemperatureControl}
 */
proto.ord.TemperatureConditions.TemperatureControl.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.TemperatureConditions.TemperatureControl;
  return proto.ord.TemperatureConditions.TemperatureControl.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.TemperatureConditions.TemperatureControl} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.TemperatureConditions.TemperatureControl}
 */
proto.ord.TemperatureConditions.TemperatureControl.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.TemperatureConditions.TemperatureControl.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.TemperatureConditions.TemperatureControl.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.TemperatureConditions.TemperatureControl} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.TemperatureConditions.TemperatureControl.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  AMBIENT: 2,
  OIL_BATH: 3,
  WATER_BATH: 4,
  SAND_BATH: 5,
  ICE_BATH: 6,
  DRY_ALUMINUM_PLATE: 7,
  MICROWAVE: 8,
  DRY_ICE_BATH: 9,
  AIR_FAN: 10,
  LIQUID_NITROGEN: 11
};

/**
 * optional TemperatureControlType type = 1;
 * @return {!proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType}
 */
proto.ord.TemperatureConditions.TemperatureControl.prototype.getType = function() {
  return /** @type {!proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType} value
 * @return {!proto.ord.TemperatureConditions.TemperatureControl} returns this
 */
proto.ord.TemperatureConditions.TemperatureControl.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.TemperatureConditions.TemperatureControl.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.TemperatureConditions.TemperatureControl} returns this
 */
proto.ord.TemperatureConditions.TemperatureControl.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.TemperatureConditions.TemperatureMeasurement.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.TemperatureConditions.TemperatureMeasurement} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    time: (f = msg.getTime()) && proto.ord.Time.toObject(includeInstance, f),
    temperature: (f = msg.getTemperature()) && proto.ord.Temperature.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.TemperatureConditions.TemperatureMeasurement;
  return proto.ord.TemperatureConditions.TemperatureMeasurement.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.TemperatureConditions.TemperatureMeasurement} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setTime(value);
      break;
    case 4:
      var value = new proto.ord.Temperature;
      reader.readMessage(value,proto.ord.Temperature.deserializeBinaryFromReader);
      msg.setTemperature(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.TemperatureConditions.TemperatureMeasurement.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.TemperatureConditions.TemperatureMeasurement} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getTime();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getTemperature();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.Temperature.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  THERMOCOUPLE_INTERNAL: 2,
  THERMOCOUPLE_EXTERNAL: 3,
  INFRARED: 4
};

/**
 * optional TemperatureMeasurementType type = 1;
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.getType = function() {
  return /** @type {!proto.ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.TemperatureConditions.TemperatureMeasurement.TemperatureMeasurementType} value
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement} returns this
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement} returns this
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Time time = 3;
 * @return {?proto.ord.Time}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.getTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 3));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement} returns this
*/
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.setTime = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement} returns this
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.clearTime = function() {
  return this.setTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.hasTime = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional Temperature temperature = 4;
 * @return {?proto.ord.Temperature}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.getTemperature = function() {
  return /** @type{?proto.ord.Temperature} */ (
    jspb.Message.getWrapperField(this, proto.ord.Temperature, 4));
};


/**
 * @param {?proto.ord.Temperature|undefined} value
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement} returns this
*/
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.setTemperature = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement} returns this
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.clearTemperature = function() {
  return this.setTemperature(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.TemperatureConditions.TemperatureMeasurement.prototype.hasTemperature = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional TemperatureControl control = 1;
 * @return {?proto.ord.TemperatureConditions.TemperatureControl}
 */
proto.ord.TemperatureConditions.prototype.getControl = function() {
  return /** @type{?proto.ord.TemperatureConditions.TemperatureControl} */ (
    jspb.Message.getWrapperField(this, proto.ord.TemperatureConditions.TemperatureControl, 1));
};


/**
 * @param {?proto.ord.TemperatureConditions.TemperatureControl|undefined} value
 * @return {!proto.ord.TemperatureConditions} returns this
*/
proto.ord.TemperatureConditions.prototype.setControl = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.TemperatureConditions} returns this
 */
proto.ord.TemperatureConditions.prototype.clearControl = function() {
  return this.setControl(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.TemperatureConditions.prototype.hasControl = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional Temperature setpoint = 2;
 * @return {?proto.ord.Temperature}
 */
proto.ord.TemperatureConditions.prototype.getSetpoint = function() {
  return /** @type{?proto.ord.Temperature} */ (
    jspb.Message.getWrapperField(this, proto.ord.Temperature, 2));
};


/**
 * @param {?proto.ord.Temperature|undefined} value
 * @return {!proto.ord.TemperatureConditions} returns this
*/
proto.ord.TemperatureConditions.prototype.setSetpoint = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.TemperatureConditions} returns this
 */
proto.ord.TemperatureConditions.prototype.clearSetpoint = function() {
  return this.setSetpoint(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.TemperatureConditions.prototype.hasSetpoint = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * repeated TemperatureMeasurement measurements = 3;
 * @return {!Array<!proto.ord.TemperatureConditions.TemperatureMeasurement>}
 */
proto.ord.TemperatureConditions.prototype.getMeasurementsList = function() {
  return /** @type{!Array<!proto.ord.TemperatureConditions.TemperatureMeasurement>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.TemperatureConditions.TemperatureMeasurement, 3));
};


/**
 * @param {!Array<!proto.ord.TemperatureConditions.TemperatureMeasurement>} value
 * @return {!proto.ord.TemperatureConditions} returns this
*/
proto.ord.TemperatureConditions.prototype.setMeasurementsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 3, value);
};


/**
 * @param {!proto.ord.TemperatureConditions.TemperatureMeasurement=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.TemperatureConditions.TemperatureMeasurement}
 */
proto.ord.TemperatureConditions.prototype.addMeasurements = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 3, opt_value, proto.ord.TemperatureConditions.TemperatureMeasurement, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.TemperatureConditions} returns this
 */
proto.ord.TemperatureConditions.prototype.clearMeasurementsList = function() {
  return this.setMeasurementsList([]);
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.PressureConditions.repeatedFields_ = [4];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.PressureConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.PressureConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.PressureConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    control: (f = msg.getControl()) && proto.ord.PressureConditions.PressureControl.toObject(includeInstance, f),
    setpoint: (f = msg.getSetpoint()) && proto.ord.Pressure.toObject(includeInstance, f),
    atmosphere: (f = msg.getAtmosphere()) && proto.ord.PressureConditions.Atmosphere.toObject(includeInstance, f),
    measurementsList: jspb.Message.toObjectList(msg.getMeasurementsList(),
    proto.ord.PressureConditions.PressureMeasurement.toObject, includeInstance)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.PressureConditions}
 */
proto.ord.PressureConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.PressureConditions;
  return proto.ord.PressureConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.PressureConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.PressureConditions}
 */
proto.ord.PressureConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.PressureConditions.PressureControl;
      reader.readMessage(value,proto.ord.PressureConditions.PressureControl.deserializeBinaryFromReader);
      msg.setControl(value);
      break;
    case 2:
      var value = new proto.ord.Pressure;
      reader.readMessage(value,proto.ord.Pressure.deserializeBinaryFromReader);
      msg.setSetpoint(value);
      break;
    case 3:
      var value = new proto.ord.PressureConditions.Atmosphere;
      reader.readMessage(value,proto.ord.PressureConditions.Atmosphere.deserializeBinaryFromReader);
      msg.setAtmosphere(value);
      break;
    case 4:
      var value = new proto.ord.PressureConditions.PressureMeasurement;
      reader.readMessage(value,proto.ord.PressureConditions.PressureMeasurement.deserializeBinaryFromReader);
      msg.addMeasurements(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.PressureConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.PressureConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.PressureConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getControl();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.PressureConditions.PressureControl.serializeBinaryToWriter
    );
  }
  f = message.getSetpoint();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Pressure.serializeBinaryToWriter
    );
  }
  f = message.getAtmosphere();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.PressureConditions.Atmosphere.serializeBinaryToWriter
    );
  }
  f = message.getMeasurementsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      4,
      f,
      proto.ord.PressureConditions.PressureMeasurement.serializeBinaryToWriter
    );
  }
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.PressureConditions.PressureControl.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.PressureConditions.PressureControl.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.PressureConditions.PressureControl} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.PressureControl.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.PressureConditions.PressureControl}
 */
proto.ord.PressureConditions.PressureControl.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.PressureConditions.PressureControl;
  return proto.ord.PressureConditions.PressureControl.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.PressureConditions.PressureControl} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.PressureConditions.PressureControl}
 */
proto.ord.PressureConditions.PressureControl.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.PressureConditions.PressureControl.PressureControlType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.PressureConditions.PressureControl.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.PressureConditions.PressureControl.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.PressureConditions.PressureControl} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.PressureControl.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.PressureConditions.PressureControl.PressureControlType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  AMBIENT: 2,
  SLIGHT_POSITIVE: 3,
  SEALED: 4,
  PRESSURIZED: 5
};

/**
 * optional PressureControlType type = 1;
 * @return {!proto.ord.PressureConditions.PressureControl.PressureControlType}
 */
proto.ord.PressureConditions.PressureControl.prototype.getType = function() {
  return /** @type {!proto.ord.PressureConditions.PressureControl.PressureControlType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.PressureConditions.PressureControl.PressureControlType} value
 * @return {!proto.ord.PressureConditions.PressureControl} returns this
 */
proto.ord.PressureConditions.PressureControl.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.PressureConditions.PressureControl.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.PressureConditions.PressureControl} returns this
 */
proto.ord.PressureConditions.PressureControl.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.PressureConditions.Atmosphere.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.PressureConditions.Atmosphere.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.PressureConditions.Atmosphere} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.Atmosphere.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.PressureConditions.Atmosphere}
 */
proto.ord.PressureConditions.Atmosphere.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.PressureConditions.Atmosphere;
  return proto.ord.PressureConditions.Atmosphere.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.PressureConditions.Atmosphere} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.PressureConditions.Atmosphere}
 */
proto.ord.PressureConditions.Atmosphere.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.PressureConditions.Atmosphere.AtmosphereType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.PressureConditions.Atmosphere.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.PressureConditions.Atmosphere.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.PressureConditions.Atmosphere} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.Atmosphere.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.PressureConditions.Atmosphere.AtmosphereType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  AIR: 2,
  NITROGEN: 3,
  ARGON: 4,
  OXYGEN: 5,
  HYDROGEN: 6,
  CARBON_MONOXIDE: 7,
  CARBON_DIOXIDE: 8,
  METHANE: 9,
  AMMONIA: 10,
  OZONE: 11,
  ETHYLENE: 12,
  ACETYLENE: 13
};

/**
 * optional AtmosphereType type = 1;
 * @return {!proto.ord.PressureConditions.Atmosphere.AtmosphereType}
 */
proto.ord.PressureConditions.Atmosphere.prototype.getType = function() {
  return /** @type {!proto.ord.PressureConditions.Atmosphere.AtmosphereType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.PressureConditions.Atmosphere.AtmosphereType} value
 * @return {!proto.ord.PressureConditions.Atmosphere} returns this
 */
proto.ord.PressureConditions.Atmosphere.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.PressureConditions.Atmosphere.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.PressureConditions.Atmosphere} returns this
 */
proto.ord.PressureConditions.Atmosphere.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.PressureConditions.PressureMeasurement.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.PressureConditions.PressureMeasurement} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.PressureMeasurement.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    time: (f = msg.getTime()) && proto.ord.Time.toObject(includeInstance, f),
    pressure: (f = msg.getPressure()) && proto.ord.Pressure.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.PressureConditions.PressureMeasurement}
 */
proto.ord.PressureConditions.PressureMeasurement.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.PressureConditions.PressureMeasurement;
  return proto.ord.PressureConditions.PressureMeasurement.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.PressureConditions.PressureMeasurement} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.PressureConditions.PressureMeasurement}
 */
proto.ord.PressureConditions.PressureMeasurement.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.PressureConditions.PressureMeasurement.PressureMeasurementType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setTime(value);
      break;
    case 4:
      var value = new proto.ord.Pressure;
      reader.readMessage(value,proto.ord.Pressure.deserializeBinaryFromReader);
      msg.setPressure(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.PressureConditions.PressureMeasurement.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.PressureConditions.PressureMeasurement} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.PressureConditions.PressureMeasurement.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getTime();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getPressure();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.Pressure.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.PressureConditions.PressureMeasurement.PressureMeasurementType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  PRESSURE_TRANSDUCER: 2
};

/**
 * optional PressureMeasurementType type = 1;
 * @return {!proto.ord.PressureConditions.PressureMeasurement.PressureMeasurementType}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.getType = function() {
  return /** @type {!proto.ord.PressureConditions.PressureMeasurement.PressureMeasurementType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.PressureConditions.PressureMeasurement.PressureMeasurementType} value
 * @return {!proto.ord.PressureConditions.PressureMeasurement} returns this
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.PressureConditions.PressureMeasurement} returns this
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Time time = 3;
 * @return {?proto.ord.Time}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.getTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 3));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.PressureConditions.PressureMeasurement} returns this
*/
proto.ord.PressureConditions.PressureMeasurement.prototype.setTime = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.PressureConditions.PressureMeasurement} returns this
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.clearTime = function() {
  return this.setTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.hasTime = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional Pressure pressure = 4;
 * @return {?proto.ord.Pressure}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.getPressure = function() {
  return /** @type{?proto.ord.Pressure} */ (
    jspb.Message.getWrapperField(this, proto.ord.Pressure, 4));
};


/**
 * @param {?proto.ord.Pressure|undefined} value
 * @return {!proto.ord.PressureConditions.PressureMeasurement} returns this
*/
proto.ord.PressureConditions.PressureMeasurement.prototype.setPressure = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.PressureConditions.PressureMeasurement} returns this
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.clearPressure = function() {
  return this.setPressure(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.PressureConditions.PressureMeasurement.prototype.hasPressure = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional PressureControl control = 1;
 * @return {?proto.ord.PressureConditions.PressureControl}
 */
proto.ord.PressureConditions.prototype.getControl = function() {
  return /** @type{?proto.ord.PressureConditions.PressureControl} */ (
    jspb.Message.getWrapperField(this, proto.ord.PressureConditions.PressureControl, 1));
};


/**
 * @param {?proto.ord.PressureConditions.PressureControl|undefined} value
 * @return {!proto.ord.PressureConditions} returns this
*/
proto.ord.PressureConditions.prototype.setControl = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.PressureConditions} returns this
 */
proto.ord.PressureConditions.prototype.clearControl = function() {
  return this.setControl(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.PressureConditions.prototype.hasControl = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional Pressure setpoint = 2;
 * @return {?proto.ord.Pressure}
 */
proto.ord.PressureConditions.prototype.getSetpoint = function() {
  return /** @type{?proto.ord.Pressure} */ (
    jspb.Message.getWrapperField(this, proto.ord.Pressure, 2));
};


/**
 * @param {?proto.ord.Pressure|undefined} value
 * @return {!proto.ord.PressureConditions} returns this
*/
proto.ord.PressureConditions.prototype.setSetpoint = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.PressureConditions} returns this
 */
proto.ord.PressureConditions.prototype.clearSetpoint = function() {
  return this.setSetpoint(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.PressureConditions.prototype.hasSetpoint = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional Atmosphere atmosphere = 3;
 * @return {?proto.ord.PressureConditions.Atmosphere}
 */
proto.ord.PressureConditions.prototype.getAtmosphere = function() {
  return /** @type{?proto.ord.PressureConditions.Atmosphere} */ (
    jspb.Message.getWrapperField(this, proto.ord.PressureConditions.Atmosphere, 3));
};


/**
 * @param {?proto.ord.PressureConditions.Atmosphere|undefined} value
 * @return {!proto.ord.PressureConditions} returns this
*/
proto.ord.PressureConditions.prototype.setAtmosphere = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.PressureConditions} returns this
 */
proto.ord.PressureConditions.prototype.clearAtmosphere = function() {
  return this.setAtmosphere(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.PressureConditions.prototype.hasAtmosphere = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * repeated PressureMeasurement measurements = 4;
 * @return {!Array<!proto.ord.PressureConditions.PressureMeasurement>}
 */
proto.ord.PressureConditions.prototype.getMeasurementsList = function() {
  return /** @type{!Array<!proto.ord.PressureConditions.PressureMeasurement>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.PressureConditions.PressureMeasurement, 4));
};


/**
 * @param {!Array<!proto.ord.PressureConditions.PressureMeasurement>} value
 * @return {!proto.ord.PressureConditions} returns this
*/
proto.ord.PressureConditions.prototype.setMeasurementsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 4, value);
};


/**
 * @param {!proto.ord.PressureConditions.PressureMeasurement=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.PressureConditions.PressureMeasurement}
 */
proto.ord.PressureConditions.prototype.addMeasurements = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 4, opt_value, proto.ord.PressureConditions.PressureMeasurement, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.PressureConditions} returns this
 */
proto.ord.PressureConditions.prototype.clearMeasurementsList = function() {
  return this.setMeasurementsList([]);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.StirringConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.StirringConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.StirringConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.StirringConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    rate: (f = msg.getRate()) && proto.ord.StirringConditions.StirringRate.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.StirringConditions}
 */
proto.ord.StirringConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.StirringConditions;
  return proto.ord.StirringConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.StirringConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.StirringConditions}
 */
proto.ord.StirringConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.StirringConditions.StirringMethodType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.StirringConditions.StirringRate;
      reader.readMessage(value,proto.ord.StirringConditions.StirringRate.deserializeBinaryFromReader);
      msg.setRate(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.StirringConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.StirringConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.StirringConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.StirringConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getRate();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.StirringConditions.StirringRate.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.StirringConditions.StirringMethodType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  NONE: 2,
  STIR_BAR: 3,
  OVERHEAD_MIXER: 4,
  AGITATION: 5,
  BALL_MILLING: 6,
  SONICATION: 7
};




if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.StirringConditions.StirringRate.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.StirringConditions.StirringRate.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.StirringConditions.StirringRate} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.StirringConditions.StirringRate.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    rpm: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.StirringConditions.StirringRate}
 */
proto.ord.StirringConditions.StirringRate.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.StirringConditions.StirringRate;
  return proto.ord.StirringConditions.StirringRate.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.StirringConditions.StirringRate} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.StirringConditions.StirringRate}
 */
proto.ord.StirringConditions.StirringRate.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.StirringConditions.StirringRate.StirringRateType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {number} */ (reader.readInt32());
      msg.setRpm(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.StirringConditions.StirringRate.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.StirringConditions.StirringRate.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.StirringConditions.StirringRate} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.StirringConditions.StirringRate.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getRpm();
  if (f !== 0) {
    writer.writeInt32(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.StirringConditions.StirringRate.StirringRateType = {
  UNSPECIFIED: 0,
  HIGH: 1,
  MEDIUM: 2,
  LOW: 3
};

/**
 * optional StirringRateType type = 1;
 * @return {!proto.ord.StirringConditions.StirringRate.StirringRateType}
 */
proto.ord.StirringConditions.StirringRate.prototype.getType = function() {
  return /** @type {!proto.ord.StirringConditions.StirringRate.StirringRateType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.StirringConditions.StirringRate.StirringRateType} value
 * @return {!proto.ord.StirringConditions.StirringRate} returns this
 */
proto.ord.StirringConditions.StirringRate.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.StirringConditions.StirringRate.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.StirringConditions.StirringRate} returns this
 */
proto.ord.StirringConditions.StirringRate.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional int32 rpm = 3;
 * @return {number}
 */
proto.ord.StirringConditions.StirringRate.prototype.getRpm = function() {
  return /** @type {number} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {number} value
 * @return {!proto.ord.StirringConditions.StirringRate} returns this
 */
proto.ord.StirringConditions.StirringRate.prototype.setRpm = function(value) {
  return jspb.Message.setProto3IntField(this, 3, value);
};


/**
 * optional StirringMethodType type = 1;
 * @return {!proto.ord.StirringConditions.StirringMethodType}
 */
proto.ord.StirringConditions.prototype.getType = function() {
  return /** @type {!proto.ord.StirringConditions.StirringMethodType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.StirringConditions.StirringMethodType} value
 * @return {!proto.ord.StirringConditions} returns this
 */
proto.ord.StirringConditions.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.StirringConditions.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.StirringConditions} returns this
 */
proto.ord.StirringConditions.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional StirringRate rate = 3;
 * @return {?proto.ord.StirringConditions.StirringRate}
 */
proto.ord.StirringConditions.prototype.getRate = function() {
  return /** @type{?proto.ord.StirringConditions.StirringRate} */ (
    jspb.Message.getWrapperField(this, proto.ord.StirringConditions.StirringRate, 3));
};


/**
 * @param {?proto.ord.StirringConditions.StirringRate|undefined} value
 * @return {!proto.ord.StirringConditions} returns this
*/
proto.ord.StirringConditions.prototype.setRate = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.StirringConditions} returns this
 */
proto.ord.StirringConditions.prototype.clearRate = function() {
  return this.setRate(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.StirringConditions.prototype.hasRate = function() {
  return jspb.Message.getField(this, 3) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.IlluminationConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.IlluminationConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.IlluminationConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.IlluminationConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    peakWavelength: (f = msg.getPeakWavelength()) && proto.ord.Wavelength.toObject(includeInstance, f),
    color: jspb.Message.getFieldWithDefault(msg, 4, ""),
    distanceToVessel: (f = msg.getDistanceToVessel()) && proto.ord.Length.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.IlluminationConditions}
 */
proto.ord.IlluminationConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.IlluminationConditions;
  return proto.ord.IlluminationConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.IlluminationConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.IlluminationConditions}
 */
proto.ord.IlluminationConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.IlluminationConditions.IlluminationType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.Wavelength;
      reader.readMessage(value,proto.ord.Wavelength.deserializeBinaryFromReader);
      msg.setPeakWavelength(value);
      break;
    case 4:
      var value = /** @type {string} */ (reader.readString());
      msg.setColor(value);
      break;
    case 5:
      var value = new proto.ord.Length;
      reader.readMessage(value,proto.ord.Length.deserializeBinaryFromReader);
      msg.setDistanceToVessel(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.IlluminationConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.IlluminationConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.IlluminationConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.IlluminationConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getPeakWavelength();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Wavelength.serializeBinaryToWriter
    );
  }
  f = message.getColor();
  if (f.length > 0) {
    writer.writeString(
      4,
      f
    );
  }
  f = message.getDistanceToVessel();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.Length.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.IlluminationConditions.IlluminationType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  AMBIENT: 2,
  DARK: 3,
  LED: 4,
  HALOGEN_LAMP: 5,
  DEUTERIUM_LAMP: 6,
  SOLAR_SIMULATOR: 7,
  BROAD_SPECTRUM: 8
};

/**
 * optional IlluminationType type = 1;
 * @return {!proto.ord.IlluminationConditions.IlluminationType}
 */
proto.ord.IlluminationConditions.prototype.getType = function() {
  return /** @type {!proto.ord.IlluminationConditions.IlluminationType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.IlluminationConditions.IlluminationType} value
 * @return {!proto.ord.IlluminationConditions} returns this
 */
proto.ord.IlluminationConditions.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.IlluminationConditions.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.IlluminationConditions} returns this
 */
proto.ord.IlluminationConditions.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Wavelength peak_wavelength = 3;
 * @return {?proto.ord.Wavelength}
 */
proto.ord.IlluminationConditions.prototype.getPeakWavelength = function() {
  return /** @type{?proto.ord.Wavelength} */ (
    jspb.Message.getWrapperField(this, proto.ord.Wavelength, 3));
};


/**
 * @param {?proto.ord.Wavelength|undefined} value
 * @return {!proto.ord.IlluminationConditions} returns this
*/
proto.ord.IlluminationConditions.prototype.setPeakWavelength = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.IlluminationConditions} returns this
 */
proto.ord.IlluminationConditions.prototype.clearPeakWavelength = function() {
  return this.setPeakWavelength(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.IlluminationConditions.prototype.hasPeakWavelength = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional string color = 4;
 * @return {string}
 */
proto.ord.IlluminationConditions.prototype.getColor = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 4, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.IlluminationConditions} returns this
 */
proto.ord.IlluminationConditions.prototype.setColor = function(value) {
  return jspb.Message.setProto3StringField(this, 4, value);
};


/**
 * optional Length distance_to_vessel = 5;
 * @return {?proto.ord.Length}
 */
proto.ord.IlluminationConditions.prototype.getDistanceToVessel = function() {
  return /** @type{?proto.ord.Length} */ (
    jspb.Message.getWrapperField(this, proto.ord.Length, 5));
};


/**
 * @param {?proto.ord.Length|undefined} value
 * @return {!proto.ord.IlluminationConditions} returns this
*/
proto.ord.IlluminationConditions.prototype.setDistanceToVessel = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.IlluminationConditions} returns this
 */
proto.ord.IlluminationConditions.prototype.clearDistanceToVessel = function() {
  return this.setDistanceToVessel(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.IlluminationConditions.prototype.hasDistanceToVessel = function() {
  return jspb.Message.getField(this, 5) != null;
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.ElectrochemistryConditions.repeatedFields_ = [8];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ElectrochemistryConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ElectrochemistryConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ElectrochemistryConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ElectrochemistryConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    current: (f = msg.getCurrent()) && proto.ord.Current.toObject(includeInstance, f),
    voltage: (f = msg.getVoltage()) && proto.ord.Voltage.toObject(includeInstance, f),
    anodeMaterial: jspb.Message.getFieldWithDefault(msg, 5, ""),
    cathodeMaterial: jspb.Message.getFieldWithDefault(msg, 6, ""),
    electrodeSeparation: (f = msg.getElectrodeSeparation()) && proto.ord.Length.toObject(includeInstance, f),
    measurementsList: jspb.Message.toObjectList(msg.getMeasurementsList(),
    proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.toObject, includeInstance),
    cell: (f = msg.getCell()) && proto.ord.ElectrochemistryConditions.ElectrochemistryCell.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ElectrochemistryConditions}
 */
proto.ord.ElectrochemistryConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ElectrochemistryConditions;
  return proto.ord.ElectrochemistryConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ElectrochemistryConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ElectrochemistryConditions}
 */
proto.ord.ElectrochemistryConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ElectrochemistryConditions.ElectrochemistryType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.Current;
      reader.readMessage(value,proto.ord.Current.deserializeBinaryFromReader);
      msg.setCurrent(value);
      break;
    case 4:
      var value = new proto.ord.Voltage;
      reader.readMessage(value,proto.ord.Voltage.deserializeBinaryFromReader);
      msg.setVoltage(value);
      break;
    case 5:
      var value = /** @type {string} */ (reader.readString());
      msg.setAnodeMaterial(value);
      break;
    case 6:
      var value = /** @type {string} */ (reader.readString());
      msg.setCathodeMaterial(value);
      break;
    case 7:
      var value = new proto.ord.Length;
      reader.readMessage(value,proto.ord.Length.deserializeBinaryFromReader);
      msg.setElectrodeSeparation(value);
      break;
    case 8:
      var value = new proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement;
      reader.readMessage(value,proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.deserializeBinaryFromReader);
      msg.addMeasurements(value);
      break;
    case 9:
      var value = new proto.ord.ElectrochemistryConditions.ElectrochemistryCell;
      reader.readMessage(value,proto.ord.ElectrochemistryConditions.ElectrochemistryCell.deserializeBinaryFromReader);
      msg.setCell(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ElectrochemistryConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ElectrochemistryConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ElectrochemistryConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ElectrochemistryConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getCurrent();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Current.serializeBinaryToWriter
    );
  }
  f = message.getVoltage();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.Voltage.serializeBinaryToWriter
    );
  }
  f = message.getAnodeMaterial();
  if (f.length > 0) {
    writer.writeString(
      5,
      f
    );
  }
  f = message.getCathodeMaterial();
  if (f.length > 0) {
    writer.writeString(
      6,
      f
    );
  }
  f = message.getElectrodeSeparation();
  if (f != null) {
    writer.writeMessage(
      7,
      f,
      proto.ord.Length.serializeBinaryToWriter
    );
  }
  f = message.getMeasurementsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      8,
      f,
      proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.serializeBinaryToWriter
    );
  }
  f = message.getCell();
  if (f != null) {
    writer.writeMessage(
      9,
      f,
      proto.ord.ElectrochemistryConditions.ElectrochemistryCell.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  CONSTANT_CURRENT: 2,
  CONSTANT_VOLTAGE: 3
};




if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.toObject = function(includeInstance, msg) {
  var f, obj = {
    time: (f = msg.getTime()) && proto.ord.Time.toObject(includeInstance, f),
    current: (f = msg.getCurrent()) && proto.ord.Current.toObject(includeInstance, f),
    voltage: (f = msg.getVoltage()) && proto.ord.Voltage.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement;
  return proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setTime(value);
      break;
    case 2:
      var value = new proto.ord.Current;
      reader.readMessage(value,proto.ord.Current.deserializeBinaryFromReader);
      msg.setCurrent(value);
      break;
    case 3:
      var value = new proto.ord.Voltage;
      reader.readMessage(value,proto.ord.Voltage.deserializeBinaryFromReader);
      msg.setVoltage(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getTime();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getCurrent();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Current.serializeBinaryToWriter
    );
  }
  f = message.getVoltage();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Voltage.serializeBinaryToWriter
    );
  }
};


/**
 * optional Time time = 1;
 * @return {?proto.ord.Time}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.getTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 1));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} returns this
*/
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.setTime = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} returns this
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.clearTime = function() {
  return this.setTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.hasTime = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional Current current = 2;
 * @return {?proto.ord.Current}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.getCurrent = function() {
  return /** @type{?proto.ord.Current} */ (
    jspb.Message.getWrapperField(this, proto.ord.Current, 2));
};


/**
 * @param {?proto.ord.Current|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} returns this
*/
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.setCurrent = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} returns this
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.clearCurrent = function() {
  return this.setCurrent(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.hasCurrent = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional Voltage voltage = 3;
 * @return {?proto.ord.Voltage}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.getVoltage = function() {
  return /** @type{?proto.ord.Voltage} */ (
    jspb.Message.getWrapperField(this, proto.ord.Voltage, 3));
};


/**
 * @param {?proto.ord.Voltage|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} returns this
*/
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.setVoltage = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement} returns this
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.clearVoltage = function() {
  return this.setVoltage(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement.prototype.hasVoltage = function() {
  return jspb.Message.getField(this, 3) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ElectrochemistryConditions.ElectrochemistryCell.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ElectrochemistryConditions.ElectrochemistryCell;
  return proto.ord.ElectrochemistryConditions.ElectrochemistryCell.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ElectrochemistryConditions.ElectrochemistryCell.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  DIVIDED_CELL: 2,
  UNDIVIDED_CELL: 3
};

/**
 * optional ElectrochemistryCellType type = 1;
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.prototype.getType = function() {
  return /** @type {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell.ElectrochemistryCellType} value
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell} returns this
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryCell} returns this
 */
proto.ord.ElectrochemistryConditions.ElectrochemistryCell.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional ElectrochemistryType type = 1;
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryType}
 */
proto.ord.ElectrochemistryConditions.prototype.getType = function() {
  return /** @type {!proto.ord.ElectrochemistryConditions.ElectrochemistryType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryType} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ElectrochemistryConditions.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Current current = 3;
 * @return {?proto.ord.Current}
 */
proto.ord.ElectrochemistryConditions.prototype.getCurrent = function() {
  return /** @type{?proto.ord.Current} */ (
    jspb.Message.getWrapperField(this, proto.ord.Current, 3));
};


/**
 * @param {?proto.ord.Current|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
*/
proto.ord.ElectrochemistryConditions.prototype.setCurrent = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.clearCurrent = function() {
  return this.setCurrent(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.prototype.hasCurrent = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional Voltage voltage = 4;
 * @return {?proto.ord.Voltage}
 */
proto.ord.ElectrochemistryConditions.prototype.getVoltage = function() {
  return /** @type{?proto.ord.Voltage} */ (
    jspb.Message.getWrapperField(this, proto.ord.Voltage, 4));
};


/**
 * @param {?proto.ord.Voltage|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
*/
proto.ord.ElectrochemistryConditions.prototype.setVoltage = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.clearVoltage = function() {
  return this.setVoltage(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.prototype.hasVoltage = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional string anode_material = 5;
 * @return {string}
 */
proto.ord.ElectrochemistryConditions.prototype.getAnodeMaterial = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 5, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.setAnodeMaterial = function(value) {
  return jspb.Message.setProto3StringField(this, 5, value);
};


/**
 * optional string cathode_material = 6;
 * @return {string}
 */
proto.ord.ElectrochemistryConditions.prototype.getCathodeMaterial = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 6, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.setCathodeMaterial = function(value) {
  return jspb.Message.setProto3StringField(this, 6, value);
};


/**
 * optional Length electrode_separation = 7;
 * @return {?proto.ord.Length}
 */
proto.ord.ElectrochemistryConditions.prototype.getElectrodeSeparation = function() {
  return /** @type{?proto.ord.Length} */ (
    jspb.Message.getWrapperField(this, proto.ord.Length, 7));
};


/**
 * @param {?proto.ord.Length|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
*/
proto.ord.ElectrochemistryConditions.prototype.setElectrodeSeparation = function(value) {
  return jspb.Message.setWrapperField(this, 7, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.clearElectrodeSeparation = function() {
  return this.setElectrodeSeparation(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.prototype.hasElectrodeSeparation = function() {
  return jspb.Message.getField(this, 7) != null;
};


/**
 * repeated ElectrochemistryMeasurement measurements = 8;
 * @return {!Array<!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement>}
 */
proto.ord.ElectrochemistryConditions.prototype.getMeasurementsList = function() {
  return /** @type{!Array<!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement, 8));
};


/**
 * @param {!Array<!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement>} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
*/
proto.ord.ElectrochemistryConditions.prototype.setMeasurementsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 8, value);
};


/**
 * @param {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement}
 */
proto.ord.ElectrochemistryConditions.prototype.addMeasurements = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 8, opt_value, proto.ord.ElectrochemistryConditions.ElectrochemistryMeasurement, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.clearMeasurementsList = function() {
  return this.setMeasurementsList([]);
};


/**
 * optional ElectrochemistryCell cell = 9;
 * @return {?proto.ord.ElectrochemistryConditions.ElectrochemistryCell}
 */
proto.ord.ElectrochemistryConditions.prototype.getCell = function() {
  return /** @type{?proto.ord.ElectrochemistryConditions.ElectrochemistryCell} */ (
    jspb.Message.getWrapperField(this, proto.ord.ElectrochemistryConditions.ElectrochemistryCell, 9));
};


/**
 * @param {?proto.ord.ElectrochemistryConditions.ElectrochemistryCell|undefined} value
 * @return {!proto.ord.ElectrochemistryConditions} returns this
*/
proto.ord.ElectrochemistryConditions.prototype.setCell = function(value) {
  return jspb.Message.setWrapperField(this, 9, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ElectrochemistryConditions} returns this
 */
proto.ord.ElectrochemistryConditions.prototype.clearCell = function() {
  return this.setCell(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ElectrochemistryConditions.prototype.hasCell = function() {
  return jspb.Message.getField(this, 9) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.FlowConditions.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.FlowConditions.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.FlowConditions} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FlowConditions.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    pumpType: jspb.Message.getFieldWithDefault(msg, 3, ""),
    tubing: (f = msg.getTubing()) && proto.ord.FlowConditions.Tubing.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.FlowConditions}
 */
proto.ord.FlowConditions.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.FlowConditions;
  return proto.ord.FlowConditions.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.FlowConditions} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.FlowConditions}
 */
proto.ord.FlowConditions.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.FlowConditions.FlowType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setPumpType(value);
      break;
    case 4:
      var value = new proto.ord.FlowConditions.Tubing;
      reader.readMessage(value,proto.ord.FlowConditions.Tubing.deserializeBinaryFromReader);
      msg.setTubing(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.FlowConditions.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.FlowConditions.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.FlowConditions} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FlowConditions.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getPumpType();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
  f = message.getTubing();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.FlowConditions.Tubing.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.FlowConditions.FlowType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  PLUG_FLOW_REACTOR: 2,
  CONTINUOUS_STIRRED_TANK_REACTOR: 3,
  PACKED_BED_REACTOR: 4
};




if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.FlowConditions.Tubing.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.FlowConditions.Tubing.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.FlowConditions.Tubing} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FlowConditions.Tubing.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    diameter: (f = msg.getDiameter()) && proto.ord.Length.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.FlowConditions.Tubing}
 */
proto.ord.FlowConditions.Tubing.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.FlowConditions.Tubing;
  return proto.ord.FlowConditions.Tubing.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.FlowConditions.Tubing} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.FlowConditions.Tubing}
 */
proto.ord.FlowConditions.Tubing.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.FlowConditions.Tubing.TubingType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.Length;
      reader.readMessage(value,proto.ord.Length.deserializeBinaryFromReader);
      msg.setDiameter(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.FlowConditions.Tubing.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.FlowConditions.Tubing.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.FlowConditions.Tubing} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FlowConditions.Tubing.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getDiameter();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Length.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.FlowConditions.Tubing.TubingType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  STEEL: 2,
  COPPER: 3,
  PFA: 4,
  FEP: 5,
  TEFLONAF: 6,
  PTFE: 7,
  GLASS: 8,
  QUARTZ: 9,
  SILICON: 10,
  PDMS: 11
};

/**
 * optional TubingType type = 1;
 * @return {!proto.ord.FlowConditions.Tubing.TubingType}
 */
proto.ord.FlowConditions.Tubing.prototype.getType = function() {
  return /** @type {!proto.ord.FlowConditions.Tubing.TubingType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.FlowConditions.Tubing.TubingType} value
 * @return {!proto.ord.FlowConditions.Tubing} returns this
 */
proto.ord.FlowConditions.Tubing.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.FlowConditions.Tubing.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.FlowConditions.Tubing} returns this
 */
proto.ord.FlowConditions.Tubing.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Length diameter = 3;
 * @return {?proto.ord.Length}
 */
proto.ord.FlowConditions.Tubing.prototype.getDiameter = function() {
  return /** @type{?proto.ord.Length} */ (
    jspb.Message.getWrapperField(this, proto.ord.Length, 3));
};


/**
 * @param {?proto.ord.Length|undefined} value
 * @return {!proto.ord.FlowConditions.Tubing} returns this
*/
proto.ord.FlowConditions.Tubing.prototype.setDiameter = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.FlowConditions.Tubing} returns this
 */
proto.ord.FlowConditions.Tubing.prototype.clearDiameter = function() {
  return this.setDiameter(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.FlowConditions.Tubing.prototype.hasDiameter = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional FlowType type = 1;
 * @return {!proto.ord.FlowConditions.FlowType}
 */
proto.ord.FlowConditions.prototype.getType = function() {
  return /** @type {!proto.ord.FlowConditions.FlowType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.FlowConditions.FlowType} value
 * @return {!proto.ord.FlowConditions} returns this
 */
proto.ord.FlowConditions.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.FlowConditions.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.FlowConditions} returns this
 */
proto.ord.FlowConditions.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string pump_type = 3;
 * @return {string}
 */
proto.ord.FlowConditions.prototype.getPumpType = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.FlowConditions} returns this
 */
proto.ord.FlowConditions.prototype.setPumpType = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};


/**
 * optional Tubing tubing = 4;
 * @return {?proto.ord.FlowConditions.Tubing}
 */
proto.ord.FlowConditions.prototype.getTubing = function() {
  return /** @type{?proto.ord.FlowConditions.Tubing} */ (
    jspb.Message.getWrapperField(this, proto.ord.FlowConditions.Tubing, 4));
};


/**
 * @param {?proto.ord.FlowConditions.Tubing|undefined} value
 * @return {!proto.ord.FlowConditions} returns this
*/
proto.ord.FlowConditions.prototype.setTubing = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.FlowConditions} returns this
 */
proto.ord.FlowConditions.prototype.clearTubing = function() {
  return this.setTubing(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.FlowConditions.prototype.hasTubing = function() {
  return jspb.Message.getField(this, 4) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionNotes.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionNotes.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionNotes} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionNotes.toObject = function(includeInstance, msg) {
  var f, obj = {
    isHeterogeneous: jspb.Message.getBooleanFieldWithDefault(msg, 1, false),
    formsPrecipitate: jspb.Message.getBooleanFieldWithDefault(msg, 2, false),
    isExothermic: jspb.Message.getBooleanFieldWithDefault(msg, 3, false),
    offgasses: jspb.Message.getBooleanFieldWithDefault(msg, 4, false),
    isSensitiveToMoisture: jspb.Message.getBooleanFieldWithDefault(msg, 5, false),
    isSensitiveToOxygen: jspb.Message.getBooleanFieldWithDefault(msg, 6, false),
    isSensitiveToLight: jspb.Message.getBooleanFieldWithDefault(msg, 7, false),
    safetyNotes: jspb.Message.getFieldWithDefault(msg, 8, ""),
    procedureDetails: jspb.Message.getFieldWithDefault(msg, 9, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionNotes}
 */
proto.ord.ReactionNotes.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionNotes;
  return proto.ord.ReactionNotes.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionNotes} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionNotes}
 */
proto.ord.ReactionNotes.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsHeterogeneous(value);
      break;
    case 2:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setFormsPrecipitate(value);
      break;
    case 3:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsExothermic(value);
      break;
    case 4:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setOffgasses(value);
      break;
    case 5:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsSensitiveToMoisture(value);
      break;
    case 6:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsSensitiveToOxygen(value);
      break;
    case 7:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsSensitiveToLight(value);
      break;
    case 8:
      var value = /** @type {string} */ (reader.readString());
      msg.setSafetyNotes(value);
      break;
    case 9:
      var value = /** @type {string} */ (reader.readString());
      msg.setProcedureDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionNotes.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionNotes.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionNotes} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionNotes.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {boolean} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeBool(
      1,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeBool(
      2,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 3));
  if (f != null) {
    writer.writeBool(
      3,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeBool(
      4,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 5));
  if (f != null) {
    writer.writeBool(
      5,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 6));
  if (f != null) {
    writer.writeBool(
      6,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 7));
  if (f != null) {
    writer.writeBool(
      7,
      f
    );
  }
  f = message.getSafetyNotes();
  if (f.length > 0) {
    writer.writeString(
      8,
      f
    );
  }
  f = message.getProcedureDetails();
  if (f.length > 0) {
    writer.writeString(
      9,
      f
    );
  }
};


/**
 * optional bool is_heterogeneous = 1;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getIsHeterogeneous = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 1, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setIsHeterogeneous = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearIsHeterogeneous = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasIsHeterogeneous = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional bool forms_precipitate = 2;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getFormsPrecipitate = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 2, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setFormsPrecipitate = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearFormsPrecipitate = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasFormsPrecipitate = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional bool is_exothermic = 3;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getIsExothermic = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 3, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setIsExothermic = function(value) {
  return jspb.Message.setField(this, 3, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearIsExothermic = function() {
  return jspb.Message.setField(this, 3, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasIsExothermic = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional bool offgasses = 4;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getOffgasses = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 4, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setOffgasses = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearOffgasses = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasOffgasses = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional bool is_sensitive_to_moisture = 5;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getIsSensitiveToMoisture = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 5, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setIsSensitiveToMoisture = function(value) {
  return jspb.Message.setField(this, 5, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearIsSensitiveToMoisture = function() {
  return jspb.Message.setField(this, 5, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasIsSensitiveToMoisture = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional bool is_sensitive_to_oxygen = 6;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getIsSensitiveToOxygen = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 6, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setIsSensitiveToOxygen = function(value) {
  return jspb.Message.setField(this, 6, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearIsSensitiveToOxygen = function() {
  return jspb.Message.setField(this, 6, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasIsSensitiveToOxygen = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * optional bool is_sensitive_to_light = 7;
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.getIsSensitiveToLight = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 7, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setIsSensitiveToLight = function(value) {
  return jspb.Message.setField(this, 7, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.clearIsSensitiveToLight = function() {
  return jspb.Message.setField(this, 7, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionNotes.prototype.hasIsSensitiveToLight = function() {
  return jspb.Message.getField(this, 7) != null;
};


/**
 * optional string safety_notes = 8;
 * @return {string}
 */
proto.ord.ReactionNotes.prototype.getSafetyNotes = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 8, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setSafetyNotes = function(value) {
  return jspb.Message.setProto3StringField(this, 8, value);
};


/**
 * optional string procedure_details = 9;
 * @return {string}
 */
proto.ord.ReactionNotes.prototype.getProcedureDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 9, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionNotes} returns this
 */
proto.ord.ReactionNotes.prototype.setProcedureDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 9, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionObservation.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionObservation.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionObservation} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionObservation.toObject = function(includeInstance, msg) {
  var f, obj = {
    time: (f = msg.getTime()) && proto.ord.Time.toObject(includeInstance, f),
    comment: jspb.Message.getFieldWithDefault(msg, 2, ""),
    image: (f = msg.getImage()) && proto.ord.Data.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionObservation}
 */
proto.ord.ReactionObservation.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionObservation;
  return proto.ord.ReactionObservation.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionObservation} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionObservation}
 */
proto.ord.ReactionObservation.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setTime(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setComment(value);
      break;
    case 3:
      var value = new proto.ord.Data;
      reader.readMessage(value,proto.ord.Data.deserializeBinaryFromReader);
      msg.setImage(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionObservation.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionObservation.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionObservation} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionObservation.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getTime();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getComment();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getImage();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Data.serializeBinaryToWriter
    );
  }
};


/**
 * optional Time time = 1;
 * @return {?proto.ord.Time}
 */
proto.ord.ReactionObservation.prototype.getTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 1));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ReactionObservation} returns this
*/
proto.ord.ReactionObservation.prototype.setTime = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionObservation} returns this
 */
proto.ord.ReactionObservation.prototype.clearTime = function() {
  return this.setTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionObservation.prototype.hasTime = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional string comment = 2;
 * @return {string}
 */
proto.ord.ReactionObservation.prototype.getComment = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionObservation} returns this
 */
proto.ord.ReactionObservation.prototype.setComment = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Data image = 3;
 * @return {?proto.ord.Data}
 */
proto.ord.ReactionObservation.prototype.getImage = function() {
  return /** @type{?proto.ord.Data} */ (
    jspb.Message.getWrapperField(this, proto.ord.Data, 3));
};


/**
 * @param {?proto.ord.Data|undefined} value
 * @return {!proto.ord.ReactionObservation} returns this
*/
proto.ord.ReactionObservation.prototype.setImage = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionObservation} returns this
 */
proto.ord.ReactionObservation.prototype.clearImage = function() {
  return this.setImage(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionObservation.prototype.hasImage = function() {
  return jspb.Message.getField(this, 3) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionWorkup.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionWorkup.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionWorkup} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionWorkup.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    duration: (f = msg.getDuration()) && proto.ord.Time.toObject(includeInstance, f),
    input: (f = msg.getInput()) && proto.ord.ReactionInput.toObject(includeInstance, f),
    amount: (f = msg.getAmount()) && proto.ord.Amount.toObject(includeInstance, f),
    temperature: (f = msg.getTemperature()) && proto.ord.TemperatureConditions.toObject(includeInstance, f),
    keepPhase: jspb.Message.getFieldWithDefault(msg, 7, ""),
    stirring: (f = msg.getStirring()) && proto.ord.StirringConditions.toObject(includeInstance, f),
    targetPh: jspb.Message.getFloatingPointFieldWithDefault(msg, 9, 0.0),
    isAutomated: jspb.Message.getBooleanFieldWithDefault(msg, 10, false)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionWorkup}
 */
proto.ord.ReactionWorkup.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionWorkup;
  return proto.ord.ReactionWorkup.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionWorkup} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionWorkup}
 */
proto.ord.ReactionWorkup.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ReactionWorkup.ReactionWorkupType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setDuration(value);
      break;
    case 4:
      var value = new proto.ord.ReactionInput;
      reader.readMessage(value,proto.ord.ReactionInput.deserializeBinaryFromReader);
      msg.setInput(value);
      break;
    case 5:
      var value = new proto.ord.Amount;
      reader.readMessage(value,proto.ord.Amount.deserializeBinaryFromReader);
      msg.setAmount(value);
      break;
    case 6:
      var value = new proto.ord.TemperatureConditions;
      reader.readMessage(value,proto.ord.TemperatureConditions.deserializeBinaryFromReader);
      msg.setTemperature(value);
      break;
    case 7:
      var value = /** @type {string} */ (reader.readString());
      msg.setKeepPhase(value);
      break;
    case 8:
      var value = new proto.ord.StirringConditions;
      reader.readMessage(value,proto.ord.StirringConditions.deserializeBinaryFromReader);
      msg.setStirring(value);
      break;
    case 9:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setTargetPh(value);
      break;
    case 10:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsAutomated(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionWorkup.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionWorkup.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionWorkup} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionWorkup.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getDuration();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getInput();
  if (f != null) {
    writer.writeMessage(
      4,
      f,
      proto.ord.ReactionInput.serializeBinaryToWriter
    );
  }
  f = message.getAmount();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.Amount.serializeBinaryToWriter
    );
  }
  f = message.getTemperature();
  if (f != null) {
    writer.writeMessage(
      6,
      f,
      proto.ord.TemperatureConditions.serializeBinaryToWriter
    );
  }
  f = message.getKeepPhase();
  if (f.length > 0) {
    writer.writeString(
      7,
      f
    );
  }
  f = message.getStirring();
  if (f != null) {
    writer.writeMessage(
      8,
      f,
      proto.ord.StirringConditions.serializeBinaryToWriter
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 9));
  if (f != null) {
    writer.writeFloat(
      9,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 10));
  if (f != null) {
    writer.writeBool(
      10,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ReactionWorkup.ReactionWorkupType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  ADDITION: 2,
  ALIQUOT: 3,
  TEMPERATURE: 4,
  CONCENTRATION: 5,
  EXTRACTION: 6,
  FILTRATION: 7,
  WASH: 8,
  DRY_IN_VACUUM: 9,
  DRY_WITH_MATERIAL: 10,
  FLASH_CHROMATOGRAPHY: 11,
  OTHER_CHROMATOGRAPHY: 12,
  SCAVENGING: 13,
  WAIT: 14,
  STIRRING: 15,
  PH_ADJUST: 16,
  DISSOLUTION: 17,
  DISTILLATION: 18
};

/**
 * optional ReactionWorkupType type = 1;
 * @return {!proto.ord.ReactionWorkup.ReactionWorkupType}
 */
proto.ord.ReactionWorkup.prototype.getType = function() {
  return /** @type {!proto.ord.ReactionWorkup.ReactionWorkupType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ReactionWorkup.ReactionWorkupType} value
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ReactionWorkup.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional Time duration = 3;
 * @return {?proto.ord.Time}
 */
proto.ord.ReactionWorkup.prototype.getDuration = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 3));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ReactionWorkup} returns this
*/
proto.ord.ReactionWorkup.prototype.setDuration = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearDuration = function() {
  return this.setDuration(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasDuration = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional ReactionInput input = 4;
 * @return {?proto.ord.ReactionInput}
 */
proto.ord.ReactionWorkup.prototype.getInput = function() {
  return /** @type{?proto.ord.ReactionInput} */ (
    jspb.Message.getWrapperField(this, proto.ord.ReactionInput, 4));
};


/**
 * @param {?proto.ord.ReactionInput|undefined} value
 * @return {!proto.ord.ReactionWorkup} returns this
*/
proto.ord.ReactionWorkup.prototype.setInput = function(value) {
  return jspb.Message.setWrapperField(this, 4, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearInput = function() {
  return this.setInput(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasInput = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional Amount amount = 5;
 * @return {?proto.ord.Amount}
 */
proto.ord.ReactionWorkup.prototype.getAmount = function() {
  return /** @type{?proto.ord.Amount} */ (
    jspb.Message.getWrapperField(this, proto.ord.Amount, 5));
};


/**
 * @param {?proto.ord.Amount|undefined} value
 * @return {!proto.ord.ReactionWorkup} returns this
*/
proto.ord.ReactionWorkup.prototype.setAmount = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearAmount = function() {
  return this.setAmount(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasAmount = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional TemperatureConditions temperature = 6;
 * @return {?proto.ord.TemperatureConditions}
 */
proto.ord.ReactionWorkup.prototype.getTemperature = function() {
  return /** @type{?proto.ord.TemperatureConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.TemperatureConditions, 6));
};


/**
 * @param {?proto.ord.TemperatureConditions|undefined} value
 * @return {!proto.ord.ReactionWorkup} returns this
*/
proto.ord.ReactionWorkup.prototype.setTemperature = function(value) {
  return jspb.Message.setWrapperField(this, 6, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearTemperature = function() {
  return this.setTemperature(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasTemperature = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * optional string keep_phase = 7;
 * @return {string}
 */
proto.ord.ReactionWorkup.prototype.getKeepPhase = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 7, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.setKeepPhase = function(value) {
  return jspb.Message.setProto3StringField(this, 7, value);
};


/**
 * optional StirringConditions stirring = 8;
 * @return {?proto.ord.StirringConditions}
 */
proto.ord.ReactionWorkup.prototype.getStirring = function() {
  return /** @type{?proto.ord.StirringConditions} */ (
    jspb.Message.getWrapperField(this, proto.ord.StirringConditions, 8));
};


/**
 * @param {?proto.ord.StirringConditions|undefined} value
 * @return {!proto.ord.ReactionWorkup} returns this
*/
proto.ord.ReactionWorkup.prototype.setStirring = function(value) {
  return jspb.Message.setWrapperField(this, 8, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearStirring = function() {
  return this.setStirring(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasStirring = function() {
  return jspb.Message.getField(this, 8) != null;
};


/**
 * optional float target_ph = 9;
 * @return {number}
 */
proto.ord.ReactionWorkup.prototype.getTargetPh = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 9, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.setTargetPh = function(value) {
  return jspb.Message.setField(this, 9, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearTargetPh = function() {
  return jspb.Message.setField(this, 9, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasTargetPh = function() {
  return jspb.Message.getField(this, 9) != null;
};


/**
 * optional bool is_automated = 10;
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.getIsAutomated = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 10, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.setIsAutomated = function(value) {
  return jspb.Message.setField(this, 10, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionWorkup} returns this
 */
proto.ord.ReactionWorkup.prototype.clearIsAutomated = function() {
  return jspb.Message.setField(this, 10, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionWorkup.prototype.hasIsAutomated = function() {
  return jspb.Message.getField(this, 10) != null;
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.ReactionOutcome.repeatedFields_ = [3];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionOutcome.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionOutcome.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionOutcome} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionOutcome.toObject = function(includeInstance, msg) {
  var f, obj = {
    reactionTime: (f = msg.getReactionTime()) && proto.ord.Time.toObject(includeInstance, f),
    conversion: (f = msg.getConversion()) && proto.ord.Percentage.toObject(includeInstance, f),
    productsList: jspb.Message.toObjectList(msg.getProductsList(),
    proto.ord.ProductCompound.toObject, includeInstance),
    analysesMap: (f = msg.getAnalysesMap()) ? f.toObject(includeInstance, proto.ord.Analysis.toObject) : []
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionOutcome}
 */
proto.ord.ReactionOutcome.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionOutcome;
  return proto.ord.ReactionOutcome.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionOutcome} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionOutcome}
 */
proto.ord.ReactionOutcome.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setReactionTime(value);
      break;
    case 2:
      var value = new proto.ord.Percentage;
      reader.readMessage(value,proto.ord.Percentage.deserializeBinaryFromReader);
      msg.setConversion(value);
      break;
    case 3:
      var value = new proto.ord.ProductCompound;
      reader.readMessage(value,proto.ord.ProductCompound.deserializeBinaryFromReader);
      msg.addProducts(value);
      break;
    case 4:
      var value = msg.getAnalysesMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Analysis.deserializeBinaryFromReader, "", new proto.ord.Analysis());
         });
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionOutcome.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionOutcome.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionOutcome} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionOutcome.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getReactionTime();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getConversion();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Percentage.serializeBinaryToWriter
    );
  }
  f = message.getProductsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      3,
      f,
      proto.ord.ProductCompound.serializeBinaryToWriter
    );
  }
  f = message.getAnalysesMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(4, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Analysis.serializeBinaryToWriter);
  }
};


/**
 * optional Time reaction_time = 1;
 * @return {?proto.ord.Time}
 */
proto.ord.ReactionOutcome.prototype.getReactionTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 1));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ReactionOutcome} returns this
*/
proto.ord.ReactionOutcome.prototype.setReactionTime = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionOutcome} returns this
 */
proto.ord.ReactionOutcome.prototype.clearReactionTime = function() {
  return this.setReactionTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionOutcome.prototype.hasReactionTime = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional Percentage conversion = 2;
 * @return {?proto.ord.Percentage}
 */
proto.ord.ReactionOutcome.prototype.getConversion = function() {
  return /** @type{?proto.ord.Percentage} */ (
    jspb.Message.getWrapperField(this, proto.ord.Percentage, 2));
};


/**
 * @param {?proto.ord.Percentage|undefined} value
 * @return {!proto.ord.ReactionOutcome} returns this
*/
proto.ord.ReactionOutcome.prototype.setConversion = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionOutcome} returns this
 */
proto.ord.ReactionOutcome.prototype.clearConversion = function() {
  return this.setConversion(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionOutcome.prototype.hasConversion = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * repeated ProductCompound products = 3;
 * @return {!Array<!proto.ord.ProductCompound>}
 */
proto.ord.ReactionOutcome.prototype.getProductsList = function() {
  return /** @type{!Array<!proto.ord.ProductCompound>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ProductCompound, 3));
};


/**
 * @param {!Array<!proto.ord.ProductCompound>} value
 * @return {!proto.ord.ReactionOutcome} returns this
*/
proto.ord.ReactionOutcome.prototype.setProductsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 3, value);
};


/**
 * @param {!proto.ord.ProductCompound=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ProductCompound}
 */
proto.ord.ReactionOutcome.prototype.addProducts = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 3, opt_value, proto.ord.ProductCompound, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ReactionOutcome} returns this
 */
proto.ord.ReactionOutcome.prototype.clearProductsList = function() {
  return this.setProductsList([]);
};


/**
 * map<string, Analysis> analyses = 4;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Analysis>}
 */
proto.ord.ReactionOutcome.prototype.getAnalysesMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Analysis>} */ (
      jspb.Message.getMapField(this, 4, opt_noLazyCreate,
      proto.ord.Analysis));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.ReactionOutcome} returns this
 */
proto.ord.ReactionOutcome.prototype.clearAnalysesMap = function() {
  this.getAnalysesMap().clear();
  return this;
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.ProductCompound.repeatedFields_ = [1,3];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ProductCompound.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ProductCompound.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ProductCompound} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductCompound.toObject = function(includeInstance, msg) {
  var f, obj = {
    identifiersList: jspb.Message.toObjectList(msg.getIdentifiersList(),
    proto.ord.CompoundIdentifier.toObject, includeInstance),
    isDesiredProduct: jspb.Message.getBooleanFieldWithDefault(msg, 2, false),
    measurementsList: jspb.Message.toObjectList(msg.getMeasurementsList(),
    proto.ord.ProductMeasurement.toObject, includeInstance),
    isolatedColor: jspb.Message.getFieldWithDefault(msg, 4, ""),
    texture: (f = msg.getTexture()) && proto.ord.ProductCompound.Texture.toObject(includeInstance, f),
    featuresMap: (f = msg.getFeaturesMap()) ? f.toObject(includeInstance, proto.ord.Data.toObject) : [],
    reactionRole: jspb.Message.getFieldWithDefault(msg, 7, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ProductCompound}
 */
proto.ord.ProductCompound.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ProductCompound;
  return proto.ord.ProductCompound.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ProductCompound} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ProductCompound}
 */
proto.ord.ProductCompound.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.CompoundIdentifier;
      reader.readMessage(value,proto.ord.CompoundIdentifier.deserializeBinaryFromReader);
      msg.addIdentifiers(value);
      break;
    case 2:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsDesiredProduct(value);
      break;
    case 3:
      var value = new proto.ord.ProductMeasurement;
      reader.readMessage(value,proto.ord.ProductMeasurement.deserializeBinaryFromReader);
      msg.addMeasurements(value);
      break;
    case 4:
      var value = /** @type {string} */ (reader.readString());
      msg.setIsolatedColor(value);
      break;
    case 5:
      var value = new proto.ord.ProductCompound.Texture;
      reader.readMessage(value,proto.ord.ProductCompound.Texture.deserializeBinaryFromReader);
      msg.setTexture(value);
      break;
    case 6:
      var value = msg.getFeaturesMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Data.deserializeBinaryFromReader, "", new proto.ord.Data());
         });
      break;
    case 7:
      var value = /** @type {!proto.ord.ReactionRole.ReactionRoleType} */ (reader.readEnum());
      msg.setReactionRole(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ProductCompound.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ProductCompound.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ProductCompound} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductCompound.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getIdentifiersList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      1,
      f,
      proto.ord.CompoundIdentifier.serializeBinaryToWriter
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeBool(
      2,
      f
    );
  }
  f = message.getMeasurementsList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      3,
      f,
      proto.ord.ProductMeasurement.serializeBinaryToWriter
    );
  }
  f = message.getIsolatedColor();
  if (f.length > 0) {
    writer.writeString(
      4,
      f
    );
  }
  f = message.getTexture();
  if (f != null) {
    writer.writeMessage(
      5,
      f,
      proto.ord.ProductCompound.Texture.serializeBinaryToWriter
    );
  }
  f = message.getFeaturesMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(6, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Data.serializeBinaryToWriter);
  }
  f = message.getReactionRole();
  if (f !== 0.0) {
    writer.writeEnum(
      7,
      f
    );
  }
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ProductCompound.Texture.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ProductCompound.Texture.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ProductCompound.Texture} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductCompound.Texture.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ProductCompound.Texture}
 */
proto.ord.ProductCompound.Texture.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ProductCompound.Texture;
  return proto.ord.ProductCompound.Texture.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ProductCompound.Texture} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ProductCompound.Texture}
 */
proto.ord.ProductCompound.Texture.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ProductCompound.Texture.TextureType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ProductCompound.Texture.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ProductCompound.Texture.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ProductCompound.Texture} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductCompound.Texture.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ProductCompound.Texture.TextureType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  POWDER: 2,
  CRYSTAL: 3,
  OIL: 4,
  AMORPHOUS_SOLID: 5,
  FOAM: 6,
  WAX: 7,
  SEMI_SOLID: 8,
  SOLID: 9,
  LIQUID: 10
};

/**
 * optional TextureType type = 1;
 * @return {!proto.ord.ProductCompound.Texture.TextureType}
 */
proto.ord.ProductCompound.Texture.prototype.getType = function() {
  return /** @type {!proto.ord.ProductCompound.Texture.TextureType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ProductCompound.Texture.TextureType} value
 * @return {!proto.ord.ProductCompound.Texture} returns this
 */
proto.ord.ProductCompound.Texture.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ProductCompound.Texture.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductCompound.Texture} returns this
 */
proto.ord.ProductCompound.Texture.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * repeated CompoundIdentifier identifiers = 1;
 * @return {!Array<!proto.ord.CompoundIdentifier>}
 */
proto.ord.ProductCompound.prototype.getIdentifiersList = function() {
  return /** @type{!Array<!proto.ord.CompoundIdentifier>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.CompoundIdentifier, 1));
};


/**
 * @param {!Array<!proto.ord.CompoundIdentifier>} value
 * @return {!proto.ord.ProductCompound} returns this
*/
proto.ord.ProductCompound.prototype.setIdentifiersList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 1, value);
};


/**
 * @param {!proto.ord.CompoundIdentifier=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.CompoundIdentifier}
 */
proto.ord.ProductCompound.prototype.addIdentifiers = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 1, opt_value, proto.ord.CompoundIdentifier, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.clearIdentifiersList = function() {
  return this.setIdentifiersList([]);
};


/**
 * optional bool is_desired_product = 2;
 * @return {boolean}
 */
proto.ord.ProductCompound.prototype.getIsDesiredProduct = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 2, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.setIsDesiredProduct = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.clearIsDesiredProduct = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductCompound.prototype.hasIsDesiredProduct = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * repeated ProductMeasurement measurements = 3;
 * @return {!Array<!proto.ord.ProductMeasurement>}
 */
proto.ord.ProductCompound.prototype.getMeasurementsList = function() {
  return /** @type{!Array<!proto.ord.ProductMeasurement>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.ProductMeasurement, 3));
};


/**
 * @param {!Array<!proto.ord.ProductMeasurement>} value
 * @return {!proto.ord.ProductCompound} returns this
*/
proto.ord.ProductCompound.prototype.setMeasurementsList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 3, value);
};


/**
 * @param {!proto.ord.ProductMeasurement=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.ProductMeasurement}
 */
proto.ord.ProductCompound.prototype.addMeasurements = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 3, opt_value, proto.ord.ProductMeasurement, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.clearMeasurementsList = function() {
  return this.setMeasurementsList([]);
};


/**
 * optional string isolated_color = 4;
 * @return {string}
 */
proto.ord.ProductCompound.prototype.getIsolatedColor = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 4, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.setIsolatedColor = function(value) {
  return jspb.Message.setProto3StringField(this, 4, value);
};


/**
 * optional Texture texture = 5;
 * @return {?proto.ord.ProductCompound.Texture}
 */
proto.ord.ProductCompound.prototype.getTexture = function() {
  return /** @type{?proto.ord.ProductCompound.Texture} */ (
    jspb.Message.getWrapperField(this, proto.ord.ProductCompound.Texture, 5));
};


/**
 * @param {?proto.ord.ProductCompound.Texture|undefined} value
 * @return {!proto.ord.ProductCompound} returns this
*/
proto.ord.ProductCompound.prototype.setTexture = function(value) {
  return jspb.Message.setWrapperField(this, 5, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.clearTexture = function() {
  return this.setTexture(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductCompound.prototype.hasTexture = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * map<string, Data> features = 6;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Data>}
 */
proto.ord.ProductCompound.prototype.getFeaturesMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Data>} */ (
      jspb.Message.getMapField(this, 6, opt_noLazyCreate,
      proto.ord.Data));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.clearFeaturesMap = function() {
  this.getFeaturesMap().clear();
  return this;
};


/**
 * optional ReactionRole.ReactionRoleType reaction_role = 7;
 * @return {!proto.ord.ReactionRole.ReactionRoleType}
 */
proto.ord.ProductCompound.prototype.getReactionRole = function() {
  return /** @type {!proto.ord.ReactionRole.ReactionRoleType} */ (jspb.Message.getFieldWithDefault(this, 7, 0));
};


/**
 * @param {!proto.ord.ReactionRole.ReactionRoleType} value
 * @return {!proto.ord.ProductCompound} returns this
 */
proto.ord.ProductCompound.prototype.setReactionRole = function(value) {
  return jspb.Message.setProto3EnumField(this, 7, value);
};



/**
 * Oneof group definitions for this message. Each group defines the field
 * numbers belonging to that group. When of these fields' value is set, all
 * other fields in the group are cleared. During deserialization, if multiple
 * fields are encountered for a group, only the last value seen will be kept.
 * @private {!Array<!Array<number>>}
 * @const
 */
proto.ord.ProductMeasurement.oneofGroups_ = [[8,9,10,11]];

/**
 * @enum {number}
 */
proto.ord.ProductMeasurement.ValueCase = {
  VALUE_NOT_SET: 0,
  PERCENTAGE: 8,
  FLOAT_VALUE: 9,
  STRING_VALUE: 10,
  AMOUNT: 11
};

/**
 * @return {proto.ord.ProductMeasurement.ValueCase}
 */
proto.ord.ProductMeasurement.prototype.getValueCase = function() {
  return /** @type {proto.ord.ProductMeasurement.ValueCase} */(jspb.Message.computeOneofCase(this, proto.ord.ProductMeasurement.oneofGroups_[0]));
};



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ProductMeasurement.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ProductMeasurement.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ProductMeasurement} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductMeasurement.toObject = function(includeInstance, msg) {
  var f, obj = {
    analysisKey: jspb.Message.getFieldWithDefault(msg, 1, ""),
    type: jspb.Message.getFieldWithDefault(msg, 2, 0),
    details: jspb.Message.getFieldWithDefault(msg, 3, ""),
    usesInternalStandard: jspb.Message.getBooleanFieldWithDefault(msg, 4, false),
    isNormalized: jspb.Message.getBooleanFieldWithDefault(msg, 5, false),
    usesAuthenticStandard: jspb.Message.getBooleanFieldWithDefault(msg, 6, false),
    authenticStandard: (f = msg.getAuthenticStandard()) && proto.ord.Compound.toObject(includeInstance, f),
    percentage: (f = msg.getPercentage()) && proto.ord.Percentage.toObject(includeInstance, f),
    floatValue: (f = msg.getFloatValue()) && proto.ord.FloatValue.toObject(includeInstance, f),
    stringValue: jspb.Message.getFieldWithDefault(msg, 10, ""),
    amount: (f = msg.getAmount()) && proto.ord.Amount.toObject(includeInstance, f),
    retentionTime: (f = msg.getRetentionTime()) && proto.ord.Time.toObject(includeInstance, f),
    massSpecDetails: (f = msg.getMassSpecDetails()) && proto.ord.ProductMeasurement.MassSpecMeasurementDetails.toObject(includeInstance, f),
    selectivity: (f = msg.getSelectivity()) && proto.ord.ProductMeasurement.Selectivity.toObject(includeInstance, f),
    wavelength: (f = msg.getWavelength()) && proto.ord.Wavelength.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ProductMeasurement}
 */
proto.ord.ProductMeasurement.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ProductMeasurement;
  return proto.ord.ProductMeasurement.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ProductMeasurement} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ProductMeasurement}
 */
proto.ord.ProductMeasurement.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {string} */ (reader.readString());
      msg.setAnalysisKey(value);
      break;
    case 2:
      var value = /** @type {!proto.ord.ProductMeasurement.ProductMeasurementType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 4:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setUsesInternalStandard(value);
      break;
    case 5:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsNormalized(value);
      break;
    case 6:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setUsesAuthenticStandard(value);
      break;
    case 7:
      var value = new proto.ord.Compound;
      reader.readMessage(value,proto.ord.Compound.deserializeBinaryFromReader);
      msg.setAuthenticStandard(value);
      break;
    case 8:
      var value = new proto.ord.Percentage;
      reader.readMessage(value,proto.ord.Percentage.deserializeBinaryFromReader);
      msg.setPercentage(value);
      break;
    case 9:
      var value = new proto.ord.FloatValue;
      reader.readMessage(value,proto.ord.FloatValue.deserializeBinaryFromReader);
      msg.setFloatValue(value);
      break;
    case 10:
      var value = /** @type {string} */ (reader.readString());
      msg.setStringValue(value);
      break;
    case 11:
      var value = new proto.ord.Amount;
      reader.readMessage(value,proto.ord.Amount.deserializeBinaryFromReader);
      msg.setAmount(value);
      break;
    case 12:
      var value = new proto.ord.Time;
      reader.readMessage(value,proto.ord.Time.deserializeBinaryFromReader);
      msg.setRetentionTime(value);
      break;
    case 13:
      var value = new proto.ord.ProductMeasurement.MassSpecMeasurementDetails;
      reader.readMessage(value,proto.ord.ProductMeasurement.MassSpecMeasurementDetails.deserializeBinaryFromReader);
      msg.setMassSpecDetails(value);
      break;
    case 14:
      var value = new proto.ord.ProductMeasurement.Selectivity;
      reader.readMessage(value,proto.ord.ProductMeasurement.Selectivity.deserializeBinaryFromReader);
      msg.setSelectivity(value);
      break;
    case 15:
      var value = new proto.ord.Wavelength;
      reader.readMessage(value,proto.ord.Wavelength.deserializeBinaryFromReader);
      msg.setWavelength(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ProductMeasurement.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ProductMeasurement.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ProductMeasurement} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductMeasurement.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getAnalysisKey();
  if (f.length > 0) {
    writer.writeString(
      1,
      f
    );
  }
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      2,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeBool(
      4,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 5));
  if (f != null) {
    writer.writeBool(
      5,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 6));
  if (f != null) {
    writer.writeBool(
      6,
      f
    );
  }
  f = message.getAuthenticStandard();
  if (f != null) {
    writer.writeMessage(
      7,
      f,
      proto.ord.Compound.serializeBinaryToWriter
    );
  }
  f = message.getPercentage();
  if (f != null) {
    writer.writeMessage(
      8,
      f,
      proto.ord.Percentage.serializeBinaryToWriter
    );
  }
  f = message.getFloatValue();
  if (f != null) {
    writer.writeMessage(
      9,
      f,
      proto.ord.FloatValue.serializeBinaryToWriter
    );
  }
  f = /** @type {string} */ (jspb.Message.getField(message, 10));
  if (f != null) {
    writer.writeString(
      10,
      f
    );
  }
  f = message.getAmount();
  if (f != null) {
    writer.writeMessage(
      11,
      f,
      proto.ord.Amount.serializeBinaryToWriter
    );
  }
  f = message.getRetentionTime();
  if (f != null) {
    writer.writeMessage(
      12,
      f,
      proto.ord.Time.serializeBinaryToWriter
    );
  }
  f = message.getMassSpecDetails();
  if (f != null) {
    writer.writeMessage(
      13,
      f,
      proto.ord.ProductMeasurement.MassSpecMeasurementDetails.serializeBinaryToWriter
    );
  }
  f = message.getSelectivity();
  if (f != null) {
    writer.writeMessage(
      14,
      f,
      proto.ord.ProductMeasurement.Selectivity.serializeBinaryToWriter
    );
  }
  f = message.getWavelength();
  if (f != null) {
    writer.writeMessage(
      15,
      f,
      proto.ord.Wavelength.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ProductMeasurement.ProductMeasurementType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  IDENTITY: 2,
  YIELD: 3,
  SELECTIVITY: 4,
  PURITY: 5,
  AREA: 6,
  COUNTS: 7,
  INTENSITY: 8,
  AMOUNT: 9
};


/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.repeatedFields_ = [5];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ProductMeasurement.MassSpecMeasurementDetails.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    ticMinimumMz: jspb.Message.getFloatingPointFieldWithDefault(msg, 3, 0.0),
    ticMaximumMz: jspb.Message.getFloatingPointFieldWithDefault(msg, 4, 0.0),
    eicMassesList: (f = jspb.Message.getRepeatedFloatingPointField(msg, 5)) == null ? undefined : f
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ProductMeasurement.MassSpecMeasurementDetails;
  return proto.ord.ProductMeasurement.MassSpecMeasurementDetails.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setTicMinimumMz(value);
      break;
    case 4:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setTicMaximumMz(value);
      break;
    case 5:
      var values = /** @type {!Array<number>} */ (reader.isDelimited() ? reader.readPackedFloat() : [reader.readFloat()]);
      for (var i = 0; i < values.length; i++) {
        msg.addEicMasses(values[i]);
      }
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ProductMeasurement.MassSpecMeasurementDetails.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 3));
  if (f != null) {
    writer.writeFloat(
      3,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeFloat(
      4,
      f
    );
  }
  f = message.getEicMassesList();
  if (f.length > 0) {
    writer.writePackedFloat(
      5,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  TIC: 2,
  TIC_POSITIVE: 3,
  TIC_NEGATIVE: 4,
  EIC: 5
};

/**
 * optional MassSpecMeasurementType type = 1;
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.getType = function() {
  return /** @type {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType} value
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional float tic_minimum_mz = 3;
 * @return {number}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.getTicMinimumMz = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 3, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.setTicMinimumMz = function(value) {
  return jspb.Message.setField(this, 3, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.clearTicMinimumMz = function() {
  return jspb.Message.setField(this, 3, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.hasTicMinimumMz = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional float tic_maximum_mz = 4;
 * @return {number}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.getTicMaximumMz = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 4, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.setTicMaximumMz = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.clearTicMaximumMz = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.hasTicMaximumMz = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * repeated float eic_masses = 5;
 * @return {!Array<number>}
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.getEicMassesList = function() {
  return /** @type {!Array<number>} */ (jspb.Message.getRepeatedFloatingPointField(this, 5));
};


/**
 * @param {!Array<number>} value
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.setEicMassesList = function(value) {
  return jspb.Message.setField(this, 5, value || []);
};


/**
 * @param {number} value
 * @param {number=} opt_index
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.addEicMasses = function(value, opt_index) {
  return jspb.Message.addToRepeatedField(this, 5, value, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ProductMeasurement.MassSpecMeasurementDetails} returns this
 */
proto.ord.ProductMeasurement.MassSpecMeasurementDetails.prototype.clearEicMassesList = function() {
  return this.setEicMassesList([]);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ProductMeasurement.Selectivity.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ProductMeasurement.Selectivity.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ProductMeasurement.Selectivity} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductMeasurement.Selectivity.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ProductMeasurement.Selectivity}
 */
proto.ord.ProductMeasurement.Selectivity.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ProductMeasurement.Selectivity;
  return proto.ord.ProductMeasurement.Selectivity.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ProductMeasurement.Selectivity} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ProductMeasurement.Selectivity}
 */
proto.ord.ProductMeasurement.Selectivity.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.ProductMeasurement.Selectivity.SelectivityType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ProductMeasurement.Selectivity.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ProductMeasurement.Selectivity.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ProductMeasurement.Selectivity} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ProductMeasurement.Selectivity.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.ProductMeasurement.Selectivity.SelectivityType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  EE: 2,
  ER: 3,
  DR: 4,
  EZ: 5,
  ZE: 6
};

/**
 * optional SelectivityType type = 1;
 * @return {!proto.ord.ProductMeasurement.Selectivity.SelectivityType}
 */
proto.ord.ProductMeasurement.Selectivity.prototype.getType = function() {
  return /** @type {!proto.ord.ProductMeasurement.Selectivity.SelectivityType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.ProductMeasurement.Selectivity.SelectivityType} value
 * @return {!proto.ord.ProductMeasurement.Selectivity} returns this
 */
proto.ord.ProductMeasurement.Selectivity.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.ProductMeasurement.Selectivity.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductMeasurement.Selectivity} returns this
 */
proto.ord.ProductMeasurement.Selectivity.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string analysis_key = 1;
 * @return {string}
 */
proto.ord.ProductMeasurement.prototype.getAnalysisKey = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 1, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setAnalysisKey = function(value) {
  return jspb.Message.setProto3StringField(this, 1, value);
};


/**
 * optional ProductMeasurementType type = 2;
 * @return {!proto.ord.ProductMeasurement.ProductMeasurementType}
 */
proto.ord.ProductMeasurement.prototype.getType = function() {
  return /** @type {!proto.ord.ProductMeasurement.ProductMeasurementType} */ (jspb.Message.getFieldWithDefault(this, 2, 0));
};


/**
 * @param {!proto.ord.ProductMeasurement.ProductMeasurementType} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 2, value);
};


/**
 * optional string details = 3;
 * @return {string}
 */
proto.ord.ProductMeasurement.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};


/**
 * optional bool uses_internal_standard = 4;
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.getUsesInternalStandard = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 4, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setUsesInternalStandard = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearUsesInternalStandard = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasUsesInternalStandard = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional bool is_normalized = 5;
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.getIsNormalized = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 5, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setIsNormalized = function(value) {
  return jspb.Message.setField(this, 5, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearIsNormalized = function() {
  return jspb.Message.setField(this, 5, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasIsNormalized = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional bool uses_authentic_standard = 6;
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.getUsesAuthenticStandard = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 6, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setUsesAuthenticStandard = function(value) {
  return jspb.Message.setField(this, 6, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearUsesAuthenticStandard = function() {
  return jspb.Message.setField(this, 6, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasUsesAuthenticStandard = function() {
  return jspb.Message.getField(this, 6) != null;
};


/**
 * optional Compound authentic_standard = 7;
 * @return {?proto.ord.Compound}
 */
proto.ord.ProductMeasurement.prototype.getAuthenticStandard = function() {
  return /** @type{?proto.ord.Compound} */ (
    jspb.Message.getWrapperField(this, proto.ord.Compound, 7));
};


/**
 * @param {?proto.ord.Compound|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setAuthenticStandard = function(value) {
  return jspb.Message.setWrapperField(this, 7, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearAuthenticStandard = function() {
  return this.setAuthenticStandard(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasAuthenticStandard = function() {
  return jspb.Message.getField(this, 7) != null;
};


/**
 * optional Percentage percentage = 8;
 * @return {?proto.ord.Percentage}
 */
proto.ord.ProductMeasurement.prototype.getPercentage = function() {
  return /** @type{?proto.ord.Percentage} */ (
    jspb.Message.getWrapperField(this, proto.ord.Percentage, 8));
};


/**
 * @param {?proto.ord.Percentage|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setPercentage = function(value) {
  return jspb.Message.setOneofWrapperField(this, 8, proto.ord.ProductMeasurement.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearPercentage = function() {
  return this.setPercentage(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasPercentage = function() {
  return jspb.Message.getField(this, 8) != null;
};


/**
 * optional FloatValue float_value = 9;
 * @return {?proto.ord.FloatValue}
 */
proto.ord.ProductMeasurement.prototype.getFloatValue = function() {
  return /** @type{?proto.ord.FloatValue} */ (
    jspb.Message.getWrapperField(this, proto.ord.FloatValue, 9));
};


/**
 * @param {?proto.ord.FloatValue|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setFloatValue = function(value) {
  return jspb.Message.setOneofWrapperField(this, 9, proto.ord.ProductMeasurement.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearFloatValue = function() {
  return this.setFloatValue(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasFloatValue = function() {
  return jspb.Message.getField(this, 9) != null;
};


/**
 * optional string string_value = 10;
 * @return {string}
 */
proto.ord.ProductMeasurement.prototype.getStringValue = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 10, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.setStringValue = function(value) {
  return jspb.Message.setOneofField(this, 10, proto.ord.ProductMeasurement.oneofGroups_[0], value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearStringValue = function() {
  return jspb.Message.setOneofField(this, 10, proto.ord.ProductMeasurement.oneofGroups_[0], undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasStringValue = function() {
  return jspb.Message.getField(this, 10) != null;
};


/**
 * optional Amount amount = 11;
 * @return {?proto.ord.Amount}
 */
proto.ord.ProductMeasurement.prototype.getAmount = function() {
  return /** @type{?proto.ord.Amount} */ (
    jspb.Message.getWrapperField(this, proto.ord.Amount, 11));
};


/**
 * @param {?proto.ord.Amount|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setAmount = function(value) {
  return jspb.Message.setOneofWrapperField(this, 11, proto.ord.ProductMeasurement.oneofGroups_[0], value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearAmount = function() {
  return this.setAmount(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasAmount = function() {
  return jspb.Message.getField(this, 11) != null;
};


/**
 * optional Time retention_time = 12;
 * @return {?proto.ord.Time}
 */
proto.ord.ProductMeasurement.prototype.getRetentionTime = function() {
  return /** @type{?proto.ord.Time} */ (
    jspb.Message.getWrapperField(this, proto.ord.Time, 12));
};


/**
 * @param {?proto.ord.Time|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setRetentionTime = function(value) {
  return jspb.Message.setWrapperField(this, 12, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearRetentionTime = function() {
  return this.setRetentionTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasRetentionTime = function() {
  return jspb.Message.getField(this, 12) != null;
};


/**
 * optional MassSpecMeasurementDetails mass_spec_details = 13;
 * @return {?proto.ord.ProductMeasurement.MassSpecMeasurementDetails}
 */
proto.ord.ProductMeasurement.prototype.getMassSpecDetails = function() {
  return /** @type{?proto.ord.ProductMeasurement.MassSpecMeasurementDetails} */ (
    jspb.Message.getWrapperField(this, proto.ord.ProductMeasurement.MassSpecMeasurementDetails, 13));
};


/**
 * @param {?proto.ord.ProductMeasurement.MassSpecMeasurementDetails|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setMassSpecDetails = function(value) {
  return jspb.Message.setWrapperField(this, 13, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearMassSpecDetails = function() {
  return this.setMassSpecDetails(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasMassSpecDetails = function() {
  return jspb.Message.getField(this, 13) != null;
};


/**
 * optional Selectivity selectivity = 14;
 * @return {?proto.ord.ProductMeasurement.Selectivity}
 */
proto.ord.ProductMeasurement.prototype.getSelectivity = function() {
  return /** @type{?proto.ord.ProductMeasurement.Selectivity} */ (
    jspb.Message.getWrapperField(this, proto.ord.ProductMeasurement.Selectivity, 14));
};


/**
 * @param {?proto.ord.ProductMeasurement.Selectivity|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setSelectivity = function(value) {
  return jspb.Message.setWrapperField(this, 14, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearSelectivity = function() {
  return this.setSelectivity(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasSelectivity = function() {
  return jspb.Message.getField(this, 14) != null;
};


/**
 * optional Wavelength wavelength = 15;
 * @return {?proto.ord.Wavelength}
 */
proto.ord.ProductMeasurement.prototype.getWavelength = function() {
  return /** @type{?proto.ord.Wavelength} */ (
    jspb.Message.getWrapperField(this, proto.ord.Wavelength, 15));
};


/**
 * @param {?proto.ord.Wavelength|undefined} value
 * @return {!proto.ord.ProductMeasurement} returns this
*/
proto.ord.ProductMeasurement.prototype.setWavelength = function(value) {
  return jspb.Message.setWrapperField(this, 15, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ProductMeasurement} returns this
 */
proto.ord.ProductMeasurement.prototype.clearWavelength = function() {
  return this.setWavelength(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ProductMeasurement.prototype.hasWavelength = function() {
  return jspb.Message.getField(this, 15) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.DateTime.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.DateTime.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.DateTime} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.DateTime.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFieldWithDefault(msg, 1, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.DateTime}
 */
proto.ord.DateTime.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.DateTime;
  return proto.ord.DateTime.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.DateTime} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.DateTime}
 */
proto.ord.DateTime.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {string} */ (reader.readString());
      msg.setValue(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.DateTime.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.DateTime.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.DateTime} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.DateTime.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getValue();
  if (f.length > 0) {
    writer.writeString(
      1,
      f
    );
  }
};


/**
 * optional string value = 1;
 * @return {string}
 */
proto.ord.DateTime.prototype.getValue = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 1, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.DateTime} returns this
 */
proto.ord.DateTime.prototype.setValue = function(value) {
  return jspb.Message.setProto3StringField(this, 1, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Analysis.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Analysis.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Analysis} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Analysis.toObject = function(includeInstance, msg) {
  var f, obj = {
    type: jspb.Message.getFieldWithDefault(msg, 1, 0),
    details: jspb.Message.getFieldWithDefault(msg, 2, ""),
    chmoId: jspb.Message.getFieldWithDefault(msg, 3, 0),
    isOfIsolatedSpecies: jspb.Message.getBooleanFieldWithDefault(msg, 4, false),
    dataMap: (f = msg.getDataMap()) ? f.toObject(includeInstance, proto.ord.Data.toObject) : [],
    instrumentManufacturer: jspb.Message.getFieldWithDefault(msg, 6, ""),
    instrumentLastCalibrated: (f = msg.getInstrumentLastCalibrated()) && proto.ord.DateTime.toObject(includeInstance, f)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Analysis}
 */
proto.ord.Analysis.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Analysis;
  return proto.ord.Analysis.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Analysis} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Analysis}
 */
proto.ord.Analysis.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {!proto.ord.Analysis.AnalysisType} */ (reader.readEnum());
      msg.setType(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    case 3:
      var value = /** @type {number} */ (reader.readInt32());
      msg.setChmoId(value);
      break;
    case 4:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsOfIsolatedSpecies(value);
      break;
    case 5:
      var value = msg.getDataMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Data.deserializeBinaryFromReader, "", new proto.ord.Data());
         });
      break;
    case 6:
      var value = /** @type {string} */ (reader.readString());
      msg.setInstrumentManufacturer(value);
      break;
    case 7:
      var value = new proto.ord.DateTime;
      reader.readMessage(value,proto.ord.DateTime.deserializeBinaryFromReader);
      msg.setInstrumentLastCalibrated(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Analysis.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Analysis.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Analysis} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Analysis.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getType();
  if (f !== 0.0) {
    writer.writeEnum(
      1,
      f
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getChmoId();
  if (f !== 0) {
    writer.writeInt32(
      3,
      f
    );
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeBool(
      4,
      f
    );
  }
  f = message.getDataMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(5, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Data.serializeBinaryToWriter);
  }
  f = message.getInstrumentManufacturer();
  if (f.length > 0) {
    writer.writeString(
      6,
      f
    );
  }
  f = message.getInstrumentLastCalibrated();
  if (f != null) {
    writer.writeMessage(
      7,
      f,
      proto.ord.DateTime.serializeBinaryToWriter
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Analysis.AnalysisType = {
  UNSPECIFIED: 0,
  CUSTOM: 1,
  LC: 2,
  GC: 3,
  IR: 4,
  NMR_1H: 5,
  NMR_13C: 6,
  NMR_OTHER: 7,
  MP: 8,
  UV: 9,
  TLC: 10,
  MS: 11,
  HRMS: 12,
  MSMS: 13,
  WEIGHT: 14,
  LCMS: 15,
  GCMS: 16,
  ELSD: 17,
  CD: 18,
  SFC: 19,
  EPR: 20,
  XRD: 21,
  RAMAN: 22,
  ED: 23
};

/**
 * optional AnalysisType type = 1;
 * @return {!proto.ord.Analysis.AnalysisType}
 */
proto.ord.Analysis.prototype.getType = function() {
  return /** @type {!proto.ord.Analysis.AnalysisType} */ (jspb.Message.getFieldWithDefault(this, 1, 0));
};


/**
 * @param {!proto.ord.Analysis.AnalysisType} value
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.setType = function(value) {
  return jspb.Message.setProto3EnumField(this, 1, value);
};


/**
 * optional string details = 2;
 * @return {string}
 */
proto.ord.Analysis.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional int32 chmo_id = 3;
 * @return {number}
 */
proto.ord.Analysis.prototype.getChmoId = function() {
  return /** @type {number} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.setChmoId = function(value) {
  return jspb.Message.setProto3IntField(this, 3, value);
};


/**
 * optional bool is_of_isolated_species = 4;
 * @return {boolean}
 */
proto.ord.Analysis.prototype.getIsOfIsolatedSpecies = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 4, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.setIsOfIsolatedSpecies = function(value) {
  return jspb.Message.setField(this, 4, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.clearIsOfIsolatedSpecies = function() {
  return jspb.Message.setField(this, 4, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Analysis.prototype.hasIsOfIsolatedSpecies = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * map<string, Data> data = 5;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Data>}
 */
proto.ord.Analysis.prototype.getDataMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Data>} */ (
      jspb.Message.getMapField(this, 5, opt_noLazyCreate,
      proto.ord.Data));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.clearDataMap = function() {
  this.getDataMap().clear();
  return this;
};


/**
 * optional string instrument_manufacturer = 6;
 * @return {string}
 */
proto.ord.Analysis.prototype.getInstrumentManufacturer = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 6, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.setInstrumentManufacturer = function(value) {
  return jspb.Message.setProto3StringField(this, 6, value);
};


/**
 * optional DateTime instrument_last_calibrated = 7;
 * @return {?proto.ord.DateTime}
 */
proto.ord.Analysis.prototype.getInstrumentLastCalibrated = function() {
  return /** @type{?proto.ord.DateTime} */ (
    jspb.Message.getWrapperField(this, proto.ord.DateTime, 7));
};


/**
 * @param {?proto.ord.DateTime|undefined} value
 * @return {!proto.ord.Analysis} returns this
*/
proto.ord.Analysis.prototype.setInstrumentLastCalibrated = function(value) {
  return jspb.Message.setWrapperField(this, 7, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.Analysis} returns this
 */
proto.ord.Analysis.prototype.clearInstrumentLastCalibrated = function() {
  return this.setInstrumentLastCalibrated(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Analysis.prototype.hasInstrumentLastCalibrated = function() {
  return jspb.Message.getField(this, 7) != null;
};



/**
 * List of repeated fields within this message type.
 * @private {!Array<number>}
 * @const
 */
proto.ord.ReactionProvenance.repeatedFields_ = [8];



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.ReactionProvenance.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.ReactionProvenance.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.ReactionProvenance} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionProvenance.toObject = function(includeInstance, msg) {
  var f, obj = {
    experimenter: (f = msg.getExperimenter()) && proto.ord.Person.toObject(includeInstance, f),
    city: jspb.Message.getFieldWithDefault(msg, 2, ""),
    experimentStart: (f = msg.getExperimentStart()) && proto.ord.DateTime.toObject(includeInstance, f),
    doi: jspb.Message.getFieldWithDefault(msg, 4, ""),
    patent: jspb.Message.getFieldWithDefault(msg, 5, ""),
    publicationUrl: jspb.Message.getFieldWithDefault(msg, 6, ""),
    recordCreated: (f = msg.getRecordCreated()) && proto.ord.RecordEvent.toObject(includeInstance, f),
    recordModifiedList: jspb.Message.toObjectList(msg.getRecordModifiedList(),
    proto.ord.RecordEvent.toObject, includeInstance),
    reactionMetadataMap: (f = msg.getReactionMetadataMap()) ? f.toObject(includeInstance, proto.ord.Data.toObject) : [],
    isMined: jspb.Message.getBooleanFieldWithDefault(msg, 10, false)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.ReactionProvenance}
 */
proto.ord.ReactionProvenance.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.ReactionProvenance;
  return proto.ord.ReactionProvenance.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.ReactionProvenance} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.ReactionProvenance}
 */
proto.ord.ReactionProvenance.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.Person;
      reader.readMessage(value,proto.ord.Person.deserializeBinaryFromReader);
      msg.setExperimenter(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setCity(value);
      break;
    case 3:
      var value = new proto.ord.DateTime;
      reader.readMessage(value,proto.ord.DateTime.deserializeBinaryFromReader);
      msg.setExperimentStart(value);
      break;
    case 4:
      var value = /** @type {string} */ (reader.readString());
      msg.setDoi(value);
      break;
    case 5:
      var value = /** @type {string} */ (reader.readString());
      msg.setPatent(value);
      break;
    case 6:
      var value = /** @type {string} */ (reader.readString());
      msg.setPublicationUrl(value);
      break;
    case 7:
      var value = new proto.ord.RecordEvent;
      reader.readMessage(value,proto.ord.RecordEvent.deserializeBinaryFromReader);
      msg.setRecordCreated(value);
      break;
    case 8:
      var value = new proto.ord.RecordEvent;
      reader.readMessage(value,proto.ord.RecordEvent.deserializeBinaryFromReader);
      msg.addRecordModified(value);
      break;
    case 9:
      var value = msg.getReactionMetadataMap();
      reader.readMessage(value, function(message, reader) {
        jspb.Map.deserializeBinary(message, reader, jspb.BinaryReader.prototype.readString, jspb.BinaryReader.prototype.readMessage, proto.ord.Data.deserializeBinaryFromReader, "", new proto.ord.Data());
         });
      break;
    case 10:
      var value = /** @type {boolean} */ (reader.readBool());
      msg.setIsMined(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.ReactionProvenance.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.ReactionProvenance.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.ReactionProvenance} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.ReactionProvenance.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getExperimenter();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.Person.serializeBinaryToWriter
    );
  }
  f = message.getCity();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getExperimentStart();
  if (f != null) {
    writer.writeMessage(
      3,
      f,
      proto.ord.DateTime.serializeBinaryToWriter
    );
  }
  f = message.getDoi();
  if (f.length > 0) {
    writer.writeString(
      4,
      f
    );
  }
  f = message.getPatent();
  if (f.length > 0) {
    writer.writeString(
      5,
      f
    );
  }
  f = message.getPublicationUrl();
  if (f.length > 0) {
    writer.writeString(
      6,
      f
    );
  }
  f = message.getRecordCreated();
  if (f != null) {
    writer.writeMessage(
      7,
      f,
      proto.ord.RecordEvent.serializeBinaryToWriter
    );
  }
  f = message.getRecordModifiedList();
  if (f.length > 0) {
    writer.writeRepeatedMessage(
      8,
      f,
      proto.ord.RecordEvent.serializeBinaryToWriter
    );
  }
  f = message.getReactionMetadataMap(true);
  if (f && f.getLength() > 0) {
    f.serializeBinary(9, writer, jspb.BinaryWriter.prototype.writeString, jspb.BinaryWriter.prototype.writeMessage, proto.ord.Data.serializeBinaryToWriter);
  }
  f = /** @type {boolean} */ (jspb.Message.getField(message, 10));
  if (f != null) {
    writer.writeBool(
      10,
      f
    );
  }
};


/**
 * optional Person experimenter = 1;
 * @return {?proto.ord.Person}
 */
proto.ord.ReactionProvenance.prototype.getExperimenter = function() {
  return /** @type{?proto.ord.Person} */ (
    jspb.Message.getWrapperField(this, proto.ord.Person, 1));
};


/**
 * @param {?proto.ord.Person|undefined} value
 * @return {!proto.ord.ReactionProvenance} returns this
*/
proto.ord.ReactionProvenance.prototype.setExperimenter = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.clearExperimenter = function() {
  return this.setExperimenter(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionProvenance.prototype.hasExperimenter = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional string city = 2;
 * @return {string}
 */
proto.ord.ReactionProvenance.prototype.getCity = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.setCity = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional DateTime experiment_start = 3;
 * @return {?proto.ord.DateTime}
 */
proto.ord.ReactionProvenance.prototype.getExperimentStart = function() {
  return /** @type{?proto.ord.DateTime} */ (
    jspb.Message.getWrapperField(this, proto.ord.DateTime, 3));
};


/**
 * @param {?proto.ord.DateTime|undefined} value
 * @return {!proto.ord.ReactionProvenance} returns this
*/
proto.ord.ReactionProvenance.prototype.setExperimentStart = function(value) {
  return jspb.Message.setWrapperField(this, 3, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.clearExperimentStart = function() {
  return this.setExperimentStart(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionProvenance.prototype.hasExperimentStart = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional string doi = 4;
 * @return {string}
 */
proto.ord.ReactionProvenance.prototype.getDoi = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 4, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.setDoi = function(value) {
  return jspb.Message.setProto3StringField(this, 4, value);
};


/**
 * optional string patent = 5;
 * @return {string}
 */
proto.ord.ReactionProvenance.prototype.getPatent = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 5, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.setPatent = function(value) {
  return jspb.Message.setProto3StringField(this, 5, value);
};


/**
 * optional string publication_url = 6;
 * @return {string}
 */
proto.ord.ReactionProvenance.prototype.getPublicationUrl = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 6, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.setPublicationUrl = function(value) {
  return jspb.Message.setProto3StringField(this, 6, value);
};


/**
 * optional RecordEvent record_created = 7;
 * @return {?proto.ord.RecordEvent}
 */
proto.ord.ReactionProvenance.prototype.getRecordCreated = function() {
  return /** @type{?proto.ord.RecordEvent} */ (
    jspb.Message.getWrapperField(this, proto.ord.RecordEvent, 7));
};


/**
 * @param {?proto.ord.RecordEvent|undefined} value
 * @return {!proto.ord.ReactionProvenance} returns this
*/
proto.ord.ReactionProvenance.prototype.setRecordCreated = function(value) {
  return jspb.Message.setWrapperField(this, 7, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.clearRecordCreated = function() {
  return this.setRecordCreated(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionProvenance.prototype.hasRecordCreated = function() {
  return jspb.Message.getField(this, 7) != null;
};


/**
 * repeated RecordEvent record_modified = 8;
 * @return {!Array<!proto.ord.RecordEvent>}
 */
proto.ord.ReactionProvenance.prototype.getRecordModifiedList = function() {
  return /** @type{!Array<!proto.ord.RecordEvent>} */ (
    jspb.Message.getRepeatedWrapperField(this, proto.ord.RecordEvent, 8));
};


/**
 * @param {!Array<!proto.ord.RecordEvent>} value
 * @return {!proto.ord.ReactionProvenance} returns this
*/
proto.ord.ReactionProvenance.prototype.setRecordModifiedList = function(value) {
  return jspb.Message.setRepeatedWrapperField(this, 8, value);
};


/**
 * @param {!proto.ord.RecordEvent=} opt_value
 * @param {number=} opt_index
 * @return {!proto.ord.RecordEvent}
 */
proto.ord.ReactionProvenance.prototype.addRecordModified = function(opt_value, opt_index) {
  return jspb.Message.addToRepeatedWrapperField(this, 8, opt_value, proto.ord.RecordEvent, opt_index);
};


/**
 * Clears the list making it empty but non-null.
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.clearRecordModifiedList = function() {
  return this.setRecordModifiedList([]);
};


/**
 * map<string, Data> reaction_metadata = 9;
 * @param {boolean=} opt_noLazyCreate Do not create the map if
 * empty, instead returning `undefined`
 * @return {!jspb.Map<string,!proto.ord.Data>}
 */
proto.ord.ReactionProvenance.prototype.getReactionMetadataMap = function(opt_noLazyCreate) {
  return /** @type {!jspb.Map<string,!proto.ord.Data>} */ (
      jspb.Message.getMapField(this, 9, opt_noLazyCreate,
      proto.ord.Data));
};


/**
 * Clears values from the map. The map will be non-null.
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.clearReactionMetadataMap = function() {
  this.getReactionMetadataMap().clear();
  return this;
};


/**
 * optional bool is_mined = 10;
 * @return {boolean}
 */
proto.ord.ReactionProvenance.prototype.getIsMined = function() {
  return /** @type {boolean} */ (jspb.Message.getBooleanFieldWithDefault(this, 10, false));
};


/**
 * @param {boolean} value
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.setIsMined = function(value) {
  return jspb.Message.setField(this, 10, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.ReactionProvenance} returns this
 */
proto.ord.ReactionProvenance.prototype.clearIsMined = function() {
  return jspb.Message.setField(this, 10, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.ReactionProvenance.prototype.hasIsMined = function() {
  return jspb.Message.getField(this, 10) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Person.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Person.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Person} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Person.toObject = function(includeInstance, msg) {
  var f, obj = {
    username: jspb.Message.getFieldWithDefault(msg, 1, ""),
    name: jspb.Message.getFieldWithDefault(msg, 2, ""),
    orcid: jspb.Message.getFieldWithDefault(msg, 3, ""),
    organization: jspb.Message.getFieldWithDefault(msg, 4, ""),
    email: jspb.Message.getFieldWithDefault(msg, 5, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Person}
 */
proto.ord.Person.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Person;
  return proto.ord.Person.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Person} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Person}
 */
proto.ord.Person.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {string} */ (reader.readString());
      msg.setUsername(value);
      break;
    case 2:
      var value = /** @type {string} */ (reader.readString());
      msg.setName(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setOrcid(value);
      break;
    case 4:
      var value = /** @type {string} */ (reader.readString());
      msg.setOrganization(value);
      break;
    case 5:
      var value = /** @type {string} */ (reader.readString());
      msg.setEmail(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Person.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Person.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Person} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Person.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getUsername();
  if (f.length > 0) {
    writer.writeString(
      1,
      f
    );
  }
  f = message.getName();
  if (f.length > 0) {
    writer.writeString(
      2,
      f
    );
  }
  f = message.getOrcid();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
  f = message.getOrganization();
  if (f.length > 0) {
    writer.writeString(
      4,
      f
    );
  }
  f = message.getEmail();
  if (f.length > 0) {
    writer.writeString(
      5,
      f
    );
  }
};


/**
 * optional string username = 1;
 * @return {string}
 */
proto.ord.Person.prototype.getUsername = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 1, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Person} returns this
 */
proto.ord.Person.prototype.setUsername = function(value) {
  return jspb.Message.setProto3StringField(this, 1, value);
};


/**
 * optional string name = 2;
 * @return {string}
 */
proto.ord.Person.prototype.getName = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 2, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Person} returns this
 */
proto.ord.Person.prototype.setName = function(value) {
  return jspb.Message.setProto3StringField(this, 2, value);
};


/**
 * optional string orcid = 3;
 * @return {string}
 */
proto.ord.Person.prototype.getOrcid = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Person} returns this
 */
proto.ord.Person.prototype.setOrcid = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};


/**
 * optional string organization = 4;
 * @return {string}
 */
proto.ord.Person.prototype.getOrganization = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 4, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Person} returns this
 */
proto.ord.Person.prototype.setOrganization = function(value) {
  return jspb.Message.setProto3StringField(this, 4, value);
};


/**
 * optional string email = 5;
 * @return {string}
 */
proto.ord.Person.prototype.getEmail = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 5, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Person} returns this
 */
proto.ord.Person.prototype.setEmail = function(value) {
  return jspb.Message.setProto3StringField(this, 5, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.RecordEvent.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.RecordEvent.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.RecordEvent} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.RecordEvent.toObject = function(includeInstance, msg) {
  var f, obj = {
    time: (f = msg.getTime()) && proto.ord.DateTime.toObject(includeInstance, f),
    person: (f = msg.getPerson()) && proto.ord.Person.toObject(includeInstance, f),
    details: jspb.Message.getFieldWithDefault(msg, 3, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.RecordEvent}
 */
proto.ord.RecordEvent.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.RecordEvent;
  return proto.ord.RecordEvent.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.RecordEvent} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.RecordEvent}
 */
proto.ord.RecordEvent.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = new proto.ord.DateTime;
      reader.readMessage(value,proto.ord.DateTime.deserializeBinaryFromReader);
      msg.setTime(value);
      break;
    case 2:
      var value = new proto.ord.Person;
      reader.readMessage(value,proto.ord.Person.deserializeBinaryFromReader);
      msg.setPerson(value);
      break;
    case 3:
      var value = /** @type {string} */ (reader.readString());
      msg.setDetails(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.RecordEvent.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.RecordEvent.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.RecordEvent} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.RecordEvent.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = message.getTime();
  if (f != null) {
    writer.writeMessage(
      1,
      f,
      proto.ord.DateTime.serializeBinaryToWriter
    );
  }
  f = message.getPerson();
  if (f != null) {
    writer.writeMessage(
      2,
      f,
      proto.ord.Person.serializeBinaryToWriter
    );
  }
  f = message.getDetails();
  if (f.length > 0) {
    writer.writeString(
      3,
      f
    );
  }
};


/**
 * optional DateTime time = 1;
 * @return {?proto.ord.DateTime}
 */
proto.ord.RecordEvent.prototype.getTime = function() {
  return /** @type{?proto.ord.DateTime} */ (
    jspb.Message.getWrapperField(this, proto.ord.DateTime, 1));
};


/**
 * @param {?proto.ord.DateTime|undefined} value
 * @return {!proto.ord.RecordEvent} returns this
*/
proto.ord.RecordEvent.prototype.setTime = function(value) {
  return jspb.Message.setWrapperField(this, 1, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.RecordEvent} returns this
 */
proto.ord.RecordEvent.prototype.clearTime = function() {
  return this.setTime(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.RecordEvent.prototype.hasTime = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional Person person = 2;
 * @return {?proto.ord.Person}
 */
proto.ord.RecordEvent.prototype.getPerson = function() {
  return /** @type{?proto.ord.Person} */ (
    jspb.Message.getWrapperField(this, proto.ord.Person, 2));
};


/**
 * @param {?proto.ord.Person|undefined} value
 * @return {!proto.ord.RecordEvent} returns this
*/
proto.ord.RecordEvent.prototype.setPerson = function(value) {
  return jspb.Message.setWrapperField(this, 2, value);
};


/**
 * Clears the message field making it undefined.
 * @return {!proto.ord.RecordEvent} returns this
 */
proto.ord.RecordEvent.prototype.clearPerson = function() {
  return this.setPerson(undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.RecordEvent.prototype.hasPerson = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional string details = 3;
 * @return {string}
 */
proto.ord.RecordEvent.prototype.getDetails = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.RecordEvent} returns this
 */
proto.ord.RecordEvent.prototype.setDetails = function(value) {
  return jspb.Message.setProto3StringField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Time.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Time.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Time} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Time.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Time}
 */
proto.ord.Time.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Time;
  return proto.ord.Time.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Time} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Time}
 */
proto.ord.Time.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Time.TimeUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Time.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Time.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Time} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Time.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Time.TimeUnit = {
  UNSPECIFIED: 0,
  DAY: 4,
  HOUR: 1,
  MINUTE: 2,
  SECOND: 3
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Time.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Time} returns this
 */
proto.ord.Time.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Time} returns this
 */
proto.ord.Time.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Time.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Time.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Time} returns this
 */
proto.ord.Time.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Time} returns this
 */
proto.ord.Time.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Time.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional TimeUnit units = 3;
 * @return {!proto.ord.Time.TimeUnit}
 */
proto.ord.Time.prototype.getUnits = function() {
  return /** @type {!proto.ord.Time.TimeUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Time.TimeUnit} value
 * @return {!proto.ord.Time} returns this
 */
proto.ord.Time.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Mass.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Mass.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Mass} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Mass.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Mass}
 */
proto.ord.Mass.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Mass;
  return proto.ord.Mass.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Mass} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Mass}
 */
proto.ord.Mass.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Mass.MassUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Mass.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Mass.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Mass} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Mass.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Mass.MassUnit = {
  UNSPECIFIED: 0,
  KILOGRAM: 1,
  GRAM: 2,
  MILLIGRAM: 3,
  MICROGRAM: 4
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Mass.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Mass} returns this
 */
proto.ord.Mass.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Mass} returns this
 */
proto.ord.Mass.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Mass.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Mass.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Mass} returns this
 */
proto.ord.Mass.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Mass} returns this
 */
proto.ord.Mass.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Mass.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional MassUnit units = 3;
 * @return {!proto.ord.Mass.MassUnit}
 */
proto.ord.Mass.prototype.getUnits = function() {
  return /** @type {!proto.ord.Mass.MassUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Mass.MassUnit} value
 * @return {!proto.ord.Mass} returns this
 */
proto.ord.Mass.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Moles.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Moles.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Moles} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Moles.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Moles}
 */
proto.ord.Moles.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Moles;
  return proto.ord.Moles.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Moles} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Moles}
 */
proto.ord.Moles.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Moles.MolesUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Moles.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Moles.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Moles} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Moles.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Moles.MolesUnit = {
  UNSPECIFIED: 0,
  MOLE: 1,
  MILLIMOLE: 2,
  MICROMOLE: 3,
  NANOMOLE: 4
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Moles.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Moles} returns this
 */
proto.ord.Moles.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Moles} returns this
 */
proto.ord.Moles.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Moles.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Moles.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Moles} returns this
 */
proto.ord.Moles.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Moles} returns this
 */
proto.ord.Moles.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Moles.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional MolesUnit units = 3;
 * @return {!proto.ord.Moles.MolesUnit}
 */
proto.ord.Moles.prototype.getUnits = function() {
  return /** @type {!proto.ord.Moles.MolesUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Moles.MolesUnit} value
 * @return {!proto.ord.Moles} returns this
 */
proto.ord.Moles.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Volume.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Volume.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Volume} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Volume.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Volume}
 */
proto.ord.Volume.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Volume;
  return proto.ord.Volume.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Volume} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Volume}
 */
proto.ord.Volume.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Volume.VolumeUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Volume.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Volume.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Volume} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Volume.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Volume.VolumeUnit = {
  UNSPECIFIED: 0,
  LITER: 1,
  MILLILITER: 2,
  MICROLITER: 3,
  NANOLITER: 4
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Volume.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Volume} returns this
 */
proto.ord.Volume.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Volume} returns this
 */
proto.ord.Volume.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Volume.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Volume.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Volume} returns this
 */
proto.ord.Volume.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Volume} returns this
 */
proto.ord.Volume.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Volume.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional VolumeUnit units = 3;
 * @return {!proto.ord.Volume.VolumeUnit}
 */
proto.ord.Volume.prototype.getUnits = function() {
  return /** @type {!proto.ord.Volume.VolumeUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Volume.VolumeUnit} value
 * @return {!proto.ord.Volume} returns this
 */
proto.ord.Volume.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Concentration.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Concentration.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Concentration} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Concentration.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Concentration}
 */
proto.ord.Concentration.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Concentration;
  return proto.ord.Concentration.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Concentration} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Concentration}
 */
proto.ord.Concentration.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Concentration.ConcentrationUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Concentration.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Concentration.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Concentration} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Concentration.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Concentration.ConcentrationUnit = {
  UNSPECIFIED: 0,
  MOLAR: 1,
  MILLIMOLAR: 2,
  MICROMOLAR: 3
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Concentration.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Concentration} returns this
 */
proto.ord.Concentration.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Concentration} returns this
 */
proto.ord.Concentration.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Concentration.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Concentration.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Concentration} returns this
 */
proto.ord.Concentration.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Concentration} returns this
 */
proto.ord.Concentration.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Concentration.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional ConcentrationUnit units = 3;
 * @return {!proto.ord.Concentration.ConcentrationUnit}
 */
proto.ord.Concentration.prototype.getUnits = function() {
  return /** @type {!proto.ord.Concentration.ConcentrationUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Concentration.ConcentrationUnit} value
 * @return {!proto.ord.Concentration} returns this
 */
proto.ord.Concentration.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Pressure.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Pressure.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Pressure} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Pressure.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Pressure}
 */
proto.ord.Pressure.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Pressure;
  return proto.ord.Pressure.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Pressure} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Pressure}
 */
proto.ord.Pressure.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Pressure.PressureUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Pressure.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Pressure.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Pressure} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Pressure.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Pressure.PressureUnit = {
  UNSPECIFIED: 0,
  BAR: 1,
  ATMOSPHERE: 2,
  PSI: 3,
  KPSI: 4,
  PASCAL: 5,
  KILOPASCAL: 6,
  TORR: 7,
  MM_HG: 8
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Pressure.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Pressure} returns this
 */
proto.ord.Pressure.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Pressure} returns this
 */
proto.ord.Pressure.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Pressure.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Pressure.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Pressure} returns this
 */
proto.ord.Pressure.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Pressure} returns this
 */
proto.ord.Pressure.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Pressure.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional PressureUnit units = 3;
 * @return {!proto.ord.Pressure.PressureUnit}
 */
proto.ord.Pressure.prototype.getUnits = function() {
  return /** @type {!proto.ord.Pressure.PressureUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Pressure.PressureUnit} value
 * @return {!proto.ord.Pressure} returns this
 */
proto.ord.Pressure.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Temperature.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Temperature.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Temperature} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Temperature.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Temperature}
 */
proto.ord.Temperature.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Temperature;
  return proto.ord.Temperature.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Temperature} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Temperature}
 */
proto.ord.Temperature.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Temperature.TemperatureUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Temperature.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Temperature.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Temperature} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Temperature.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Temperature.TemperatureUnit = {
  UNSPECIFIED: 0,
  CELSIUS: 1,
  FAHRENHEIT: 2,
  KELVIN: 3
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Temperature.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Temperature} returns this
 */
proto.ord.Temperature.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Temperature} returns this
 */
proto.ord.Temperature.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Temperature.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Temperature.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Temperature} returns this
 */
proto.ord.Temperature.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Temperature} returns this
 */
proto.ord.Temperature.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Temperature.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional TemperatureUnit units = 3;
 * @return {!proto.ord.Temperature.TemperatureUnit}
 */
proto.ord.Temperature.prototype.getUnits = function() {
  return /** @type {!proto.ord.Temperature.TemperatureUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Temperature.TemperatureUnit} value
 * @return {!proto.ord.Temperature} returns this
 */
proto.ord.Temperature.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Current.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Current.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Current} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Current.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Current}
 */
proto.ord.Current.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Current;
  return proto.ord.Current.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Current} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Current}
 */
proto.ord.Current.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Current.CurrentUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Current.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Current.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Current} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Current.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Current.CurrentUnit = {
  UNSPECIFIED: 0,
  AMPERE: 1,
  MILLIAMPERE: 2
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Current.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Current} returns this
 */
proto.ord.Current.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Current} returns this
 */
proto.ord.Current.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Current.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Current.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Current} returns this
 */
proto.ord.Current.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Current} returns this
 */
proto.ord.Current.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Current.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional CurrentUnit units = 3;
 * @return {!proto.ord.Current.CurrentUnit}
 */
proto.ord.Current.prototype.getUnits = function() {
  return /** @type {!proto.ord.Current.CurrentUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Current.CurrentUnit} value
 * @return {!proto.ord.Current} returns this
 */
proto.ord.Current.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Voltage.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Voltage.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Voltage} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Voltage.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Voltage}
 */
proto.ord.Voltage.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Voltage;
  return proto.ord.Voltage.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Voltage} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Voltage}
 */
proto.ord.Voltage.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Voltage.VoltageUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Voltage.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Voltage.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Voltage} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Voltage.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Voltage.VoltageUnit = {
  UNSPECIFIED: 0,
  VOLT: 1,
  MILLIVOLT: 2
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Voltage.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Voltage} returns this
 */
proto.ord.Voltage.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Voltage} returns this
 */
proto.ord.Voltage.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Voltage.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Voltage.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Voltage} returns this
 */
proto.ord.Voltage.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Voltage} returns this
 */
proto.ord.Voltage.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Voltage.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional VoltageUnit units = 3;
 * @return {!proto.ord.Voltage.VoltageUnit}
 */
proto.ord.Voltage.prototype.getUnits = function() {
  return /** @type {!proto.ord.Voltage.VoltageUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Voltage.VoltageUnit} value
 * @return {!proto.ord.Voltage} returns this
 */
proto.ord.Voltage.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Length.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Length.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Length} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Length.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Length}
 */
proto.ord.Length.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Length;
  return proto.ord.Length.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Length} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Length}
 */
proto.ord.Length.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Length.LengthUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Length.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Length.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Length} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Length.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Length.LengthUnit = {
  UNSPECIFIED: 0,
  CENTIMETER: 1,
  MILLIMETER: 2,
  METER: 3,
  INCH: 4,
  FOOT: 5
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Length.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Length} returns this
 */
proto.ord.Length.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Length} returns this
 */
proto.ord.Length.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Length.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Length.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Length} returns this
 */
proto.ord.Length.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Length} returns this
 */
proto.ord.Length.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Length.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional LengthUnit units = 3;
 * @return {!proto.ord.Length.LengthUnit}
 */
proto.ord.Length.prototype.getUnits = function() {
  return /** @type {!proto.ord.Length.LengthUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Length.LengthUnit} value
 * @return {!proto.ord.Length} returns this
 */
proto.ord.Length.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Wavelength.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Wavelength.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Wavelength} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Wavelength.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Wavelength}
 */
proto.ord.Wavelength.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Wavelength;
  return proto.ord.Wavelength.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Wavelength} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Wavelength}
 */
proto.ord.Wavelength.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.Wavelength.WavelengthUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Wavelength.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Wavelength.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Wavelength} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Wavelength.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.Wavelength.WavelengthUnit = {
  UNSPECIFIED: 0,
  NANOMETER: 1,
  WAVENUMBER: 2
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Wavelength.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Wavelength} returns this
 */
proto.ord.Wavelength.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Wavelength} returns this
 */
proto.ord.Wavelength.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Wavelength.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Wavelength.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Wavelength} returns this
 */
proto.ord.Wavelength.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Wavelength} returns this
 */
proto.ord.Wavelength.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Wavelength.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional WavelengthUnit units = 3;
 * @return {!proto.ord.Wavelength.WavelengthUnit}
 */
proto.ord.Wavelength.prototype.getUnits = function() {
  return /** @type {!proto.ord.Wavelength.WavelengthUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.Wavelength.WavelengthUnit} value
 * @return {!proto.ord.Wavelength} returns this
 */
proto.ord.Wavelength.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.FlowRate.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.FlowRate.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.FlowRate} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FlowRate.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0),
    units: jspb.Message.getFieldWithDefault(msg, 3, 0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.FlowRate}
 */
proto.ord.FlowRate.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.FlowRate;
  return proto.ord.FlowRate.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.FlowRate} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.FlowRate}
 */
proto.ord.FlowRate.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    case 3:
      var value = /** @type {!proto.ord.FlowRate.FlowRateUnit} */ (reader.readEnum());
      msg.setUnits(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.FlowRate.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.FlowRate.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.FlowRate} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FlowRate.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
  f = message.getUnits();
  if (f !== 0.0) {
    writer.writeEnum(
      3,
      f
    );
  }
};


/**
 * @enum {number}
 */
proto.ord.FlowRate.FlowRateUnit = {
  UNSPECIFIED: 0,
  MICROLITER_PER_MINUTE: 1,
  MICROLITER_PER_SECOND: 2,
  MILLILITER_PER_MINUTE: 3,
  MILLILITER_PER_SECOND: 4,
  MICROLITER_PER_HOUR: 5
};

/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.FlowRate.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.FlowRate} returns this
 */
proto.ord.FlowRate.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.FlowRate} returns this
 */
proto.ord.FlowRate.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.FlowRate.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.FlowRate.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.FlowRate} returns this
 */
proto.ord.FlowRate.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.FlowRate} returns this
 */
proto.ord.FlowRate.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.FlowRate.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional FlowRateUnit units = 3;
 * @return {!proto.ord.FlowRate.FlowRateUnit}
 */
proto.ord.FlowRate.prototype.getUnits = function() {
  return /** @type {!proto.ord.FlowRate.FlowRateUnit} */ (jspb.Message.getFieldWithDefault(this, 3, 0));
};


/**
 * @param {!proto.ord.FlowRate.FlowRateUnit} value
 * @return {!proto.ord.FlowRate} returns this
 */
proto.ord.FlowRate.prototype.setUnits = function(value) {
  return jspb.Message.setProto3EnumField(this, 3, value);
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Percentage.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Percentage.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Percentage} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Percentage.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Percentage}
 */
proto.ord.Percentage.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Percentage;
  return proto.ord.Percentage.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Percentage} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Percentage}
 */
proto.ord.Percentage.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Percentage.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Percentage.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Percentage} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Percentage.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
};


/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.Percentage.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Percentage} returns this
 */
proto.ord.Percentage.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Percentage} returns this
 */
proto.ord.Percentage.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Percentage.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.Percentage.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Percentage} returns this
 */
proto.ord.Percentage.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Percentage} returns this
 */
proto.ord.Percentage.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Percentage.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};





if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.FloatValue.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.FloatValue.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.FloatValue} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FloatValue.toObject = function(includeInstance, msg) {
  var f, obj = {
    value: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    precision: jspb.Message.getFloatingPointFieldWithDefault(msg, 2, 0.0)
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.FloatValue}
 */
proto.ord.FloatValue.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.FloatValue;
  return proto.ord.FloatValue.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.FloatValue} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.FloatValue}
 */
proto.ord.FloatValue.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setPrecision(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.FloatValue.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.FloatValue.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.FloatValue} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.FloatValue.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeFloat(
      2,
      f
    );
  }
};


/**
 * optional float value = 1;
 * @return {number}
 */
proto.ord.FloatValue.prototype.getValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.FloatValue} returns this
 */
proto.ord.FloatValue.prototype.setValue = function(value) {
  return jspb.Message.setField(this, 1, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.FloatValue} returns this
 */
proto.ord.FloatValue.prototype.clearValue = function() {
  return jspb.Message.setField(this, 1, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.FloatValue.prototype.hasValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional float precision = 2;
 * @return {number}
 */
proto.ord.FloatValue.prototype.getPrecision = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 2, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.FloatValue} returns this
 */
proto.ord.FloatValue.prototype.setPrecision = function(value) {
  return jspb.Message.setField(this, 2, value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.FloatValue} returns this
 */
proto.ord.FloatValue.prototype.clearPrecision = function() {
  return jspb.Message.setField(this, 2, undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.FloatValue.prototype.hasPrecision = function() {
  return jspb.Message.getField(this, 2) != null;
};



/**
 * Oneof group definitions for this message. Each group defines the field
 * numbers belonging to that group. When of these fields' value is set, all
 * other fields in the group are cleared. During deserialization, if multiple
 * fields are encountered for a group, only the last value seen will be kept.
 * @private {!Array<!Array<number>>}
 * @const
 */
proto.ord.Data.oneofGroups_ = [[1,2,3,4,5]];

/**
 * @enum {number}
 */
proto.ord.Data.KindCase = {
  KIND_NOT_SET: 0,
  FLOAT_VALUE: 1,
  INTEGER_VALUE: 2,
  BYTES_VALUE: 3,
  STRING_VALUE: 4,
  URL: 5
};

/**
 * @return {proto.ord.Data.KindCase}
 */
proto.ord.Data.prototype.getKindCase = function() {
  return /** @type {proto.ord.Data.KindCase} */(jspb.Message.computeOneofCase(this, proto.ord.Data.oneofGroups_[0]));
};



if (jspb.Message.GENERATE_TO_OBJECT) {
/**
 * Creates an object representation of this proto.
 * Field names that are reserved in JavaScript and will be renamed to pb_name.
 * Optional fields that are not set will be set to undefined.
 * To access a reserved field use, foo.pb_<name>, eg, foo.pb_default.
 * For the list of reserved names please see:
 *     net/proto2/compiler/js/internal/generator.cc#kKeyword.
 * @param {boolean=} opt_includeInstance Deprecated. whether to include the
 *     JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @return {!Object}
 */
proto.ord.Data.prototype.toObject = function(opt_includeInstance) {
  return proto.ord.Data.toObject(opt_includeInstance, this);
};


/**
 * Static version of the {@see toObject} method.
 * @param {boolean|undefined} includeInstance Deprecated. Whether to include
 *     the JSPB instance for transitional soy proto support:
 *     http://goto/soy-param-migration
 * @param {!proto.ord.Data} msg The msg instance to transform.
 * @return {!Object}
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Data.toObject = function(includeInstance, msg) {
  var f, obj = {
    floatValue: jspb.Message.getFloatingPointFieldWithDefault(msg, 1, 0.0),
    integerValue: jspb.Message.getFieldWithDefault(msg, 2, 0),
    bytesValue: msg.getBytesValue_asB64(),
    stringValue: jspb.Message.getFieldWithDefault(msg, 4, ""),
    url: jspb.Message.getFieldWithDefault(msg, 5, ""),
    description: jspb.Message.getFieldWithDefault(msg, 6, ""),
    format: jspb.Message.getFieldWithDefault(msg, 7, "")
  };

  if (includeInstance) {
    obj.$jspbMessageInstance = msg;
  }
  return obj;
};
}


/**
 * Deserializes binary data (in protobuf wire format).
 * @param {jspb.ByteSource} bytes The bytes to deserialize.
 * @return {!proto.ord.Data}
 */
proto.ord.Data.deserializeBinary = function(bytes) {
  var reader = new jspb.BinaryReader(bytes);
  var msg = new proto.ord.Data;
  return proto.ord.Data.deserializeBinaryFromReader(msg, reader);
};


/**
 * Deserializes binary data (in protobuf wire format) from the
 * given reader into the given message object.
 * @param {!proto.ord.Data} msg The message object to deserialize into.
 * @param {!jspb.BinaryReader} reader The BinaryReader to use.
 * @return {!proto.ord.Data}
 */
proto.ord.Data.deserializeBinaryFromReader = function(msg, reader) {
  while (reader.nextField()) {
    if (reader.isEndGroup()) {
      break;
    }
    var field = reader.getFieldNumber();
    switch (field) {
    case 1:
      var value = /** @type {number} */ (reader.readFloat());
      msg.setFloatValue(value);
      break;
    case 2:
      var value = /** @type {number} */ (reader.readInt32());
      msg.setIntegerValue(value);
      break;
    case 3:
      var value = /** @type {!Uint8Array} */ (reader.readBytes());
      msg.setBytesValue(value);
      break;
    case 4:
      var value = /** @type {string} */ (reader.readString());
      msg.setStringValue(value);
      break;
    case 5:
      var value = /** @type {string} */ (reader.readString());
      msg.setUrl(value);
      break;
    case 6:
      var value = /** @type {string} */ (reader.readString());
      msg.setDescription(value);
      break;
    case 7:
      var value = /** @type {string} */ (reader.readString());
      msg.setFormat(value);
      break;
    default:
      reader.skipField();
      break;
    }
  }
  return msg;
};


/**
 * Serializes the message to binary data (in protobuf wire format).
 * @return {!Uint8Array}
 */
proto.ord.Data.prototype.serializeBinary = function() {
  var writer = new jspb.BinaryWriter();
  proto.ord.Data.serializeBinaryToWriter(this, writer);
  return writer.getResultBuffer();
};


/**
 * Serializes the given message to binary data (in protobuf wire
 * format), writing to the given BinaryWriter.
 * @param {!proto.ord.Data} message
 * @param {!jspb.BinaryWriter} writer
 * @suppress {unusedLocalVariables} f is only used for nested messages
 */
proto.ord.Data.serializeBinaryToWriter = function(message, writer) {
  var f = undefined;
  f = /** @type {number} */ (jspb.Message.getField(message, 1));
  if (f != null) {
    writer.writeFloat(
      1,
      f
    );
  }
  f = /** @type {number} */ (jspb.Message.getField(message, 2));
  if (f != null) {
    writer.writeInt32(
      2,
      f
    );
  }
  f = /** @type {!(string|Uint8Array)} */ (jspb.Message.getField(message, 3));
  if (f != null) {
    writer.writeBytes(
      3,
      f
    );
  }
  f = /** @type {string} */ (jspb.Message.getField(message, 4));
  if (f != null) {
    writer.writeString(
      4,
      f
    );
  }
  f = /** @type {string} */ (jspb.Message.getField(message, 5));
  if (f != null) {
    writer.writeString(
      5,
      f
    );
  }
  f = message.getDescription();
  if (f.length > 0) {
    writer.writeString(
      6,
      f
    );
  }
  f = message.getFormat();
  if (f.length > 0) {
    writer.writeString(
      7,
      f
    );
  }
};


/**
 * optional float float_value = 1;
 * @return {number}
 */
proto.ord.Data.prototype.getFloatValue = function() {
  return /** @type {number} */ (jspb.Message.getFloatingPointFieldWithDefault(this, 1, 0.0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setFloatValue = function(value) {
  return jspb.Message.setOneofField(this, 1, proto.ord.Data.oneofGroups_[0], value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.clearFloatValue = function() {
  return jspb.Message.setOneofField(this, 1, proto.ord.Data.oneofGroups_[0], undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Data.prototype.hasFloatValue = function() {
  return jspb.Message.getField(this, 1) != null;
};


/**
 * optional int32 integer_value = 2;
 * @return {number}
 */
proto.ord.Data.prototype.getIntegerValue = function() {
  return /** @type {number} */ (jspb.Message.getFieldWithDefault(this, 2, 0));
};


/**
 * @param {number} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setIntegerValue = function(value) {
  return jspb.Message.setOneofField(this, 2, proto.ord.Data.oneofGroups_[0], value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.clearIntegerValue = function() {
  return jspb.Message.setOneofField(this, 2, proto.ord.Data.oneofGroups_[0], undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Data.prototype.hasIntegerValue = function() {
  return jspb.Message.getField(this, 2) != null;
};


/**
 * optional bytes bytes_value = 3;
 * @return {!(string|Uint8Array)}
 */
proto.ord.Data.prototype.getBytesValue = function() {
  return /** @type {!(string|Uint8Array)} */ (jspb.Message.getFieldWithDefault(this, 3, ""));
};


/**
 * optional bytes bytes_value = 3;
 * This is a type-conversion wrapper around `getBytesValue()`
 * @return {string}
 */
proto.ord.Data.prototype.getBytesValue_asB64 = function() {
  return /** @type {string} */ (jspb.Message.bytesAsB64(
      this.getBytesValue()));
};


/**
 * optional bytes bytes_value = 3;
 * Note that Uint8Array is not supported on all browsers.
 * @see http://caniuse.com/Uint8Array
 * This is a type-conversion wrapper around `getBytesValue()`
 * @return {!Uint8Array}
 */
proto.ord.Data.prototype.getBytesValue_asU8 = function() {
  return /** @type {!Uint8Array} */ (jspb.Message.bytesAsU8(
      this.getBytesValue()));
};


/**
 * @param {!(string|Uint8Array)} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setBytesValue = function(value) {
  return jspb.Message.setOneofField(this, 3, proto.ord.Data.oneofGroups_[0], value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.clearBytesValue = function() {
  return jspb.Message.setOneofField(this, 3, proto.ord.Data.oneofGroups_[0], undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Data.prototype.hasBytesValue = function() {
  return jspb.Message.getField(this, 3) != null;
};


/**
 * optional string string_value = 4;
 * @return {string}
 */
proto.ord.Data.prototype.getStringValue = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 4, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setStringValue = function(value) {
  return jspb.Message.setOneofField(this, 4, proto.ord.Data.oneofGroups_[0], value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.clearStringValue = function() {
  return jspb.Message.setOneofField(this, 4, proto.ord.Data.oneofGroups_[0], undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Data.prototype.hasStringValue = function() {
  return jspb.Message.getField(this, 4) != null;
};


/**
 * optional string url = 5;
 * @return {string}
 */
proto.ord.Data.prototype.getUrl = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 5, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setUrl = function(value) {
  return jspb.Message.setOneofField(this, 5, proto.ord.Data.oneofGroups_[0], value);
};


/**
 * Clears the field making it undefined.
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.clearUrl = function() {
  return jspb.Message.setOneofField(this, 5, proto.ord.Data.oneofGroups_[0], undefined);
};


/**
 * Returns whether this field is set.
 * @return {boolean}
 */
proto.ord.Data.prototype.hasUrl = function() {
  return jspb.Message.getField(this, 5) != null;
};


/**
 * optional string description = 6;
 * @return {string}
 */
proto.ord.Data.prototype.getDescription = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 6, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setDescription = function(value) {
  return jspb.Message.setProto3StringField(this, 6, value);
};


/**
 * optional string format = 7;
 * @return {string}
 */
proto.ord.Data.prototype.getFormat = function() {
  return /** @type {string} */ (jspb.Message.getFieldWithDefault(this, 7, ""));
};


/**
 * @param {string} value
 * @return {!proto.ord.Data} returns this
 */
proto.ord.Data.prototype.setFormat = function(value) {
  return jspb.Message.setProto3StringField(this, 7, value);
};


goog.object.extend(exports, proto.ord);
