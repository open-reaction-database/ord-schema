/**
 * Copyright 2020 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

goog.module('ord.pressure');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  addMeasurement,
  validatePressure
};

goog.require('proto.ord.Pressure');
goog.require('proto.ord.PressureConditions');
goog.require('proto.ord.PressureConditions.Measurement');
goog.require('proto.ord.Time');

/**
 * Populates the pressure conditions section in the form.
 * @param {!proto.ord.PressureConditions} pressure
 */
function load(pressure) {
  const control = pressure.getControl();
  if (control) {
    ord.reaction.setSelector($('#pressure_control_type'), control.getType());
    $('#pressure_control_details').text(control.getDetails());
  }
  const measurements = pressure.getMeasurementsList();
  measurements.forEach(function(measurement) {
    const node = addMeasurement();
    loadMeasurement(measurement, node);
  });
  const setpoint = pressure.getSetpoint();
  ord.reaction.writeMetric('#pressure_setpoint', setpoint);

  const atmosphere = pressure.getAtmosphere();
  if (atmosphere) {
    ord.reaction.setSelector(
        $('#pressure_atmosphere_type'), atmosphere.getType());
    $('#pressure_atmosphere_details').text(atmosphere.getDetails());
  }
}

/**
 * Populates a pressure measurement section in the form.
 * @param {!proto.ord.PressureConditions.Measurement} measurement
 * @param {!Node} node The target div.
 */
function loadMeasurement(measurement, node) {
  const type = measurement.getType();
  ord.reaction.setSelector($('.pressure_measurement_type', node), type);
  $('.pressure_measurement_details', node).text(measurement.getDetails());

  const pressure = measurement.getPressure();
  ord.reaction.writeMetric('.pressure_measurement_pressure', pressure, node);

  const time = measurement.getTime();
  ord.reaction.writeMetric('.pressure_measurement_time', time, node);
}

/**
 * Fetches pressure conditions from the form.
 * @return {!proto.ord.PressureConditions}
 */
function unload() {
  const pressure = new proto.ord.PressureConditions();

  const control = new proto.ord.PressureConditions.PressureControl();
  control.setType(ord.reaction.getSelector($('#pressure_control_type')));
  control.setDetails($('#pressure_control_details').text());
  if (!ord.reaction.isEmptyMessage(control)) {
    pressure.setControl(control);
  }

  const setpoint =
      ord.reaction.readMetric('#pressure_setpoint', new proto.ord.Pressure());
  if (!ord.reaction.isEmptyMessage(setpoint)) {
    pressure.setSetpoint(setpoint);
  }

  const atmosphere = new proto.ord.PressureConditions.Atmosphere();
  atmosphere.setType(ord.reaction.getSelector('#pressure_atmosphere_type'));
  atmosphere.setDetails($('#pressure_atmosphere_details').text());
  if (!ord.reaction.isEmptyMessage(atmosphere)) {
    pressure.setAtmosphere(atmosphere);
  }

  const measurements = [];
  $('.pressure_measurement').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const measurement = unloadMeasurement(node);
      if (!ord.reaction.isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  pressure.setMeasurementsList(measurements);

  return pressure;
}

/**
 * Fetches a pressure measurement from the form.
 * @param {!Node} node The div of the measurement to fetch.
 * @return {!proto.ord.PressureConditions.Measurement}
 */
function unloadMeasurement(node) {
  const measurement = new proto.ord.PressureConditions.Measurement();
  const type = ord.reaction.getSelector($('.pressure_measurement_type', node));
  measurement.setType(type);
  const details = $('.pressure_measurement_details', node).text();
  measurement.setDetails(details);
  const pressure = ord.reaction.readMetric(
      '.pressure_measurement_pressure', new proto.ord.Pressure(), node);
  if (!ord.reaction.isEmptyMessage(pressure)) {
    measurement.setPressure(pressure);
  }
  const time = ord.reaction.readMetric(
      '.pressure_measurement_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(time)) {
    measurement.setTime(time);
  }

  return measurement;
}

/**
 * Adds a pressure measurement section to the form.
 * @return {!Node} The node of the new measurement div.
 */
function addMeasurement() {
  return ord.reaction.addSlowly(
      '#pressure_measurement_template', '#pressure_measurements');
}

/**
 * Validates pressure conditions defined in the form.
 * @param {!Node} node The node containing the pressure conditions div.
 * @param {!Node} validateNode The target div for validation results.
 */
function validatePressure(node, validateNode) {
  const pressure = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(pressure, 'PressureConditions', validateNode);
}
