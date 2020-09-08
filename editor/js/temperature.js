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

goog.module('ord.temperature');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  addMeasurement,
  validateTemperature
};

goog.require('proto.ord.Temperature');
goog.require('proto.ord.TemperatureConditions');
goog.require('proto.ord.TemperatureConditions.Measurement');
goog.require('proto.ord.TemperatureConditions.TemperatureControl');
goog.require('proto.ord.Time');

/**
 * Populates the temperature conditions section in the form.
 * @param {!proto.ord.TemperatureConditions} temperature
 */
function load(temperature) {
  const control = temperature.getControl();
  if (control) {
    ord.reaction.setSelector($('#temperature_control'), control.getType());
    $('#temperature_control_details').text(control.getDetails());
  }
  const measurements = temperature.getMeasurementsList();
  measurements.forEach(function(measurement) {
    const node = addMeasurement();
    loadMeasurement(measurement, node);
  });
  const setpoint = temperature.getSetpoint();
  ord.reaction.writeMetric('#temperature_setpoint', setpoint);
}

/**
 * Populates a temperature measurement section in the form.
 * @param {!proto.ord.TemperatureConditions.Measurement} measurement
 * @param {!Node} node The target div.
 */
function loadMeasurement(measurement, node) {
  const type = measurement.getType();
  ord.reaction.setSelector($('.temperature_measurement_type', node), type);
  $('.temperature_measurement_details', node).text(measurement.getDetails());

  const temperature = measurement.getTemperature();
  ord.reaction.writeMetric(
      '.temperature_measurement_temperature', temperature, node);

  const time = measurement.getTime();
  ord.reaction.writeMetric('.temperature_measurement_time', time, node);
}

/**
 * Fetches temperature conditions from the form.
 * @return {!proto.ord.TemperatureConditions}
 */
function unload() {
  const temperature = new proto.ord.TemperatureConditions();

  const control = new proto.ord.TemperatureConditions.TemperatureControl();
  control.setType(ord.reaction.getSelector($('#temperature_control')));
  control.setDetails($('#temperature_control_details').text());
  if (!ord.reaction.isEmptyMessage(control)) {
    temperature.setControl(control);
  }

  const setpoint = ord.reaction.readMetric(
      '#temperature_setpoint', new proto.ord.Temperature());
  if (!ord.reaction.isEmptyMessage(setpoint)) {
    temperature.setSetpoint(setpoint);
  }

  const measurements = [];
  $('.temperature_measurement').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const measurement = unloadMeasurement(node);
      if (!ord.reaction.isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  temperature.setMeasurementsList(measurements);
  return temperature;
}

/**
 * Fetches a temperature measurement from the form.
 * @param {!Node} node The div of the measurement to fetch.
 * @return {!proto.ord.TemperatureConditions.Measurement}
 */
function unloadMeasurement(node) {
  const measurement = new proto.ord.TemperatureConditions.Measurement();
  const type =
      ord.reaction.getSelector($('.temperature_measurement_type', node));
  measurement.setType(type);
  const details = $('.temperature_measurement_details', node).text();
  measurement.setDetails(details);
  const temperature = ord.reaction.readMetric(
      '.temperature_measurement_temperature', new proto.ord.Temperature(),
      node);
  if (!ord.reaction.isEmptyMessage(temperature)) {
    measurement.setTemperature(temperature);
  }
  const time = ord.reaction.readMetric(
      '.temperature_measurement_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(time)) {
    measurement.setTime(time);
  }
  return measurement;
}

/**
 * Adds a temperature measurment section to the form.
 * @return {!Node} The node of the new measurement div.
 */
function addMeasurement() {
  return ord.reaction.addSlowly(
      '#temperature_measurement_template', '#temperature_measurements');
}

/**
 * Validates temperature conditions defined in the form.
 * @param {!Node} node The node containing the temperature conditions div.
 * @param {!Node} validateNode The target div for validation results.
 */
function validateTemperature(node, validateNode) {
  const temperature = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(temperature, 'TemperatureConditions', validateNode);
}
