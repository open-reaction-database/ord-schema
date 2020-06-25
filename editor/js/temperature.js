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

goog.provide('ord.temperature');

goog.require('proto.ord.Temperature');
goog.require('proto.ord.TemperatureConditions');
goog.require('proto.ord.TemperatureConditions.Measurement');
goog.require('proto.ord.TemperatureConditions.TemperatureControl');
goog.require('proto.ord.Time');

ord.temperature.load = function (temperature) {
  const control = temperature.getControl();
  if (control) {
    setSelector($('#temperature_control'), control.getType());
    $('#temperature_control_details').text(control.getDetails());
  }
  const measurements = temperature.getMeasurementsList();
  measurements.forEach(function (measurement) {
    const node = ord.temperature.addMeasurement();
    ord.temperature.loadMeasurement(measurement, node);
  });
  const setpoint = temperature.getSetpoint();
  writeMetric('#temperature_setpoint', setpoint);
};

ord.temperature.loadMeasurement = function (measurement, node) {
  const type = measurement.getType();
  setSelector($('.temperature_measurement_type', node), type);
  $('.temperature_measurement_details', node).text(measurement.getDetails());

  const temperature = measurement.getTemperature();
  writeMetric('.temperature_measurement_temperature', temperature, node);

  const time = measurement.getTime();
  writeMetric('.temperature_measurement_time', time, node);
};

ord.temperature.unload = function () {
  const temperature = new proto.ord.TemperatureConditions();

  const control = new proto.ord.TemperatureConditions.TemperatureControl();
  control.setType(getSelector($('#temperature_control')));
  control.setDetails($('#temperature_control_details').text());
  temperature.setControl(control);

  const setpoint =
      readMetric('#temperature_setpoint', new proto.ord.Temperature());
  temperature.setSetpoint(setpoint);

  const measurements = [];
  $('.temperature_measurement').each(function (index, node) {
    node = $(node)
    if (!node.attr('id')) {
      const measurement = ord.temperature.unloadMeasurement(node);
      measurements.push(measurement);
    }
  });
  temperature.setMeasurementsList(measurements);
  return temperature;
};

ord.temperature.unloadMeasurement = function (node) {
  const measurement = new proto.ord.TemperatureConditions.Measurement();
  const type = getSelector($('.temperature_measurement_type', node));
  measurement.setType(type);
  const details = $('.temperature_measurement_details', node).text();
  measurement.setDetails(details);
  const temperature = readMetric(
      '.temperature_measurement_temperature',
      new proto.ord.Temperature(),
      node);
  measurement.setTemperature(temperature);
  const time =
      readMetric('.temperature_measurement_time', new proto.ord.Time(), node);
  measurement.setTime(time);
  return measurement;
};

ord.temperature.addMeasurement = function () {
  return addSlowly('#temperature_measurement_template', '#temperature_measurements');
};

ord.temperature.validateTemperature = function(node, validateNode) {
  const temperature = ord.temperature.unload();
  if (typeof validateNode === 'undefined') {
    validateNode = $('.validate', node).first();
  }
  validate(temperature, "TemperatureConditions", validateNode);
};