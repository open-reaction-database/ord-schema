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

goog.provide('ord.pressure');

goog.require('proto.ord.Pressure');
goog.require('proto.ord.PressureConditions');
goog.require('proto.ord.PressureConditions.Measurement');
goog.require('proto.ord.Time');

ord.pressure.load = function (pressure) {
  const control = pressure.getControl();
  if (control) {
    setSelector($('#pressure_control_type'), control.getType());
    $('#pressure_control_details').text(control.getDetails());
  }
  const measurements = pressure.getMeasurementsList();
  measurements.forEach(function (measurement) {
    const node = ord.pressure.addMeasurement();
    ord.pressure.loadMeasurement(measurement, node);
  });
  const setpoint = pressure.getSetpoint();
  writeMetric('#pressure_setpoint', setpoint);

  const atmosphere = pressure.getAtmosphere();
  if (atmosphere) {
    setSelector($('#pressure_atmosphere_type'), atmosphere.getType());
    $('#pressure_atmosphere_details').text(atmosphere.getDetails());
  }
};

ord.pressure.loadMeasurement = function (measurement, node) {
  const type = measurement.getType();
  setSelector($('.pressure_measurement_type', node), type);
  $('.pressure_measurement_details', node).text(measurement.getDetails());

  const pressure = measurement.getPressure();
  writeMetric('.pressure_measurement_pressure', pressure, node);

  const time = measurement.getTime();
  writeMetric('.pressure_measurement_time', time, node);
};

ord.pressure.unload = function () {
  const pressure = new proto.ord.PressureConditions();

  const control = new proto.ord.PressureConditions.PressureControl();
  control.setType(getSelector($('#pressure_control_type')));
  control.setDetails($('#pressure_control_details').text());
  if (!isEmptyMessage(control)) {
    pressure.setControl(control);
  }

  const setpoint = readMetric('#pressure_setpoint', new proto.ord.Pressure());
  if (!isEmptyMessage(setpoint)) {
    pressure.setSetpoint(setpoint);
  }

  const atmosphere = new proto.ord.PressureConditions.Atmosphere();
  atmosphere.setType(getSelector('#pressure_atmosphere_type'));
  atmosphere.setDetails($('#pressure_atmosphere_details').text());
  if (!isEmptyMessage(atmosphere)) {
    pressure.setAtmosphere(atmosphere);
  }

  const measurements = [];
  $('.pressure_measurement').each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const measurement = ord.pressure.unloadMeasurement(node);
      if (!isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  pressure.setMeasurementsList(measurements);

  return pressure;
};

ord.pressure.unloadMeasurement = function (node) {
  const measurement = new proto.ord.PressureConditions.Measurement();
  const type = getSelector($('.pressure_measurement_type', node));
  measurement.setType(type);
  const details = $('.pressure_measurement_details', node).text();
  measurement.setDetails(details);
  const pressure = readMetric(
      '.pressure_measurement_pressure', new proto.ord.Pressure(), node);
  if (!isEmptyMessage(pressure)) {
    measurement.setPressure(pressure);
  }
  const time =
      readMetric('.pressure_measurement_time', new proto.ord.Time(), node);
  if (!isEmptyMessage(time)) {
    measurement.setTime(time);
  }

  return measurement;
};

ord.pressure.addMeasurement = function () {
  return addSlowly('#pressure_measurement_template', '#pressure_measurements');
};

ord.pressure.validatePressure = function(node, validateNode) {
  const pressure = ord.pressure.unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(pressure, 'PressureConditions', validateNode);
};