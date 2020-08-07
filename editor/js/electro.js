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

goog.provide('ord.electro');

goog.require('proto.ord.ElectrochemistryConditions');
goog.require('proto.ord.ElectrochemistryConditions.Measurement');

// Freely create radio button groups by generating new input names.
ord.electro.radioGroupCounter = 0;

ord.electro.load = function(electro) {
  const type = electro.getElectrochemistryType();
  if (type) {
    setSelector($('#electro_type'), type.getType());
    $('#electro_details').text(type.getDetails());
  }
  writeMetric('#electro_current', electro.getCurrent());
  writeMetric('#electro_voltage', electro.getVoltage());
  $('#electro_anode').text(electro.getAnodeMaterial());
  $('#electro_cathode').text(electro.getCathodeMaterial());
  writeMetric('#electro_separation', electro.getElectrodeSeparation());

  const cell = electro.getCell();
  if (cell) {
    setSelector($('#electro_cell_type'), cell.getType());
    $('#electro_cell_details').text(cell.getDetails());
  }
  electro.getMeasurementsList().forEach(function(measurement) {
    const node = ord.electro.addMeasurement();
    ord.electro.loadMeasurement(node, measurement);
  });
};

ord.electro.loadMeasurement = function(node, measurement) {
  const time = measurement.getTime();
  if (time) {
    writeMetric('.electro_measurement_time', time, node);
  }
  const current = measurement.getCurrent();
  const voltage = measurement.getVoltage();
  if (current) {
    writeMetric('.electro_measurement_current', current, node);
    $('input[value=\'current\']', node).prop('checked', true);
    $('.electro_measurement_current_fields', node).show();
    $('.electro_measurement_voltage_fields', node).hide();
  }
  if (voltage) {
    $('input[value=\'voltage\']', node).prop('checked', true);
    writeMetric('.electro_measurement_voltage', voltage, node);
    $('.electro_measurement_current_fields', node).hide();
    $('.electro_measurement_voltage_fields', node).show();
  }
};

ord.electro.unload = function() {
  const electro = new proto.ord.ElectrochemistryConditions();

  const type = new proto.ord.ElectrochemistryConditions.ElectrochemistryType();
  type.setType(getSelector($('#electro_type')));
  type.setDetails($('#electro_details').text());
  if (!isEmptyMessage(type)) {
    electro.setElectrochemistryType(type);
  }

  const current = readMetric('#electro_current', new proto.ord.Current());
  if (!isEmptyMessage(current)) {
    electro.setCurrent(current);
  }
  const voltage = readMetric('#electro_voltage', new proto.ord.Voltage());
  if (!isEmptyMessage(voltage)) {
    electro.setVoltage(voltage);
  }
  electro.setAnodeMaterial($('#electro_anode').text());
  electro.setCathodeMaterial($('#electro_cathode').text());
  const electrodeSeparation =
      readMetric('#electro_separation', new proto.ord.Length());
  if (!isEmptyMessage(electrodeSeparation)) {
    electro.setElectrodeSeparation(electrodeSeparation);
  }

  const cell = new proto.ord.ElectrochemistryConditions.ElectrochemistryCell();
  cell.setType(getSelector($('#electro_cell_type')));
  cell.setDetails($('#electro_cell_details').text());
  if (!isEmptyMessage(cell)) {
    electro.setCell(cell);
  }

  const measurements = []
  $('.electro_measurement').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const measurement = ord.electro.unloadMeasurement(node);
      if (!isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  electro.setMeasurementsList(measurements);
  return electro;
};

ord.electro.unloadMeasurement = function(node) {
  const measurement = new proto.ord.ElectrochemistryConditions.Measurement();
  const time =
      readMetric('.electro_measurement_time', new proto.ord.Time(), node);
  if (!isEmptyMessage(time)) {
    measurement.setTime(time);
  }

  if ($('.electro_measurement_current', node).is(':checked')) {
    const current = readMetric(
        '.electro_measurement_current', new proto.ord.Current(), node);
    if (!isEmptyMessage(current)) {
      measurement.setCurrent(current);
    }
  }
  if ($('.electro_measurement_voltage', node).is(':checked')) {
    const voltage = readMetric(
        '.electro_measurement_voltage', new proto.ord.Voltage(), node);
    if (!isEmptyMessage(voltage)) {
      measurement.setVoltage(voltage);
    }
  }
  return measurement;
};

ord.electro.addMeasurement = function() {
  const node =
      addSlowly('#electro_measurement_template', '#electro_measurements');

  const metricButtons = $('input', node);
  metricButtons.attr('name', 'electro_' + ord.electro.radioGroupCounter++);
  metricButtons.change(function() {
    if (this.value == 'current') {
      $('.electro_measurement_current_fields', node).show();
      $('.electro_measurement_voltage_fields', node).hide();
    }
    if (this.value == 'voltage') {
      $('.electro_measurement_current_fields', node).hide();
      $('.electro_measurement_voltage_fields', node).show();
    }
  });

  return node;
};

ord.electro.validateElectro = function(node, validateNode) {
  const electro = ord.electro.unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(electro, 'ElectrochemistryConditions', validateNode);
};