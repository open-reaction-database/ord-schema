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

goog.module('ord.electro');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  addMeasurement,
  validateElectro
};

goog.require('proto.ord.ElectrochemistryConditions');
goog.require('proto.ord.ElectrochemistryConditions.Measurement');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the electrochemistry conditions section in the form.
 * @param {!proto.ord.ElectrochemistryConditions} electro
 */
function load(electro) {
  const type = electro.getElectrochemistryType();
  if (type) {
    ord.reaction.setSelector($('#electro_type'), type.getType());
    $('#electro_details').text(type.getDetails());
  }
  ord.reaction.writeMetric('#electro_current', electro.getCurrent());
  ord.reaction.writeMetric('#electro_voltage', electro.getVoltage());
  $('#electro_anode').text(electro.getAnodeMaterial());
  $('#electro_cathode').text(electro.getCathodeMaterial());
  ord.reaction.writeMetric(
      '#electro_separation', electro.getElectrodeSeparation());

  const cell = electro.getCell();
  if (cell) {
    ord.reaction.setSelector($('#electro_cell_type'), cell.getType());
    $('#electro_cell_details').text(cell.getDetails());
  }
  electro.getMeasurementsList().forEach(function(measurement) {
    const node = addMeasurement();
    loadMeasurement(node, measurement);
  });
}

/**
 * Adds and populates an electrochemistry measurement section in the form.
 * @param {!Node} node The target div.
 * @param {!proto.ord.ElectrochemistryConditions.Measurement} measurement
 */
function loadMeasurement(node, measurement) {
  const time = measurement.getTime();
  if (time) {
    ord.reaction.writeMetric('.electro_measurement_time', time, node);
  }
  const current = measurement.getCurrent();
  const voltage = measurement.getVoltage();
  if (current) {
    ord.reaction.writeMetric('.electro_measurement_current', current, node);
    $('input[value=\'current\']', node).prop('checked', true);
    $('.electro_measurement_current_fields', node).show();
    $('.electro_measurement_voltage_fields', node).hide();
  }
  if (voltage) {
    $('input[value=\'voltage\']', node).prop('checked', true);
    ord.reaction.writeMetric('.electro_measurement_voltage', voltage, node);
    $('.electro_measurement_current_fields', node).hide();
    $('.electro_measurement_voltage_fields', node).show();
  }
}

/**
 * Fetches the electrochemistry conditions defined in the form.
 * @return {!proto.ord.ElectrochemistryConditions}
 */
function unload() {
  const electro = new proto.ord.ElectrochemistryConditions();

  const type = new proto.ord.ElectrochemistryConditions.ElectrochemistryType();
  type.setType(ord.reaction.getSelector($('#electro_type')));
  type.setDetails($('#electro_details').text());
  if (!ord.reaction.isEmptyMessage(type)) {
    electro.setElectrochemistryType(type);
  }

  const current =
      ord.reaction.readMetric('#electro_current', new proto.ord.Current());
  if (!ord.reaction.isEmptyMessage(current)) {
    electro.setCurrent(current);
  }
  const voltage =
      ord.reaction.readMetric('#electro_voltage', new proto.ord.Voltage());
  if (!ord.reaction.isEmptyMessage(voltage)) {
    electro.setVoltage(voltage);
  }
  electro.setAnodeMaterial($('#electro_anode').text());
  electro.setCathodeMaterial($('#electro_cathode').text());
  const electrodeSeparation =
      ord.reaction.readMetric('#electro_separation', new proto.ord.Length());
  if (!ord.reaction.isEmptyMessage(electrodeSeparation)) {
    electro.setElectrodeSeparation(electrodeSeparation);
  }

  const cell = new proto.ord.ElectrochemistryConditions.ElectrochemistryCell();
  cell.setType(ord.reaction.getSelector($('#electro_cell_type')));
  cell.setDetails($('#electro_cell_details').text());
  if (!ord.reaction.isEmptyMessage(cell)) {
    electro.setCell(cell);
  }

  const measurements = [];
  $('.electro_measurement').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const measurement = unloadMeasurement(node);
      if (!ord.reaction.isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  electro.setMeasurementsList(measurements);
  return electro;
}

/**
 * Fetches an electrochemistry measurement from the form.
 * @param {!Node} node Root node of the measurement.
 * @return {!proto.ord.ElectrochemistryConditions.Measurement}
 */
function unloadMeasurement(node) {
  const measurement = new proto.ord.ElectrochemistryConditions.Measurement();
  const time = ord.reaction.readMetric(
      '.electro_measurement_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(time)) {
    measurement.setTime(time);
  }

  if ($('.electro_measurement_current', node).is(':checked')) {
    const current = ord.reaction.readMetric(
        '.electro_measurement_current', new proto.ord.Current(), node);
    if (!ord.reaction.isEmptyMessage(current)) {
      measurement.setCurrent(current);
    }
  }
  if ($('.electro_measurement_voltage', node).is(':checked')) {
    const voltage = ord.reaction.readMetric(
        '.electro_measurement_voltage', new proto.ord.Voltage(), node);
    if (!ord.reaction.isEmptyMessage(voltage)) {
      measurement.setVoltage(voltage);
    }
  }
  return measurement;
}

/**
 * Adds an electrochemistry measurement section to the form.
 * @return {!Node} The newly added parent node for the measurement.
 */
function addMeasurement() {
  const node = ord.reaction.addSlowly(
      '#electro_measurement_template', '#electro_measurements');

  const metricButtons = $('input', node);
  metricButtons.attr('name', 'electro_' + radioGroupCounter++);
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
}

/**
 * Validates the electrochemistry conditions defined in the form.
 * @param {!Node} node Root node for the electrochemistry conditions.
 * @param {?Node} validateNode Target node for validation results.
 */
function validateElectro(node, validateNode) {
  const electro = unload();
  ord.reaction.validate(
      electro, 'ElectrochemistryConditions', node, validateNode);
}
