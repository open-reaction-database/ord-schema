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

goog.module('ord.workups');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  add,
  addMeasurement,
  validateWorkup
};

goog.require('ord.inputs');
goog.require('proto.ord.ReactionWorkup');

/**
 * Adds and populates the reaction workup sections in the form.
 * @param {!Array<!proto.ord.ReactionWorkup>} workups
 */
function load(workups) {
  workups.forEach(workup => loadWorkup(workup));
}

/**
 * Adds and populates a reaction workup section in the form.
 * @param {!proto.ord.ReactionWorkup} workup
 */
function loadWorkup(workup) {
  const node = add();
  ord.reaction.setSelector($('.workup_type', node), workup.getType());
  $('.workup_details', node).text(workup.getDetails());
  const duration = workup.getDuration();
  if (duration) {
    ord.reaction.writeMetric('.workup_duration', duration, node);
  }

  const input = workup.getInput();
  if (input) {
    ord.inputs.loadInputUnnamed($('.workup_input', node), input);
  }

  const temperature = workup.getTemperature();
  if (temperature) {
    const control = temperature.getControl();
    if (control) {
      ord.reaction.setSelector(
          $('.workup_temperature_control_type', node), control.getType());
      $('.workup_temperature_details', node).text(control.getDetails());
    }
    const setpoint = temperature.getSetpoint();
    if (setpoint) {
      ord.reaction.writeMetric('.workup_temperature_setpoint', setpoint, node);
    }

    temperature.getMeasurementsList().forEach(
        measurement => loadMeasurement(node, measurement));
  }

  $('.workup_keep_phase', node).text(workup.getKeepPhase());

  const stirring = workup.getStirring();
  if (stirring) {
    const method = stirring.getMethod();
    if (method) {
      ord.reaction.setSelector(
          $('.workup_stirring_method_type', node), method.getType());
      $('.workup_stirring_method_details', node).text(method.getDetails());
    }
    const rate = stirring.getRate();
    if (rate) {
      ord.reaction.setSelector(
          $('.workup_stirring_rate_type', node), rate.getType());
      $('.workup_stirring_rate_details', node).text(rate.getDetails());
      const rpm = rate.getRpm();
      if (rpm != 0) {
        $('.workup_stirring_rate_rpm', node).text(rpm);
      }
    }
  }
  if (workup.hasTargetPh()) {
    $('.workup_target_ph', node).text(workup.getTargetPh());
  }
  ord.reaction.setOptionalBool(
      $('.workup_automated', node),
      workup.hasIsAutomated() ? workup.getIsAutomated() : null);
}

/**
 * Loads a measurement into the given node in a workup.
 * @param {!Node} workupNode The div corresponding to the workup whose fields
 *     should be updated.
 * @param {!proto.ord.TemperatureConditions.Measurement} measurement
 */
function loadMeasurement(workupNode, measurement) {
  const node = addMeasurement(workupNode);
  ord.reaction.setSelector(
      $('.workup_temperature_measurement_type', node), measurement.getType());
  $('.workup_temperature_measurement_details', node)
      .text(measurement.getDetails());
  const time = measurement.getTime();
  if (time) {
    ord.reaction.writeMetric(
        '.workup_temperature_measurement_time', time, node);
  }
  const temperature = measurement.getTemperature();
  if (temperature) {
    ord.reaction.writeMetric(
        '.workup_temperature_measurement_temperature', temperature, node);
  }
}

/**
 * Fetches a list of workups defined in the form.
 * @return {!Array<!proto.ord.ReactionWorkup>} workups
 */
function unload() {
  const workups = [];
  $('.workup').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const workup = unloadWorkup(node);
      if (!ord.reaction.isEmptyMessage(workup)) {
        workups.push(workup);
      }
    }
  });
  return workups;
}

/**
 * Fetches a single workup from the form.
 * @param {!Node} node The div corresponding to the workup to fetch.
 * @return {!proto.ord.ReactionWorkup}
 */
function unloadWorkup(node) {
  const workup = new proto.ord.ReactionWorkup();

  workup.setType(ord.reaction.getSelector($('.workup_type', node)));

  workup.setDetails($('.workup_details', node).text());

  const duration =
      ord.reaction.readMetric('.workup_duration', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(duration)) {
    workup.setDuration(duration);
  }

  const input = ord.inputs.unloadInputUnnamed(node);
  if (!ord.reaction.isEmptyMessage(input)) {
    workup.setInput(input);
  }

  const control = new proto.ord.TemperatureConditions.TemperatureControl();
  control.setType(
      ord.reaction.getSelector($('.workup_temperature_control_type', node)));
  control.setDetails($('.workup_temperature_details', node).text());

  const temperature = new proto.ord.TemperatureConditions();
  if (!ord.reaction.isEmptyMessage(control)) {
    temperature.setControl(control);
  }

  const setpoint = ord.reaction.readMetric(
      '.workup_temperature_setpoint', new proto.ord.Temperature(), node);
  if (!ord.reaction.isEmptyMessage(setpoint)) {
    temperature.setSetpoint(setpoint);
  }

  const measurements = [];
  const measurementNodes = $('.workup_temperature_measurement', node);
  measurementNodes.each(function(index, measurementNode) {
    measurementNode = $(measurementNode);
    if (!measurementNode.attr('id')) {
      // Not a template.
      const measurement = unloadMeasurement(measurementNode);
      if (!ord.reaction.isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  temperature.setMeasurementsList(measurements);
  if (!ord.reaction.isEmptyMessage(temperature)) {
    workup.setTemperature(temperature);
  }

  workup.setKeepPhase($('.workup_keep_phase', node).text());

  const stirring = new proto.ord.StirringConditions();

  const method = new proto.ord.StirringConditions.StirringMethod();
  method.setType(
      ord.reaction.getSelector($('.workup_stirring_method_type', node)));
  method.setDetails($('.workup_stirring_method_details').text());
  if (!ord.reaction.isEmptyMessage(method)) {
    stirring.setMethod(method);
  }

  const rate = new proto.ord.StirringConditions.StirringRate();
  rate.setType(ord.reaction.getSelector($('.workup_stirring_rate_type', node)));
  rate.setDetails($('.workup_stirring_rate_details').text());
  const rpm = parseFloat($('.workup_stirring_rate_rpm', node).text());
  if (!isNaN(rpm)) {
    rate.setRpm(rpm);
  }
  if (!ord.reaction.isEmptyMessage(rate)) {
    stirring.setRate(rate);
  }

  if (!ord.reaction.isEmptyMessage(stirring)) {
    workup.setStirring(stirring);
  }

  const targetPh = parseFloat($('.workup_target_ph', node).text());
  if (!isNaN(targetPh)) {
    workup.setTargetPh(targetPh);
  }
  workup.setIsAutomated(
      ord.reaction.getOptionalBool($('.workup_automated', node)));
  return workup;
}

/**
 * Fetches a single workup temperature measurement from the form.
 * @param {!Node} node The div corresponding to the measurement to fetch.
 * @return {!proto.ord.TemperatureConditions.Measurement}
 */
function unloadMeasurement(node) {
  const measurement = new proto.ord.TemperatureConditions.Measurement();
  measurement.setType(ord.reaction.getSelector(
      $('.workup_temperature_measurement_type', node)));
  measurement.setDetails(
      $('.workup_temperature_measurement_details', node).text());
  const time = ord.reaction.readMetric(
      '.workup_temperature_measurement_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(time)) {
    measurement.setTime(time);
  }
  const temperature = ord.reaction.readMetric(
      '.workup_temperature_measurement_temperature',
      new proto.ord.Temperature(), node);
  if (!ord.reaction.isEmptyMessage(temperature)) {
    measurement.setTemperature(temperature);
  }
  return measurement;
}

/**
 * Adds a new reaction workup section to the form.
 * @return {!Node} The newly added parent node for the reaction workup.
 */
function add() {
  const workupNode = ord.reaction.addSlowly('#workup_template', '#workups');
  const inputNode = $('.workup_input', workupNode);
  // The template for ReactionWorkup.input is taken from Reaction.inputs.
  const workupInputNode = ord.inputs.add(inputNode);
  // Workup inputs start collapsed by default.
  workupInputNode.find('.collapse').trigger('click');
  // Temperature conditions and stirring fields also start collapsed.
  workupNode.find('.workup_temperature').trigger('click');
  workupNode.find('.workup_temperature_measurements_wrap').trigger('click');
  workupNode.find('.workup_stirring').trigger('click');
  // Unlike Reaction.inputs, this ReactionInput has no name.
  $('.input_name_label', inputNode).hide();
  $('.input_name', inputNode).hide();
  // Unlike Reaction.inputs, this ReactionInput is not repeated.
  $('.remove', inputNode).hide();

  // Add live validation handling.
  ord.reaction.addChangeHandler(workupNode, () => {
    validateWorkup(workupNode);
  });

  return workupNode;
}

/**
 * Adds a new measurement section to the current workup in the form.
 * @param {!Node} node The workup div where the new measurement should be added.
 * @return {!Node} The node of the new measurement div.
 */
function addMeasurement(node) {
  return ord.reaction.addSlowly(
      '#workup_temperature_measurement_template',
      $('.workup_temperature_measurements', node));
}

/**
 * Validates a workup as defined in the form.
 * @param {!Node} node The div containing to the workup in the form.
 * @param {?Node} validateNode The target div for validation results.
 */
function validateWorkup(node, validateNode) {
  const workup = unloadWorkup(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(workup, 'ReactionWorkup', validateNode);
}
