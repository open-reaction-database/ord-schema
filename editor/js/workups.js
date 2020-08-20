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

goog.provide('ord.workups');

goog.require('ord.inputs');
goog.require('proto.ord.ReactionWorkup');

ord.workups.load = function(workups) {
  workups.forEach(workup => ord.workups.loadWorkup(workup));
};

ord.workups.loadWorkup = function(workup) {
  const node = ord.workups.add();
  setSelector($('.workup_type', node), workup.getType());
  $('.workup_details', node).text(workup.getDetails());
  const duration = workup.getDuration();
  if (duration) {
    writeMetric('.workup_duration', duration, node);
  }

  const input = workup.getInput();
  if (input) {
    ord.inputs.loadInputUnnamed($('.workup_input', node), input);
  }

  const temperature = workup.getTemperature();
  if (temperature) {
    const control = temperature.getControl();
    if (control) {
      setSelector(
          $('.workup_temperature_control_type', node), control.getType());
      $('.workup_temperature_details', node).text(control.getDetails());
    }
    const setpoint = temperature.getSetpoint();
    if (setpoint) {
      writeMetric('.workup_temperature_setpoint', setpoint, node);
    }

    temperature.getMeasurementsList().forEach(
        measurement => ord.workups.loadMeasurement(node, measurement));
  }

  $('.workup_keep_phase', node).text(workup.getKeepPhase());

  const stirring = workup.getStirring();
  if (stirring) {
    const method = stirring.getMethod();
    if (method) {
      setSelector($('.workup_stirring_method_type', node), method.getType());
      $('.workup_stirring_method_details', node).text(method.getDetails());
    }
    const rate = stirring.getRate();
    if (rate) {
      setSelector($('.workup_stirring_rate_type', node), rate.getType());
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
  setOptionalBool(
      $('.workup_automated', node),
      workup.hasIsAutomated() ? workup.getIsAutomated() : null);
};

ord.workups.loadMeasurement = function(workupNode, measurement) {
  const node = ord.workups.addMeasurement(workupNode);
  setSelector(
      $('.workup_temperature_measurement_type', node), measurement.getType());
  $('.workup_temperature_measurement_details', node)
      .text(measurement.getDetails());
  const time = measurement.getTime();
  if (time) {
    writeMetric('.workup_temperature_measurement_time', time, node);
  }
  const temperature = measurement.getTemperature();
  if (temperature) {
    writeMetric(
        '.workup_temperature_measurement_temperature', temperature, node);
  }
};

ord.workups.unload = function() {
  const workups = [];
  $('.workup').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const workup = ord.workups.unloadWorkup(node);
      if (!isEmptyMessage(workup)) {
        workups.push(workup);
      }
    }
  });
  return workups;
};

ord.workups.unloadWorkup = function(node) {
  const workup = new proto.ord.ReactionWorkup();

  workup.setType(getSelector($('.workup_type', node)));

  workup.setDetails($('.workup_details', node).text());

  const duration = readMetric('.workup_duration', new proto.ord.Time(), node);
  if (!isEmptyMessage(duration)) {
    workup.setDuration(duration);
  }

  const input = ord.inputs.unloadInputUnnamed(node);
  if (!isEmptyMessage(input)) {
    workup.setInput(input);
  }

  const control = new proto.ord.TemperatureConditions.TemperatureControl();
  control.setType(getSelector($('.workup_temperature_control_type', node)));
  control.setDetails($('.workup_temperature_details', node).text());

  const temperature = new proto.ord.TemperatureConditions();
  if (!isEmptyMessage(control)) {
    temperature.setControl(control);
  }

  const setpoint = readMetric(
      '.workup_temperature_setpoint', new proto.ord.Temperature(), node);
  if (!isEmptyMessage(setpoint)) {
    temperature.setSetpoint(setpoint);
  }

  const measurements = [];
  const measurementNodes = $('.workup_temperature_measurement', node);
  measurementNodes.each(function(index, measurementNode) {
    measurementNode = $(measurementNode);
    if (!measurementNode.attr('id')) {
      // Not a template.
      const measurement = ord.workups.unloadMeasurement(measurementNode);
      if (!isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  temperature.setMeasurementsList(measurements);
  if (!isEmptyMessage(temperature)) {
    workup.setTemperature(temperature);
  }

  workup.setKeepPhase($('.workup_keep_phase', node).text());

  const stirring = new proto.ord.StirringConditions();

  const method = new proto.ord.StirringConditions.StirringMethod();
  method.setType(getSelector($('.workup_stirring_method_type', node)));
  method.setDetails($('.workup_stirring_method_details').text());
  if (!isEmptyMessage(method)) {
    stirring.setMethod(method);
  }

  const rate = new proto.ord.StirringConditions.StirringRate();
  rate.setType(getSelector($('.workup_stirring_rate_type', node)));
  rate.setDetails($('.workup_stirring_rate_details').text());
  const rpm = parseFloat($('.workup_stirring_rate_rpm', node).text());
  if (!isNaN(rpm)) {
    rate.setRpm(rpm);
  }
  if (!isEmptyMessage(rate)) {
    stirring.setRate(rate);
  }

  if (!isEmptyMessage(stirring)) {
    workup.setStirring(stirring);
  }

  const targetPh = parseFloat($('.workup_target_ph', node).text());
  if (!isNaN(targetPh)) {
    workup.setTargetPh(targetPh);
  }
  workup.setIsAutomated(getOptionalBool($('.workup_automated', node)));
  return workup;
};

ord.workups.unloadMeasurement = function(node) {
  const measurement = new proto.ord.TemperatureConditions.Measurement();
  measurement.setType(
      getSelector($('.workup_temperature_measurement_type', node)));
  measurement.setDetails(
      $('.workup_temperature_measurement_details', node).text());
  const time = readMetric(
      '.workup_temperature_measurement_time', new proto.ord.Time(), node);
  if (!isEmptyMessage(time)) {
    measurement.setTime(time);
  }
  const temperature = readMetric(
      '.workup_temperature_measurement_temperature',
      new proto.ord.Temperature(), node);
  if (!isEmptyMessage(temperature)) {
    measurement.setTemperature(temperature);
  }
  return measurement;
};

ord.workups.add = function() {
  const workupNode = addSlowly('#workup_template', '#workups');
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
  addChangeHandler(workupNode, () => {
    ord.workups.validateWorkup(workupNode);
  });

  return workupNode;
};

ord.workups.addMeasurement = function(node) {
  return addSlowly(
      '#workup_temperature_measurement_template',
      $('.workup_temperature_measurements', node));
};

ord.workups.validateWorkup = function(node, validateNode) {
  const workup = ord.workups.unloadWorkup(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(workup, 'ReactionWorkup', validateNode);
};