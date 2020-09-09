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

goog.module('ord.inputs');
goog.module.declareLegacyNamespace();
exports = {
  load,
  loadInputUnnamed,
  unload,
  unloadInputUnnamed,
  add,
  validateInput
};

goog.require('ord.compounds');
goog.require('ord.crudes');
goog.require('proto.ord.FlowRate');
goog.require('proto.ord.ReactionInput');
goog.require('proto.ord.Time');

/**
 * Populates the reaction input sections in the form.
 * @param {!jspb.Map<string, !proto.ord.ReactionInput>} inputs
 */
function load(inputs) {
  const names = inputs.stringKeys_();
  names.forEach(function(name) {
    const input = inputs.get(name);
    loadInput('#inputs', name, input);
  });
}

/**
 * Populates a single reaction input section in the form.
 * @param {!Node} root Root node for the reaction input.
 * @param {string} name The name of this input.
 * @param {!proto.ord.ReactionInput} input
 */
function loadInput(root, name, input) {
  const node = add(root);
  loadInputUnnamed(node, input);
  $('.input_name', node).text(name);
}

/**
 * Populates a single reaction input section in the form without assigning it a
 * name.
 * @param {!Node} node Root node for the reaction input.
 * @param {!proto.ord.ReactionInput} input
 * @return {!Node} The original root node.
 */
function loadInputUnnamed(node, input) {
  const compounds = input.getComponentsList();
  ord.compounds.load(node, compounds);

  const crudes = input.getCrudeComponentsList();
  ord.crudes.load(node, crudes);

  const additionOrder = input.getAdditionOrder();
  if (additionOrder != 0) {
    $('.input_addition_order', node).text(additionOrder);
  }

  const additionTime = input.getAdditionTime();
  if (additionTime) {
    ord.reaction.writeMetric('.input_addition_time', additionTime, node);
  }
  const additionSpeed = input.getAdditionSpeed();
  if (additionSpeed) {
    ord.reaction.setSelector(
        $('.input_addition_speed_type', node), additionSpeed.getType());
    $('.input_addition_speed_details', node).text(additionSpeed.getDetails());
  }
  const additionDevice = input.getAdditionDevice();
  if (additionDevice) {
    ord.reaction.setSelector(
        $('.input_addition_device_type', node), additionDevice.getType());
    $('.input_addition_device_details', node).text(additionDevice.getDetails());
  }
  const duration = input.getAdditionDuration();
  if (duration) {
    ord.reaction.writeMetric('.input_addition_duration', duration, node);
  }
  const temperature = input.getAdditionTemperature();
  if (temperature) {
    ord.reaction.writeMetric('.input_addition_temperature', temperature, node);
  }
  const flowRate = input.getFlowRate();
  if (flowRate) {
    ord.reaction.writeMetric('.input_flow_rate', flowRate, node);
  }
  return node;
}

/**
 * Fetches the reaction inputs defined in the form and adds them to `inputs`.
 * @param {!jspb.Map<string, !proto.ord.ReactionInput>} inputs
 */
function unload(inputs) {
  $('#inputs > div.input').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      unloadInput(inputs, node);
    }
  });
}

/**
 * Fetches a single reaction input defined in the form and adds it to `inputs`.
 * @param {!jspb.Map<string, !proto.ord.ReactionInput>} inputs
 * @param {!Node} node Root node for the reaction input.
 */
function unloadInput(inputs, node) {
  const name = $('.input_name', node).text();
  const input = unloadInputUnnamed(node);
  if (!ord.reaction.isEmptyMessage(input) ||
      !ord.reaction.isEmptyMessage(name)) {
    inputs.set(name, input);
  }
}

/**
 * Fetches a single reaction input defined in the form.
 * @param {!Node} node Root node for the reaction input.
 * @return {!proto.ord.ReactionInput}
 */
function unloadInputUnnamed(node) {
  const input = new proto.ord.ReactionInput();

  const compounds = ord.compounds.unload(node);
  if (!ord.reaction.isEmptyMessage(compounds)) {
    input.setComponentsList(compounds);
  }

  const crudes = ord.crudes.unload(node);
  if (!ord.reaction.isEmptyMessage(crudes)) {
    input.setCrudeComponentsList(crudes);
  }

  const additionOrder = parseInt($('.input_addition_order', node).text());
  if (!isNaN(additionOrder)) {
    input.setAdditionOrder(additionOrder);
  }
  const additionTime = ord.reaction.readMetric(
      '.input_addition_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(additionTime)) {
    input.setAdditionTime(additionTime);
  }

  const additionSpeed = new proto.ord.ReactionInput.AdditionSpeed();
  additionSpeed.setType(
      ord.reaction.getSelector($('.input_addition_speed_type', node)));
  additionSpeed.setDetails($('.input_addition_speed_details', node).text());
  if (!ord.reaction.isEmptyMessage(additionSpeed)) {
    input.setAdditionSpeed(additionSpeed);
  }

  const additionDevice = new proto.ord.ReactionInput.AdditionDevice();
  additionDevice.setType(
      ord.reaction.getSelector($('.input_addition_device_type', node)));
  additionDevice.setDetails($('.input_addition_device_details', node).text());
  if (!ord.reaction.isEmptyMessage(additionDevice)) {
    input.setAdditionDevice(additionDevice);
  }

  const additionDuration = ord.reaction.readMetric(
      '.input_addition_duration', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(additionDuration)) {
    input.setAdditionDuration(additionDuration);
  }

  const additionTemperature = ord.reaction.readMetric(
      '.input_addition_temperature', new proto.ord.Temperature(), node);
  if (!ord.reaction.isEmptyMessage(additionTemperature)) {
    input.setAdditionTemperature(additionTemperature);
  }

  const flowRate = ord.reaction.readMetric(
      '.input_flow_rate', new proto.ord.FlowRate(), node);
  if (!ord.reaction.isEmptyMessage(flowRate)) {
    input.setFlowRate(flowRate);
  }

  return input;
}

/**
 * Adds a reaction input section to the form.
 * @param {!Node} root Parent node for reaction inputs.
 * @return {!Node} The newly added parent node for the reaction input.
 */
function add(root) {
  const node = ord.reaction.addSlowly('#input_template', root);
  // Add live validation handling.
  ord.reaction.addChangeHandler(node, () => {
    validateInput(node);
  });
  return node;
}

/**
 * Validates a reaction input defined in the form.
 * @param {!Node} node Root node for the reaction input.
 * @param {!Node} validateNode Target node for validation results.
 */
function validateInput(node, validateNode) {
  const input = unloadInputUnnamed(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(input, 'ReactionInput', validateNode);
}
