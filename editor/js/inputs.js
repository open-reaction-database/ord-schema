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

goog.provide('ord.inputs');

goog.require('ord.compounds');
goog.require('ord.crudes');
goog.require('proto.ord.FlowRate');
goog.require('proto.ord.ReactionInput');
goog.require('proto.ord.Time');

ord.inputs.load = function(inputs) {
  const names = inputs.stringKeys_();
  names.forEach(function(name) {
    const input = inputs.get(name);
    ord.inputs.loadInput('#inputs', name, input);
  });
};

ord.inputs.loadInput = function(root, name, input) {
  const node = ord.inputs.add(root);
  ord.inputs.loadInputUnnamed(node, input);
  $('.input_name', node).text(name);
};

ord.inputs.loadInputUnnamed = function(node, input) {
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
    writeMetric('.input_addition_time', additionTime, node);
  }
  const additionSpeed = input.getAdditionSpeed();
  if (additionSpeed) {
    setSelector($('.input_addition_speed_type', node), additionSpeed.getType());
    $('.input_addition_speed_details', node).text(additionSpeed.getDetails());
  }
  const additionDevice = input.getAdditionDevice();
  if (additionDevice) {
    setSelector(
        $('.input_addition_device_type', node), additionDevice.getType());
    $('.input_addition_device_details', node).text(additionDevice.getDetails());
  }
  const duration = input.getAdditionDuration();
  if (duration) {
    writeMetric('.input_addition_duration', duration, node);
  }
  const flowRate = input.getFlowRate();
  if (flowRate) {
    writeMetric('.input_flow_rate', flowRate, node);
  }
  return node;
};

ord.inputs.unload = function(inputs) {
  $('#inputs > div.input').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      ord.inputs.unloadInput(inputs, node);
    }
  });
};

ord.inputs.unloadInput = function(inputs, node) {
  const name = $('.input_name', node).text();
  const input = ord.inputs.unloadInputUnnamed(node);
  if (!isEmptyMessage(input) || !isEmptyMessage(name)) {
    inputs.set(name, input);
  }
};

ord.inputs.unloadInputUnnamed = function(node) {
  const input = new proto.ord.ReactionInput();

  const compounds = ord.compounds.unload(node);
  if (!isEmptyMessage(compounds)) {
    input.setComponentsList(compounds);
  }

  const crudes = ord.crudes.unload(node);
  if (!isEmptyMessage(crudes)) {
    input.setCrudeComponentsList(crudes);
  }

  const additionOrder = parseInt($('.input_addition_order', node).text());
  if (!isNaN(additionOrder)) {
    input.setAdditionOrder(additionOrder);
  }
  const additionTime =
      readMetric('.input_addition_time', new proto.ord.Time(), node);
  if (!isEmptyMessage(additionTime)) {
    input.setAdditionTime(additionTime);
  }

  const additionSpeed = new proto.ord.ReactionInput.AdditionSpeed();
  additionSpeed.setType(getSelector($('.input_addition_speed_type', node)));
  additionSpeed.setDetails($('.input_addition_speed_details', node).text());
  if (!isEmptyMessage(additionSpeed)) {
    input.setAdditionSpeed(additionSpeed);
  }

  const additionDevice = new proto.ord.ReactionInput.AdditionDevice();
  additionDevice.setType(getSelector($('.input_addition_device_type', node)));
  additionDevice.setDetails($('.input_addition_device_details', node).text());
  if (!isEmptyMessage(additionDevice)) {
    input.setAdditionDevice(additionDevice);
  }

  const additionDuration =
      readMetric('.input_addition_duration', new proto.ord.Time(), node);
  if (!isEmptyMessage(additionDuration)) {
    input.setAdditionDuration(additionDuration);
  }

  const flowRate =
      readMetric('.input_flow_rate', new proto.ord.FlowRate(), node);
  if (!isEmptyMessage(flowRate)) {
    input.setFlowRate(flowRate);
  }

  return input;
};

ord.inputs.add = function(root) {
  const node = addSlowly('#input_template', root);
  // Add live validation handling.
  addChangeHandler(node, () => {
    ord.inputs.validateInput(node);
  });
  return node;
};

ord.inputs.validateInput = function(node, validateNode) {
  const input = ord.inputs.unloadInputUnnamed(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(input, 'ReactionInput', validateNode);
};