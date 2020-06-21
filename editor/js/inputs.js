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
goog.require('proto.ord.FlowRate');
goog.require('proto.ord.ReactionInput');
goog.require('proto.ord.Time');

ord.inputs.load = function (inputs) {
  const names = inputs.stringKeys_();
  names.forEach(function (name) {
    const input = inputs.get(name);
    ord.inputs.loadInput('#inputs', name, input);
  });
};

ord.inputs.loadInput = function (root, name, input) {
  const node = ord.inputs.add(root); 
  ord.inputs.loadInputUnnamed(node, input);
  $('.input_name', node).text(name);
};

ord.inputs.loadInputUnnamed = function (node, input) {
  const compounds = input.getComponentsList();
  ord.compounds.load(node, compounds);

  const additionOrder = input.getAdditionOrder();
  $('.input_addition_order', node).text(additionOrder);

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

ord.inputs.unload = function (inputs) {
  $('.input').each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      ord.inputs.unloadInput(inputs, node);
    }
  });
};

ord.inputs.unloadInput = function (inputs, node) {
  const name = $('.input_name', node).text();
  const input = ord.inputs.unloadInputUnnamed(node);
  inputs.set(name, input);
};

ord.inputs.unloadInputUnnamed = function (node) {
  const input = new proto.ord.ReactionInput();

  const compounds = ord.compounds.unload(node);
  input.setComponentsList(compounds);

  const additionOrder = parseInt($('.input_addition_order', node).text());
  if (!isNaN(additionOrder)) {
    input.setAdditionOrder(additionOrder);
  }
  const additionTime =
      readMetric('.input_addition_time', new proto.ord.Time(), node);
  input.setAdditionTime(additionTime);

  const additionSpeed = new proto.ord.ReactionInput.AdditionSpeed();
  additionSpeed.setType(getSelector($('.input_addition_speed_type', node)));
  additionSpeed.setDetails($('.input_addition_speed_details', node).text());
  input.setAdditionSpeed(additionSpeed);

  const additionDevice = new proto.ord.ReactionInput.AdditionDevice();
  additionDevice.setType(getSelector($('.input_addition_device_type', node)));
  additionDevice.setDetails($('.input_addition_device_details', node).text());
  input.setAdditionDevice(additionDevice);

  const additionDuration =
      readMetric('.input_addition_duration', new proto.ord.Time(), node);
  input.setAdditionDuration(additionDuration);

  const flowRate =
      readMetric('.input_flow_rate', new proto.ord.FlowRate(), node);
  input.setFlowRate(flowRate);

  return input;
};

ord.inputs.add = function (root) {
  const node = addSlowly('#input_template', root);
  // node.css('backgroundColor', 'green')

  // statusNode = $('.validate_status', node); 
  // statusNode.css('backgroundColor', 'red')
  // handler = function () {ord.inputs.validateInput(node, statusNode)};
  // ^ unfortunately, the statusNode is not queried properly in the first line
  // (bizarrely, however, the correct statusNode does turn red...
  // perhaps the first line occcurs, then the addSlowly occurs, then the second line occurs?

  
  $('.validate_status', node).css('backgroundColor', 'yellow')
  statusNode = $('.validate_status', node); 
    // handler = function () {ord.inputs.validateInput(node, statusNode)};
    handler = function () {ord.inputs.validateInput(node, $('.validate_status', node))};


  // handler = function () {ord.inputs.validateInput(node, $('.validate_status', node))};
  // ^this works !! ig statusNode

  addChangeHandler(node, handler);
  return node;
};

ord.inputs.validateInput = function(node, statusNode) {
  console.log("attempting to validate input");
  console.log("node")
  console.log(node)
  node.css('backgroundColor', 'green')

  statusNode = $(statusNode)
  console.log("statusNode");
  console.log(statusNode);
  statusNode.css('backgroundColor', 'blue')
  // when triggered by the handler, nothing turns blue

  // // what if we try to find it on the fly?
  // statusNode = $('.validate_status', node)
  // // this works! it's not ideal tho, bc if there are multiple validate status, we'd modify all;
  // // using .first() isn't great either bc the statusNode we're trying to find might not always be first;
  // // could use more descriptive class name (eg input_validate_status) but thats less copy-pastable
  // console.log("statusNode");
  // console.log(statusNode);
  // statusNode.css('backgroundColor', 'blue')

  console.log('attempt to find messageNode:')
  console.log($('.validate_message', statusNode))
  const input = ord.inputs.unloadInputUnnamed(node);
  validate(input, "ReactionInput", statusNode);
};