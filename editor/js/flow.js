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

goog.provide('ord.flows');

goog.require('proto.ord.FlowConditions');
goog.require('proto.ord.FlowConditions.Tubing');

ord.flows.load = function(flow) {
  const type = flow.getFlowType();
  if (type) {
    setSelector($('#flow_type'), type.getType());
    $('#flow_details').text(type.getDetails());
  }
  $('#flow_pump').text(flow.getPumpType());

  const tubing = flow.getTubing();
  setSelector($('#flow_tubing_type'), tubing.getType());
  $('#flow_tubing_details').text(tubing.getDetails());
  writeMetric('#flow_tubing', tubing.getDiameter());
};

ord.flows.unload = function() {
  const flow = new proto.ord.FlowConditions();

  const type = new proto.ord.FlowConditions.FlowType();
  type.setType(getSelector('#flow_type'));
  type.setDetails($('#flow_details').text());
  if (!isEmptyMessage(type)) {
    flow.setFlowType(type);
  }

  flow.setPumpType($('#flow_pump').text());

  const tubing = new proto.ord.FlowConditions.Tubing();
  tubing.setType(getSelector('#flow_tubing_type'));
  tubing.setDetails($('#flow_tubing_details').text());
  const diameter = readMetric('#flow_tubing', new proto.ord.Length());
  if (!isEmptyMessage(diameter)) {
    tubing.setDiameter(diameter);
  }

  if (!isEmptyMessage(tubing)) {
    flow.setTubing(tubing);
  }
  return flow;
};

ord.flows.validateFlow = function(node, validateNode) {
  const flow = ord.flows.unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(flow, 'FlowConditions', validateNode);
};