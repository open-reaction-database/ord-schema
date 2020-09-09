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

goog.module('ord.flows');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  validateFlow
};

goog.require('proto.ord.FlowConditions');
goog.require('proto.ord.FlowConditions.Tubing');

/**
 * Populates the flow conditions section in the form.
 * @param {!proto.ord.FlowConditions} flow
 */
function load(flow) {
  const type = flow.getFlowType();
  if (type) {
    ord.reaction.setSelector($('#flow_type'), type.getType());
    $('#flow_details').text(type.getDetails());
  }
  $('#flow_pump').text(flow.getPumpType());

  const tubing = flow.getTubing();
  ord.reaction.setSelector($('#flow_tubing_type'), tubing.getType());
  $('#flow_tubing_details').text(tubing.getDetails());
  ord.reaction.writeMetric('#flow_tubing', tubing.getDiameter());
}

/**
 * Fetches the flow conditions defined in the form.
 * @return {!proto.ord.FlowConditions}
 */
function unload() {
  const flow = new proto.ord.FlowConditions();

  const type = new proto.ord.FlowConditions.FlowType();
  type.setType(ord.reaction.getSelector('#flow_type'));
  type.setDetails($('#flow_details').text());
  if (!ord.reaction.isEmptyMessage(type)) {
    flow.setFlowType(type);
  }

  flow.setPumpType($('#flow_pump').text());

  const tubing = new proto.ord.FlowConditions.Tubing();
  tubing.setType(ord.reaction.getSelector('#flow_tubing_type'));
  tubing.setDetails($('#flow_tubing_details').text());
  const diameter =
      ord.reaction.readMetric('#flow_tubing', new proto.ord.Length());
  if (!ord.reaction.isEmptyMessage(diameter)) {
    tubing.setDiameter(diameter);
  }

  if (!ord.reaction.isEmptyMessage(tubing)) {
    flow.setTubing(tubing);
  }
  return flow;
}

/**
 * Validates the flow conditions defined in the form.
 * @param {!Node} node Root node for the flow conditions.
 * @param {!Node} validateNode Target node for validation results.
 */
function validateFlow(node, validateNode) {
  const flow = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(flow, 'FlowConditions', validateNode);
}
