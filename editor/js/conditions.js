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

goog.module('ord.conditions');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  validateConditions
};

goog.require('ord.electro');
goog.require('ord.flows');
goog.require('ord.illumination');
goog.require('ord.pressure');
goog.require('ord.stirring');
goog.require('ord.temperature');
goog.require('proto.ord.ReactionConditions');

function load(conditions) {
  const temperature = conditions.getTemperature();
  if (temperature) {
    ord.temperature.load(temperature);
  }
  const pressure = conditions.getPressure();
  if (pressure) {
    ord.pressure.load(pressure);
  }
  const stirring = conditions.getStirring();
  if (stirring) {
    ord.stirring.load(stirring);
  }
  const illumination = conditions.getIllumination();
  if (illumination) {
    ord.illumination.load(illumination);
  }
  const electro = conditions.getElectrochemistry();
  if (electro) {
    ord.electro.load(electro);
  }
  const flow = conditions.getFlow();
  if (flow) {
    ord.flows.load(flow);
  }
  const reflux = conditions.hasReflux() ? conditions.getReflux() : null;
  ord.reaction.setOptionalBool($('#condition_reflux'), reflux);
  if (conditions.hasPh()) {
    $('#condition_ph').text(conditions.getPh());
  }
  const dynamic = conditions.hasConditionsAreDynamic() ?
      conditions.getConditionsAreDynamic() :
      null;
  ord.reaction.setOptionalBool($('#condition_dynamic'), dynamic);
  $('#condition_details').text(conditions.getDetails());
}

function unload() {
  const conditions = new proto.ord.ReactionConditions();
  const temperature = ord.temperature.unload();
  if (!ord.reaction.isEmptyMessage(temperature)) {
    conditions.setTemperature(temperature);
  }
  const pressure = ord.pressure.unload();
  if (!ord.reaction.isEmptyMessage(pressure)) {
    conditions.setPressure(pressure);
  }
  const stirring = ord.stirring.unload();
  if (!ord.reaction.isEmptyMessage(stirring)) {
    conditions.setStirring(stirring);
  }
  const illumination = ord.illumination.unload();
  if (!ord.reaction.isEmptyMessage(illumination)) {
    conditions.setIllumination(illumination);
  }
  const electro = ord.electro.unload();
  if (!ord.reaction.isEmptyMessage(electro)) {
    conditions.setElectrochemistry(electro);
  }
  const flow = ord.flows.unload();
  if (!ord.reaction.isEmptyMessage(flow)) {
    conditions.setFlow(flow);
  }

  const reflux = ord.reaction.getOptionalBool($('#condition_reflux'));
  conditions.setReflux(reflux);
  const ph = parseFloat($('#condition_ph').text());
  if (!isNaN(ph)) {
    conditions.setPh(ph);
  }
  const dynamic = ord.reaction.getOptionalBool($('#condition_dynamic'));
  conditions.setConditionsAreDynamic(dynamic);
  const details = $('#condition_details').text();
  conditions.setDetails(details);
  return conditions;
}

function validateConditions(node, validateNode) {
  const condition = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(condition, 'ReactionConditions', validateNode);
}