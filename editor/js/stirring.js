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

goog.module('ord.stirring');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  validateStirring
};

goog.require('proto.ord.StirringConditions');

/**
 * Populates the stirring conditions section in the form.
 * @param {!proto.ord.StirringConditions} stirring
 */
function load(stirring) {
  const method = stirring.getMethod();
  if (method) {
    ord.reaction.setSelector($('#stirring_method_type'), method.getType());
    $('#stirring_method_details').text(method.getDetails());
  }
  const rate = stirring.getRate();
  if (rate) {
    ord.reaction.setSelector($('#stirring_rate_type'), rate.getType());
    $('#stirring_rate_details').text(rate.getDetails());
    const rpm = rate.getRpm();
    if (rpm != 0) {
      $('#stirring_rpm').text(rpm);
    }
  }
}

/**
 * Fetches the stirring conditions from the form.
 * @return {!proto.ord.StirringConditions}
 */
function unload() {
  const stirring = new proto.ord.StirringConditions();

  const method = new proto.ord.StirringConditions.StirringMethod();
  method.setType(ord.reaction.getSelector($('#stirring_method_type')));
  method.setDetails($('#stirring_method_details').text());
  if (!ord.reaction.isEmptyMessage(method)) {
    stirring.setMethod(method);
  }

  const rate = new proto.ord.StirringConditions.StirringRate();
  rate.setType(ord.reaction.getSelector($('#stirring_rate_type')));
  rate.setDetails($('#stirring_rate_details').text());
  const rpm = parseFloat($('#stirring_rpm').text());
  if (!isNaN(rpm)) {
    rate.setRpm(rpm);
  }
  if (!ord.reaction.isEmptyMessage(rate)) {
    stirring.setRate(rate);
  }
  return stirring;
}

/**
 * Validates the stirring conditions defined in the form.
 * @param {!Node} node The div containing the stirring conditions.
 * @param {!Node} validateNode The target div for validation results.
 */
function validateStirring(node, validateNode) {
  const stirring = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(stirring, 'StirringConditions', validateNode);
}
