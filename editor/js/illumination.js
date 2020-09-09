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

goog.module('ord.illumination');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  validateIllumination
};

goog.require('proto.ord.IlluminationConditions');
goog.require('proto.ord.Length');
goog.require('proto.ord.Wavelength');

/**
 * Populates the illumination conditions section in the form.
 * @param {!proto.ord.IlluminationConditions} illumination
 */
function load(illumination) {
  const type = illumination.getType();
  if (type) {
    ord.reaction.setSelector($('#illumination_type'), type.getType());
    $('#illumination_details').text(type.getDetails());
  }
  const wavelength = illumination.getPeakWavelength();
  ord.reaction.writeMetric('#illumination_wavelength', wavelength);
  $('#illumination_color').text(illumination.getColor());
  const distance = illumination.getDistanceToVessel();
  ord.reaction.writeMetric('#illumination_distance', distance);
}

/**
 * Fetches the illumination conditions defined in the form.
 * @return {!proto.ord.IlluminationConditions}
 */
function unload() {
  const illumination = new proto.ord.IlluminationConditions();

  const type = new proto.ord.IlluminationConditions.IlluminationType();
  type.setType(ord.reaction.getSelector($('#illumination_type')));
  type.setDetails($('#illumination_details').text());
  if (!ord.reaction.isEmptyMessage(type)) {
    illumination.setType(type);
  }

  const wavelength = ord.reaction.readMetric(
      '#illumination_wavelength', new proto.ord.Wavelength());
  if (!ord.reaction.isEmptyMessage(wavelength)) {
    illumination.setPeakWavelength(wavelength);
  }
  illumination.setColor($('#illumination_color').text());
  const distance =
      ord.reaction.readMetric('#illumination_distance', new proto.ord.Length());
  if (!ord.reaction.isEmptyMessage(distance)) {
    illumination.setDistanceToVessel(distance);
  }
  return illumination;
}

/**
 * Validates the illumination conditions defined in the form.
 * @param {!Node} node Root node for the illumination conditions.
 * @param {?Node} validateNode Target node for validation results.
 */
function validateIllumination(node, validateNode) {
  const illumination = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(illumination, 'IlluminationConditions', validateNode);
}
