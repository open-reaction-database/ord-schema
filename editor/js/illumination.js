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

goog.provide('ord.illumination');

goog.require('proto.ord.IlluminationConditions');
goog.require('proto.ord.Length');
goog.require('proto.ord.Wavelength');

ord.illumination.load = function (illumination) {
  const type = illumination.getType();
  if (type) {
    setSelector($('#illumination_type'), type.getType());
    $('#illumination_details').text(type.getDetails());
  }
  const wavelength = illumination.getPeakWavelength();
  writeMetric('#illumination_wavelength', wavelength);
  $('#illumination_color').text(illumination.getColor());
  const distance = illumination.getDistanceToVessel();
  writeMetric('#illumination_distance', distance);
};

ord.illumination.unload = function () {
  const illumination = new proto.ord.IlluminationConditions();

  const type = new proto.ord.IlluminationConditions.IlluminationType();
  type.setType(getSelector($('#illumination_type')));
  type.setDetails($('#illumination_details').text());
  illumination.setType(type);

  const wavelength =
      readMetric('#illumination_wavelength', new proto.ord.Wavelength());
  illumination.setPeakWavelength(wavelength);
  illumination.setColor($('#illumination_color').text());
  const distance =
      readMetric('#illumination_distance', new proto.ord.Length());
  illumination.setDistanceToVessel(distance);
  return illumination;
};

ord.illumination.validateIllumination = function(node, validateNode) {
  const illumination = ord.illumination.unload();
  if (typeof validateNode === 'undefined') {
    validateNode = $('.validate', node).first();
  }
  validate(illumination, "IlluminationConditions", validateNode);
};