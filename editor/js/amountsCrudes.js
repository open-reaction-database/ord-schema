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

goog.module('ord.amountsCrudes');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload
};

goog.require('proto.ord.Mass');
goog.require('proto.ord.Volume');

/**
 * Populates the form's fields describing the amount of a crude compound.
 * @param {!Node} node The div corresponding to the crude compound whose amount
 *     fields on the form should be updated.
 * @param {?proto.ord.Mass} mass
 * @param {?proto.ord.Volume} volume
 */
function load(node, mass, volume) {
  const amount = $('.amount', node);
  $('.crude_amount_units_mass', node).hide();
  $('.crude_amount_units_volume', node).hide();
  if (mass) {
    $('input[value=\'mass\']', amount).prop('checked', true);
    if (mass.hasValue()) {
      $('.crude_amount_value', node).text(mass.getValue());
    }
    if (mass.hasPrecision()) {
      $('.crude_amount_precision', node).text(mass.getPrecision());
    }
    $('.crude_amount_units_mass', node).show();
    ord.reaction.setSelector(
        $('.crude_amount_units_mass', amount), mass.getUnits());
  }
  if (volume) {
    $('input[value=\'volume\']', amount).prop('checked', true);
    if (volume.hasValue()) {
      $('.crude_amount_value', node).text(volume.getValue());
    }
    if (volume.hasPrecision()) {
      $('.crude_amount_precision', node).text(volume.getPrecision());
    }
    $('.crude_amount_units_volume', node).show();
    ord.reaction.setSelector(
        $('.crude_amount_units_volume', amount), volume.getUnits());
  }
}

/**
 * Sets the amount fields of a crude component message according to the form.
 * @param {!Node} node The div corresponding to the crude component whose amount
 *     fields should be read from the form.
 * @param {!proto.ord.CrudeComponent} crude
 */
function unload(node, crude) {
  const mass = unloadMass(node);
  const volume = unloadVolume(node);
  if (mass) {
    if (!ord.reaction.isEmptyMessage(mass)) {
      crude.setMass(mass);
    }
  }
  if (volume) {
    if (!ord.reaction.isEmptyMessage(volume)) {
      crude.setVolume(volume);
    }
  }
}

/**
 * Reads and returns a mass amount as defined in the form for a crude
 * component.
 * @param {!Node} node The div corresponding to the crude component whose mass
 *     fields should be read from the form.
 * @return {?proto.ord.Mass} mass
 */
function unloadMass(node) {
  if (!$('.crude_amount_mass', node).is(':checked')) {
    return null;
  }
  const mass = new proto.ord.Mass();
  const value = parseFloat($('.crude_amount_value', node).text());
  if (!isNaN(value)) {
    mass.setValue(value);
  }
  const units = ord.reaction.getSelector($('.crude_amount_units_mass', node));
  mass.setUnits(units);
  const precision = parseFloat($('.crude_amount_precision', node).text());
  if (!isNaN(precision)) {
    mass.setPrecision(precision);
  }
  return mass;
}

/**
 * Reads and returns a volumetric amount as defined in the form for a crude
 * component.
 * @param {!Node} node The div corresponding to the crude component whose
 *     volume fields should be read from the form.
 * @return {?proto.ord.Volume} volume
 */
function unloadVolume(node) {
  if (!$('.crude_amount_volume', node).is(':checked')) {
    return null;
  }
  const volume = new proto.ord.Volume();
  const value = parseFloat($('.crude_amount_value', node).text());
  if (!isNaN(value)) {
    volume.setValue(value);
  }
  const units = ord.reaction.getSelector($('.crude_amount_units_volume', node));
  volume.setUnits(units);
  const precision = parseFloat($('.crude_amount_precision', node).text());
  if (!isNaN(precision)) {
    volume.setPrecision(precision);
  }
  return volume;
}
