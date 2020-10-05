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

goog.module('ord.amountsWorkups');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload
};

goog.require('proto.ord.Mass');
goog.require('proto.ord.Volume');

/**
 * Adds and populates the form's fields describing the amount of an aliquot.
 * @param {!Node} node The div corresponding to the workup step whose amount
 *     fields on the form should be updated.
 * @param {?proto.ord.Mass} mass
 * @param {?proto.ord.Volume} volume
 */
function load(node, mass, volume) {
  const amount = $('.amount', node);
  if (mass) {
    $('input[value=\'mass\']', amount).prop('checked', true);
    if (mass.hasValue()) {
      $('.workup_amount_value', node).text(mass.getValue());
    }
    if (mass.hasPrecision()) {
      $('.workup_amount_precision', node).text(mass.getPrecision());
    }
    $('.workup_amount_units_mass', node).show();
    $('.workup_amount_units_volume', node).hide();
    ord.reaction.setSelector(
        $('.workup_amount_units_mass', amount), mass.getUnits());
  }
  if (volume) {
    $('input[value=\'volume\']', amount).prop('checked', true);
    if (volume.hasValue()) {
      $('.workup_amount_value', node).text(volume.getValue());
    }
    if (volume.hasPrecision()) {
      $('.workup_amount_precision', node).text(volume.getPrecision());
    }
    $('.workup_amount_units_volume', node).show();
    $('.workup_amount_units_mass', node).hide();
    ord.reaction.setSelector(
        $('.workup_amount_units_volume', amount), volume.getUnits());
  }
}

/**
 * Sets the amount fields of a workup message according to the form.
 * @param {!Node} node The div corresponding to the workup step whose amount
 *     fields should be read from the form.
 * @param {!proto.ord.ReactionWorkup} workup
 */
function unload(node, workup) {
  const mass = unloadMass(node);
  const volume = unloadVolume(node);
  if (mass) {
    if (!ord.reaction.isEmptyMessage(mass)) {
      workup.setMass(mass);
    }
  }
  if (volume) {
    if (!ord.reaction.isEmptyMessage(volume)) {
      workup.setVolume(volume);
    }
  }
}

/**
 * Reads and returns a mass amount as defined in the form for a workup step.
 * @param {!Node} node The div corresponding to the workup step whose mass
 *     fields should be read from the form.
 * @return {?proto.ord.Mass}
 */
function unloadMass(node) {
  if (!$('.workup_amount_mass', node).is(':checked')) {
    return null;
  }
  const mass = new proto.ord.Mass();
  const value = parseFloat($('.workup_amount_value', node).text());
  if (!isNaN(value)) {
    mass.setValue(value);
  }
  const units = ord.reaction.getSelector($('.workup_amount_units_mass', node));
  mass.setUnits(units);
  const precision = parseFloat($('.workup_amount_precision', node).text());
  if (!isNaN(precision)) {
    mass.setPrecision(precision);
  }
  return mass;
}

/**
 * Reads and returns a volumetric amount as defined in the form for a workup
 * step.
 * @param {!Node} node The div corresponding to the workups tep whose
 *     volume fields should be read from the form.
 * @return {?proto.ord.Volume}
 */
function unloadVolume(node) {
  if (!$('.workup_amount_volume', node).is(':checked')) {
    return null;
  }
  const volume = new proto.ord.Volume();
  const value = parseFloat($('.workup_amount_value', node).text());
  if (!isNaN(value)) {
    volume.setValue(value);
  }
  const units =
      ord.reaction.getSelector($('.workup_amount_units_volume', node));
  volume.setUnits(units);
  const precision = parseFloat($('.workup_amount_precision', node).text());
  if (!isNaN(precision)) {
    volume.setPrecision(precision);
  }
  return volume;
}
