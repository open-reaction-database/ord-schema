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

goog.module('ord.amounts');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  unloadVolume
};

goog.require('proto.ord.Mass');
goog.require('proto.ord.Moles');
goog.require('proto.ord.Volume');

/**
 * Populates the form's fields describing the amount of a compound.
 * @param {Node} node The div corresponding to the compound whose amount fields
 *     on the form should be updated.
 * @param {?proto.ord.Mass} mass
 * @param {?proto.ord.Moles} moles
 * @param {?proto.ord.Volume} volume
 */
function load(node, mass, moles, volume) {
  const amount = $('.amount', node);
  $('.component_amount_units_mass', node).hide();
  $('.component_amount_units_moles', node).hide();
  $('.component_amount_units_volume', node).hide();
  $('.includes_solutes', node).hide();
  if (mass) {
    $('input[value=\'mass\']', amount).prop('checked', true);
    if (mass.hasValue()) {
      $('.component_amount_value', node).text(mass.getValue());
    }
    if (mass.hasPrecision()) {
      $('.component_amount_precision', node).text(mass.getPrecision());
    }
    $('.component_amount_units_mass', node).show();
    ord.reaction.setSelector(
        $('.component_amount_units_mass', amount), mass.getUnits());
  }
  if (moles) {
    $('input[value=\'moles\']', amount).prop('checked', true);
    if (moles.hasValue()) {
      $('.component_amount_value', node).text(moles.getValue());
    }
    if (moles.hasPrecision()) {
      $('.component_amount_precision', node).text(moles.getPrecision());
    }
    $('.component_amount_units_moles', node).show();
    ord.reaction.setSelector(
        $('.component_amount_units_moles', amount), moles.getUnits());
  }
  if (volume) {
    $('input[value=\'volume\']', amount).prop('checked', true);
    if (volume.hasValue()) {
      $('.component_amount_value', node).text(volume.getValue());
    }
    if (volume.hasPrecision()) {
      $('.component_amount_precision', node).text(volume.getPrecision());
    }
    $('.component_amount_units_volume', node).show();
    $('.includes_solutes', node).show().css('display', 'inline-block');
    ord.reaction.setSelector(
        $('.component_amount_units_volume', amount), volume.getUnits());
  }
}

/**
 * Sets the amount fields of a compound message according to the form.
 * @param {Node} node The div corresponding to the compound whose amount fields
 *     should be read from the form.
 * @param {proto.ord.Compound} compound
 */
function unload(node, compound) {
  const mass = unloadMass(node);
  const moles = unloadMoles(node);
  const volume = unloadVolume(node);
  if (mass) {
    if (!ord.reaction.isEmptyMessage(mass)) {
      compound.setMass(mass);
    }
  }
  if (moles) {
    if (!ord.reaction.isEmptyMessage(moles)) {
      compound.setMoles(moles);
    }
  }
  if (volume) {
    if (!ord.reaction.isEmptyMessage(volume)) {
      compound.setVolume(volume);
    }
  }
}

/**
 * Reads and returns a mass amount as defined in the form.
 * @param {Node} node The div corresponding to the compound whose mass fields
 *     should be read from the form.
 * @return {?proto.ord.Mass} mass
 */
function unloadMass(node) {
  if (!$('.component_amount_mass', node).is(':checked')) {
    return null;
  }
  const mass = new proto.ord.Mass();
  const value = parseFloat($('.component_amount_value', node).text());
  if (!isNaN(value)) {
    mass.setValue(value);
  }
  const units =
      ord.reaction.getSelector($('.component_amount_units_mass', node));
  mass.setUnits(units);
  const precision = parseFloat($('.component_amount_precision', node).text());
  if (!isNaN(precision)) {
    mass.setPrecision(precision);
  }
  return mass;
}

/**
 * Reads and returns a molar amount as defined in the form.
 * @param {Node} node The div corresponding to the compound whose moles fields
 *     should be read from the form.
 * @return {?proto.ord.Moles} moles
 */
function unloadMoles(node) {
  if (!$('.component_amount_moles', node).is(':checked')) {
    return null;
  }
  const moles = new proto.ord.Moles();
  const value = parseFloat($('.component_amount_value', node).text());
  if (!isNaN(value)) {
    moles.setValue(value);
  }
  const units =
      ord.reaction.getSelector($('.component_amount_units_moles', node));
  moles.setUnits(units);
  const precision = parseFloat($('.component_amount_precision', node).text());
  if (!isNaN(precision)) {
    moles.setPrecision(precision);
  }
  return moles;
}

/**
 * Reads and returns a volumetric amount as defined in the form.
 * @param {Node} node The div corresponding to the compound whose volume fields
 *     should be read from the form.
 * @return {?proto.ord.Volume} volume
 */
function unloadVolume(node) {
  if (!$('.component_amount_volume', node).is(':checked')) {
    return null;
  }
  const volume = new proto.ord.Volume();
  const value = parseFloat($('.component_amount_value', node).text());
  if (!isNaN(value)) {
    volume.setValue(value);
  }
  const units =
      ord.reaction.getSelector($('.component_amount_units_volume', node));
  volume.setUnits(units);
  const precision = parseFloat($('.component_amount_precision', node).text());
  if (!isNaN(precision)) {
    volume.setPrecision(precision);
  }
  return volume;
}
