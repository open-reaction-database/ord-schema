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

goog.provide('ord.amounts');

goog.require('proto.ord.Mass');
goog.require('proto.ord.Moles');
goog.require('proto.ord.Volume');

ord.amounts.load = function(node, mass, moles, volume) {
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
};

ord.amounts.unload = function(node, compound) {
  const mass = ord.amounts.unloadMass(node);
  const moles = ord.amounts.unloadMoles(node);
  const volume = ord.amounts.unloadVolume(node);
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
};

ord.amounts.unloadMass = function(node) {
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
};

ord.amounts.unloadMoles = function(node) {
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
};

ord.amounts.unloadVolume = function(node) {
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
};
