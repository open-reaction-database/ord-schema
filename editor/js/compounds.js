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

goog.provide('ord.compounds');

goog.require('ord.amounts');
goog.require('ord.features');
goog.require('proto.ord.Compound');
goog.require('proto.ord.CompoundIdentifier');

// Freely create radio button groups by generating new input names.
ord.compounds.radioGroupCounter = 0;

ord.compounds.load = function (node, compounds) {
  compounds.forEach(compound => ord.compounds.loadCompound(node, compound));
};

ord.compounds.loadCompound = function (root, compound) {
  const node = ord.compounds.add(root);
  ord.compounds.loadIntoCompound(node, compound);
};

ord.compounds.loadIntoCompound = function (node, compound) {
  const reactionRole = compound.getReactionRole();
  setSelector($('.component_reaction_role', node), reactionRole);

  const isLimiting = compound.hasIsLimiting() ? compound.getIsLimiting() : null;
  setOptionalBool($('.component_limiting', node), isLimiting);

  const solutes = compound.hasVolumeIncludesSolutes() ?
      compound.getVolumeIncludesSolutes() : null;
  setOptionalBool($('.component_includes_solutes', node), solutes);

  const identifiers = compound.getIdentifiersList();
  identifiers.forEach(
      identifier => ord.compounds.loadIdentifier(node, identifier));

  const mass = compound.getMass();
  const moles = compound.getMoles();
  const volume = compound.getVolume();
  ord.amounts.load(node, mass, moles, volume);

  const preparations = compound.getPreparationsList();
  preparations.forEach(preparation => {
    const preparationNode = ord.compounds.addPreparation(node);
    ord.compounds.loadPreparation(preparationNode, preparation);
  });
  const vendorSource = compound.getVendorSource();
  const vendorLot = compound.getVendorLot();
  const vendorId = compound.getVendorId();
  ord.compounds.loadVendor(node, vendorSource, vendorLot, vendorId);

  const features = compound.getFeaturesList();
  ord.features.load(node, features);
};

ord.compounds.loadIdentifier = function (compoundNode, identifier) {
  const node = ord.compounds.addIdentifier(compoundNode);
  const bytesValue = identifier.getBytesValue();
  if (bytesValue) {
    $('.component_identifier_upload', node).prop('checked', true);
    $('.component_identifier_value', node).hide();
    ord.uploads.load(node, bytesValue);
  } else {
    const value = identifier.getValue();
    $('.component_identifier_value', node).text(value);
  }
  setSelector(node, identifier.getType());
  $('.component_identifier_details', node).text(identifier.getDetails());
};

ord.compounds.loadPreparation = function (node, preparation) {
  const type = preparation.getType();
  setSelector($('.component_compound_preparation_type', node), type);
  const details = preparation.getDetails();
  $('.component_compound_preparation_details', node).text(details);
  const reaction = preparation.getReactionId();
  $('.component_compound_preparation_reaction', node).text(reaction);
};

ord.compounds.loadVendor =
    function (compoundNode, vendorSource, vendorLot, vendorId) {
  const node = $('fieldset.vendor', compoundNode);
  $('.component_vendor_source', node).text(vendorSource);
  $('.component_vendor_lot', node).text(vendorLot);
  $('.component_vendor_id', node).text(vendorId);
};

ord.compounds.unload = function (node) {
  const compounds = [];
  $('.component', node).each(function (index, compoundNode) {
    compoundNode = $(compoundNode);
    if (!compoundNode.attr('id')) {
      // Not a template.
      const compound = ord.compounds.unloadCompound(compoundNode);
      if (!isEmptyMessage(compound)) {
        compounds.push(compound);
      }
    }
  });
  return compounds;
};

ord.compounds.unloadCompound = function (node) {
  const compound = new proto.ord.Compound();

  const reactionRole = getSelector($('.component_reaction_role', node));
  compound.setReactionRole(reactionRole);

  // Only call setIsLimiting if this is a reactant Compound.
  if (getSelectorText($('.component_reaction_role', node)[0]) == 'REACTANT') {
    const isLimiting = getOptionalBool($('.component_limiting', node));
    compound.setIsLimiting(isLimiting);
  }

  // Only call setVolumeIncludesSolutes if the amount is defined as a volume.
  if (ord.amounts.unloadVolume(node)) {
    const solutes = getOptionalBool($('.component_includes_solutes', node));
    compound.setVolumeIncludesSolutes(solutes);
  }

  const identifiers = ord.compounds.unloadIdentifiers(node);
  if (!isEmptyMessage(identifiers)) {
    compound.setIdentifiersList(identifiers);
  }

  ord.amounts.unload(node, compound);

  const preparations = [];
  $('.component_preparation', node).each(function (index, preparationNode) {
    const preparation = ord.compounds.unloadPreparation(preparationNode);
    if (!isEmptyMessage(preparation)) {
      preparations.push(preparation);
    }
  });
  compound.setPreparationsList(preparations);

  ord.compounds.unloadVendor(node, compound);

  const features = ord.features.unload(node);
  compound.setFeaturesList(features);

  return compound;
};

ord.compounds.unloadIdentifiers = function (node) {
  const identifiers = [];
  $('.component_identifier', node).each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      const identifier = ord.compounds.unloadIdentifier(node);
      if (!isEmptyMessage(identifier)) {
        identifiers.push(identifier);
      }
    }
  });
  return identifiers;
};

ord.compounds.unloadIdentifier = function (node) {
  const identifier = new proto.ord.CompoundIdentifier();

  if ($('.component_identifier_upload', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!isEmptyMessage(bytesValue)) {
      identifier.setBytesValue(bytesValue);
    }
  } else {
    const value = $('.component_identifier_value', node).text();
    if (!isEmptyMessage(value)) {
      identifier.setValue(value);
    }
  }
  const type = getSelector(node);
  identifier.setType(type);
  const details = $('.component_identifier_details', node).text();
  identifier.setDetails(details);
  return identifier;
};

ord.compounds.unloadPreparation = function (node) {
  const preparation = new proto.ord.CompoundPreparation();
  const type =
      getSelector($('.component_compound_preparation_type', node));
  preparation.setType(type);
  const details = $('.component_compound_preparation_details', node).text();
  preparation.setDetails(details);
  const reaction = $('.component_compound_preparation_reaction', node).text();
  preparation.setReactionId(reaction);
  return preparation;
};

ord.compounds.unloadVendor = function (node, compound) {
  const vendorSource = $('.component_vendor_source', node).text();
  compound.setVendorSource(vendorSource);
  const vendorLot = $('.component_vendor_lot', node).text();
  compound.setVendorLot(vendorLot);
  const vendorId = $('.component_vendor_id', node).text();
  compound.setVendorId(vendorId);
};

ord.compounds.add = function (root) {
  const node = addSlowly(
      '#component_template', $('.components', root));

  // Connect reaction role selection to limiting reactant field.
  const roleSelector = $('.component_reaction_role', node);
  roleSelector.change(function() {
    $('.limiting_reactant').css('display', 'none');
    if (getSelectorText(this) == 'REACTANT') {
      $('.limiting_reactant').css('display', '');
    }
  });

  // Create an "amount" radio button group and connect it to the unit selectors.
  const amountButtons = $('.amount input', node);
  amountButtons.attr('name', 'compounds_' + ord.compounds.radioGroupCounter++);
  amountButtons.change(function () {
    $('.amount .selector', node).hide();
    if (this.value == 'mass') {
      $('.component_amount_units_mass', root).show();
      $('.includes_solutes', root).hide();
    }
    if (this.value == 'moles') {
      $('.component_amount_units_moles', root).show();
      $('.includes_solutes', root).hide();
    }
    if (this.value == 'volume') {
      $('.component_amount_units_volume', root).show();
      $('.includes_solutes', root).show();
    }
  });

  handler = function () {ord.compounds.validateCompound(node)};
  addChangeHandler(node, handler);

  return node;
};

ord.compounds.addIdentifier = function (node) {
  const identifierNode =
      addSlowly('#component_identifier_template', $('.identifiers', node));

  const uploadButton = $('.component_identifier_upload', identifierNode);
  uploadButton.change(function () {
    if ($(this).is(':checked')) {
      $('.uploader', identifierNode).show();
      $('.component_identifier_value', identifierNode).hide();
    } else {
      $('.uploader', identifierNode).hide();
      $('.component_identifier_value', identifierNode).show();
    }
  });
  ord.uploads.initialize(identifierNode);
  return identifierNode;
};

ord.compounds.addPreparation = function (node) {
  return addSlowly('#component_preparation_template', $('.preparations', node));
};

ord.compounds.validateCompound = function(node, validateNode) {
  const compound = ord.compounds.unloadCompound(node);
  if (typeof validateNode === 'undefined') {
    validateNode = $('.validate', node).first();
  }
  validate(compound, "Compound", validateNode);
};