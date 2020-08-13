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

ord.compounds.load = function(node, compounds) {
  compounds.forEach(compound => ord.compounds.loadCompound(node, compound));
};

ord.compounds.loadCompound = function(root, compound) {
  const node = ord.compounds.add(root);
  ord.compounds.loadIntoCompound(node, compound);
};

ord.compounds.loadIntoCompound = function(node, compound) {
  const reactionRole = compound.getReactionRole();
  setSelector($('.component_reaction_role', node), reactionRole);
  $('.component_reaction_role', node).trigger('change')

  const isLimiting = compound.hasIsLimiting() ? compound.getIsLimiting() : null;
  setOptionalBool($('.component_limiting', node), isLimiting);

  const solutes = compound.hasVolumeIncludesSolutes() ?
      compound.getVolumeIncludesSolutes() :
      null;
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

ord.compounds.loadIdentifier = function(compoundNode, identifier) {
  const node = ord.compounds.addIdentifier(compoundNode);
  const value = identifier.getValue();
  $('.component_identifier_value', node).text(value);
  setSelector(node, identifier.getType());
  $('.component_identifier_details', node).text(identifier.getDetails());
};

ord.compounds.loadPreparation = function(node, preparation) {
  const type = preparation.getType();
  setSelector($('.component_compound_preparation_type', node), type);
  const details = preparation.getDetails();
  $('.component_compound_preparation_details', node).text(details);
  const reaction = preparation.getReactionId();
  $('.component_compound_preparation_reaction', node).text(reaction);
};

ord.compounds.loadVendor = function(
    compoundNode, vendorSource, vendorLot, vendorId) {
  const node = $('fieldset.vendor', compoundNode);
  $('.component_vendor_source', node).text(vendorSource);
  $('.component_vendor_lot', node).text(vendorLot);
  $('.component_vendor_id', node).text(vendorId);
};

ord.compounds.unload = function(node) {
  const compounds = [];
  $('.component', node).each(function(index, compoundNode) {
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

ord.compounds.unloadCompound = function(node) {
  const compound = new proto.ord.Compound();

  const reactionRole = getSelector($('.component_reaction_role', node));
  compound.setReactionRole(reactionRole);

  // Only call setIsLimiting if this is a reactant Compound.
  if (getSelectorText($('.component_reaction_role', node)[0]) === 'REACTANT') {
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
  $('.component_preparation', node).each(function(index, preparationNode) {
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

ord.compounds.unloadIdentifiers = function(node) {
  const identifiers = [];
  $('.component_identifier', node).each(function(index, node) {
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

ord.compounds.unloadIdentifier = function(node) {
  const identifier = new proto.ord.CompoundIdentifier();

  const value = $('.component_identifier_value', node).text();
  if (!isEmptyMessage(value)) {
    identifier.setValue(value);
  }
  const type = getSelector(node);
  identifier.setType(type);
  const details = $('.component_identifier_details', node).text();
  identifier.setDetails(details);
  return identifier;
};

ord.compounds.unloadPreparation = function(node) {
  const preparation = new proto.ord.CompoundPreparation();
  const type = getSelector($('.component_compound_preparation_type', node));
  preparation.setType(type);
  const details = $('.component_compound_preparation_details', node).text();
  preparation.setDetails(details);
  const reaction = $('.component_compound_preparation_reaction', node).text();
  preparation.setReactionId(reaction);
  return preparation;
};

ord.compounds.unloadVendor = function(node, compound) {
  const vendorSource = $('.component_vendor_source', node).text();
  compound.setVendorSource(vendorSource);
  const vendorLot = $('.component_vendor_lot', node).text();
  compound.setVendorLot(vendorLot);
  const vendorId = $('.component_vendor_id', node).text();
  compound.setVendorId(vendorId);
};

ord.compounds.add = function(root) {
  const node = addSlowly('#component_template', $('.components', root));

  // Connect reaction role selection to limiting reactant field.
  const roleSelector = $('.component_reaction_role', node);
  roleSelector.change(function() {
    if (getSelectorText(this) === 'REACTANT') {
      $('.limiting_reactant', node).show();
    } else {
      $('.limiting_reactant', node).hide();
    }
  });

  // Create an "amount" radio button group and connect it to the unit selectors.
  const amountButtons = $('.amount input', node);
  amountButtons.attr('name', 'compounds_' + ord.compounds.radioGroupCounter++);
  amountButtons.change(function() {
    $('.amount .selector', node).hide();
    if (this.value == 'mass') {
      $('.component_amount_units_mass', node).show();
      $('.includes_solutes', node).hide();
    }
    if (this.value == 'moles') {
      $('.component_amount_units_moles', node).show();
      $('.includes_solutes', node).hide();
    }
    if (this.value == 'volume') {
      $('.component_amount_units_volume', node).show();
      $('.includes_solutes', node).show().css('display', 'inline-block');
    }
  });

  // Add live validation handling.
  addChangeHandler(node, () => {ord.compounds.validateCompound(node)});

  return node;
};

ord.compounds.addIdentifier = function(node) {
  const identifierNode =
      addSlowly('#component_identifier_template', $('.identifiers', node));

  const uploadButton = $('.component_identifier_upload', identifierNode);
  uploadButton.change(function() {
    if ($(this).is(':checked')) {
      $('.uploader', identifierNode).show();
      $('.component_identifier_value', identifierNode).hide();
      $('.text_upload', identifierNode).hide();
    } else {
      $('.uploader', identifierNode).hide();
      $('.component_identifier_value', identifierNode).show();
    }
  });
  ord.uploads.initialize(identifierNode);
  return identifierNode;
};

// Shortcut to add an identifier based on name.
ord.compounds.addNameIdentifier = function(node) {
  var name = prompt('Compound name: ');
  if (!(name)) {
    return;
  }
  const identifier = new proto.ord.CompoundIdentifier();
  identifier.setValue(name);
  identifier.setType(proto.ord.CompoundIdentifier.IdentifierType.NAME);
  ord.compounds.loadIdentifier(node, identifier);

  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/resolve/name');
  xhr.responseType = 'json';
  xhr.onload = function() {
    if (xhr.response) {
      const smiles = xhr.response[0];
      const resolver = xhr.response[1];
      const identifier = new proto.ord.CompoundIdentifier();
      identifier.setValue(smiles);
      identifier.setType(proto.ord.CompoundIdentifier.IdentifierType.SMILES);
      identifier.setDetails('NAME resolved by the ' + resolver);
      ord.compounds.loadIdentifier(node, identifier);
    };
    ord.compounds.validateCompound(node);
  };
  xhr.send(name);
};

// Shortcut to add an identifier by drawing.
ord.compounds.drawIdentifier = function(node) {
  // Get a reference to Ketcher, and to look nice, clear any old drawings.
  const ketcher =
      document.getElementById('ketcher-iframe').contentWindow.ketcher;
  ketcher.editor.struct(null);
  // Start the loading spinner.
  $('#ketcher-spinner').fadeIn(0);

  // First, pack the current Compound into a message.
  const compound = new proto.ord.Compound();
  const identifiers = ord.compounds.unloadIdentifiers(node);
  if (!isEmptyMessage(identifiers)) {
    compound.setIdentifiersList(identifiers);
  }
  // Then, try to resolve compound into a MolBlock.
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/ketcher/molfile');
  const binary = compound.serializeBinary();
  xhr.responseType = 'json';
  xhr.onload = function() {
    if ((xhr.status == 200)) {
      const molblock = xhr.response;
      // Set the molecule in ketcher.
      // Note: In case async / callback issues prove difficult,
      // a cleaner fix may be to put this entire xhr in a modal callback, then
      // toggle the modal.

      // If the modal is already open, we can simply set the molecule.
      const ketcherModal = $('#ketcher_modal');
      if (ketcherModal.hasClass('show')) {
        ketcher.setMolecule(molblock);
      }
      // Otherwise, we need to set up a callback, so that the molecule is set
      // only when Ketcher is open. (to prevent a graphical glitch)
      else {
        ketcherModal.on('shown.bs.modal', function() {
          // This callback should only be ever run once, so make sure to remove
          // it.
          ketcherModal.off('shown.bs.modal');
          ketcher.setMolecule(molblock);
        })
      }
    }
    // Now that we're done with (trying to) loading the molecule, hide the
    // spinner.
    $('#ketcher-spinner').fadeOut();
  };
  xhr.send(binary);
  // Finally, open the ketcher modal.
  $('#ketcher_modal').modal('toggle');
  // Define a callback so that when a user is done drawing, the new SMILES
  // string gets saved.
  ketcher.successCallback = function() {
    // Check if an existing SMILES/MolBlock identifier exists. If yes, remove.
    $('.component_identifier', node).each(function(index, node) {
      node = $(node);
      if (!node.attr('id')) {
        // Not a template.
        const identifier = ord.compounds.unloadIdentifier(node);
        if ((identifier.getType() ===
             proto.ord.CompoundIdentifier.IdentifierType.SMILES) ||
            (identifier.getType() ===
             proto.ord.CompoundIdentifier.IdentifierType.MOLBLOCK)) {
          removeSlowly(node, '.component_identifier');
        }
      }
    });
    // Create new identifiers.
    if (ketcher.getSmiles()) {
      const identifier = new proto.ord.CompoundIdentifier();
      identifier.setType(proto.ord.CompoundIdentifier.IdentifierType.SMILES);
      identifier.setValue(ketcher.getSmiles());
      identifier.setDetails('Drawn with Ketcher');
      ord.compounds.loadIdentifier(node, identifier);
      identifier.setType(proto.ord.CompoundIdentifier.IdentifierType.MOLBLOCK);
      identifier.setValue(ketcher.getMolfile());
      ord.compounds.loadIdentifier(node, identifier);
    }
    ord.compounds.validateCompound(node);
  }
};

ord.compounds.addPreparation = function(node) {
  const PreparationNode =
      addSlowly('#component_preparation_template', $('.preparations', node));

  const typeSelector =
      $('.component_compound_preparation_type', PreparationNode);
  typeSelector.change(function() {
    if (getSelectorText(this) == 'SYNTHESIZED') {
      $('.component_compound_preparation_reaction_id', PreparationNode)
          .css('display', 'inline-block');
    } else {
      $('.component_compound_preparation_reaction_id', PreparationNode)
          .css('display', 'none');
    }
  });

  return PreparationNode
};

// Update the image tag with a drawing of this component.
ord.compounds.renderCompound = function(node, compound) {
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/render/compound');
  const binary = compound.serializeBinary();
  xhr.responseType = 'json';
  xhr.onload = function() {
    const png_data = xhr.response;
    if (png_data) {
      $('.component_rendering', node)[0].src = 'data:image/png;base64,' + png_data;
    } else {
      $('.component_rendering', node)[0].src = '';
    }
  };
  xhr.send(binary);
};

ord.compounds.validateCompound = function(node, validateNode) {
  const compound = ord.compounds.unloadCompound(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(compound, 'Compound', validateNode);

  // Try to resolve compound structural identifiers. This is tied to
  // validation so the same trigger is used and we only have to unload the
  // compound once per update.
  ord.compounds.renderCompound(node, compound);
};
