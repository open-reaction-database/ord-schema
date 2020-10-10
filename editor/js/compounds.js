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

goog.module('ord.compounds');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  unloadCompound,
  loadIntoCompound,
  add,
  validateCompound,
  drawIdentifier,
  addNameIdentifier,
  addIdentifier,
  addPreparation
};

goog.require('ord.amounts');
goog.require('ord.features');
goog.require('proto.ord.Compound');
goog.require('proto.ord.CompoundIdentifier');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the form's fields describing multiple compounds for a
 * single reaction input.
 * @param {!Node} node The div corresponding to the reaction input to which
 *     compound definitions should be added.
 * @param {!Array<!proto.ord.Compound>} compounds
 */
function load(node, compounds) {
  compounds.forEach(compound => loadCompound(node, compound));
}

/**
 * Adds fields describing a new component to an existing reaction input in the
 * form and populates them according to the provided compound.
 * @param {!Node} root The div corresponding to the reaction input to which
 *     a new compound definition should be added.
 * @param {!proto.ord.Compound} compound
 */
function loadCompound(root, compound) {
  const node = add(root);
  loadIntoCompound(node, compound);
}

/**
 * Adds and populates the form's fields describing a compound.
 * @param {!Node} node The div corresponding to the compound whose fields
 *     should be updated.
 * @param {!proto.ord.Compound} compound
 */
function loadIntoCompound(node, compound) {
  const reactionRole = compound.getReactionRole();
  ord.reaction.setSelector($('.component_reaction_role', node), reactionRole);
  $('.component_reaction_role', node).trigger('change');

  const isLimiting = compound.hasIsLimiting() ? compound.getIsLimiting() : null;
  ord.reaction.setOptionalBool($('.component_limiting', node), isLimiting);

  const solutes = compound.hasVolumeIncludesSolutes() ?
      compound.getVolumeIncludesSolutes() :
      null;
  ord.reaction.setOptionalBool($('.component_includes_solutes', node), solutes);

  const identifiers = compound.getIdentifiersList();
  identifiers.forEach(identifier => loadIdentifier(node, identifier));

  const mass = compound.getMass();
  const moles = compound.getMoles();
  const volume = compound.getVolume();
  ord.amounts.load(node, mass, moles, volume);

  const preparations = compound.getPreparationsList();
  preparations.forEach(preparation => {
    const preparationNode = addPreparation(node);
    loadPreparation(preparationNode, preparation);
  });
  const vendorSource = compound.getVendorSource();
  const vendorLot = compound.getVendorLot();
  const vendorId = compound.getVendorId();
  loadVendor(node, vendorSource, vendorLot, vendorId);

  const features = compound.getFeaturesList();
  ord.features.load(node, features);
}

/**
 * Adds fields describing a new identifier to an existing compound in the form
 * and populates them according to the provided identifier.
 * @param {!Node} compoundNode The div corresponding to the compound to which a
 *     new compound definition should be added.
 * @param {!proto.ord.CompoundIdentifier} identifier
 */
function loadIdentifier(compoundNode, identifier) {
  const node = addIdentifier(compoundNode);
  const value = identifier.getValue();
  $('.component_identifier_value', node).text(value);
  ord.reaction.setSelector(node, identifier.getType());
  $('.component_identifier_details', node).text(identifier.getDetails());
}

/**
 * Adds and populates the form's fields describing a compound preparation.
 * @param {!Node} node The div corresponding to the preparation that should be
 *     updated on the form.
 * @param {!proto.ord.CompoundPreparation} preparation
 */
function loadPreparation(node, preparation) {
  const type = preparation.getType();
  ord.reaction.setSelector(
      $('.component_compound_preparation_type', node), type);
  const details = preparation.getDetails();
  $('.component_compound_preparation_details', node).text(details);
  const reaction = preparation.getReactionId();
  $('.component_compound_preparation_reaction', node).text(reaction);
}

/**
 * Adds and populates the form's fields describing a compound's source.
 * @param {!Node} compoundNode The div corresponding to the compound whose
 *     vendor information should be updated on the form.
 * @param {string} vendorSource
 * @param {string} vendorLot
 * @param {string} vendorId
 */
function loadVendor(compoundNode, vendorSource, vendorLot, vendorId) {
  const node = $('fieldset.vendor', compoundNode);
  $('.component_vendor_source', node).text(vendorSource);
  $('.component_vendor_lot', node).text(vendorLot);
  $('.component_vendor_id', node).text(vendorId);
}

/**
 * Reads and returns a list of compounds defined within part of the form.
 * @param {!Node} node The div corresponding to the reaction inputs whose
 *     compounds should be read from the form.
 * @return {!Array<!proto.ord.Compound>}
 */
function unload(node) {
  const compounds = [];
  $('.component', node).each(function(index, compoundNode) {
    compoundNode = $(compoundNode);
    if (!compoundNode.attr('id')) {
      // Not a template.
      const compound = unloadCompound(compoundNode);
      if (!ord.reaction.isEmptyMessage(compound)) {
        compounds.push(compound);
      }
    }
  });
  return compounds;
}

/**
 * Reads and returns a single compound as defined on the form.
 * @param {!Node} node The div corresponding to the compound whose definition
 *     should be read from the form.
 * @return {!proto.ord.Compound}
 */
function unloadCompound(node) {
  const compound = new proto.ord.Compound();

  const reactionRole =
      ord.reaction.getSelector($('.component_reaction_role', node));
  compound.setReactionRole(reactionRole);

  // Only call setIsLimiting if this is a reactant Compound.
  if (ord.reaction.getSelectorText($('.component_reaction_role', node)[0]) ===
      'REACTANT') {
    const isLimiting =
        ord.reaction.getOptionalBool($('.component_limiting', node));
    compound.setIsLimiting(isLimiting);
  }

  // Only call setVolumeIncludesSolutes if the amount is defined as a volume.
  if (ord.amounts.unloadVolume(node)) {
    const solutes =
        ord.reaction.getOptionalBool($('.component_includes_solutes', node));
    compound.setVolumeIncludesSolutes(solutes);
  }

  const identifiers = unloadIdentifiers(node);
  if (!ord.reaction.isEmptyMessage(identifiers)) {
    compound.setIdentifiersList(identifiers);
  }

  ord.amounts.unload(node, compound);

  const preparations = [];
  $('.component_preparation', node).each(function(index, preparationNode) {
    const preparation = unloadPreparation(preparationNode);
    if (!ord.reaction.isEmptyMessage(preparation)) {
      preparations.push(preparation);
    }
  });
  compound.setPreparationsList(preparations);

  unloadVendor(node, compound);

  const features = ord.features.unload(node);
  compound.setFeaturesList(features);

  return compound;
}

/**
 * Reads and returns a list of compound identifiers for a single compound as
 * defined on the form.
 * @param {!Node} node The div corresponding to the compound whose identifiers
 *     should be read from the form.
 * @return {!Array<!proto.ord.CompoundIdentifier>}
 */
function unloadIdentifiers(node) {
  const identifiers = [];
  $('.component_identifier', node).each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      const identifier = unloadIdentifier(node);
      if (!ord.reaction.isEmptyMessage(identifier)) {
        identifiers.push(identifier);
      }
    }
  });
  return identifiers;
}

/**
 * Reads and returns a single compound identifier as defined on the form.
 * @param {!Node} node The div corresponding to the compound identifier that
 *     should be read from the form.
 * @return {!proto.ord.CompoundIdentifier}
 */
function unloadIdentifier(node) {
  const identifier = new proto.ord.CompoundIdentifier();

  const value = $('.component_identifier_value', node).text();
  if (!ord.reaction.isEmptyMessage(value)) {
    identifier.setValue(value);
  }
  const type = ord.reaction.getSelector(node);
  identifier.setType(type);
  const details = $('.component_identifier_details', node).text();
  identifier.setDetails(details);
  return identifier;
}

/**
 * Reads and returns a single compound preparation as defined on the form.
 * @param {!Node} node The div corresponding to a compound preparation that
 *     should be read from the form.
 * @return {!proto.ord.CompoundPreparation}
 */
function unloadPreparation(node) {
  const preparation = new proto.ord.CompoundPreparation();
  const type =
      ord.reaction.getSelector($('.component_compound_preparation_type', node));
  preparation.setType(type);
  const details = $('.component_compound_preparation_details', node).text();
  preparation.setDetails(details);
  const reaction = $('.component_compound_preparation_reaction', node).text();
  preparation.setReactionId(reaction);
  return preparation;
}

/**
 * Sets the vendor information fields of a compound according to the form.
 * @param {!Node} node The div corresponding to the compound whose vendor
 *     information should be read from the form.
 * @param {!proto.ord.Compound} compound
 */
function unloadVendor(node, compound) {
  const vendorSource = $('.component_vendor_source', node).text();
  compound.setVendorSource(vendorSource);
  const vendorLot = $('.component_vendor_lot', node).text();
  compound.setVendorLot(vendorLot);
  const vendorId = $('.component_vendor_id', node).text();
  compound.setVendorId(vendorId);
}

/**
 * Adds fields to the form corresponding to a new, empty compound definition as
 * specified by the component template with ID "component_template".
 * @param {!Node} root The div within which the new compound should be added.
 * @return {!Node} The node of the new component div.
 */
function add(root) {
  const node =
      ord.reaction.addSlowly('#component_template', $('.components', root));

  // Connect reaction role selection to limiting reactant field.
  const roleSelector = $('.component_reaction_role', node);
  roleSelector.change(function() {
    if (ord.reaction.getSelectorText(this) === 'REACTANT') {
      $('.limiting_reactant', node).show();
    } else {
      $('.limiting_reactant', node).hide();
    }
  });

  // Create an "amount" radio button group and connect it to the unit selectors.
  const amountButtons = $('.amount input', node);
  amountButtons.attr('name', 'compounds_' + radioGroupCounter++);
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
  ord.reaction.addChangeHandler(node, () => {
    validateCompound(node);
  });

  return node;
}

/**
 * Adds fields to the form corresponding to a new, empty compound identifier as
 * specified by the component identifier template with ID
 * "component_identifier_template".
 * @param {!Node} node The div within which the new identifier should be added.
 * @return {!Node} The node of the new compound identifier div.
 */
function addIdentifier(node) {
  const identifierNode = ord.reaction.addSlowly(
      '#component_identifier_template', $('.identifiers', node));

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
}

/**
 * Adds new compound identifier(s) after prompting the user for the name of
 * a compound to add. A NAME-type identifier is always included. A SMILES-type
 * identifier is included when the name can be parsed.
 * @param {!Node} node The div corresponding to the compound to which the new
 *     identifiers should be added.
 */
function addNameIdentifier(node) {
  const name = prompt('Compound name: ');
  if (!(name)) {
    return;
  }
  const identifier = new proto.ord.CompoundIdentifier();
  identifier.setValue(name);
  identifier.setType(proto.ord.CompoundIdentifier.IdentifierType.NAME);
  loadIdentifier(node, identifier);

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
      loadIdentifier(node, identifier);
    }
    validateCompound(node);
  };
  xhr.send(name);
}

/**
 * Displays the Ketcher drawing tool for defining compound identifiers. A
 * callback is defined to create SMILES and MOLBLOCK-type identifiers; this
 * callback is triggered upon submission of a molecular drawing from the
 * Ketcher window.
 * @param {!Node} node The div corresponding to the compound to which the new
 *     identifiers should be added.
 */
function drawIdentifier(node) {
  // Get a reference to Ketcher, and to look nice, clear any old drawings.
  const ketcher =
      document.getElementById('ketcher-iframe').contentWindow.ketcher;
  ketcher.editor.struct(null);
  // Start the loading spinner.
  $('#ketcher-spinner').fadeIn(0);

  // First, pack the current Compound into a message.
  const compound = new proto.ord.Compound();
  const identifiers = unloadIdentifiers(node);
  if (!ord.reaction.isEmptyMessage(identifiers)) {
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
        });
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
        const identifier = unloadIdentifier(node);
        if ((identifier.getType() ===
             proto.ord.CompoundIdentifier.IdentifierType.SMILES) ||
            (identifier.getType() ===
             proto.ord.CompoundIdentifier.IdentifierType.MOLBLOCK)) {
          ord.reaction.removeSlowly(node, '.component_identifier');
        }
      }
    });
    // Create new identifiers.
    if (ketcher.getSmiles()) {
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/canonicalize');
      xhr.responseType = 'json';
      xhr.onload = function() {
        const smilesIdentifier = new proto.ord.CompoundIdentifier();
        smilesIdentifier.setType(
            proto.ord.CompoundIdentifier.IdentifierType.SMILES);
        const smiles = xhr.response;
        smilesIdentifier.setValue(smiles);
        smilesIdentifier.setDetails('Drawn with Ketcher');
        loadIdentifier(node, smilesIdentifier);
      };
      xhr.send(ketcher.getSmiles());
      const molfileIdentifier = new proto.ord.CompoundIdentifier();
      molfileIdentifier.setType(
          proto.ord.CompoundIdentifier.IdentifierType.MOLBLOCK);
      molfileIdentifier.setValue(ketcher.getMolfile());
      molfileIdentifier.setDetails('Drawn with Ketcher');
      loadIdentifier(node, molfileIdentifier);
    }
    validateCompound(node);
  };
}

/**
 * Adds fields to the form corresponding to a new, empty compound preparation
 * as specified by the component preparation template with ID
 * "component_preparation_template".
 * @param {!Node} node The div corresponding to the compound to which the new
 *     preparation should be added.
 * @return {!Node} The div corresponding to the new compound preparation.
 */
function addPreparation(node) {
  const PreparationNode = ord.reaction.addSlowly(
      '#component_preparation_template', $('.preparations', node));

  const typeSelector =
      $('.component_compound_preparation_type', PreparationNode);
  typeSelector.change(function() {
    if (ord.reaction.getSelectorText(this) == 'SYNTHESIZED') {
      $('.component_compound_preparation_reaction_id', PreparationNode)
          .css('display', 'inline-block');
    } else {
      $('.component_compound_preparation_reaction_id', PreparationNode)
          .css('display', 'none');
    }
  });

  return PreparationNode;
}

/**
 * Updates a png rendering of a compound as defined by its identifiers.
 * @param {!Node} node The div corresponding to the compound whose rendering
 *     should be updated.
 * @param {!proto.ord.Compound} compound
 */
function renderCompound(node, compound) {
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/render/compound');
  const binary = compound.serializeBinary();
  xhr.responseType = 'json';
  xhr.onload = function() {
    const svg_data = xhr.response;
    if (svg_data) {
      $('.component_rendering', node).html(svg_data);
    } else {
      $('.component_rendering', node).html('');
    }
  };
  xhr.send(binary);
}

/**
 * Validates the definition of a compound and updates the validation error
 * display node.
 * @param {!Node} node The div corresponding to the compound that should be
 *     read from the form and validated.
 * @param {?Node} validateNode The div that is used to show the results of
 *     validation (i.e., success or errors).
 */
function validateCompound(node, validateNode) {
  const compound = unloadCompound(node);
  ord.reaction.validate(compound, 'Compound', node, validateNode);

  // Try to resolve compound structural identifiers. This is tied to
  // validation so the same trigger is used and we only have to unload the
  // compound once per update.
  renderCompound(node, compound);
}
