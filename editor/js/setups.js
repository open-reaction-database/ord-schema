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

goog.module('ord.setups');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  addVesselPreparation,
  validateSetup
};

goog.require('ord.codes');
goog.require('proto.ord.ReactionSetup');
goog.require('proto.ord.Vessel');
goog.require('proto.ord.Volume');

/**
 * Adds and populates the reaction setup section in the form.
 * @param {!proto.ord.ReactionSetup} setup
 */
function load(setup) {
  const vessel = setup.getVessel();
  if (vessel) {
    loadVessel(vessel);
  }
  const isAutomated = setup.hasIsAutomated() ? setup.getIsAutomated() : null;
  ord.reaction.setOptionalBool($('#setup_automated'), isAutomated);
  if (isAutomated) {
    $('#automation_platform').show();
  }

  $('#setup_automated').change(function() {
    if (ord.reaction.getSelectorText(this) == 'TRUE') {
      $('#automation_platform').show();
    } else {
      $('#automation_platform').hide();
    }
  });

  const platform = setup.getAutomationPlatform();
  $('#setup_platform').text(platform);

  const codes = setup.getAutomationCodeMap();
  ord.codes.load(codes);

  const environment = setup.getEnvironment();
  if (environment != null) {
    ord.reaction.setSelector(
        $('#setup_environment_type'), environment.getType());
    $('#setup_environment_details').text(environment.getDetails());
  }
}

/**
 * Adds and populates the reaction vessel section of the form.
 * @param {!proto.ord.Vessel} vessel
 */
function loadVessel(vessel) {
  const type = vessel.getType();
  if (type) {
    ord.reaction.setSelector($('#setup_vessel_type'), type.getType());
    $('#setup_vessel_details').text(type.getDetails());
  }
  const material = vessel.getMaterial();
  if (material) {
    ord.reaction.setSelector($('#setup_vessel_material'), material.getType());
    $('#setup_vessel_material_details').text(material.getDetails());
  }
  const preparations = vessel.getPreparationsList();
  preparations.forEach(preparation => {
    const node = addVesselPreparation();
    ord.reaction.setSelector(
        $('.setup_vessel_preparation_type', node), preparation.getType());
    $('.setup_vessel_preparation_details', node).text(preparation.getDetails());
  });
  const attachments = vessel.getAttachmentsList();
  attachments.forEach(attachment => {
    const node = addVesselAttachment();
    ord.reaction.setSelector(
        $('.setup_vessel_attachment_type', node), attachment.getType());
    $('.setup_vessel_attachment_details', node).text(attachment.getDetails());
  });
  if (vessel.hasVolume()) {
    const volume = vessel.getVolume();
    ord.reaction.writeMetric('#setup_vessel_volume', volume);
  }
}

/**
 * Fetches the reaction setup from the form.
 * @return {!proto.ord.ReactionSetup}
 */
function unload() {
  const setup = new proto.ord.ReactionSetup();

  const vessel = unloadVessel();
  if (!ord.reaction.isEmptyMessage(vessel)) {
    setup.setVessel(vessel);
  }

  const isAutomated = ord.reaction.getOptionalBool($('#setup_automated'));
  setup.setIsAutomated(isAutomated);

  const platform = $('#setup_platform').text();
  setup.setAutomationPlatform(platform);

  const codes = setup.getAutomationCodeMap();
  ord.codes.unload(codes);

  const environment = new proto.ord.ReactionSetup.ReactionEnvironment();
  environment.setType(ord.reaction.getSelector($('#setup_environment_type')));
  environment.setDetails($('#setup_environment_details').text());
  if (!ord.reaction.isEmptyMessage(environment)) {
    setup.setEnvironment(environment);
  }

  return setup;
}

/**
 * Fetches the reaction vessel information from the form.
 * @return {!proto.ord.Vessel}
 */
function unloadVessel() {
  const vessel = new proto.ord.Vessel();

  const type = new proto.ord.VesselType();
  type.setType(ord.reaction.getSelector($('#setup_vessel_type')));
  type.setDetails($('#setup_vessel_details').text());
  if (!ord.reaction.isEmptyMessage(type)) {
    vessel.setType(type);
  }

  const material = new proto.ord.VesselMaterial();
  material.setType(ord.reaction.getSelector('#setup_vessel_material'));
  material.setDetails($('#setup_vessel_material_details').text());
  if (!ord.reaction.isEmptyMessage(material)) {
    vessel.setMaterial(material);
  }

  const preparations = [];
  $('.setup_vessel_preparation').each(function(index, node) {
    node = $(node);
    if (node.attr('id')) {
      // The template.
      return;
    }
    const preparation = new proto.ord.VesselPreparation();
    preparation.setType(
        ord.reaction.getSelector($('.setup_vessel_preparation_type', node)));
    preparation.setDetails($('.setup_vessel_preparation_details', node).text());
    if (!ord.reaction.isEmptyMessage(preparation)) {
      preparations.push(preparation);
    }
  });
  vessel.setPreparationsList(preparations);

  const attachments = [];
  $('.setup_vessel_attachment').each(function(index, node) {
    node = $(node);
    if (node.attr('id')) {
      // The template.
      return;
    }
    const attachment = new proto.ord.VesselAttachment();
    attachment.setType(
        ord.reaction.getSelector($('.setup_vessel_attachment_type', node)));
    attachment.setDetails($('.setup_vessel_attachment_details', node).text());
    if (!ord.reaction.isEmptyMessage(attachment)) {
      attachments.push(attachment);
    }
  });
  vessel.setAttachmentsList(attachments);

  const volume =
      ord.reaction.readMetric('#setup_vessel_volume', new proto.ord.Volume());
  if (!ord.reaction.isEmptyMessage(volume)) {
    vessel.setVolume(volume);
  }

  return vessel;
}

/**
 * Adds a new vessel preparation section to the form.
 * @return {!Node} The node of the newly added div.
 */
function addVesselPreparation() {
  return ord.reaction.addSlowly(
      '#setup_vessel_preparation_template', '#setup_vessel_preparations');
}

/**
 * Adds a new vessel attachment section to the form.
 * @return {!Node} The node of the newly added div.
 */
function addVesselAttachment() {
  return ord.reaction.addSlowly(
      '#setup_vessel_attachment_template', '#setup_vessel_attachments');
}

/**
 * Validates the reaction setup defined in the form.
 * @param {!Node} node The node containing the reaction setup div.
 * @param {?Node} validateNode The target div for validation results.
 */
function validateSetup(node, validateNode) {
  const setup = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(setup, 'ReactionSetup', validateNode);
}
