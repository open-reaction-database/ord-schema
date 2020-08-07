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

goog.provide('ord.setups');

goog.require('ord.codes');
goog.require('proto.ord.ReactionSetup');
goog.require('proto.ord.Vessel');
goog.require('proto.ord.Volume');

ord.setups.load = function(setup) {
  const vessel = setup.getVessel();
  if (vessel) {
    ord.setups.loadVessel(vessel);
  }
  const isAutomated = setup.hasIsAutomated() ? setup.getIsAutomated() : null;
  setOptionalBool($('#setup_automated'), isAutomated);
  if (isAutomated) {
    $('#automation_platform').show();
  }

  $('#setup_automated').change(function() {
    if (getSelectorText(this) == 'TRUE') {
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
    setSelector($('#setup_environment_type'), environment.getType());
    $('#setup_environment_details').text(environment.getDetails());
  }
};

ord.setups.loadVessel = function(vessel) {
  const type = vessel.getType();
  if (type) {
    setSelector($('#setup_vessel_type'), type.getType());
    $('#setup_vessel_details').text(type.getDetails());
  }
  const material = vessel.getMaterial();
  if (material) {
    setSelector($('#setup_vessel_material'), material.getType());
    $('#setup_vessel_material_details').text(material.getDetails());
  }
  const preparations = vessel.getPreparationsList();
  preparations.forEach(preparation => {
    const node = ord.setups.addVesselPreparation();
    setSelector(
        $('.setup_vessel_preparation_type', node), preparation.getType());
    $('.setup_vessel_preparation_details', node).text(preparation.getDetails());
  });
  const attachments = vessel.getAttachmentsList();
  attachments.forEach(attachment => {
    const node = ord.setups.addVesselAttachment();
    setSelector($('.setup_vessel_attachment_type', node), attachment.getType());
    $('.setup_vessel_attachment_details', node).text(attachment.getDetails());
  });
  if (vessel.hasVolume()) {
    const volume = vessel.getVolume();
    writeMetric('#setup_vessel_volume', volume);
  }
};

ord.setups.unload = function() {
  const setup = new proto.ord.ReactionSetup();

  const vessel = ord.setups.unloadVessel();
  if (!isEmptyMessage(vessel)) {
    setup.setVessel(vessel);
  }

  const isAutomated = getOptionalBool($('#setup_automated'));
  setup.setIsAutomated(isAutomated);

  const platform = $('#setup_platform').text();
  setup.setAutomationPlatform(platform);

  const codes = setup.getAutomationCodeMap();
  ord.codes.unload(codes);

  const environment = new proto.ord.ReactionSetup.ReactionEnvironment();
  environment.setType(getSelector($('#setup_environment_type')));
  environment.setDetails($('#setup_environment_details').text());
  if (!isEmptyMessage(environment)) {
    setup.setEnvironment(environment);
  }

  return setup;
};

ord.setups.unloadVessel = function() {
  const vessel = new proto.ord.Vessel();

  type = new proto.ord.VesselType();
  type.setType(getSelector($('#setup_vessel_type')));
  type.setDetails($('#setup_vessel_details').text());
  if (!isEmptyMessage(type)) {
    vessel.setType(type);
  }

  const material = new proto.ord.VesselMaterial();
  material.setType(getSelector('#setup_vessel_material'));
  material.setDetails($('#setup_vessel_material_details').text());
  if (!isEmptyMessage(material)) {
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
    preparation.setType(getSelector($('.setup_vessel_preparation_type', node)));
    preparation.setDetails($('.setup_vessel_preparation_details', node).text());
    if (!isEmptyMessage(preparation)) {
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
    attachment.setType(getSelector($('.setup_vessel_attachment_type', node)));
    attachment.setDetails($('.setup_vessel_attachment_details', node).text());
    if (!isEmptyMessage(attachment)) {
      attachments.push(attachment);
    }
  });
  vessel.setAttachmentsList(attachments);

  const volume = readMetric('#setup_vessel_volume', new proto.ord.Volume());
  if (!isEmptyMessage(volume)) {
    vessel.setVolume(volume);
  }

  return vessel;
};

ord.setups.addVesselPreparation = function() {
  return addSlowly(
      '#setup_vessel_preparation_template', '#setup_vessel_preparations');
};

ord.setups.addVesselAttachment = function() {
  return addSlowly(
      '#setup_vessel_attachment_template', '#setup_vessel_attachments');
};

ord.setups.validateSetup = function(node, validateNode) {
  const setup = ord.setups.unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  validate(setup, 'ReactionSetup', validateNode);
};