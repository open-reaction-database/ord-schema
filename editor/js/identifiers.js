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

goog.provide('ord.identifiers');

goog.require('ord.uploads');
goog.require('proto.ord.ReactionIdentifier');

ord.identifiers.load = function (identifiers) {
  identifiers.forEach(identifier => ord.identifiers.loadIdentifier(identifier));
};

ord.identifiers.loadIdentifier = function (identifier) {
  const node = ord.identifiers.add();
  const bytesValue = identifier.getBytesValue();
  if (bytesValue) {
    $('.reaction_identifier_upload', node).prop('checked', true);
    $('.reaction_identifier_value', node).hide();
    ord.uploads.load(node, bytesValue);
  } else {
    const value = identifier.getValue();
    $('.reaction_identifier_value', node).text(value);
  }
  setSelector(node, identifier.getType());
  $('.reaction_identifier_details', node).text(identifier.getDetails());
};

ord.identifiers.unload = function () {
  const identifiers = [];
  $('.reaction_identifier').each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      const identifier = ord.identifiers.unloadIdentifier(node);
      if (!isEmptyMessage(identifier)) {
        identifiers.push(identifier);
      }
    }
  });
  return identifiers;
};

ord.identifiers.unloadIdentifier = function (node) {
  const identifier = new proto.ord.ReactionIdentifier();

  if ($('.reaction_identifier_upload', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!isEmptyMessage(bytesValue)) {
      identifier.setBytesValue(bytesValue);
    }
  }
  else {
    const value = $('.reaction_identifier_value', node).text();
    if (!isEmptyMessage(value)) {
      identifier.setValue(value);
    }
  }
  const type = getSelector(node);
  if (!isEmptyMessage(type)) {
    identifier.setType(type);
  }
  const details = $('.reaction_identifier_details', node).text();
  if (!isEmptyMessage(details)) {
    identifier.setDetails(details);
  }
  return identifier;
};

ord.identifiers.add = function () {
  const node = addSlowly('#reaction_identifier_template', '#identifiers');

  const uploadButton = $('.reaction_identifier_upload', node);
  uploadButton.change(function () {
    if ($(this).is(':checked')) {
      $('.uploader', node).show();
      $('.reaction_identifier_value', node).hide();
    } else {
      $('.uploader', node).hide();
      $('.reaction_identifier_value', node).show();
    }
  });
  ord.uploads.initialize(node);
  return node;
};
