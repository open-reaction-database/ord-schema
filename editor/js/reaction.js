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

goog.provide('ord.reaction');

goog.require('ord.conditions');
goog.require('ord.enums');
goog.require('ord.identifiers');
goog.require('ord.inputs');
goog.require('ord.notes');
goog.require('ord.observations');
goog.require('ord.outcomes');
goog.require('ord.provenance');
goog.require('ord.setups');
goog.require('ord.workups');
goog.require('proto.ord.Dataset');

// Remember the dataset and reaction we are editing.
const session = {
  fileName: null,
  dataset: null,
  index: null // Ordinal position of the Reaction in its Dataset.
};

async function init(fileName, index) {
  session.fileName = fileName;
  session.index = index;
  // Initialize all the template popup menus.
  $('.selector').each((index, node) => initSelector($(node)));
  $('.optional_bool').each((index, node) => initOptionalBool($(node)));
  // Enable all the editable text fields.
  $('.edittext').attr('contentEditable', 'true');
  // Show "save" on modifications.
  listen('body');
  // Fetch the Dataset containing the Reaction proto.
  session.dataset = await getDataset(fileName);
  // Initialize the UI with the Reaction.
  const reaction = session.dataset.getReactionsList()[index];
  loadReaction(reaction);
  clean();
  // Signal to tests that the DOM is initialized.
  ready();
}

function ready() {
  $('body').attr('ready', true);
}

function listen(node) {
  addChangeHandler($(node), dirty);
  $('.edittext').on('focus', event => selectText(event.target));
}

function dirty() {
  $('#save').css('visibility', 'visible');
  $('#reaction_validate').css('visibility', 'visible');
}

function clean() {
  $('#save').css('visibility', 'hidden');
  $('#save').text('save');
}

function selectText(node) {
  const range = document.createRange();
  range.selectNodeContents(node);
  const selection = window.getSelection();
  selection.removeAllRanges();
  selection.addRange(range);
}

// Adds a jQuery handler to a node such that the handler is run once 
// whenever data entry within that node is changed.
function addChangeHandler (node, handler) {
  // For textboxes
  node.on('blur', '.edittext', handler);
  // For selectors, optional bool selectors,
  // and checkboxes/radio buttons/file upload, respectively
  node.on('change', '.selector, .optional_bool, input', handler);
  // For add and remove buttons
  node.on('click', '.add, .remove', handler);
}

// Generic validator for many message types, not just reaction
// note: does not commit or save anything!
function validate(message, messageTypeString, statusNode) {
  // eg message is a type of reaction, messageTypeString = "Reaction"
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/validate/' + messageTypeString);
  const binary = message.serializeBinary();

  xhr.responseType = 'json';
  xhr.onload = function () {
    const errors = xhr.response;
    resultNode = $('.validate_result', statusNode);
    messageNode = $('.validate_message', statusNode);
    if (errors.length) {
      resultNode.css('backgroundColor', 'pink');
      resultNode.text('invalid (click)');
      messageNode.text(errors);
    }
    else {
      resultNode.css('backgroundColor', 'lightgreen');
      resultNode.text('valid');
      messageNode.text('');
      messageNode.css('visibility', 'hidden');
    }
  };
  xhr.send(binary);
}

toggleValidateMessage = function(node) {
  messageNode = $('.validate_message', node);
  switch (messageNode.css('visibility')) {
    case 'visible':
      messageNode.css('visibility', 'hidden');
      break;
    case 'hidden':
      messageNode.css('visibility', 'visible');
      break;
  }
};

function validateReaction() {
  $('#reaction_validate').css('visibility', 'hidden');
  const reaction = unloadReaction();
  statusNode = $('#reaction_validate_status');
  messageNode = $('#reaction_validate_message');
  validate(reaction, "Reaction", statusNode, messageNode);
}

function commit() {
  const reaction = unloadReaction();
  const reactions = session.dataset.getReactionsList();
  reactions[session.index] = reaction;
  putDataset(session.fileName, session.dataset);
  ord.uploads.putAll(session.fileName);
}

async function getDataset(fileName) {
  return new Promise(resolve => {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', '/dataset/proto/read/' + fileName);
    xhr.responseType = 'arraybuffer';
    xhr.onload = function (event) {
      const bytes = new Uint8Array(xhr.response);
      const dataset = proto.ord.Dataset.deserializeBinary(bytes);
      resolve(dataset);
    };
    xhr.send();
  });
}

function putDataset(fileName, dataset) {
  $('#save').text('saving');
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/write/' + fileName);
  const binary = dataset.serializeBinary();
  xhr.onload = clean;
  xhr.send(binary);
}

// For testing, compare a local Dataset to a Dataset on the server.
async function compareDataset(fileName, dataset) {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open('POST', '/dataset/proto/compare/' + fileName);
    const binary = dataset.serializeBinary();
    xhr.onload = () => {
      if (xhr.status == 200) {
        resolve();
      } else {
        reject();
      }
    };
    xhr.onerror = reject;
    xhr.send(binary);
  });
}

function loadReaction(reaction) {
  const identifiers = reaction.getIdentifiersList();
  ord.identifiers.load(identifiers);
  const inputs = reaction.getInputsMap();
  ord.inputs.load(inputs);
  const setup = reaction.getSetup();
  if (setup) {
    ord.setups.load(setup);
  }
  const conditions = reaction.getConditions();
  if (conditions) {
    ord.conditions.load(conditions);
  }
  const notes = reaction.getNotes();
  if (notes) {
    ord.notes.load(notes);
  }
  const observations = reaction.getObservationsList();
  ord.observations.load(observations);

  const workups = reaction.getWorkupList();
  ord.workups.load(workups);

  const outcomes = reaction.getOutcomesList();
  ord.outcomes.load(outcomes);

  const provenance = reaction.getProvenance();
  if (provenance) {
    ord.provenance.load(provenance);
  }
  $('#reaction_id').text(reaction.getReactionId());
}

function unloadReaction() {
  const reaction = new proto.ord.Reaction();
  const identifiers = ord.identifiers.unload();
  reaction.setIdentifiersList(identifiers);
  const inputs = reaction.getInputsMap();
  ord.inputs.unload(inputs);
  const setup = ord.setups.unload();
  reaction.setSetup(setup);
  const conditions = ord.conditions.unload();
  reaction.setConditions(conditions);
  const notes = ord.notes.unload();
  reaction.setNotes(notes);
  const observations = ord.observations.unload();
  reaction.setObservationsList(observations);
  const workups = ord.workups.unload();
  reaction.setWorkupList(workups);
  const outcomes = ord.outcomes.unload();
  reaction.setOutcomesList(outcomes);
  const provenance = ord.provenance.unload();
  reaction.setProvenance(provenance);
  reaction.setReactionId($('#reaction_id').text());
  return reaction;
}

// The template is a jQuery selector. The root may be a jQuery object.
function addSlowly(template, root) {
  const node = $(template).clone();
  node.removeAttr('id');
  $(root).append(node);
  node.show('slow');
  dirty();
  listen(node);
  return node;
}

// Removes from the DOM the nearest ancestor element matching the pattern.
function removeSlowly(button, pattern) {
  const node = $(button).closest(pattern);
  node.hide('slow', () => node.remove());
  dirty();
}

// Unpack a value/units/precision triple into the given type.
function readMetric(prefix, proto, node) {
  const value = parseFloat($(prefix + '_value', node).text());
  if (!isNaN(value)) {
    proto.setValue(value);
  }
  if (proto.setUnits) {
    // proto.ord.Percentage doesn't have units.
    proto.setUnits(getSelector($(prefix + '_units', node)));
  }
  const precision = parseFloat($(prefix + '_precision', node).text());
  if (!isNaN(precision)) {
    proto.setPrecision(precision);
  }
  return proto;
}

// Pack a value/units/precision triple into the elements specified.
function writeMetric(prefix, proto, node) {
  if (proto.hasValue()) {
    $(prefix + '_value', node).text(proto.getValue());
  }
  if (proto.getUnits) {
    // proto.ord.Percentage doesn't have units.
    setSelector($(prefix + '_units', node), proto.getUnits());
  }
  if (proto.hasPrecision()) {
    $(prefix + '_precision', node).text(proto.getPrecision());
  }
}

// Populate a <select/> node according to its data-proto type declaration.
function initSelector(node) {
  const protoName = node.attr('data-proto');
  const protoEnum = nameToProto(protoName);
  if (!protoEnum) {
    console.log('missing require: "' + protoName + '"');
  }
  const options = enumToStrings(protoEnum);
  const select = $('<select>');
  for (let i = 0; i < options.length; i++) {
    const option = $('<option>').text(options[i]);
    option.attr('value', i);
    if (options[i] == 'UNSPECIFIED') {
      option.attr('selected', 'selected');
    }
    select.append(option);
  }
  node.append(select);
}

// Select an <option/> under a <select/>. The "value" is an integer.
function setSelector(node, value) {
  $('option', node).removeAttr('selected');
  $('option[value=' + value + ']', node).attr('selected', 'selected');
}

// Find the selected <option/> and map its text onto a proto Enum.
function getSelector(node) {
  return parseInt($('select', node).val());
}

// Set up the three-way popup, "true"/"false"/"unspecified".
function initOptionalBool(node) {
  const select = $('<select>');
  const options = ['UNSPECIFIED', 'TRUE', 'FALSE'];
  for (let i = 0; i < options.length; i++) {
    const option = $('<option>').text(options[i]);
    option.attr('value', options[i]);
    if (options[i] == 'UNSPECIFIED') {
      option.attr('selected', 'selected');
    }
    select.append(option);
  }
  node.append(select);
}

function setOptionalBool(node, value) {
  $('option', node).removeAttr('selected');
  if (value == true) {
    $('option[value=TRUE]', node).attr('selected', 'selected');
  }
  if (value == false) {
    $('option[value=FALSE]', node).attr('selected', 'selected');
  }
  if (value == null) {
    $('option[value=UNSPECIFIED]', node).attr('selected', 'selected');
  }
}

function getOptionalBool(node) {
  const value = $('select', node).val();
  if (value == 'TRUE') {
    return true;
  }
  if (value == 'FALSE') {
    return false;
  }
  return null;
}

// Convert a Message_Field name from a data-proto attribute into a proto class.
function nameToProto(protoName) {
  let clazz = proto.ord;
  protoName.split('_').forEach(function (name) {
    clazz = clazz[name];
    if (!clazz) {
      return null;
    }
  });
  return clazz;
}

// Convert an Enum protobuf class to an arrray of strings.
function enumToStrings(protoEnum) {
  const types = Object.entries(protoEnum);
  const strings = [];
  for (let i = 0; i < types.length; i++) {
    strings.push(types[i][0]);
  }
  return strings;
}

// Convert an Enum name string to its protobuf member int.
function stringToEnum(name, protoEnum) {
  return protoEnum[name];
}
