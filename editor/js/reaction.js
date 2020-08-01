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
  // TODO respond to start_collapsed
  $('.collapse').each((index, node) => initCollapse($(node)));
  // Initialize all the validators.
  $('.validate').each((index, node) => initValidateNode($(node)));
  // Initialize validation handlers that don't go in "add" methods.
  initValidateHandlers();
  // Initailize tooltips.
  $("[data-toggle='tooltip']").tooltip();
  // Prevent tooltip pop-ups from blurring. 
  // (see github.com/twbs/bootstrap/issues/22610)
  Popper.Defaults.modifiers.computeStyle.gpuAcceleration = false;
  // Show "save" on modifications.
  listen('body');
  // Load Ketcher content into an element with attribute role="application".
  // TODO there's a better way to do this; will fix later
  $("#ketcher-iframe")[0].contentWindow.ketcher.initKetcher();

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
// whenever data entry within that node is changed,
// *except through remove* -- this must be handled manually.
// (This prevents inconsistent timing in the ordering of 
// the element being removed and the handler being called)
function addChangeHandler (node, handler) {
  // For textboxes
  node.on('blur', '.edittext', handler);
  // For selectors, optional bool selectors,
  // and checkboxes/radio buttons/file upload, respectively
  node.on('change', '.selector, .optional_bool, input', handler);
  // For add buttons
  node.on('click', '.add', handler);
}

// Generic validator for many message types, not just reaction
// note: does not commit or save anything!
function validate(message, messageTypeString, validateNode) {
  // eg message is a type of reaction, messageTypeString = 'Reaction'
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/validate/' + messageTypeString);
  const binary = message.serializeBinary();

  xhr.responseType = 'json';
  xhr.onload = function () {
    const errors = xhr.response;
    statusNode = $('.validate_status', validateNode);
    messageNode = $('.validate_message', validateNode);
    statusNode.removeClass('fa-check');
    statusNode.removeClass('fa-exclamation-triangle');
    statusNode.css('backgroundColor', null);
    statusNode.text('');
    if (errors.length) {
      statusNode.addClass('fa fa-exclamation-triangle')
      statusNode.css('color', 'red');
      statusNode.text(' ' + errors.length);

      messageNode.empty();
      for (index = 0; index < errors.length; index++) { 
        error = errors[index];
        errorNode = $('<div></div>');
        errorNode.text('\u2022 ' + error);
        messageNode.append(errorNode);
      } 
      messageNode.css('backgroundColor', 'pink');
    }
    else {
      statusNode.addClass('fa fa-check')
      statusNode.css('color', 'green');

      messageNode.html('');
      messageNode.css('backgroundColor', '');
      messageNode.css('visibility', 'hidden');
    }
  };
  xhr.send(binary);
}

function toggleValidateMessage(node) {
  let messageNode = $('.validate_message', node);
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
  var validateNode = $('#reaction_validate');
  const reaction = unloadReaction();
  validate(reaction, 'Reaction', validateNode);
  // Trigger all submessages to validate.
  $('.validate_button:visible:not(#reaction_validate_button)').trigger('click');
}

function commit() {
  const reaction = unloadReaction();
  const reactions = session.dataset.getReactionsList();
  reactions[session.index] = reaction;
  putDataset(session.fileName, session.dataset);
  ord.uploads.putAll(session.fileName);
}

function downloadReaction() {
  const reaction = unloadReaction();
  const binary = reaction.serializeBinary();
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/reaction/download');
  xhr.onload = () => {
    // Make the browser write the file.
    const url = URL.createObjectURL(new Blob([xhr.response]));
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', 'reaction.pbtxt');
    document.body.appendChild(link);
    link.click();
  };
  xhr.send(binary);
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
  // Reactions start with an input by default.
  if (inputs.arr_.length) {
    ord.inputs.load(inputs);
  } else {
    ord.inputs.add('#inputs');
  }
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
  // Reactions start with an outcome by default.
  if (outcomes.length) {
    ord.outcomes.load(outcomes);
  } else {
    ord.outcomes.add();
  }

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
  // isEmptyMessage check occurs in inputs.unload.
  ord.inputs.unload(inputs);

  const setup = ord.setups.unload();
  if (!isEmptyMessage(setup)) {
    reaction.setSetup(setup);
  }

  const conditions = ord.conditions.unload();
  if (!isEmptyMessage(conditions)) {
    reaction.setConditions(conditions);
  }

  const notes = ord.notes.unload();
  if (!isEmptyMessage(notes)) {
    reaction.setNotes(notes);
  }

  const observations = ord.observations.unload();
  reaction.setObservationsList(observations);

  const workups = ord.workups.unload();
  reaction.setWorkupList(workups);

  const outcomes = ord.outcomes.unload();
  reaction.setOutcomesList(outcomes);

  const provenance = ord.provenance.unload();
  if (!isEmptyMessage(provenance)) {
    reaction.setProvenance(provenance);
  }

  // Setter does nothing when passed an empty string.
  reaction.setReactionId($('#reaction_id').text());
  return reaction;
}

// Checks if the argument represents an empty protobuf message 
// (that is, the argument's nested arrays only contains null or empty values),
// or is null or undefined.
// We use this check on both primitives and arrays/messages.
// Note: Unlike other primitive types, using a setter to set a oneof string field to “” 
// causes the message to include the field and “”, which would be unwanted. 
// So we instead claim that empty strings are empty messages. 
// (Hence we don’t set _any_ empty string)
// Note: In a submessage, setting a meaningful value (e.g. optional float to 0)
// will result in a non-null/undefined value in the submessage array. 
// So, if the array of a submessage only contains null and undefined vals, 
// we can assume that the message is truly “empty” (that is, 
// doesn’t have anything meaningful that is set) 
// and can be omitted when constructing the surrounding message.
function isEmptyMessage(obj) {
  if ([null, undefined, ''].includes(obj)) {
    return true;
  }
  array = obj.array;
  if (array !== undefined) {
    // message is a protobuf message, test underlying array
    return isEmptyMessage(array);
  }
  if (Array.isArray(obj)) {
    // message arg is an array, test as-is
    return obj.every(e => isEmptyMessage(e));
  }
  if (obj.byteLength !== undefined) {
    // message arg is a byteArray; has no meaning only if empty
    // (it's possible that a byteArray filled with 0's may be meaningful)
    return (obj.length == 0);
  }
  return false;
}

// The template is a jQuery selector. The root may be a jQuery object.
function addSlowly(template, root) {
  const node = $(template).clone();
  node.removeAttr('id');
  $(root).append(node);
  node.show('slow');
  dirty();
  listen(node);
  $("[data-toggle='tooltip']", node).tooltip();
  return node;
}

// Removes from the DOM the nearest ancestor element matching the pattern.
function removeSlowly(button, pattern) {
  const node = $(button).closest(pattern);
  // Must call necessary validators only after the node is removed,
  // but we can only figure out which validators these are before removal.
  // We do so, and after removal, click the corresponding buttons to trigger validation.
  let buttonsToClick = $();
  node.parents('fieldset').each(function () {
    buttonsToClick = buttonsToClick.add($(this).children('legend').find('.validate_button'));
  })
  node.hide('slow', function () {
    node.remove();
    buttonsToClick.trigger('click');
  });
  dirty();
}

// Toggle visibility of all siblings of an element, 
// or if a pattern is provided, toggle visibility of all siblings of  
// the nearest ancestor element matching the pattern.
function toggleSlowly(node, pattern) {
  node = $(node);
  if (pattern) {
    node = node.closest(pattern);
  }
  // 'collapsed' tag is used to hold previously collapsed siblings, 
  // and would be stored as node's next sibling;
  // the following line checks whether a collapse has occured.
  if (node.next('collapsed').length !== 0) {
    // Need to uncollapse.
    const collapsedNode = node.next('collapsed');
    collapsedNode.toggle('slow', () => {
      collapsedNode.children().unwrap();
    });
  }
  else {
    // Need to collapse.
    node.siblings().wrapAll('<collapsed>');
    node.next('collapsed').toggle('slow');
  }
}

function collapseToggle(button) {
  $(button).toggleClass('fa-chevron-down fa-chevron-right');
  toggleSlowly(button, 'legend');
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

// Prompt user to upload a file and return its content as a string.
function setTextFromFile(identifierNode, valueClass) {
  var input = document.createElement('input');
  input.type = 'file';
  input.onchange = event => {
    var file = event.target.files[0];
    var reader = new FileReader();
    reader.readAsText(file);
    reader.onload = readerEvent => {
      const contents = readerEvent.target.result;
      $('.' + valueClass, identifierNode).text(contents);
     }
  }
  input.click();
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

// Find the selected <option/> and return its text.
function getSelectorText(node) {
  const selectorElement = node.getElementsByTagName('select')[0]
  return selectorElement.options[selectorElement.selectedIndex].text;
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

// Set up and initialize a collapse button,
// by adding attributes into a div in reaction.html
function initCollapse (node) {
  node.addClass('fa');
  node.addClass('fa-chevron-down');
  node.attr('onclick', 'collapseToggle(this)');
  if (node.hasClass('starts_collapsed')) {
    node.trigger('click');
  }
}

// Set up a validator div (button, status indicator, error list, etc.),
// inserting contents into a div in reaction.html
function initValidateNode (oldNode) {
  let newNode = $('#validate_template').clone();
  // Add attributes necessary for validation functions.
  $('.validate_button', newNode).attr('onclick', oldNode.attr('button-onclick'));
  if (oldNode.attr('id')) {
    $('.validate_button', newNode).attr('id', oldNode.attr('id') + '_button');
  }
  oldNode.append(newNode.children());
}

// Some nodes are dynamically added / removed; 
// we add their validation handlers when the nodes themselves are added.
// However, other nodes are always present in the html, and aren't dynamically added nor removed.
// To add live validation to these nodes, we do so here.
function initValidateHandlers () {
  // For setup
  var setupNode = $('#section_setup');
  addChangeHandler(setupNode, () => {ord.setups.validateSetup(setupNode)});

  // For conditions
  var conditionNode = $('#section_conditions');
  addChangeHandler(conditionNode, () => {ord.conditions.validateConditions(conditionNode)});

  // For temperature
  var temperatureNode = $('#section_conditions_temperature');
  addChangeHandler(temperatureNode, () => {ord.temperature.validateTemperature(temperatureNode)});

  // For pressure
  var pressureNode = $('#section_conditions_pressure');
  addChangeHandler(pressureNode, () => {ord.pressure.validatePressure(pressureNode)});

  // For stirring
  var stirringNode = $('#section_conditions_stirring');
  addChangeHandler(stirringNode, () => {ord.stirring.validateStirring(stirringNode)});

  // For illumination
  var illuminationNode = $('#section_conditions_illumination');
  addChangeHandler(illuminationNode, () => {ord.illumination.validateIllumination(illuminationNode)});

  // For electro
  var electroNode = $('#section_conditions_electro');
  addChangeHandler(electroNode, () => {ord.electro.validateElectro(electroNode)});

  // For flow
  var flowNode = $('#section_conditions_flow');
  addChangeHandler(flowNode, () => {ord.flows.validateFlow(flowNode)});

  // For notes
  var notesNode = $('#section_notes');
  addChangeHandler(notesNode, () => {ord.notes.validateNotes(notesNode)});

  // For provenance
  var provenanceNode = $('#section_provenance');
  addChangeHandler(provenanceNode, () => {ord.provenance.validateProvenance(provenanceNode)});
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
