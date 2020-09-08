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

goog.module('ord.reaction');
goog.module.declareLegacyNamespace();
exports = {
  init,
  commit,
  downloadReaction,
  validateReaction,
  setTextFromFile,
  removeSlowly,
  collapseToggle,
  addSlowly,
  addChangeHandler,
  validate,
  toggleValidateMessage,
  compareDataset,
  unloadReaction,
  isEmptyMessage,
  readMetric,
  writeMetric,
  setSelector,
  getSelector,
  getSelectorText,
  setOptionalBool,
  getOptionalBool
};

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
  index: null  // Ordinal position of the Reaction in its Dataset.
};
// Export session, because it's used by test.js.
exports.session = session;

/**
 * Initializes the form.
 * @param {string} fileName Path to the current Dataset proto.
 * @param {number} index The index of this Reaction in the current Dataset.
 */
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
  $('[data-toggle=\'tooltip\']').tooltip();
  // Prevent tooltip pop-ups from blurring.
  // (see github.com/twbs/bootstrap/issues/22610)
  Popper.Defaults.modifiers.computeStyle.gpuAcceleration = false;
  // Show "save" on modifications.
  listen('body');
  // Load Ketcher content into an element with attribute role="application".
  document.getElementById('ketcher-iframe').contentWindow.ketcher.initKetcher();
  // Fetch the Dataset containing the Reaction proto.
  session.dataset = await getDataset(fileName);
  // Initialize the UI with the Reaction.
  const reaction = session.dataset.getReactionsList()[index];
  loadReaction(reaction);
  clean();
  // Trigger reaction-level validation.
  validateReaction();
  // Signal to tests that the DOM is initialized.
  ready();
}

/**
 * Sets the `ready` value to true.
 */
function ready() {
  $('body').attr('ready', true);
}

/**
 * Adds a change handler to the given node that shows the 'save' button when
 * the node text is edited.
 * @param {!Node} node
 */
function listen(node) {
  addChangeHandler($(node), dirty);
  $('.edittext').on('focus', event => selectText(event.target));
}

/**
 * Shows the 'save' button.
 */
function dirty() {
  $('#save').css('visibility', 'visible');
}

/**
 * Hides the 'save' button.
 */
function clean() {
  $('#save').css('visibility', 'hidden');
  $('#save').text('save');
}

/**
 * Selects the contents of the given node.
 * @param {!Node} node
 */
function selectText(node) {
  const range = document.createRange();
  range.selectNodeContents(node);
  const selection = window.getSelection();
  selection.removeAllRanges();
  selection.addRange(range);
}

/**
 * Adds a jQuery handler to a node such that the handler is run once whenever
 * data entry within that node is changed, *except through remove* -- this must
 * be handled manually. (This prevents inconsistent timing in the ordering of
 * the element being removed and the handler being called.)
 * @param {!Node} node
 * @param {!Function} handler
 */
function addChangeHandler(node, handler) {
  // For textboxes
  node.on('blur', '.edittext', handler);
  // For selectors, optional bool selectors,
  // and checkboxes/radio buttons/file upload, respectively
  node.on('change', '.selector, .optional_bool, input', handler);
  // For add buttons
  node.on('click', '.add', handler);
}

/**
 * Generic validator for many message types, not just reaction.
 * NOTE: This function does not commit or save anything!
 * @param {!jspb.Message} message The proto to validate.
 * @param {string} messageTypeString The message type.
 * @param {!Node} validateNode Target div for validation output.
 */
function validate(message, messageTypeString, validateNode) {
  // eg message is a type of reaction, messageTypeString = 'Reaction'
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/validate/' + messageTypeString);
  const binary = message.serializeBinary();

  xhr.responseType = 'json';
  xhr.onload = function() {
    const errors = xhr.response;
    const statusNode = $('.validate_status', validateNode);
    const messageNode = $('.validate_message', validateNode);
    statusNode.removeClass('fa-check');
    statusNode.removeClass('fa-exclamation-triangle');
    statusNode.css('backgroundColor', null);
    statusNode.text('');
    if (errors.length) {
      statusNode.addClass('fa fa-exclamation-triangle');
      statusNode.css('color', 'red');
      statusNode.text(' ' + errors.length);

      messageNode.empty();
      for (let index = 0; index < errors.length; index++) {
        const error = errors[index];
        const errorNode = $('<div></div>');
        errorNode.text('\u2022 ' + error);
        messageNode.append(errorNode);
      }
      messageNode.css('backgroundColor', 'pink');
    } else {
      statusNode.addClass('fa fa-check');
      statusNode.css('color', 'green');

      messageNode.html('');
      messageNode.css('backgroundColor', '');
      messageNode.css('visibility', 'hidden');
    }
  };
  xhr.send(binary);
}

/**
 * Toggles the visibility of the 'validate' button for a given node.
 * @param {!Node} node
 */
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
}

/**
 * Updates the visual summary of the current reaction.
 * @param {!proto.ord.Reaction} reaction
 */
function renderReaction(reaction) {
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/render/reaction');
  const binary = reaction.serializeBinary();
  xhr.responseType = 'json';
  xhr.onload = function() {
    const html_block = xhr.response;
    $('#reaction_render').html(html_block);
  };
  xhr.send(binary);
}

/**
 * Validates the current reaction.
 */
function validateReaction() {
  const validateNode = $('#reaction_validate');
  const reaction = unloadReaction();
  validate(reaction, 'Reaction', validateNode);
  // Trigger all submessages to validate.
  $('.validate_button:visible:not(#reaction_validate_button)').trigger('click');
  // Render reaction as an HTML block.
  renderReaction(reaction);
}

/**
 * Writes the current reaction to disk.
 */
function commit() {
  const reaction = unloadReaction();
  const reactions = session.dataset.getReactionsList();
  reactions[session.index] = reaction;
  putDataset(session.fileName, session.dataset);
  ord.uploads.putAll(session.fileName);
}

/**
 * Downloads the current reaction as a serialized Reaction proto.
 */
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

/**
 * Downloads a dataset as a serialized Dataset proto.
 * @param {string} fileName The name of the dataset to fetch.
 * @return {!Uint8Array}
 */
async function getDataset(fileName) {
  return new Promise(resolve => {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', '/dataset/proto/read/' + fileName);
    xhr.responseType = 'arraybuffer';
    xhr.onload = function(event) {
      const bytes = new Uint8Array(xhr.response);
      const dataset = proto.ord.Dataset.deserializeBinary(bytes);
      resolve(dataset);
    };
    xhr.send();
  });
}

/**
 * Uploads a serialized Dataset proto.
 * @param {string} fileName The name of the new dataset.
 * @param {!proto.ord.Dataset} dataset
 */
function putDataset(fileName, dataset) {
  $('#save').text('saving');
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/write/' + fileName);
  const binary = dataset.serializeBinary();
  xhr.onload = clean;
  xhr.send(binary);
}

/**
 * Compares a local Dataset to a Dataset on the server (used for testing).
 * @param {string} fileName The name of a dataset on the server.
 * @param {!proto.ord.Dataset} dataset A local Dataset.
 * @return {!Uint8Array}
 */
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

/**
 * Populates the form with the given reaction.
 * @param {!proto.ord.Reaction} reaction
 */
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

/**
 * Fetches the current reaction from the form.
 * @return {!proto.ord.Reaction}
 */
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

/**
 * Checks if the argument represents an empty protobuf message (that is, the
 * argument's nested arrays only contains null or empty values), or is null or
 * undefined. We use this check on both primitives and arrays/messages.
 * NOTE: Unlike other primitive types, using a setter to set a oneof string
 * field to “” causes the message to include the field and “”, which would be
 * unwanted. So we instead claim that empty strings are empty messages. (Hence
 * we don’t set _any_ empty string.)
 * NOTE: In a submessage, setting a meaningful value (e.g. optional float to 0)
 * will result in a non-null/undefined value in the submessage array. So, if
 * the array of a submessage only contains null and undefined vals, we can
 * assume that the message is truly “empty” (that is, doesn’t have anything
 * meaningful that is set) and can be omitted when constructing the surrounding
 * message.
 * @param {?jspb.Message|?Array|?Uint8Array} obj The object to test.
 * @return {boolean} Whether the message is empty.
 */
function isEmptyMessage(obj) {
  if ([null, undefined, ''].includes(obj)) {
    return true;
  }
  const array = obj.array;
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

/**
 * Adds an instance of `template` to the root node.
 * @param {string} template A jQuery selector.
 * @param {!Node} root A jQuery object.
 * @return {!Node}
 */
function addSlowly(template, root) {
  const node = $(template).clone();
  node.removeAttr('id');
  $(root).append(node);
  node.show('slow');
  dirty();
  listen(node);
  $('[data-toggle=\'tooltip\']', node).tooltip();
  return node;
}

/**
 * Removes from the DOM the nearest ancestor element matching the pattern.
 * @param {string} button The element from which to start the search.
 * @param {string} pattern The pattern for the element to remove.
 */
function removeSlowly(button, pattern) {
  const node = $(button).closest(pattern);
  // Must call necessary validators only after the node is removed,
  // but we can only figure out which validators these are before removal.
  // We do so, and after removal, click the corresponding buttons to trigger
  // validation.
  let buttonsToClick = $();
  node.parents('fieldset').each(function() {
    buttonsToClick =
        buttonsToClick.add($(this).children('legend').find('.validate_button'));
  });
  node.hide('slow', function() {
    node.remove();
    buttonsToClick.trigger('click');
  });
  dirty();
}

/**
 * Toggles the visibility of all siblings of an element, or if a pattern is
 * provided, toggles the visibility of all siblings of the nearest ancestor
 * element matching the pattern.
 * @param {!Node} node The element to toggle or use as the search root.
 * @param {string} pattern The pattern to match for finding siblings to toggle.
 */
//
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
  } else {
    // Need to collapse.
    node.siblings().wrapAll('<collapsed>');
    node.next('collapsed').toggle('slow');
  }
}

/**
 * Toggles the collapse of a section in the form.
 * @param {string} button The element to toggle.
 */
function collapseToggle(button) {
  $(button).toggleClass('fa-chevron-down fa-chevron-right');
  toggleSlowly(button, 'legend');
}

/**
 * Unpacks a (value, units, precision) tuple into the given type.
 * @param {string} prefix The prefix for element attributes.
 * @param {!jspb.Message} proto A protocol buffer with `value`, `precision`,
 *     and `units` fields.
 * @param {!Node} node The node containing the tuple.
 * @return {!jspb.Message} The updated protocol buffer. Note that the message
 *     is modified in-place.
 */
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

/**
 * Packs a (value, units, precision) tuple into form elements.
 * @param {string} prefix The prefix for element attributes.
 * @param {!jspb.Message} proto A protocol buffer with `value`, `precision`,
 *     and`units` fields.
 * @param {!Node} node The target node for the tuple.
 */
function writeMetric(prefix, proto, node) {
  if (!(proto)) {
    return;
  }
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

/**
 * Prompts the user to upload a file and sets the target node text with its
 * contents.
 * @param {!Node} identifierNode The node to update with the file contents.
 * @param {string} valueClass The class containing `identifierNode`.
 */
function setTextFromFile(identifierNode, valueClass) {
  let input = document.createElement('input');
  input.type = 'file';
  input.onchange = (event => {
    const file = event.target.files[0];
    let reader = new FileReader();
    reader.readAsText(file);
    reader.onload = readerEvent => {
      const contents = readerEvent.target.result;
      $('.' + valueClass, identifierNode).text(contents);
    };
  });
  input.click();
}

/**
 * Populates a <select/> node according to its data-proto type declaration.
 * @param {!Node} node A node containing a `data-proto` attribute.
 */
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

/**
 * Selects an <option/> under a <select/>.
 * @param {!Node} node A <select/> element.
 * @param {number} value
 */
function setSelector(node, value) {
  $('option', node).removeAttr('selected');
  $('option[value=' + value + ']', node).attr('selected', 'selected');
}

/**
 * Finds the selected <option/> and maps its text onto a proto Enum.
 * @param {!Node} node A <select/> element.
 * @return {number}
 */
function getSelector(node) {
  return parseInt($('select', node).val());
}

/**
 * Finds the selected <option/> and returns its text.
 * @param {!Node} node A node containing one or more <select/> elements.
 * @return {string}
 */
function getSelectorText(node) {
  const selectorElement = node.getElementsByTagName('select')[0];
  return selectorElement.options[selectorElement.selectedIndex].text;
}

/**
 * Sets up a three-way popup (true/false/unspecified).
 * @param {!Node} node Target node for the new <select/> element.
 */
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

/**
 * Sets the value of a three-way popup (true/false/unspecified).
 * @param {!Node} node A node containing a three-way selector.
 * @param {boolean|null} value The value to select.
 */
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

/**
 * Fetches the value of a three-way popup (true/false/unspecified).
 * @param {!Node} node A node containing a three-way selector.
 * @return {boolean|null}
 */
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

/**
 * Sets up and initializes a collapse button by adding attributes into a div in
 * reaction.html.
 * @param {!Node} node Target node for the new button.
 */
function initCollapse(node) {
  node.addClass('fa');
  node.addClass('fa-chevron-down');
  node.attr('onclick', 'ord.reaction.collapseToggle(this)');
  if (node.hasClass('starts_collapsed')) {
    node.trigger('click');
  }
}

/**
 * Sets up a validator div (button, status indicator, error list, etc.) by
 * inserting contents into a div in reaction.html.
 * @param {!Node} oldNode Target node for the new validation elements.
 */
function initValidateNode(oldNode) {
  let newNode = $('#validate_template').clone();
  // Add attributes necessary for validation functions:
  // Convert the placeholder onclick method into the button's onclick method.
  $('.validate_button', newNode).attr('onclick', oldNode.attr('onclick'));
  oldNode.removeAttr('onclick');
  // Add an id to the button.
  if (oldNode.attr('id')) {
    $('.validate_button', newNode).attr('id', oldNode.attr('id') + '_button');
  }
  oldNode.append(newNode.children());
}

/**
 * Initializes the validation handlers. Some nodes are dynamically added or
 * removed; we add their validation handlers when the nodes themselves are
 * added. However, other nodes are always present in the HTML, and aren't
 * dynamically added nor removed. We add live validation to these nodes here.
 */
function initValidateHandlers() {
  // For setup
  let setupNode = $('#section_setup');
  addChangeHandler(setupNode, () => {
    ord.setups.validateSetup(setupNode);
  });

  // For conditions
  let conditionNode = $('#section_conditions');
  addChangeHandler(conditionNode, () => {
    ord.conditions.validateConditions(conditionNode);
  });

  // For temperature
  let temperatureNode = $('#section_conditions_temperature');
  addChangeHandler(temperatureNode, () => {
    ord.temperature.validateTemperature(temperatureNode);
  });

  // For pressure
  let pressureNode = $('#section_conditions_pressure');
  addChangeHandler(pressureNode, () => {
    ord.pressure.validatePressure(pressureNode);
  });

  // For stirring
  let stirringNode = $('#section_conditions_stirring');
  addChangeHandler(stirringNode, () => {
    ord.stirring.validateStirring(stirringNode);
  });

  // For illumination
  let illuminationNode = $('#section_conditions_illumination');
  addChangeHandler(illuminationNode, () => {
    ord.illumination.validateIllumination(illuminationNode);
  });

  // For electro
  let electroNode = $('#section_conditions_electro');
  addChangeHandler(electroNode, () => {
    ord.electro.validateElectro(electroNode);
  });

  // For flow
  let flowNode = $('#section_conditions_flow');
  addChangeHandler(flowNode, () => {
    ord.flows.validateFlow(flowNode);
  });

  // For notes
  let notesNode = $('#section_notes');
  addChangeHandler(notesNode, () => {
    ord.notes.validateNotes(notesNode);
  });

  // For provenance
  let provenanceNode = $('#section_provenance');
  addChangeHandler(provenanceNode, () => {
    ord.provenance.validateProvenance(provenanceNode);
  });
}

/**
 * Converts a Message_Field name from a data-proto attribute into a proto class.
 * @param {string} protoName Underscore-delimited protocol buffer field name,
 *     such as Reaction_provenance.
 * @return {!typeof jspb.Message|null}
 */
function nameToProto(protoName) {
  let clazz = proto.ord;
  protoName.split('_').forEach(function(name) {
    clazz = clazz[name];
    if (!clazz) {
      return null;
    }
  });
  return clazz;
}

/**
 * Converts an Enum protobuf class to an arrray of strings.
 * @param {!enum} protoEnum
 * @return {!Array<string>}
 */
function enumToStrings(protoEnum) {
  const types = Object.entries(protoEnum);
  const strings = [];
  for (let i = 0; i < types.length; i++) {
    strings.push(types[i][0]);
  }
  return strings;
}

/**
 * Converts an Enum name string to its protobuf member value.
 * @param {string} name A text representation of an enum member.
 * @param {!enum} protoEnum The protocol buffer enum to search.
 * @return {number}
 */
function stringToEnum(name, protoEnum) {
  return protoEnum[name];
}
