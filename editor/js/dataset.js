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

goog.module('ord.dataset');
goog.module.declareLegacyNamespace();
exports = {
  init,
  download,
  commit,
  deleteReaction,
  newReaction,
  removeReactionId,
  addReactionId,
  freeze
};

goog.require('proto.ord.Dataset');

const session = {
  fileName: null,
  dataset: null
};

/**
 * Initializes the dataset landing page.
 * @param {string} fileName The filename of the dataset to load.
 */
function init(fileName) {
  session.fileName = fileName;
  $('.edittext').attr('contentEditable', 'true');
  getDataset(fileName, loadDataset);
  listenDirty($('#text_fields'));
}

/**
 * Adds change handlers to a newly added reaction or reaction ID node.
 * @param {!Node} node Root node for the reaction or reaction ID.
 */
function listenDirty(node) {
  $('.edittext', node).on('input', dirty);
  $('.selector', node).on('input', dirty);
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
 * Writes the current dataset to disk.
 */
function commit() {
  const dataset = unloadDataset();
  $('#save').text('saving');
  const xhr = new XMLHttpRequest();
  xhr.open(
      'POST', '/dataset/proto/write/' + session.fileName, true /* async */);
  const binary = dataset.serializeBinary();
  xhr.onload = clean;
  xhr.send(binary);
}

/**
 * Downloads the current dataset.
 */
function download() {
  const xhr = new XMLHttpRequest();
  xhr.open('GET', '/dataset/' + session.fileName + '/download');
  xhr.onload = () => {
    // Make the browser write the file.
    const url = URL.createObjectURL(new Blob([xhr.response]));
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', session.fileName + '.pbtxt');
    document.body.appendChild(link);
    link.click();
  };
  xhr.send();
}

/**
 * Fetches a dataset from the server.
 * @param {string} fileName The filename of the dataset to fetch.
 * @param {?Function} listener Function used to load the dataset into the
 *     editor.
 */
function getDataset(fileName, listener) {
  if (!listener) {
    return;
  }
  const xhr = new XMLHttpRequest();
  xhr.open('GET', '/dataset/proto/read/' + session.fileName, true /* async */);
  xhr.responseType = 'arraybuffer';
  xhr.onload = () => {
    const bytes = new Uint8Array(xhr.response);
    const dataset = proto.ord.Dataset.deserializeBinary(bytes);
    session.dataset = dataset;
    listener(dataset);
  };
  xhr.send();
}

/**
 * Loads a dataset into the editor.
 * @param {!proto.ord.Dataset} dataset
 */
function loadDataset(dataset) {
  $('#name').text(dataset.getName());
  $('#description').text(dataset.getDescription());
  $('#dataset_id').text(dataset.getDatasetId());

  const reactions = dataset.getReactionsList();
  loadReactions(reactions);

  const reactionIds = dataset.getReactionIdsList();
  loadReactionIds(reactionIds);

  clean();
}

/**
 * Loads a list of reactions into the editor.
 * @param {!Array<!proto.ord.Reaction>} reactions
 */
function loadReactions(reactions) {
  for (var i = 0; i < reactions.length; i++) {
    const reaction = reactions[i];
    loadReaction(i, reaction);
  }
}

/**
 * Loads a single reaction into the editor.
 * @param {number} index The index of the new reaction.
 * @param {!proto.ord.Reaction} reaction
 */
function loadReaction(index, reaction) {
  const node = addReaction(index);
  const id = reaction.getReactionId();
  $('.reaction_id', node).text(id);
}

/**
 * Loads a list of reaction IDs into the editor.
 * @param {!Array<string>} reactionIds
 */
function loadReactionIds(reactionIds) {
  reactionIds.forEach(reactionId => loadReactionId(reactionId));
}

/**
 * Loads a single reaction ID into the editor.
 * @param {string} reactionId
 */
function loadReactionId(reactionId) {
  const node = addReactionId();
  $('.other_reaction_id_text', node).text(reactionId);
}

/**
 * Fetches the current dataset.
 * @return {!proto.ord.Dataset}
 */
function unloadDataset() {
  const dataset = session.dataset;
  dataset.setName($('#name').text());
  dataset.setDescription($('#description').text());
  dataset.setDatasetId($('#dataset_id').text());
  const reactionIds = [];
  $('.other_reaction_id').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      reactionIds.push($('.other_reaction_id_text', node).text());
    }
  });
  dataset.setReactionIdsList(reactionIds);
  // Do not mutate Reactions. They are edited separately.
  return dataset;
}

/**
 * Adds a new reaction to the current dataset.
 * @param {number} index The index of the new reaction.
 * @return {!Node} The newly added root node for the reaction.
 */
function addReaction(index) {
  const node = $('#reaction_template').clone();
  node.removeAttr('id');
  const anchor = $('.reaction_index', node);
  anchor.text(index);
  anchor.attr('href', '/dataset/' + session.fileName + '/reaction/' + index);
  const root = $('#reactions');
  root.append(node);
  node.show('slow');
  listenDirty(node);
  dirty();
  return node;
}

/**
 * Adds a new reaction ID to the current dataset.
 * @return {!Node} The newly added root node for the reaction ID.
 */
function addReactionId() {
  const node = $('#other_reaction_id_template').clone();
  node.removeAttr('id');
  const root = $('#other_reaction_ids');
  root.append(node);
  node.show('slow');
  listenDirty(node);
  dirty();
  return node;
}

/**
 * Loads the Reaction editor immediately without waiting for 'save'.
 */
function newReaction() {
  window.location.href = '/dataset/' + session.fileName + '/new/reaction';
}

/**
 * Deletes a Reaction immediately without waiting for 'save'.
 * @param {!Node} button The node of the 'remove' button.
 */
function deleteReaction(button) {
  const node = $(button).closest('.reaction');
  const index = parseInt($('a', node).text());
  window.location.href =
      '/dataset/' + session.fileName + '/delete/reaction/' + index;
}

/**
 * Deletes a Reaction ID.
 * @param {!Node} button The node of the 'remove' button.
 */
function removeReactionId(button) {
  removeSlowly(button, '.other_reaction_id');
}

/**
 * Deletes an element matching `pattern`.
 * @param {!Node} button The node of the 'remove' button.
 * @param {string} pattern The element pattern to match.
 */
function removeSlowly(button, pattern) {
  const node = $(button).closest(pattern);
  node.hide('slow', () => node.remove());
  dirty();
}

/**
 * Switches the UI into a read-only mode. This is irreversible.
 */
function freeze() {
  $('.edittext').attr('contenteditable', 'false');
  $('.remove').hide();
  $('#save').hide();
  $('.edittext').each((i, x) => {
    // Ensure non-editable divs have a text node to preserve vertical alignment.
    if (!$(x).text()) {
      $(x).text(' ');
    }
  });
}
