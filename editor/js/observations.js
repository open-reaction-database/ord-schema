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

goog.module('ord.observations');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  add,
  validateObservation
};

goog.require('proto.ord.ReactionObservation');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Populates the reaction observation sections in the form.
 * @param {!Array<!proto.ord.ReactionObservation>} observations
 */
function load(observations) {
  observations.forEach(observation => loadObservation(observation));
}

/**
 * Populates a single reaction observation section in the form.
 * @param {!proto.ord.ReactionObservation} observation
 */
function loadObservation(observation) {
  const node = add();
  ord.reaction.writeMetric('.observation_time', observation.getTime(), node);

  $('.observation_comment', node).text(observation.getComment());

  const image = observation.getImage();
  $('.observation_image_description', node).text(image.getDescription());
  $('.observation_image_format', node).text(image.getFormat());

  const stringValue = image.getStringValue();
  const floatValue = image.getFloatValue();
  const bytesValue = image.getBytesValue();
  const url = image.getUrl();
  if (stringValue) {
    $('.observation_image_text', node).show();
    $('.uploader', node).hide();
    $('.observation_image_text', node).text(stringValue);
    $('input[value=\'text\']', node).prop('checked', true);
  }
  if (floatValue) {
    $('.observation_image_text', node).show();
    $('.uploader', node).hide();
    $('.observation_image_text', node).text(floatValue);
    $('input[value=\'number\']', node).prop('checked', true);
  }
  if (bytesValue) {
    $('.observation_image_text', node).hide();
    $('.uploader', node).show();
    ord.uploads.load(node, bytesValue);
    $('input[value=\'upload\']', node).prop('checked', true);
  }
  if (url) {
    $('.observation_image_text', node).show();
    $('.uploader', node).hide();
    $('.observation_image_text', node).text(url);
    $('input[value=\'url\']', node).prop('checked', true);
  }
}

/**
 * Fetches the reaction observations defined in the form.
 * @return {!Array<!proto.ord.ReactionObservation>}
 */
function unload() {
  const observations = [];
  $('.observation').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template
      const observation = unloadObservation(node);
      if (!ord.reaction.isEmptyMessage(observation)) {
        observations.push(observation);
      }
    }
  });
  return observations;
}

/**
 * Fetches a single reaction observation defined in the form.
 * @param {!Node} node Root node for the reaction observation.
 * @return {!proto.ord.ReactionObservation}
 */
function unloadObservation(node) {
  const observation = new proto.ord.ReactionObservation();
  const time =
      ord.reaction.readMetric('.observation_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(time)) {
    observation.setTime(time);
  }

  observation.setComment($('.observation_comment', node).text());

  const image = new proto.ord.Data();
  image.setDescription($('.observation_image_description', node).text());
  image.setFormat($('.observation_image_format', node).text());

  if ($('input[value=\'text\']', node).is(':checked')) {
    const stringValue = $('.observation_image_text', node).text();
    if (!ord.reaction.isEmptyMessage(stringValue)) {
      image.setStringValue(stringValue);
    }
  }
  if ($('input[value=\'number\']', node).is(':checked')) {
    const floatValue = parseFloat($('.setup_code_text', node).text());
    if (!isNaN(floatValue)) {
      image.setFloatValue(floatValue);
    }
  }
  if ($('input[value=\'upload\']', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!ord.reaction.isEmptyMessage(bytesValue)) {
      image.setBytesValue(bytesValue);
    }
  }
  if ($('input[value=\'url\']', node).is(':checked')) {
    const url = $('.observation_image_text', node).text();
    if (!ord.reaction.isEmptyMessage(url)) {
      image.setUrl(url);
    }
  }
  if (!ord.reaction.isEmptyMessage(image)) {
    observation.setImage(image);
  }
  return observation;
}

/**
 * Adds a reaction observation section to the form.
 * @return {!Node} The newly added parent node for the reaction observation.
 */
function add() {
  const node = ord.reaction.addSlowly('#observation_template', '#observations');

  const typeButtons = $('input[type=\'radio\']', node);
  typeButtons.attr('name', 'observations_' + radioGroupCounter++);
  typeButtons.change(function() {
    if ((this.value == 'text') || (this.value == 'number') ||
        (this.value == 'url')) {
      $('.observation_image_text', node).show();
      $('.uploader', node).hide();
    } else {
      $('.observation_image_text', node).hide();
      $('.uploader', node).show();
    }
  });
  ord.uploads.initialize(node);

  // Add live validation handling.
  ord.reaction.addChangeHandler(node, () => {
    validateObservation(node);
  });

  return node;
}

/**
 * Validates a single reaction observation defined in the form.
 * @param {!Node} node Root node for the reaction observation.
 * @param {!Node} validateNode Target node for validation results.
 */
function validateObservation(node, validateNode) {
  const observation = unloadObservation(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(observation, 'ReactionObservation', validateNode);
}
