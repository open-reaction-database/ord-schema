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

goog.require('ord.data');
goog.require('proto.ord.ReactionObservation');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the reaction observation sections in the form.
 * @param {!Array<!proto.ord.ReactionObservation>} observations
 */
function load(observations) {
  observations.forEach(observation => loadObservation(observation));
}

/**
 * Adds and populates a single reaction observation section in the form.
 * @param {!proto.ord.ReactionObservation} observation
 */
function loadObservation(observation) {
  const node = add();
  ord.reaction.writeMetric('.observation_time', observation.getTime(), node);
  $('.observation_comment', node).text(observation.getComment());
  ord.data.loadData(node, observation.getImage());
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
  const image = ord.data.unloadData(node);
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
  ord.data.addData(node);
  // Add live validation handling.
  ord.reaction.addChangeHandler(node, () => {
    validateObservation(node);
  });
  return node;
}

/**
 * Validates a single reaction observation defined in the form.
 * @param {!Node} node Root node for the reaction observation.
 * @param {?Node} validateNode Target node for validation results.
 */
function validateObservation(node, validateNode) {
  const observation = unloadObservation(node);
  ord.reaction.validate(observation, 'ReactionObservation', node, validateNode);
}
