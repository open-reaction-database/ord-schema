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

goog.module('ord.features');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  add
};

goog.require('proto.ord.Compound.Feature');

/**
 * Adds and populates the feature sections for a Compound.
 * @param {!Node} node Root node for the parent Compound.
 * @param {!Array<!proto.ord.Compound.Feature>} features
 */
function load(node, features) {
  features.forEach(feature => loadFeature(node, feature));
}

/**
 * Adds and populates a single feature section in the form.
 * @param {!Node} compoundNode Root node for the parent Compound.
 * @param {!proto.ord.Compound.Feature} feature
 */
function loadFeature(compoundNode, feature) {
  const node = add(compoundNode);
  const name = feature.getName();
  $('.component_feature_name', node).text(name);
  const valueText = feature.getStringValue();
  const valueFloat = feature.getFloatValue();
  if (valueText) {
    $('.component_feature_value', node).text(valueText);
  }
  if (valueFloat) {
    $('.component_feature_value', node).text(valueFloat);
  }
  const how = feature.getHowComputed();
  $('.component_feature_how', node).text(how);
}

/**
 * Fetches the features for a compound from the form.
 * @param {!Node} compoundNode Root node for the parent Compound.
 * @return {!Array<!proto.ord.Compound.Feature>}
 */
function unload(compoundNode) {
  const features = [];
  $('.component_feature', compoundNode).each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const feature = unloadFeature(node);
      if (!ord.reaction.isEmptyMessage(feature)) {
        features.push(feature);
      }
    }
  });
  return features;
}

/**
 * Fetches a single feature defined in the form.
 * @param {!Node} node Root node for the compound feature.
 * @return {!proto.ord.Compound.Feature}
 */
function unloadFeature(node) {
  const feature = new proto.ord.Compound.Feature();
  const name = $('.component_feature_name', node).text();
  feature.setName(name);

  const valueText = $('.component_feature_value', node).text();
  const valueFloat = parseFloat(valueText);
  if (isNaN(valueFloat)) {
    if (!ord.reaction.isEmptyMessage(valueText)) {
      feature.setStringValue(valueText);
    }
  } else {
    if (!ord.reaction.isEmptyMessage(valueFloat)) {
      feature.setFloatValue(valueFloat);
    }
  }
  const how = $('.component_feature_how', node).text();
  feature.setHowComputed(how);
  return feature;
}

/**
 * Adds a feature section to the given compound in the form.
 * @param {!Node} compoundNode Root node for the parent Compound.
 * @return {!Node} The newly added parent node for the compound feature.
 */
function add(compoundNode) {
  return ord.reaction.addSlowly(
      '#component_feature_template', $('.features', compoundNode));
}
