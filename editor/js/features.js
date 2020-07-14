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

goog.provide('ord.features');

goog.require('proto.ord.Compound.Feature');

ord.features.load = function (node, features) {
  features.forEach(feature => ord.features.loadFeature(node, feature));
};

ord.features.loadFeature = function (compoundNode, feature) {
  const node = ord.features.add(compoundNode);
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
};

ord.features.unload = function (compoundNode) {
  const features = [];
  $('.component_feature', compoundNode).each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const feature = ord.features.unloadFeature(node);
      if (!isEmptyMessage(feature)) {
        features.push(feature);
      }
    }
  });
  return features;
};

ord.features.unloadFeature = function (node) {
  const feature = new proto.ord.Compound.Feature();
  const name = $('.component_feature_name', node).text();
  feature.setName(name);

  const valueText = $('.component_feature_value', node).text();
  const valueFloat = parseFloat(valueText);
  if (isNaN(valueFloat)) {
    if (!isEmptyMessage(valueText)) {
      feature.setStringValue(valueText);
    }
  } else {
    if (!isEmptyMessage(valueFloat)) {
      feature.setFloatValue(valueFloat);
    }
  }
  const how = $('.component_feature_how', node).text();
  feature.setHowComputed(how);
  return feature;
};

ord.features.add = function (compoundNode) {
  return addSlowly('#component_feature_template', $('.features', compoundNode));
};
