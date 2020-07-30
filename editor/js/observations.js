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

goog.provide('ord.observations');

goog.require('proto.ord.ReactionObservation');

// Freely create radio button groups by generating new input names.
ord.observations.radioGroupCounter = 0;

ord.observations.load = function (observations) {
  observations.forEach(
      observation => ord.observations.loadObservation(observation));
};

ord.observations.loadObservation = function (observation) {
  const node = ord.observations.add();
  writeMetric('.observation_time', observation.getTime(), node);

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
    $("input[value='text']", node).prop('checked', true);
  }
  if (floatValue) {
    $('.observation_image_text', node).show();
    $('.uploader', node).hide();
    $('.observation_image_text', node).text(floatValue);
    $("input[value='number']", node).prop('checked', true);
  }
  if (bytesValue) {
    $('.observation_image_text', node).hide();
    $('.uploader', node).show();
    ord.uploads.load(node, bytesValue);
    $("input[value='upload']", node).prop('checked', true);
  }
  if (url) {
    $('.observation_image_text', node).show();
    $('.uploader', node).hide();
    $('.observation_image_text', node).text(url);
    $("input[value='url']", node).prop('checked', true);
  }
};

ord.observations.unload = function () {
  const observations = [];
  $('.observation').each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template
      const observation = ord.observations.unloadObservation(node);
      if (!isEmptyMessage(observation)) {
        observations.push(observation);
      }
    }
  });
  return observations;
};

ord.observations.unloadObservation = function (node) {
  const observation = new proto.ord.ReactionObservation();
  const time =
      readMetric('.observation_time', new proto.ord.Time(), node);
  if (!isEmptyMessage(time)) {
    observation.setTime(time);
  }

  observation.setComment($('.observation_comment', node).text());

  const image = new proto.ord.Data();
  image.setDescription($('.observation_image_description', node).text());
  image.setFormat($('.observation_image_format', node).text());

  if ($("input[value='text']", node).is(':checked')) {
    const stringValue = $('.observation_image_text', node).text();
    if (!isEmptyMessage(stringValue)) {
      image.setStringValue(stringValue);
    }
  }
  if ($("input[value='number']", node).is(':checked')) {
    const floatValue = parseFloat($('.setup_code_text', node).text());
    if (!isNaN(floatValue)) {
      image.setFloatValue(floatValue);
    }
  }
  if ($("input[value='upload']", node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!isEmptyMessage(bytesValue)) {
      image.setBytesValue(bytesValue);
    }
  }
  if ($("input[value='url']", node).is(':checked')) {
    const url = $('.observation_image_text', node).text();
    if (!isEmptyMessage(url)) {
      image.setUrl(url);
    }
  }
  if (!isEmptyMessage(image)) {
    observation.setImage(image);
  }
  return observation;
};

ord.observations.add = function () {
  const node = addSlowly('#observation_template', '#observations');

  const typeButtons = $("input[type='radio']", node);
  typeButtons.attr(
      'name', 'observations_' + ord.observations.radioGroupCounter++);
  typeButtons.change(function () {
    if ((this.value == 'text')
        || (this.value == 'number')
        || (this.value == 'url')) {
      $('.observation_image_text', node).show();
      $('.uploader', node).hide();
    } else {
      $('.observation_image_text', node).hide();
      $('.uploader', node).show();
    }
  });
  ord.uploads.initialize(node);

  // Add live validation handling.
  addChangeHandler(node, () => {ord.observations.validateObservation(node)});

  return node;
};

ord.observations.validateObservation = function(node, validateNode) {
  const observation = ord.observations.unloadObservation(node);
  if (typeof validateNode === 'undefined') {
    validateNode = $('.validate', node).first();
  }
  validate(observation, 'ReactionObservation', validateNode);
};