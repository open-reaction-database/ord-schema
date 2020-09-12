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

goog.module('ord.data');
goog.module.declareLegacyNamespace();
exports = {
  addData,
  loadData,
  unloadData,
};

goog.require('proto.ord.Data');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds a new Data section to the form.
 * @param {!Node} parentNode Parent node.
 * @return {!Node} The newly added node for the Data record.
 */
function addData(parentNode) {
  const target = parentNode.children('fieldset').first();
  const node = ord.reaction.addSlowly('#data_template', target);
  const typeButtons = $('input[type=\'radio\']', node);
  typeButtons.attr('name', 'data_' + radioGroupCounter++);
  typeButtons.change(function() {
    if ((this.value === 'text') || (this.value === 'number') ||
        (this.value === 'url')) {
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
    } else {
      $('.data_text', node).hide();
      $('.data_uploader', node).show();
    }
  });
  ord.uploads.initialize(node);
  return node;
}

/**
 * Populates an existing Data section in the form.
 * @param {!Node} node Root node.
 * @param {!proto.ord.Data} data
 */
function loadData(node, data) {
  $('.data_description', node).text(data.getDescription());
  $('.data_format', node).text(data.getFormat());
  let value;
  switch (data.getKindCase()) {
    case proto.ord.Data.KindCase.FLOAT_VALUE:
      value = data.getFloatValue().toString();
      if (value.indexOf('.') === -1) {
        value = value.concat('.');
      }
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'number\']', node).prop('checked', true);
      break;
    case proto.ord.Data.KindCase.INTEGER_VALUE:
      value = data.getIntegerValue();
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'number\']', node).prop('checked', true);
      break;
    case proto.ord.Data.KindCase.BYTES_VALUE:
      value = data.getBytesValue();
      $('.data_text', node).hide();
      $('.data_uploader', node).show();
      ord.uploads.load(node, value);
      $('input[value=\'upload\']', node).prop('checked', true);
      break;
    case proto.ord.Data.KindCase.STRING_VALUE:
      value = data.getStringValue();
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'text\']', node).prop('checked', true);
      break;
    case proto.ord.Data.KindCase.URL:
      value = data.getUrl();
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'url\']', node).prop('checked', true);
      break;
    default:
      break;
  }
}

/**
 * Fetches a Data section from the form.
 * @param {!Node} node Root node of the Data section to fetch.
 * @return {!proto.ord.Data}
 */
function unloadData(node) {
  const data = new proto.ord.Data();
  const description = $('.data_description', node).text();
  data.setDescription(description);
  const format = $('.data_format', node).text();
  data.setFormat(format);
  if ($('input[value=\'text\']', node).is(':checked')) {
    const stringValue = $('.data_text', node).text();
    if (!ord.reaction.isEmptyMessage(stringValue)) {
      data.setStringValue(stringValue);
    }
  }
  if ($('input[value=\'number\']', node).is(':checked')) {
    const stringValue = $('.setup_code_text', node).text();
    const value = parseFloat(stringValue);
    if (Number.isInteger(value) && stringValue.indexOf('.') === -1) {
      data.setIntegerValue(value);
    } else if (!Number.isNaN(value)) {
      data.setFloatValue(value);
    }
  }
  if ($('input[value=\'upload\']', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!ord.reaction.isEmptyMessage(bytesValue)) {
      data.setBytesValue(bytesValue);
    }
  }
  if ($('input[value=\'url\']', node).is(':checked')) {
    const url = $('.setup_code_text', node).text();
    if (!ord.reaction.isEmptyMessage(url)) {
      data.setUrl(url);
    }
  }
  return data;
}
