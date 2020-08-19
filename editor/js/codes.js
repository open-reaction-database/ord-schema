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

goog.provide('ord.codes');

goog.require('proto.ord.Data');

// Freely create radio button groups by generating new input names.
ord.codes.radioGroupCounter = 0;

ord.codes.load = function(codes) {
  const names = codes.stringKeys_();
  names.forEach(function(name) {
    const code = codes.get(name);
    ord.codes.loadCode(name, code);
  });
};

ord.codes.loadCode = function(name, code) {
  const node = ord.codes.addCode();
  $('.setup_code_name', node).text(name);
  $('.setup_code_description', node).text(code.getDescription());
  $('.setup_code_format', node).text(code.getFormat());

  const stringValue = code.getStringValue();
  const floatValue = code.getFloatValue();
  const bytesValue = code.getBytesValue();
  const url = code.getUrl();
  if (stringValue) {
    $('.setup_code_text', node).show();
    $('.uploader', node).hide();
    $('.setup_code_text', node).text(stringValue);
    $('input[value=\'text\']', node).prop('checked', true);
  }
  if (floatValue) {
    $('.setup_code_text', node).show();
    $('.uploader', node).hide();
    $('.setup_code_text', node).text(floatValue);
    $('input[value=\'number\']', node).prop('checked', true);
  }
  if (bytesValue) {
    $('.setup_code_text', node).hide();
    $('.uploader', node).show();
    ord.uploads.load(node, bytesValue);
    $('input[value=\'upload\']', node).prop('checked', true);
  }
  if (url) {
    $('.setup_code_text', node).show();
    $('.uploader', node).hide();
    $('.setup_code_text', node).text(url);
    $('input[value=\'url\']', node).prop('checked', true);
  }
};

ord.codes.unload = function(codes) {
  $('.setup_code').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      ord.codes.unloadCode(codes, node);
    }
  });
};

ord.codes.unloadCode = function(codes, node) {
  const name = $('.setup_code_name', node).text();

  const code = new proto.ord.Data();

  const description = $('.setup_code_description', node).text();
  code.setDescription(description);
  const format = $('.setup_code_format', node).text();
  code.setFormat(format);

  if ($('input[value=\'text\']', node).is(':checked')) {
    const stringValue = $('.setup_code_text', node).text();
    if (!ord.reaction.isEmptyMessage(stringValue)) {
      code.setStringValue(stringValue);
    }
  }
  if ($('input[value=\'number\']', node).is(':checked')) {
    const floatValue = parseFloat($('.setup_code_text', node).text());
    if (!isNaN(floatValue)) {
      code.setFloatValue(floatValue);
    }
  }
  if ($('input[value=\'upload\']', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!ord.reaction.isEmptyMessage(bytesValue)) {
      code.setBytesValue(bytesValue);
    }
  }
  if ($('input[value=\'url\']', node).is(':checked')) {
    const url = $('.setup_code_text', node).text();
    if (!ord.reaction.isEmptyMessage(url)) {
      code.setUrl(url);
    }
  }
  if (!ord.reaction.isEmptyMessage(name) ||
      !ord.reaction.isEmptyMessage(code)) {
    codes.set(name, code);
  }
};

ord.codes.addCode = function() {
  const node = ord.reaction.addSlowly('#setup_code_template', '#setup_codes');

  const typeButtons = $('input[type=\'radio\']', node);
  typeButtons.attr('name', 'codes_' + ord.codes.radioGroupCounter++);
  typeButtons.change(function() {
    if ((this.value == 'text') || (this.value == 'number') ||
        (this.value == 'url')) {
      $('.setup_code_text', node).show();
      $('.uploader', node).hide();
    } else {
      $('.setup_code_text', node).hide();
      $('.uploader', node).show();
    }
  });
  ord.uploads.initialize(node);
  return node;
};
