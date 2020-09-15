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

goog.module('ord.codes');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  addCode
};

goog.require('ord.data');
goog.require('proto.ord.Data');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the automation_code sections in the form.
 * @param {!jspb.Map<string, !proto.ord.Data>} codes
 */
function load(codes) {
  const names = codes.stringKeys_();
  names.forEach(function(name) {
    const code = codes.get(name);
    loadCode(name, code);
  });
}

/**
 * Adds and populates a single automation_code section in the form.
 * @param {string} name The name of this automation code.
 * @param {!proto.ord.Data} code
 */
function loadCode(name, code) {
  const node = addCode();
  $('.setup_code_name', node).text(name);
  ord.data.loadData(node, code);
}

/**
 * Fetches the automation_code sections from the form and adds them to `codes`.
 * @param {!jspb.Map<string, !proto.ord.Data>} codes
 */
function unload(codes) {
  $('.setup_code').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      unloadCode(codes, node);
    }
  });
}

/**
 * Fetches a single automation_code section from the form and adds it to
 * `codes`.
 * @param {!jspb.Map<string, !proto.ord.Data>} codes
 * @param {!Node} node The root node of the automation_code section to fetch.
 */
function unloadCode(codes, node) {
  const name = $('.setup_code_name', node).text();
  const code = ord.data.unloadData(node);
  if (!ord.reaction.isEmptyMessage(name) ||
      !ord.reaction.isEmptyMessage(code)) {
    codes.set(name, code);
  }
}

/**
 * Adds an automation_code section to the form.
 * @return {!Node} The newly added root node for the automation_code section.
 */
function addCode() {
  const node = ord.reaction.addSlowly('#setup_code_template', '#setup_codes');
  ord.data.addData(node);
  return node;
}
