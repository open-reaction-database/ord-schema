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

goog.module('ord.crudes');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  add
};

goog.require('ord.amountsCrudes');
goog.require('proto.ord.CrudeComponent');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the crude components of a reaction input.
 * @param {!Node} node Root node for the parent reaction input.
 * @param {!Array<!proto.ord.CrudeComponent>} crudes
 */
function load(node, crudes) {
  crudes.forEach(crude => loadCrude(node, crude));
}

/**
 * Adds and populates a single crude component section in the form.
 * @param {!Node} root Root node for the parent reaction input.
 * @param {!proto.ord.CrudeComponent} crude
 */
function loadCrude(root, crude) {
  const node = add(root);

  const reactionId = crude.getReactionId();
  $('.crude_reaction', node).text(reactionId);

  const workup = crude.hasIncludesWorkup() ? crude.getIncludesWorkup() : null;
  ord.reaction.setOptionalBool($('.crude_includes_workup', node), workup);

  const derived =
      crude.hasHasDerivedAmount() ? crude.getHasDerivedAmount() : null;
  ord.reaction.setOptionalBool($('.crude_has_derived', node), derived);

  const mass = crude.getMass();
  const volume = crude.getVolume();
  ord.amountsCrudes.load(node, mass, volume);
}

/**
 * Fetches the crude components defined for a reaction input in the form.
 * @param {!Node} node Root node for the parent reaction input.
 * @return {!Array<!proto.ord.CrudeComponent>}
 */
function unload(node) {
  const crudes = [];
  $('.crude', node).each(function(index, crudeNode) {
    crudeNode = $(crudeNode);
    if (!crudeNode.attr('id')) {
      // Not a template.
      const crude = unloadCrude(crudeNode);
      if (!ord.reaction.isEmptyMessage(crude)) {
        crudes.push(crude);
      }
    }
  });
  return crudes;
}

/**
 * Fetches a single crude component defined in the form.
 * @param {!Node} node Root node for the crude component.
 * @return {!proto.ord.CrudeComponent}
 */
function unloadCrude(node) {
  const crude = new proto.ord.CrudeComponent();

  const reactionId = $('.crude_reaction', node).text();
  crude.setReactionId(reactionId);

  const workup =
      ord.reaction.getOptionalBool($('.crude_includes_workup', node));
  crude.setIncludesWorkup(workup);

  const derived = ord.reaction.getOptionalBool($('.crude_has_derived', node));
  crude.setHasDerivedAmount(derived);

  ord.amountsCrudes.unload(node, crude);

  return crude;
}

/**
 * Adds a crude component section to the given reaction input.
 * @param {!Node} root Root node for the parent reaction input.
 * @return {!Node} The newly added root node for the crude component.
 */
function add(root) {
  const node = ord.reaction.addSlowly('#crude_template', $('.crudes', root));

  // Create an "amount" radio button group and connect it to the unit selectors.
  const amountButtons = $('.amount input', node);
  amountButtons.attr('name', 'crudes_' + radioGroupCounter++);
  amountButtons.change(function() {
    $('.amount .selector', node).hide();
    if (this.value == 'mass') {
      $('.crude_amount_units_mass', node).show();
    }
    if (this.value == 'moles') {
      $('.crude_amount_units_moles', node).show();
    }
    if (this.value == 'volume') {
      $('.crude_amount_units_volume', node).show();
    }
  });
  return node;
}
