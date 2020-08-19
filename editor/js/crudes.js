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

goog.provide('ord.crudes');

goog.require('ord.amountsCrudes');
goog.require('proto.ord.CrudeComponent');

// Freely create radio button groups by generating new input names.
ord.crudes.radioGroupCounter = 0;

ord.crudes.load = function(node, crudes) {
  crudes.forEach(crude => ord.crudes.loadCrude(node, crude));
};

ord.crudes.loadCrude = function(root, crude) {
  const node = ord.crudes.add(root);

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
};

ord.crudes.unload = function(node) {
  const crudes = [];
  $('.crude', node).each(function(index, crudeNode) {
    crudeNode = $(crudeNode);
    if (!crudeNode.attr('id')) {
      // Not a template.
      const crude = ord.crudes.unloadCrude(crudeNode);
      if (!ord.reaction.isEmptyMessage(crude)) {
        crudes.push(crude);
      }
    }
  });
  return crudes;
};

ord.crudes.unloadCrude = function(node) {
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
};

ord.crudes.add = function(root) {
  const node = ord.reaction.addSlowly('#crude_template', $('.crudes', root));

  // Create an "amount" radio button group and connect it to the unit selectors.
  const amountButtons = $('.amount input', node);
  amountButtons.attr('name', 'crudes_' + ord.crudes.radioGroupCounter++);
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
};
