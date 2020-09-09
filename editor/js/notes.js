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

goog.module('ord.notes');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  validateNotes
};

goog.require('proto.ord.ReactionNotes');

/**
 * Populates the reaction nodes section in the form.
 * @param {!proto.ord.ReactionNotes} notes
 */
function load(notes) {
  ord.reaction.setOptionalBool(
      $('#notes_heterogeneous'),
      notes.hasIsHeterogeneous() ? notes.getIsHeterogeneous() : null);
  ord.reaction.setOptionalBool(
      $('#notes_precipitate'),
      notes.hasFormsPrecipitate() ? notes.getFormsPrecipitate() : null);
  ord.reaction.setOptionalBool(
      $('#notes_exothermic'),
      notes.hasIsExothermic() ? notes.getIsExothermic() : null);
  ord.reaction.setOptionalBool(
      $('#notes_offgas'), notes.hasOffgasses() ? notes.getOffgasses() : null);
  ord.reaction.setOptionalBool(
      $('#notes_moisture'),
      notes.hasIsSensitiveToMoisture() ? notes.getIsSensitiveToMoisture() :
                                         null);
  ord.reaction.setOptionalBool(
      $('#notes_oxygen'),
      notes.hasIsSensitiveToOxygen() ? notes.getIsSensitiveToOxygen() : null);
  ord.reaction.setOptionalBool(
      $('#notes_light'),
      notes.hasIsSensitiveToLight() ? notes.getIsSensitiveToLight() : null);
  $('#notes_safety').text(notes.getSafetyNotes());
  $('#notes_details').text(notes.getProcedureDetails());
}

/**
 * Fetches the reaction notes defined in the form.
 * @return {!proto.ord.ReactionNotes}
 */
function unload() {
  const notes = new proto.ord.ReactionNotes();
  notes.setIsHeterogeneous(
      ord.reaction.getOptionalBool($('#notes_heterogeneous')));
  notes.setFormsPrecipitate(
      ord.reaction.getOptionalBool($('#notes_precipitate')));
  notes.setIsExothermic(ord.reaction.getOptionalBool($('#notes_exothermic')));
  notes.setOffgasses(ord.reaction.getOptionalBool($('#notes_offgas')));
  notes.setIsSensitiveToMoisture(
      ord.reaction.getOptionalBool($('#notes_moisture')));
  notes.setIsSensitiveToOxygen(
      ord.reaction.getOptionalBool($('#notes_oxygen')));
  notes.setIsSensitiveToLight(ord.reaction.getOptionalBool($('#notes_light')));
  notes.setSafetyNotes($('#notes_safety').text());
  notes.setProcedureDetails($('#notes_details').text());
  return notes;
}

/**
 * Validates the reaction notes defined in the form.
 * @param {!Node} node Root node for the reaction notes.
 * @param {?Node} validateNode Target node for validation results.
 */
function validateNotes(node, validateNode) {
  const notes = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(notes, 'ReactionNotes', validateNode);
}
