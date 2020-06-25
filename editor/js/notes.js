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

goog.provide('ord.notes');

goog.require('proto.ord.ReactionNotes');

ord.notes.load = function (notes) {
  setOptionalBool($('#notes_heterogeneous'), 
      notes.hasIsHeterogeneous() ? notes.getIsHeterogeneous() : null);
  setOptionalBool($('#notes_precipitate'),
      notes.hasFormsPrecipitate() ? notes.getFormsPrecipitate() : null);
  setOptionalBool($('#notes_exothermic'),
      notes.hasIsExothermic() ? notes.getIsExothermic() : null);
  setOptionalBool($('#notes_offgas'),
      notes.hasOffgasses() ? notes.getOffgasses() : null);
  setOptionalBool($('#notes_moisture'),
      notes.hasIsSensitiveToMoisture() ?
          notes.getIsSensitiveToMoisture() : null);
  setOptionalBool($('#notes_oxygen'),
      notes.hasIsSensitiveToOxygen() ? notes.getIsSensitiveToOxygen() : null);
  setOptionalBool($('#notes_light'),
      notes.hasIsSensitiveToLight() ? notes.getIsSensitiveToLight() : null);
  $('#notes_safety').text(notes.getSafetyNotes());
  $('#notes_details').text(notes.getProcedureDetails());
};

ord.notes.unload = function () {
  const notes = new proto.ord.ReactionNotes();
  notes.setIsHeterogeneous(getOptionalBool($('#notes_heterogeneous')));
  notes.setFormsPrecipitate(getOptionalBool($('#notes_precipitate')));
  notes.setIsExothermic(getOptionalBool($('#notes_exothermic')));
  notes.setOffgasses(getOptionalBool($('#notes_offgas')));
  notes.setIsSensitiveToMoisture(getOptionalBool($('#notes_moisture')));
  notes.setIsSensitiveToOxygen(getOptionalBool($('#notes_oxygen')));
  notes.setIsSensitiveToLight(getOptionalBool($('#notes_light')));
  notes.setSafetyNotes($('#notes_safety').text());
  notes.setProcedureDetails($('#notes_details').text());
  return notes;
};

ord.notes.validateNotes = function(node, validateNode) {
    const notes = ord.notes.unload();
    if (typeof validateNode === 'undefined') {
      validateNode = $('.validate', node).first();
    }
    validate(notes, "ReactionNotes", validateNode);
  };