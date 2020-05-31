/**
 * Copyright 2020 The Open Reaction Database Authors
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
  setSelector($('#notes_heterogeneous'), notes.getIsHeterogeneous());
  setSelector($('#notes_precipitate'), notes.getFormsPrecipitate());
  setSelector($('#notes_exothermic'), notes.getIsExothermic());
  setSelector($('#notes_offgas'), notes.getOffgasses());
  setSelector($('#notes_moisture'), notes.getIsSensitiveToMoisture());
  setSelector($('#notes_oxygen'), notes.getIsSensitiveToOxygen());
  setSelector($('#notes_light'), notes.getIsSensitiveToLight());
  $('#notes_safety').text(notes.getSafetyNotes());
  $('#notes_details').text(notes.getProcedureDetails());
};

ord.notes.unload = function () {
  const notes = new proto.ord.ReactionNotes();
  notes.setIsHeterogeneous(getSelector($('#notes_heterogeneous')));
  notes.setFormsPrecipitate(getSelector($('#notes_precipitate')));
  notes.setIsExothermic(getSelector($('#notes_exothermic')));
  notes.setOffgasses(getSelector($('#notes_offgas')));
  notes.setIsSensitiveToMoisture(getSelector($('#notes_moisture')));
  notes.setIsSensitiveToOxygen(getSelector($('#notes_oxygen')));
  notes.setIsSensitiveToLight(getSelector($('#notes_light')));
  notes.setSafetyNotes($('#notes_safety').text());
  notes.setProcedureDetails($('#notes_details').text());
  return notes;
};
