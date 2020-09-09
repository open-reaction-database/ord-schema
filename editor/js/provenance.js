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

goog.module('ord.provenance');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  addModification,
  validateProvenance
};

goog.require('proto.ord.DateTime');
goog.require('proto.ord.Person');
goog.require('proto.ord.ReactionProvenance');
goog.require('proto.ord.RecordEvent');

/**
 * Adds and populates the provenance section in the form.
 * @param {!proto.ord.ReactionProvenance} provenance
 */
function load(provenance) {
  const experimenter = provenance.getExperimenter();
  if (experimenter) {
    loadPerson($('#provenance_experimenter'), experimenter);
  }
  $('#provenance_city').text(provenance.getCity());

  const start = provenance.getExperimentStart();
  if (start) {
    $('#provenance_start').text(start.getValue());
  }

  $('#provenance_doi').text(provenance.getDoi());
  $('#provenance_patent').text(provenance.getPatent());
  $('#provenance_url').text(provenance.getPublicationUrl());

  const created = provenance.getRecordCreated();
  if (created) {
    loadRecordEvent($('#provenance_created'), created);
  }
  provenance.getRecordModifiedList().forEach(modified => {
    const node = addModification();
    loadRecordEvent(node, modified);
  });
}

/**
 * Adds and populates a record event section in the form.
 * @param {!Node} node The target div.
 * @param {!proto.ord.RecordEvent} record
 */
function loadRecordEvent(node, record) {
  const time = record.getTime();
  if (time) {
    $('.provenance_time', node).text(time.getValue());
  }
  const person = record.getPerson();
  if (person) {
    loadPerson(node, record.getPerson());
  }
  $('.provenance_details', node).text(record.getDetails());
}

/**
 * Adds and populates a person section in the form.
 * @param {!Node} node The target div.
 * @param {!proto.ord.Person} person
 */
function loadPerson(node, person) {
  $('.provenance_username', node).text(person.getUsername());
  $('.provenance_name', node).text(person.getName());
  $('.provenance_orcid', node).text(person.getOrcid());
  $('.provenance_organization', node).text(person.getOrganization());
  $('.provenance_email', node).text(person.getEmail());
}

/**
 * Fetches reaction provenance as defined in the form.
 * @return {!proto.ord.ReactionProvenance}
 */
function unload() {
  const provenance = new proto.ord.ReactionProvenance();

  const experimenter = unloadPerson($('#provenance_experimenter'));
  if (!ord.reaction.isEmptyMessage(experimenter)) {
    provenance.setExperimenter(experimenter);
  }

  provenance.setCity($('#provenance_city').text());

  const start = new proto.ord.DateTime();
  start.setValue($('#provenance_start').text());
  if (!ord.reaction.isEmptyMessage(start)) {
    provenance.setExperimentStart(start);
  }

  provenance.setDoi($('#provenance_doi').text());
  provenance.setPatent($('#provenance_patent').text());
  provenance.setPublicationUrl($('#provenance_url').text());

  const created = unloadRecordEvent($('#provenance_created'));
  if (!ord.reaction.isEmptyMessage(created)) {
    provenance.setRecordCreated(created);
  }

  const modifieds = [];
  $('.provenance_modified', '#provenance_modifieds')
      .each(function(index, node) {
        node = $(node);
        if (!node.attr('id')) {
          // Not a template.
          const modified = unloadRecordEvent(node);
          if (!ord.reaction.isEmptyMessage(modified)) {
            modifieds.push(modified);
          }
        }
      });
  provenance.setRecordModifiedList(modifieds);
  return provenance;
}

/**
 * Fetches a record event as defined in the form.
 * @param {!Node} node Parent node containing the record event.
 * @return {!proto.ord.RecordEvent}
 */
function unloadRecordEvent(node) {
  const created = new proto.ord.RecordEvent();
  const createdTime = new proto.ord.DateTime();
  createdTime.setValue($('.provenance_time', node).text());
  if (!ord.reaction.isEmptyMessage(createdTime)) {
    created.setTime(createdTime);
  }
  const createdPerson = unloadPerson(node);
  if (!ord.reaction.isEmptyMessage(createdPerson)) {
    created.setPerson(createdPerson);
  }
  const createdDetails = $('.provenance_details', node).text();
  created.setDetails(createdDetails);
  return created;
}

/**
 * Fetches a person message as defined in the form.
 * @param {!Node} node Parent node containing the person message.
 * @return {!proto.ord.Person}
 */
function unloadPerson(node) {
  const person = new proto.ord.Person();
  person.setUsername($('.provenance_username', node).text());
  person.setName($('.provenance_name', node).text());
  person.setOrcid($('.provenance_orcid', node).text());
  person.setOrganization($('.provenance_organization', node).text());
  person.setEmail($('.provenance_email', node).text());
  return person;
}

/**
 * Adds a record_modified section to the form.
 * @return {!Node} The div containing the new event.
 */
function addModification() {
  return ord.reaction.addSlowly(
      '#provenance_modified_template', '#provenance_modifieds');
}

/**
 * Validates the reaction provenence defined in the form.
 * @param {!Node} node The node containing reaction provenance information.
 * @param {?Node} validateNode The target div for validation results.
 */
function validateProvenance(node, validateNode) {
  const provenance = unload();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(provenance, 'ReactionProvenance', validateNode);
}
