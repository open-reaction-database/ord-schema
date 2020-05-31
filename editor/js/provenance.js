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

goog.provide('ord.provenance');

goog.require('proto.ord.DateTime');
goog.require('proto.ord.Person');
goog.require('proto.ord.ReactionProvenance');
goog.require('proto.ord.RecordEvent');

ord.provenance.load = function (provenance) {
  const experimenter = provenance.getExperimenter();
  if (experimenter) {
    ord.provenance.loadPerson($('#provenance_experimenter'), experimenter);
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
    ord.provenance.loadRecordEvent($('#provenance_created'), created);
  }
  provenance.getRecordModifiedList().forEach(modified => {
    const node = ord.provenance.addModification();
    ord.provenance.loadRecordEvent(node, modified);
  });
};

ord.provenance.loadRecordEvent = function(node, record) {
  const time = record.getTime();
  if (time) {
    $('.provenance_time', node).text(time.getValue());
  }
  const person = record.getPerson();
  if (person) {
    ord.provenance.loadPerson(node, record.getPerson());
  }
  $('.provenance_details', node).text(record.getDetails());
};

ord.provenance.loadPerson = function (node, person) {
  $('.provenance_username', node).text(person.getUsername());
  $('.provenance_name', node).text(person.getName());
  $('.provenance_orcid', node).text(person.getOrcid());
  $('.provenance_organization', node).text(person.getOrganization());
};

ord.provenance.unload = function () {
  const provenance = new proto.ord.ReactionProvenance();

  const experimenter =
      ord.provenance.unloadPerson($('#provenance_experimenter'));
  provenance.setExperimenter(experimenter);

  provenance.setCity($('#provenance_city').text());

  const start = new proto.ord.DateTime();
  start.setValue($('#provenance_start').text());
  provenance.setExperimentStart(start);

  provenance.setDoi($('#provenance_doi').text());
  provenance.setPatent($('#provenance_patent').text());
  provenance.setPublicationUrl($('#provenance_url').text());

  const created = ord.provenance.unloadRecordEvent($('#provenance_created'));
  provenance.setRecordCreated(created);

  const modifieds = [];
  $('.provenance_modified', '#provenance_modifieds').each(
    function (index, node) {
      node = $(node);
      if (!node.attr('id')) {
        // Not a template.
        const modified = ord.provenance.unloadRecordEvent(node);
        modifieds.push(modified);
      }
    }
  );
  provenance.setRecordModifiedList(modifieds);
  return provenance;
};

ord.provenance.unloadRecordEvent = function (node) {
  const created = new proto.ord.RecordEvent();
  const createdTime = new proto.ord.DateTime();
  createdTime.setValue($('.provenance_time', node).text());
  created.setTime(createdTime);
  const createdPerson = ord.provenance.unloadPerson(node);
  created.setPerson(createdPerson);
  const createdDetails = $('.provenance_details', node).text();
  created.setDetails(createdDetails);
  return created;
};

ord.provenance.unloadPerson = function (node) {
  const person = new proto.ord.Person();
  person.setUsername($('.provenance_username', node).text());
  person.setName($('.provenance_name', node).text());
  person.setOrcid($('.provenance_orcid', node).text());
  person.setOrganization($('.provenance_organization', node).text());
  return person;
};

ord.provenance.addModification = function () {
  return addSlowly('#provenance_modified_template', '#provenance_modifieds');
};
