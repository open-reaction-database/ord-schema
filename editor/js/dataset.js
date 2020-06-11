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

goog.provide('ord.dataset');

goog.require('proto.ord.Dataset');

const session = {
  fileName: null,
  dataset: null
};

function init(fileName) {
  session.fileName = fileName;
  $('.edittext').attr('contentEditable', 'true');
  getDataset(fileName, loadDataset);
  listenDirty($('#text_fields'));
}

function listenDirty(node) {
  $('.edittext', node).on('input', dirty);
  $('.selector', node).on('input', dirty);
}

function dirty() {
  $('#save').css('visibility', 'visible');
}

function clean() {
  $('#save').css('visibility', 'hidden');
  $('#save').text('save');
}

function commit() {
  const dataset = unloadDataset();
  $('#save').text('saving');
  const xhr = new XMLHttpRequest();
  xhr.open(
      'POST', '/dataset/proto/write/' + session.fileName, true /* async */);
  const binary = dataset.serializeBinary();
  xhr.onload = clean;
  xhr.send(binary);
}

function download() {
  const xhr = new XMLHttpRequest();
  xhr.open('GET', '/dataset/' + session.fileName + '/download');
  xhr.onload = () => {
    // Make the browser write the file.
    const url = URL.createObjectURL(new Blob([xhr.response]));
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', session.fileName + '.pbtxt');
    document.body.appendChild(link);
    link.click();
  };
  xhr.send();
}

function getDataset(fileName, listener) {
  if (!listener) {
    return;
  }
  const xhr = new XMLHttpRequest();
  xhr.open('GET', '/dataset/proto/read/' + session.fileName, true /* async */);
  xhr.responseType = 'arraybuffer';
  xhr.onload = () => {
    const bytes = new Uint8Array(xhr.response);
    const dataset = proto.ord.Dataset.deserializeBinary(bytes);
    session.dataset = dataset;
    listener(dataset);
  };
  xhr.send();
}

function loadDataset(dataset) {
  $('#name').text(dataset.getName());
  $('#description').text(dataset.getDescription());
  $('#dataset_id').text(dataset.getDatasetId());

  const reactions = dataset.getReactionsList();
  loadReactions(reactions);

  const reactionIds = dataset.getReactionIdsList();
  loadReactionIds(reactionIds);

  const examples = dataset.getExamplesList();
  loadExamples(examples);

  clean();
}

function loadReactions(reactions) {
  for (var i = 0; i < reactions.length; i++) {
    const reaction = reactions[i];
    loadReaction(i, reaction);
  }
}

function loadReaction(index, reaction) {
  const node = addReaction(index);
  const id = reaction.getReactionId();
  $('.reaction_id', node).text(id);
}

function loadReactionIds(reactionIds) {
  reactionIds.forEach(reactionId => loadReactionId(reactionId));
}

function loadReactionId(reactionId) {
  const node = addReactionId();
  $('.other_reaction_id_text', node).text(reactionId);
}

function loadExamples(examples) {
  examples.forEach(example => loadExample(example));
}

function loadExample(example) {
  const node = addExample();
  $('.example_url', node).text(example.getUrl());
  $('.example_description', node).text(example.getDescription());
  const created = example.getCreated();
  if (created) {
    $('.example_created_time', node).text(created.getTime().getValue());
    $('.example_created_details', node).text(created.getDetails());
    const person = created.getPerson();
    if (person) {
      $('.example_created_person_username', node).text(person.getUsername());
      $('.example_created_person_name', node).text(person.getName());
      $('.example_created_person_orcid', node).text(person.getOrcid());
      $('.example_created_person_org', node).text(person.getOrganization());
    }
  }
}

function unloadDataset() {
  const dataset = session.dataset;
  dataset.setName($('#name').text());
  dataset.setDescription($('#description').text());
  dataset.setDatasetId($('#dataset_id').text());
  const reactionIds = [];
  $('.other_reaction_id').each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      reactionIds.push($('.other_reaction_id_text', node).text());
    }
  });
  dataset.setReactionIdsList(reactionIds);
  const examples = [];
  $('.example').each(function (index, node) {
    node = $(node);
    if (!node.attr('id')) {
      const example = unloadExample(node);
      examples.push(example);
    }
  });
  dataset.setExamplesList(examples);
  // Do not mutate Reactions. They are edited separately.
  return dataset;
}

function unloadExample(node) {
  const example = new proto.ord.DatasetExample();
  const url = $('.example_url', node).text();
  example.setUrl(url);
  const description = $('.example_description', node).text();
  example.setDescription(description);

  const created = new proto.ord.RecordEvent();
  const timeValue = $('.example_created_time', node).text();
  time = new proto.ord.DateTime();
  time.setValue(timeValue);
  created.setTime(time);
  const details = $('.example_created_details', node).text();
  created.setDetails(details);

  const person = new proto.ord.Person();
  const username = $('.example_created_person_username', node).text();
  person.setUsername(username);
  const name = $('.example_created_person_name', node).text();
  person.setName(name);
  const orcid = $('.example_created_person_orcid', node).text();
  person.setOrcid(orcid);
  const org = $('.example_created_person_org', node).text();
  person.setOrganization(org);

  created.setPerson(person);
  example.setCreated(created);

  return example;
}

function addReaction(index) {
  const node = $('#reaction_template').clone();
  node.removeAttr('id');
  const anchor = $('.reaction_index', node);
  anchor.text(index);
  anchor.attr('href', '/dataset/' + session.fileName + '/reaction/' + index);
  const root = $('#reactions');
  root.append(node);
  node.show('slow');
  listenDirty(node);
  dirty();
  return node;
}

function addReactionId() {
  const node = $('#other_reaction_id_template').clone();
  node.removeAttr('id');
  const root = $('#other_reaction_ids');
  root.append(node);
  node.show('slow');
  listenDirty(node);
  dirty();
  return node;
}

function addExample() {
  const node = $('#example_template').clone();
  node.removeAttr('id');
  const root = $('#examples');
  root.append(node);
  node.show('slow');
  listenDirty(node);
  dirty();
  return node;
}

function newReaction() {
  // Load the Reaction editor immediately without waiting for "save".
  window.location.href = '/dataset/' + session.fileName + '/new/reaction';
}

function deleteReaction(button) {
  // Delete the Reaction immediately without waiting for "save".
  const node = $(button).closest('.reaction');
  const index = parseInt($('a', node).text());
  window.location.href =
      '/dataset/' + session.fileName + '/delete/reaction/' + index;
}

function removeReactionId(button) {
  removeSlowly(button, '.other_reaction_id');
}

function removeExample(button) {
  removeSlowly(button, '.example');
}

function removeSlowly(button, pattern) {
  const node = $(button).closest(pattern);
  node.hide('slow', () => node.remove());
  dirty();
}
