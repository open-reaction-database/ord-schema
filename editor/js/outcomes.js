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

goog.module('ord.outcomes');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  add,
  addAnalysis,
  addProcess,
  addRaw,
  validateOutcome,
  validateAnalysis
};

goog.require('ord.products');
goog.require('proto.ord.ReactionOutcome');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Populates the reaction outcome sections in the form.
 * @param {!Array<!proto.ord.ReactionOutcome>} outcomes
 */
function load(outcomes) {
  outcomes.forEach(outcome => loadOutcome(outcome));
}

/**
 * Populates a reaction outcome section in the form.
 * @param outcome
 */
function loadOutcome(outcome) {
  const node = add();

  const time = outcome.getReactionTime();
  if (time != null) {
    ord.reaction.writeMetric('.outcome_time', time, node);
  }
  const conversion = outcome.getConversion();
  if (conversion) {
    ord.reaction.writeMetric(
        '.outcome_conversion', outcome.getConversion(), node);
  }

  const analyses = outcome.getAnalysesMap();
  const names = analyses.stringKeys_();
  names.forEach(function(name) {
    const analysis = analyses.get(name);
    loadAnalysis(node, name, analysis);
  });

  const products = outcome.getProductsList();
  ord.products.load(node, products);
}

/**
 * Populates a reaction analysis section in the form.
 * @param {!Node} outcomeNode Parent reaction outcome node.
 * @param {string} name The name of this analysis.
 * @param {!proto.ord.ReactionAnalysis} analysis
 */
function loadAnalysis(outcomeNode, name, analysis) {
  const node = addAnalysis(outcomeNode);

  $('.outcome_analysis_name', node).text(name).trigger('input');

  ord.reaction.setSelector(
      $('.outcome_analysis_type', node), analysis.getType());
  const chmoId = analysis.getChmoId();
  if (chmoId != 0) {
    $('.outcome_analysis_chmo_id', node).text(analysis.getChmoId());
  }
  $('.outcome_analysis_details', node).text(analysis.getDetails());

  const processes = analysis.getProcessedDataMap();
  const processNames = processes.stringKeys_();
  processNames.forEach(function(name) {
    const process = processes.get(name);
    const processNode = addProcess(node);
    loadProcess(processNode, name, process);
  });

  const raws = analysis.getRawDataMap();
  const rawNames = raws.stringKeys_();
  rawNames.forEach(function(name) {
    const raw = raws.get(name);
    const rawNode = addRaw(node);
    loadRaw(rawNode, name, raw);
  });
  $('.outcome_analysis_manufacturer', node)
      .text(analysis.getInstrumentManufacturer());
  const calibrated = analysis.getInstrumentLastCalibrated();
  if (calibrated) {
    $('.outcome_analysis_calibrated', node).text(calibrated.getValue());
  }
  ord.reaction.setOptionalBool(
      $('.outcome_analysis_internal_standard', node),
      analysis.hasUsesInternalStandard() ? analysis.getUsesInternalStandard() :
                                           null);
  ord.reaction.setOptionalBool(
      $('.outcome_analysis_authentic_standard', node),
      analysis.hasUsesAuthenticStandard() ?
          analysis.getUsesAuthenticStandard() :
          null);
}

/**
 * Populates a processed_data section in a reaction analysis.
 * @param {!Node} node Parent reaction analysis node.
 * @param {string} name The name of this Data record.
 * @param {!proto.ord.Data} process
 */
function loadProcess(node, name, process) {
  $('.outcome_process_name', node).text(name);
  $('.outcome_process_description', node).text(process.getDescription());
  $('.outcome_process_format', node).text(process.getFormat());

  const stringValue = process.getStringValue();
  const floatValue = process.getFloatValue();
  const bytesValue = process.getBytesValue();
  const url = process.getUrl();
  if (stringValue) {
    $('.outcome_process_text', node).show();
    $('.uploader', node).hide();
    $('.outcome_process_text', node).text(stringValue);
    $('input[value=\'text\']', node).prop('checked', true);
  }
  if (floatValue) {
    $('.outcome_process_text', node).show();
    $('.uploader', node).hide();
    $('.outcome_process_text', node).text(floatValue);
    $('input[value=\'number\']', node).prop('checked', true);
  }
  if (bytesValue) {
    $('.outcome_process_text', node).hide();
    $('.uploader', node).show();
    ord.uploads.load(node, bytesValue);
    $('input[value=\'upload\']', node).prop('checked', true);
  }
  if (url) {
    $('.outcome_process_text', node).show();
    $('.uploader', node).hide();
    $('.outcome_process_text', node).text(url);
    $('input[value=\'url\']', node).prop('checked', true);
  }
}

/**
 * Populates a raw_data section in a reaction analysis.
 * @param {!Node} node Parent reaction analysis node.
 * @param {string} name The name of this Data record.
 * @param {!proto.ord.Data} raw
 */
function loadRaw(node, name, raw) {
  $('.outcome_raw_name', node).text(name);
  $('.outcome_raw_description', node).text(raw.getDescription());
  $('.outcome_raw_format', node).text(raw.getFormat());

  const stringValue = raw.getStringValue();
  const floatValue = raw.getFloatValue();
  const bytesValue = raw.getBytesValue();
  const url = raw.getUrl();
  if (stringValue) {
    $('.outcome_raw_text', node).show();
    $('.uploader', node).hide();
    $('.outcome_raw_text', node).text(stringValue);
    $('input[value=\'text\']', node).prop('checked', true);
  }
  if (floatValue) {
    $('.outcome_raw_text', node).show();
    $('.uploader', node).hide();
    $('.outcome_raw_text', node).text(floatValue);
    $('input[value=\'number\']', node).prop('checked', true);
  }
  if (bytesValue) {
    $('.outcome_raw_text', node).hide();
    $('.uploader', node).show();
    ord.uploads.load(node, bytesValue);
    $('input[value=\'upload\']', node).prop('checked', true);
  }
  if (url) {
    $('.outcome_raw_text', node).show();
    $('.uploader', node).hide();
    $('.outcome_raw_text', node).text(url);
    $('input[value=\'url\']', node).prop('checked', true);
  }
}

/**
 * Fetches the reaction outcomes defined in the form.
 * @return {!Array<!proto.ord.ReactionOutcome>}
 */
function unload() {
  const outcomes = [];
  $('.outcome').each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      const outcome = unloadOutcome(node);
      if (!ord.reaction.isEmptyMessage(outcome)) {
        outcomes.push(outcome);
      }
    }
  });
  return outcomes;
}

/**
 * Fetches a reaction outcome defined in the form.
 * @param {!Node} node Root node for the reaction outcome.
 * @return {!proto.ord.ReactionOutcome}
 */
function unloadOutcome(node) {
  const outcome = new proto.ord.ReactionOutcome();

  const time =
      ord.reaction.readMetric('.outcome_time', new proto.ord.Time(), node);
  if (!ord.reaction.isEmptyMessage(time)) {
    outcome.setReactionTime(time);
  }

  const conversion = ord.reaction.readMetric(
      '.outcome_conversion', new proto.ord.Percentage(), node);
  if (!ord.reaction.isEmptyMessage(conversion)) {
    outcome.setConversion(conversion);
  }

  const products = ord.products.unload(node);
  outcome.setProductsList(products);

  const analyses = outcome.getAnalysesMap();
  $('.outcome_analysis', node).each(function(index, node) {
    node = $(node);
    if (!node.attr('id')) {
      // Not a template.
      unloadAnalysis(node, analyses);
    }
  });
  return outcome;
}

/**
 * Fetches a reaction analysis defined in the form.
 * @param {!Node} analysisNode Root node for the reaction analysis.
 * @return {!proto.ord.ReactionAnalysis}
 */
function unloadAnalysisSingle(analysisNode) {
  const analysis = new proto.ord.ReactionAnalysis();
  analysis.setType(
      ord.reaction.getSelector($('.outcome_analysis_type', analysisNode)));
  const chmoId = $('.outcome_analysis_chmo_id', analysisNode).text();
  if (!isNaN(chmoId)) {
    analysis.setChmoId(chmoId);
  }
  analysis.setDetails($('.outcome_analysis_details', analysisNode).text());

  const processes = analysis.getProcessedDataMap();
  $('.outcome_process', analysisNode).each(function(index, processNode) {
    processNode = $(processNode);
    if (!processNode.attr('id')) {
      unloadProcess(processNode, processes);
    }
  });
  const raws = analysis.getRawDataMap();
  $('.outcome_raw', analysisNode).each(function(index, rawNode) {
    rawNode = $(rawNode);
    if (!rawNode.attr('id')) {
      unloadRaw(rawNode, raws);
    }
  });
  analysis.setInstrumentManufacturer(
      $('.outcome_analysis_manufacturer', analysisNode).text());
  const calibrated = new proto.ord.DateTime();
  calibrated.setValue($('.outcome_analysis_calibrated', analysisNode).text());
  if (!ord.reaction.isEmptyMessage(calibrated)) {
    analysis.setInstrumentLastCalibrated(calibrated);
  }
  analysis.setUsesInternalStandard(ord.reaction.getOptionalBool(
      $('.outcome_analysis_internal_standard', analysisNode)));
  analysis.setUsesAuthenticStandard(ord.reaction.getOptionalBool(
      $('.outcome_analysis_authentic_standard', analysisNode)));

  return analysis;
}

/**
 * Fetches a reaction analysis defined in the form and adds it to `analyses`.
 * @param {!Node} analysisNode Root node for the reaction analysis.
 * @param {!jspb.Map<string, !proto.ord.ReactionAnalysis>} analyses
 */
function unloadAnalysis(analysisNode, analyses) {
  const analysis = unloadAnalysisSingle(analysisNode);
  const name = $('.outcome_analysis_name', analysisNode).text();
  if (!ord.reaction.isEmptyMessage(name) ||
      !ord.reaction.isEmptyMessage(analysis)) {
    analyses.set(name, analysis);
  }
}

/**
 * Fetches a processed_data record defined in the form and adds it to
 * `processes`.
 * @param {!Node} node Root node for the Data record.
 * @param {!jspb.Map<string, !proto.ord.Data>} processes
 */
function unloadProcess(node, processes) {
  const name = $('.outcome_process_name', node).text();

  const process = new proto.ord.Data();
  process.setDescription($('.outcome_process_description').text());
  process.setFormat($('.outcome_process_format').text());

  if ($('input[value=\'text\']', node).is(':checked')) {
    const stringValue = $('.outcome_process_text', node).text();
    if (!ord.reaction.isEmptyMessage(stringValue)) {
      process.setStringValue(stringValue);
    }
  }
  if ($('input[value=\'number\']', node).is(':checked')) {
    const floatValue = parseFloat($('.outcome_process_text', node).text());
    if (!isNaN(floatValue)) {
      process.setFloatValue(floatValue);
    }
  }
  if ($('input[value=\'upload\']', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!ord.reaction.isEmptyMessage(bytesValue)) {
      process.setBytesValue(bytesValue);
    }
  }
  if ($('input[value=\'url\']', node).is(':checked')) {
    const url = $('.outcome_process_text', node).text();
    if (!ord.reaction.isEmptyMessage(url)) {
      process.setUrl(url);
    }
  }
  if (!ord.reaction.isEmptyMessage(name) ||
      !ord.reaction.isEmptyMessage(process)) {
    processes.set(name, process);
  }
}

/**
 * Fetches a raw_data record defined in the form and adds it to
 * `processes`.
 * @param {!Node} node Root node for the Data record.
 * @param {!jspb.Map<string, !proto.ord.Data>} raws
 */
function unloadRaw(node, raws) {
  const name = $('.outcome_raw_name', node).text();

  const raw = new proto.ord.Data();
  raw.setDescription($('.outcome_raw_description', node).text());
  raw.setFormat($('.outcome_raw_format', node).text());

  if ($('input[value=\'text\']', node).is(':checked')) {
    const stringValue = $('.outcome_raw_text', node).text();
    if (!ord.reaction.isEmptyMessage(stringValue)) {
      raw.setStringValue(stringValue);
    }
  }
  if ($('input[value=\'number\']', node).is(':checked')) {
    const floatValue = parseFloat($('.outcome_raw_text', node).text());
    if (!isNaN(floatValue)) {
      raw.setFloatValue(floatValue);
    }
  }
  if ($('input[value=\'upload\']', node).is(':checked')) {
    const bytesValue = ord.uploads.unload(node);
    if (!ord.reaction.isEmptyMessage(bytesValue)) {
      raw.setBytesValue(bytesValue);
    }
  }
  if ($('input[value=\'url\']', node).is(':checked')) {
    const url = $('.outcome_raw_text', node).text();
    if (!ord.reaction.isEmptyMessage(url)) {
      raw.setUrl(url);
    }
  }
  if (!ord.reaction.isEmptyMessage(name) || !ord.reaction.isEmptyMessage(raw)) {
    raws.set(name, raw);
  }
}

/**
 * Adds a reaction outcome section to the form.
 * @return {!Node} The newly added parent node for the reaction outcome.
 */
function add() {
  const node = ord.reaction.addSlowly('#outcome_template', '#outcomes');
  // Add live validation handling.
  ord.reaction.addChangeHandler(node, () => {
    validateOutcome(node);
  });
  return node;
}

/**
 * Adds a reaction analysis section to the form.
 * @param {!Node} node Parent reaction outcome node.
 * @return {!Node} The newly added parent node for the reaction analysis.
 */
function addAnalysis(node) {
  const analysisNode = ord.reaction.addSlowly(
      '#outcome_analysis_template', $('.outcome_analyses', node));

  // Handle name changes.
  const nameNode = $('.outcome_analysis_name', analysisNode);
  nameNode.on('focusin', function() {
    // Store old value in val attribute.
    nameNode.data('val', nameNode.text());
  });
  nameNode.on('input', function() {
    var old_name = nameNode.data('val');
    var name = nameNode.text();
    // Remove old key.
    if (old_name) {
      // If any selector had this value selected, reset it.
      $('.analysis_key_selector', node).each(function() {
        if ($(this).val() == old_name) {
          $(this).val('');
        }
      });
      $('.analysis_key_selector option[value="' + old_name + '"]', node)
          .remove();
    }
    // Add new key.
    if (name) {
      $('.analysis_key_selector', node)
          .append('<option value="' + name + '">' + name + '</option>');
      // Ensure old value stored (necessary if focus does not change).
      nameNode.data('val', name);
    }
  });

  // Add live validation handling.
  ord.reaction.addChangeHandler(analysisNode, () => {
    validateAnalysis(analysisNode);
  });
  return analysisNode;
}

/**
 * Adds a new processed_data section to the form.
 * @param {!Node} node Parent reaction outcome node.
 * @return {!Node} The newly added parent node for the Data record.
 */
function addProcess(node) {
  const processNode = ord.reaction.addSlowly(
      '#outcome_process_template', $('.outcome_processes', node));

  const typeButtons = $('input[type=\'radio\']', processNode);
  typeButtons.attr('name', 'outcomes_' + radioGroupCounter++);
  typeButtons.change(function() {
    if ((this.value == 'text') || (this.value == 'number') ||
        (this.value == 'url')) {
      $('.outcome_process_text', processNode).show();
      $('.uploader', processNode).hide();
    } else {
      $('.outcome_process_text', processNode).hide();
      $('.uploader', processNode).show();
    }
  });
  ord.uploads.initialize(processNode);
  return processNode;
}

/**
 * Adds a new raw_data section to the form.
 * @param {!Node} node Parent reaction outcome node.
 * @return {!Node} The newly added parent node for the Data record.
 */
function addRaw(node) {
  const rawNode =
      ord.reaction.addSlowly('#outcome_raw_template', $('.outcome_raws', node));

  const typeButtons = $('input[type=\'radio\']', rawNode);
  typeButtons.attr('name', 'outcomes_' + radioGroupCounter++);
  typeButtons.change(function() {
    if ((this.value == 'text') || (this.value == 'number') ||
        (this.value == 'url')) {
      $('.outcome_raw_text', rawNode).show();
      $('.uploader', rawNode).hide();
    } else {
      $('.outcome_raw_text', rawNode).hide();
      $('.uploader', rawNode).show();
    }
  });
  ord.uploads.initialize(rawNode);
  return rawNode;
}

/**
 * Validates a reaction outcome defined in the form.
 * @param {!Node} node Root node for the reaction outcome.
 * @param {!Node} validateNode The target node for validation results.
 */
function validateOutcome(node, validateNode) {
  const outcome = unloadOutcome(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(outcome, 'ReactionOutcome', validateNode);
}

/**
 * Validates a reaction analysis defined in the form.
 * @param {!Node} node Root node for the reaction analysis.
 * @param {!Node} validateNode The target node for validation results.
 */
function validateAnalysis(node, validateNode) {
  const analysis = unloadAnalysisSingle(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(analysis, 'ReactionAnalysis', validateNode);
}
