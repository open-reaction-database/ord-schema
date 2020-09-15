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
  addProcessedData,
  addRawData,
  validateOutcome,
  validateAnalysis
};

goog.require('ord.data');
goog.require('ord.products');
goog.require('proto.ord.ReactionOutcome');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the reaction outcome sections in the form.
 * @param {!Array<!proto.ord.ReactionOutcome>} outcomes
 */
function load(outcomes) {
  outcomes.forEach(outcome => loadOutcome(outcome));
}

/**
 * Adds and populates a reaction outcome section in the form.
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
 * Adds and populates a reaction analysis section in the form.
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
    const processNode = addProcessedData(node);
    loadProcessedData(processNode, name, process);
  });

  const raws = analysis.getRawDataMap();
  const rawNames = raws.stringKeys_();
  rawNames.forEach(function(name) {
    const raw = raws.get(name);
    const rawNode = addRawData(node);
    loadRawData(rawNode, name, raw);
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
 * Adds and populates a processed_data section in a reaction analysis.
 * @param {!Node} node Parent reaction analysis node.
 * @param {string} name The name of this Data record.
 * @param {!proto.ord.Data} processedData
 */
function loadProcessedData(node, name, processedData) {
  $('.outcome_processed_data_name', node).text(name);
  ord.data.loadData(node, processedData);
}

/**
 * Adds and populates a raw_data section in a reaction analysis.
 * @param {!Node} node Parent reaction analysis node.
 * @param {string} name The name of this Data record.
 * @param {!proto.ord.Data} rawData
 */
function loadRawData(node, name, rawData) {
  $('.outcome_raw_data_name', node).text(name);
  ord.data.loadData(node, rawData);
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

  const processedDataMap = analysis.getProcessedDataMap();
  $('.outcome_processed_data', analysisNode)
      .each(function(index, processedDataNode) {
        processedDataNode = $(processedDataNode);
        if (!processedDataNode.attr('id')) {
          unloadProcessedData(processedDataNode, processedDataMap);
        }
      });
  const rawDataMap = analysis.getRawDataMap();
  $('.outcome_raw_data', analysisNode).each(function(index, rawDataNode) {
    rawDataNode = $(rawDataNode);
    if (!rawDataNode.attr('id')) {
      unloadRawData(rawDataNode, rawDataMap);
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
 * @param {!jspb.Map<string, !proto.ord.Data>} processedDataMap
 */
function unloadProcessedData(node, processedDataMap) {
  const name = $('.outcome_processed_data_name', node).text();
  const processedData = ord.data.unloadData(node);
  if (!ord.reaction.isEmptyMessage(name) ||
      !ord.reaction.isEmptyMessage(processedData)) {
    processedDataMap.set(name, processedData);
  }
}

/**
 * Fetches a raw_data record defined in the form and adds it to
 * `processes`.
 * @param {!Node} node Root node for the Data record.
 * @param {!jspb.Map<string, !proto.ord.Data>} rawDataMap
 */
function unloadRawData(node, rawDataMap) {
  const name = $('.outcome_raw_data_name', node).text();
  const rawData = ord.data.unloadData(node);
  if (!ord.reaction.isEmptyMessage(name) || !ord.reaction.isEmptyMessage(raw)) {
    rawDataMap.set(name, rawData);
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
function addProcessedData(node) {
  const processNode = ord.reaction.addSlowly(
      '#outcome_processed_data_template',
      $('.outcome_processed_data_repeated', node));
  ord.data.addData(processNode);
  return processNode;
}

/**
 * Adds a new raw_data section to the form.
 * @param {!Node} node Parent reaction outcome node.
 * @return {!Node} The newly added parent node for the Data record.
 */
function addRawData(node) {
  const rawNode = ord.reaction.addSlowly(
      '#outcome_raw_data_template', $('.outcome_raw_data_repeated', node));
  ord.data.addData(rawNode);
  return rawNode;
}

/**
 * Validates a reaction outcome defined in the form.
 * @param {!Node} node Root node for the reaction outcome.
 * @param {?Node} validateNode The target node for validation results.
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
 * @param {?Node} validateNode The target node for validation results.
 */
function validateAnalysis(node, validateNode) {
  const analysis = unloadAnalysisSingle(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(analysis, 'ReactionAnalysis', validateNode);
}
