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

goog.module('ord.products');
goog.module.declareLegacyNamespace();
exports = {load, unload, add, addIdentity, addYield, addPurity, addSelectivity, validateProduct};

goog.require('ord.compounds');
goog.require('proto.ord.ReactionProduct');

function load (node, products) {
  products.forEach(product => loadProduct(node, product));
};

function loadProduct (outcomeNode, product) {
  const node = add(outcomeNode);

  const compound = product.getCompound();
  if (compound) {
    ord.compounds.loadIntoCompound(node, compound);
  }

  ord.reaction.setOptionalBool(
      $('.outcome_product_desired', node),
      product.hasIsDesiredProduct() ? product.getIsDesiredProduct() : null);

  const compoundYield = product.getCompoundYield();
  if (compoundYield) {
    ord.reaction.writeMetric('.outcome_product_yield', compoundYield, node);
  }

  const purity = product.getPurity();
  if (purity) {
    ord.reaction.writeMetric('.outcome_product_purity', purity, node);
  }

  const selectivity = product.getSelectivity();
  if (selectivity) {
    ord.reaction.setSelector(
        $('.outcome_product_selectivity_type', node), selectivity.getType());
    $('.outcome_product_selectivity_details', node)
        .text(selectivity.getDetails());
    if (selectivity.hasValue()) {
      $('.outcome_product_selectivity_value', node)
          .text(selectivity.getValue());
    }
    if (selectivity.hasPrecision()) {
      $('.outcome_product_selectivity_precision', node)
          .text(selectivity.getPrecision());
    }
  }

  const identities = product.getAnalysisIdentityList();
  identities.forEach(identity => {
    const analysisNode = addIdentity(node);
    $('.analysis_key_selector', analysisNode).val(identity);
  });
  const yields = product.getAnalysisYieldList();
  yields.forEach(yeild => {
    const analysisNode = addYield(node);
    $('.analysis_key_selector', analysisNode).val(yeild);
  });
  const purities = product.getAnalysisPurityList();
  purities.forEach(purity => {
    const analysisNode = addPurity(node);
    $('.analysis_key_selector', analysisNode).val(purity);
  });
  const selectivities = product.getAnalysisSelectivityList();
  selectivities.forEach(selectivity => {
    const analysisNode = addSelectivity(node);
    $('.analysis_key_selector', analysisNode).val(selectivity);
  });
  $('.outcome_product_color', node).text(product.getIsolatedColor());

  const texture = product.getTexture();
  if (texture) {
    ord.reaction.setSelector(
        $('.outcome_product_texture_type', node), texture.getType());
    $('.outcome_product_texture_details', node).text(texture.getDetails());
  }
};

function unload (node) {
  const products = [];
  $('.outcome_product', node).each(function(index, productNode) {
    productNode = $(productNode);
    if (!productNode.attr('id')) {
      // Not a template.
      const product = unloadProduct(productNode);
      if (!ord.reaction.isEmptyMessage(product)) {
        products.push(product);
      }
    }
  });
  return products;
};

function unloadProduct (node) {
  const product = new proto.ord.ReactionProduct();

  const compoundNode = $('.outcome_product_compound');
  const compound = ord.compounds.unloadCompound(compoundNode);
  if (!ord.reaction.isEmptyMessage(compound)) {
    product.setCompound(compound);
  }

  product.setIsDesiredProduct(
      ord.reaction.getOptionalBool($('.outcome_product_desired', node)));

  const yeild = ord.reaction.readMetric(
      '.outcome_product_yield', new proto.ord.Percentage(), node);
  if (!ord.reaction.isEmptyMessage(yeild)) {
    product.setCompoundYield(yeild);
  }

  const purity = ord.reaction.readMetric(
      '.outcome_product_purity', new proto.ord.Percentage(), node);
  if (!ord.reaction.isEmptyMessage(purity)) {
    product.setPurity(purity);
  }

  const selectivity = new proto.ord.Selectivity();
  selectivity.setType(
      ord.reaction.getSelector($('.outcome_product_selectivity_type', node)));
  selectivity.setDetails(
      $('.outcome_product_selectivity_details', node).text());
  const selectivityValue =
      parseFloat($('.outcome_product_selectivity_value', node).text());
  if (!isNaN(selectivityValue)) {
    selectivity.setValue(selectivityValue);
  }
  const selectivityPrecision =
      parseFloat($('.outcome_product_selectivity_precision', node).text());
  if (!isNaN(selectivityPrecision)) {
    selectivity.setPrecision(selectivityPrecision);
  }
  if (!ord.reaction.isEmptyMessage(selectivity)) {
    product.setSelectivity(selectivity);
  }

  const identities = unloadAnalysisKeys(node, 'identity');
  product.setAnalysisIdentityList(identities);

  const yields = unloadAnalysisKeys(node, 'yield');
  product.setAnalysisYieldList(yields);

  const purities = unloadAnalysisKeys(node, 'purity');
  product.setAnalysisPurityList(purities);

  const selectivities = unloadAnalysisKeys(node, 'selectivity');
  product.setAnalysisSelectivityList(selectivities);

  const color = $('.outcome_product_color', node).text();
  product.setIsolatedColor(color);

  const texture = new proto.ord.ReactionProduct.Texture();
  texture.setType(
      ord.reaction.getSelector($('.outcome_product_texture_type', node)));
  texture.setDetails($('.outcome_product_texture_details', node).text());
  if (!ord.reaction.isEmptyMessage(texture)) {
    product.setTexture(texture);
  }

  return product;
};

function unloadAnalysisKeys (node, tag) {
  const values = [];
  $('.outcome_product_analysis_' + tag, node).each(function(index, tagNode) {
    tagNode = $(tagNode);
    if (!tagNode.attr('id')) {
      // Not a template.
      const value = $('.analysis_key_selector', tagNode).val();
      if (value != '') {
        values.push(value);
      }
    }
  });
  return values;
};

function add (node) {
  const productNode = ord.reaction.addSlowly(
      '#outcome_product_template', $('.outcome_products', node));

  // Add an empty compound node.
  ord.compounds.add(productNode);
  // The "compound" field is not repeated in ReactionProduct and so
  // ReactionComponents should not be added or removed.
  $('.component fieldset .remove', productNode).hide();
  // Product components implicitly have role Product.
  $('.component .component_role_limiting', productNode).hide();
  // Volume measurements of product components do not include solutes.
  $('.component .includes_solutes', productNode).remove();
  // Product components do not have Preparations nor Vendor information.
  $('.component .preparations_fieldset', productNode).hide();
  $('.component .vendor', productNode).hide();

  // Add live validation handling.
  ord.reaction.addChangeHandler(productNode, () => {
    validateProduct(productNode);
  });
  return productNode;
};

function addIdentity (node) {
  return ord.reaction.addSlowly(
      '#outcome_product_analysis_identity_template',
      $('.outcome_product_analysis_identities', node));
};

function addYield (node) {
  return ord.reaction.addSlowly(
      '#outcome_product_analysis_yield_template',
      $('.outcome_product_analysis_yields', node));
};

function addPurity (node) {
  return ord.reaction.addSlowly(
      '#outcome_product_analysis_purity_template',
      $('.outcome_product_analysis_purities', node));
};

function addSelectivity (node) {
  return ord.reaction.addSlowly(
      '#outcome_product_analysis_selectivity_template',
      $('.outcome_product_analysis_selectivities', node));
};

function validateProduct (node, validateNode) {
  const product = unloadProduct(node);
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  ord.reaction.validate(product, 'ReactionProduct', validateNode);
};