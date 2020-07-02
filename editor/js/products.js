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

goog.provide('ord.products');

goog.require('ord.compounds');
goog.require('proto.ord.ReactionProduct');

ord.products.load = function (node, products) {
  products.forEach(product => ord.products.loadProduct(node, product));
};

ord.products.loadProduct = function (outcomeNode, product) {
  const node = ord.products.add(outcomeNode);

  const compound = product.getCompound();
  if (compound) {
    // Creates an empty compound node, and loads it.
    ord.compounds.loadCompound(node, compound);
  }
  else {
    // Add an empty compound node.
    ord.compounds.add(node);
  }
  // The "compound" field is not repeated in ReactionProduct and so
  // ReactionComponents should not be added or removed.
  $('.component fieldset .remove', node).hide();
  // Product components implicitly have role Product.
  $('.component .component_role_limiting', node).hide();
  // Volume measurements of product components do not include solutes. 
  $('.component .includes_solutes', node).remove();
  // Product components do not have Preparations nor Vendor information.
  $('.component .preparations_fieldset', node).hide();
  $('.component .vendor', node).hide();

  setOptionalBool(
      $('.outcome_product_desired', node),
      product.hasIsDesiredProduct() ? product.getIsDesiredProduct() : null);

  writeMetric('.outcome_product_yield', product.getCompoundYield(), node);

  writeMetric('.outcome_product_purity', product.getPurity(), node);

  const selectivity = product.getSelectivity();
  setSelector(
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
  const identities = product.getAnalysisIdentityList();
  identities.forEach(identity => {
    const analysisNode = ord.products.addIdentity(node);
    $('.outcome_product_analysis_identity_value', analysisNode).text(identity);
  });
  const yields = product.getAnalysisYieldList();
  yields.forEach(yeild => {
    const analysisNode = ord.products.addYield(node);
    $('.outcome_product_analysis_yield_value', analysisNode).text(yeild);
  });
  const purities = product.getAnalysisPurityList();
  purities.forEach(purity => {
    const analysisNode = ord.products.addPurity(node);
    $('.outcome_product_analysis_purity_value', analysisNode).text(purity);
  });
  const selectivities = product.getAnalysisSelectivityList();
  selectivities.forEach(selectivity => {
    const analysisNode = ord.products.addSelectivity(node);
    $('.outcome_product_analysis_selectivity_value', analysisNode)
        .text(selectivity);
  });
  $('.outcome_product_color', node).text(product.getIsolatedColor());

  const texture = product.getTexture();
  if (texture) {
    setSelector($('.outcome_product_texture_type', node), texture.getType());
    $('.outcome_product_texture_details', node).text(texture.getDetails());
  }
};

ord.products.unload = function (node) {
  const products = [];
  $('.outcome_product', node).each(function (index, productNode) {
    productNode = $(productNode);
    if (!productNode.attr('id')) {
      // Not a template.
      const product = ord.products.unloadProduct(productNode);
      products.push(product);
    }
  });
  return products;
};

ord.products.unloadProduct = function (node) {
  const product = new proto.ord.ReactionProduct();

  const compounds = ord.compounds.unload(node);
  if (compounds) {
    product.setCompound(compounds[0]);
  }
  product.setIsDesiredProduct(
      getOptionalBool($('.outcome_product_desired', node)));

  const yeild =
      readMetric('.outcome_product_yield', new proto.ord.Percentage(), node);
  product.setCompoundYield(yeild);

  const purity =
      readMetric('.outcome_product_purity', new proto.ord.Percentage(), node);
  product.setPurity(purity);

  const selectivity = new proto.ord.Selectivity();
  selectivity.setType(
      getSelector($('.outcome_product_selectivity_type', node)));
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
  product.setSelectivity(selectivity);

  const identities = ord.products.unloadAnalysisKeys(node, 'identity');
  product.setAnalysisIdentityList(identities);

  const yields = ord.products.unloadAnalysisKeys(node, 'yield');
  product.setAnalysisYieldList(yields);

  const purities = ord.products.unloadAnalysisKeys(node, 'purity');
  product.setAnalysisPurityList(purities);

  const selectivities = ord.products.unloadAnalysisKeys(node, 'selectivity');
  product.setAnalysisSelectivityList(selectivities);

  const color = $('.outcome_product_color', node).text();
  product.setIsolatedColor(color);

  const texture = new proto.ord.ReactionProduct.Texture();
  texture.setType(getSelector($('.outcome_product_texture_type', node)));
  texture.setDetails($('.outcome_product_texture_details', node).text());
  product.setTexture(texture);

  return product;
};

ord.products.unloadAnalysisKeys = function (node, tag) {
  const values = [];
  $('.outcome_product_analysis_' + tag, node).each(
    function (index, tagNode) {
      tagNode = $(tagNode);
      if (!tagNode.attr('id')) {
        // Not a template.
        const value =
            $('.outcome_product_analysis_' + tag + '_value', tagNode).text();
        values.push(value);
      }
    }
  );
  return values;
};

ord.products.add = function (node) {
  const productNode = addSlowly('#outcome_product_template', $('.outcome_products', node));
  handler = function () {ord.products.validateProduct(productNode)};
  addChangeHandler(productNode, handler);
  return productNode;
};

ord.products.addIdentity = function (node) {
  return addSlowly(
      '#outcome_product_analysis_identity_template',
      $('.outcome_product_analysis_identities', node));
};

ord.products.addYield = function (node) {
  return addSlowly(
      '#outcome_product_analysis_yield_template',
      $('.outcome_product_analysis_yields', node));
};

ord.products.addPurity = function (node) {
  return addSlowly(
      '#outcome_product_analysis_purity_template',
      $('.outcome_product_analysis_purities', node));
};

ord.products.addSelectivity = function (node) {
  return addSlowly(
      '#outcome_product_analysis_selectivity_template',
      $('.outcome_product_analysis_selectivities', node));
};

ord.products.validateProduct = function(node, validateNode) {
  const product = ord.products.unloadProduct(node);
  if (typeof validateNode === 'undefined') {
    validateNode = $('.validate', node).first();
  }
  validate(product, "ReactionProduct", validateNode);
};