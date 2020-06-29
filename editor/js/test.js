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

/**
 * Automated tests of Reaction proto round-trips: proto -> DOM -> proto.
 *
 *   $ node js/test.js
 */

const puppeteer = require('puppeteer');

(async () => {
  const browser = await puppeteer.launch();
  const [page] = await browser.pages();

  // Relay console messages.
  page.on('console', msg => console.log('console>', msg.text()));

  // Round-trip these reactions through the DOM and compare at the server.
  const tests = [
    'http://localhost:5000/dataset/empty/reaction/0?user=test',
    'http://localhost:5000/dataset/full/reaction/0?user=test',
  ];

  for (let i = 0; i < tests.length; i++) {
    const url = tests[i];
    await page.goto(url);
    await page.waitFor('body[ready=true]');
    await page.evaluate(async url => {
      const reaction = unloadReaction();
      const reactions = session.dataset.getReactionsList();
      reactions[session.index] = reaction;
      await compareDataset(session.fileName, session.dataset)
          .then(() => console.log('PASS', url))
          .catch(() => console.log('FAIL', url));
    }, url);
  }
  await browser.close();
})();
