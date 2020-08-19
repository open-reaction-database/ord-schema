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
    const test = page.evaluate(async function(url) {
      const reaction = ord.reaction.unloadReaction();
      const reactions = session.dataset.getReactionsList();
      reactions[session.index] = reaction;
      return await ord.reaction
          .compareDataset(session.fileName, session.dataset)
          .then(() => {
            console.log('PASS', url);
            return 0;
          })
          .catch(() => {
            console.log('FAIL', url);
            return 1;
          })
    }, url);
    const testResult = await test;

    // Report results of testing to environment (shell, Git CI, etc.)
    // If _any_ test fails (i.e. testResult 1), then the entire process must
    // fail too. We still run all tests though, for convenience's sake.
    if (testResult == 1) {
      process.exitCode = 1;
    }
  }

  await browser.close();
})();