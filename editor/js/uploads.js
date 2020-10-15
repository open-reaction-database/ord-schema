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

goog.module('ord.uploads');
goog.module.declareLegacyNamespace();
exports = {
  putAll,
  retrieve,
  initialize,
  load,
  unload
};

// Map token strings to byte arrays to round-trip uploaded byte_values.
let tokenBytes = {};

// Map token strings to file system paths for new uploads.
let tokenFiles = {};

/**
 * Returns a short random string of ASCII characters.
 * @return {string}
 */
function newToken() {
  return 'upload_' + Math.random().toString(36).substring(2);
}

/**
 * Tags this file for upload and returns a random token to access it later.
 * @param {string} file The filename to be uploaded.
 * @return {string} A new key into `tokenFiles`.
 */
function newFile(file) {
  const token = newToken();
  tokenFiles[token] = file;
  return token;
}

/**
 * Returns the filename of a file that was previously recorded at newFile().
 * @param {string} token A key into `tokenFiles`.
 * @return {string}
 */
function getFile(token) {
  return tokenFiles[token];
}

/**
 * Stores bytes at load() so they can be restored at unload().
 * @param {!Uint8Array} bytesValue The bytes to store.
 * @returns {string} A new key into `tokenBytes`.
 */
function stashUpload(bytesValue) {
  const token = newToken();
  tokenBytes[token] = bytesValue;
  return token;
}

/**
 * Recalls bytes that were stashed at load().
 * @param {string} token A key into `tokenBytes`.
 * @return {!Uint8Array} The stored bytes.
 */
function unstashUpload(token) {
  return tokenBytes[token];
}

/**
 * Sends all files referenced in tokenFiles to the server.
 * @param {string} dirName Server directory in which to store files.
 */
function putAll(dirName) {
  const tokens = Object.getOwnPropertyNames(tokenFiles);
  tokens.forEach(token => {
    const file = tokenFiles[token];
    const reader = new FileReader();
    reader.readAsArrayBuffer(file);
    reader.onload = (event) => {
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/dataset/proto/upload/' + dirName + '/' + token);
      const payload = event.target.result;
      xhr.send(payload);
    };
  });
}

/**
 * Looks up the bytesValue of the given uploader and sends back it as a
 * download.
 * @param {!Node} uploader A `.data_uploader` div.
 */
function retrieve(uploader) {
  const token = uploader.attr('data-token');
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/download/' + token);
  xhr.onload = () => {
    // Make the browser write the file.
    const url = URL.createObjectURL(new Blob([xhr.response]));
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', token);
    document.body.appendChild(link);
    link.click();
  };
  const bytesValue = tokenBytes[token];
  xhr.send(bytesValue);
}

/**
 * Configures the behavior of an uploader.
 * @param {!Node} node A `.data_uploader` div.
 */
function initialize(node) {
  $('.data_uploader_chooser_file', node).on('input', (event) => {
    const file = event.target.files[0];
    $('.data_uploader_file_name', node).show();
    $('.data_uploader_file_name', node).text(file.name);
    const token = newFile(file);
    $('.data_uploader', node).attr('data-token', token);
    $('.data_uploader_file_retrieve', node).hide();
  });
}

/**
 * Loads a token and filename into an uploader.
 * @param {!Node} node A `.data_uploader` div.
 * @param {!Uint8Array} bytesValue File content as bytes.
 */
function load(node, bytesValue) {
  const token = stashUpload(bytesValue);
  $('.data_uploader', node).show();
  $('.data_uploader', node).attr('data-token', token);
  $('.data_uploader_chooser_button', node).text('Replace...');
  $('.data_uploader_file_name', node).hide();
  $('.data_uploader_file_retrieve', node).show();
}

/**
 * Retrieves the stored bytes from an uploader.
 * @param {!Node} node A `.data_uploader` div.
 * @returns {!Uint8Array}
 */
function unload(node) {
  const token = $('.data_uploader', node).attr('data-token');
  const bytesValue = unstashUpload(token);
  if (bytesValue) {
    // This is just a round-trip for bytesValue.
    return bytesValue;
  } else {
    // A new file has been chosen for upload.
    return (new TextEncoder()).encode(token);
  }
}
