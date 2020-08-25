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

// Return a short random string of ASCII characters.
function newToken() {
  return 'upload_' + Math.random().toString(36).substring(2);
}

// Tag this file for upload and return a random token to access it later.
function newFile(file) {
  const token = newToken();
  tokenFiles[token] = file;
  return token;
}

// Get back a file that was previously recorded at newFile().
function getFile(token) {
  return tokenFiles[token];
}

// Remember bytes at load() so they can be restored at unload().
function stashUpload(bytesValue) {
  const token = newToken();
  tokenBytes[token] = bytesValue;
  return token;
}

// At unload, recall bytes that were stashed at load().
function unstashUpload(token) {
  return tokenBytes[token];
}

// Sends all files referenced in tokenFiles to the server.
function putAll(fileName) {
  const tokens = Object.getOwnPropertyNames(tokenFiles);
  tokens.forEach(token => {
    const file = tokenFiles[token];
    const reader = new FileReader();
    reader.readAsBinaryString(file);
    reader.onload = (event) => {
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/dataset/proto/upload/' + fileName + '/' + token);
      const payload = event.target.result;
      xhr.send(payload);
    };
  });
}

// Look up the bytesValue of the given uploader and send back it as a download.
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

// Configure behaviors of a ".uploader" DOM element.
function initialize(node) {
  $('.uploader_chooser_file', node).on('input', (event) => {
    const file = event.target.files[0];
    $('.uploader_file_name', node).show();
    $('.uploader_file_name', node).text(file.name);
    const token = newFile(file);
    $('.uploader', node).attr('data-token', token);
    $('.uploader_file_retrieve', node).hide();
  });
}

function load(node, bytesValue) {
  const token = stashUpload(bytesValue);
  $('.uploader', node).show();
  $('.uploader', node).attr('data-token', token);
  $('.uploader_chooser_button', node).text('Replace...');
  $('.uploader_file_name', node).hide();
  $('.uploader_file_retrieve', node).show();
}

function unload(node) {
  const token = $('.uploader', node).attr('data-token');
  const bytesValue = unstashUpload(token);
  if (bytesValue) {
    // This is just a round-trip for bytesValue.
    return bytesValue;
  } else {
    // A new file has been chosen for upload.
    const bytes = (new TextEncoder()).encode(token);
    return bytes;
  }
}
