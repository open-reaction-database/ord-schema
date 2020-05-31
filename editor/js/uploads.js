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

goog.provide('ord.uploads');

// Map token strings to byte arrays to round-trip uploaded byte_values.
ord.uploads.tokenBytes = {};

// Map token strings to file system paths for new uploads.
ord.uploads.tokenFiles = {};

// Return a short random string of ASCII characters.
ord.uploads.newToken = function () {
  return 'upload_' + Math.random().toString(36).substring(2);
};

// Tag this file for upload and return a random token to access it later.
ord.uploads.newFile = function (file) {
  const token = ord.uploads.newToken();
  ord.uploads.tokenFiles[token] = file;
  return token;
};

// Get back a file that was previously recorded at newFile().
ord.uploads.getFile = function (token) {
  return ord.uploads.tokenFiles[token];
};

// Remember bytes at load() so they can be restored at unload().
ord.uploads.stashUpload = function (bytesValue) {
  const token = ord.uploads.newToken();
  ord.uploads.tokenBytes[token] = bytesValue;
  return token;
};

// At unload, recall bytes that were stashed at load().
ord.uploads.unstashUpload = function (token) {
  return ord.uploads.tokenBytes[token];
};

// Sends all files referenced in tokenFiles to the server.
ord.uploads.putAll = function (fileName) {
  const tokens = Object.getOwnPropertyNames(ord.uploads.tokenFiles);
  tokens.forEach(token => {
    const file = ord.uploads.tokenFiles[token];
    const reader = new FileReader();
    reader.readAsBinaryString(file);
    reader.onload = (event) => {
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/dataset/proto/upload/' + fileName + '/' + token);
      const payload = event.target.result;
      xhr.send(payload);
    };
  });
};

// Look up the bytesValue of the given uploader and send back it as a download.
ord.uploads.retrieve = function (uploader) {
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
  const bytesValue = ord.uploads.tokenBytes[token];
  xhr.send(bytesValue);
};

// Configure behaviors of a ".uploader" DOM element.
ord.uploads.initialize = function (node) {
  $('.uploader_chooser_file', node).on('input', (event) => {
    const file = event.target.files[0];
    $('.uploader_file_name', node).show();
    $('.uploader_file_name', node).text(file.name);
    const token = ord.uploads.newFile(file);
    $('.uploader', node).attr('data-token', token);
    $('.uploader_file_retrieve', node).hide();
  });
};

ord.uploads.load = function (node, bytesValue) {
  const token = ord.uploads.stashUpload(bytesValue);
  $('.uploader', node).show();
  $('.uploader', node).attr('data-token', token);
  $('.uploader_chooser_button', node).text('Replace...');
  $('.uploader_file_name', node).hide();
  $('.uploader_file_retrieve', node).show();
};

ord.uploads.unload = function (node) {
  const token = $('.uploader', node).attr('data-token');
  const bytesValue = ord.uploads.unstashUpload(token);
  if (bytesValue) {
    // This is just a round-trip for bytesValue.
    return bytesValue;
  } else {
    // A new file has been chosen for upload.
    const bytes = (new TextEncoder()).encode(token);
    return bytes;
  }
};
