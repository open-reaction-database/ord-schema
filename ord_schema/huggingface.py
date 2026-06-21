# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Fetches datasets from the Hugging Face ord-data mirror.

Importing this module requires the optional ``huggingface`` extra
(``pip install 'ord-schema[huggingface]'``); ``huggingface_hub`` is not a default
dependency because it pulls in a heavy chain (typer, fsspec, the hf-xet native
binary). Only this download helper needs it.
"""

from huggingface_hub import hf_hub_download
from huggingface_hub.errors import EntryNotFoundError

from ord_schema import validations
from ord_schema.message_helpers import id_filename

# Hugging Face dataset that mirrors the ord-data repository; see
# https://huggingface.co/datasets/open-reaction-database/ord-data.
ORD_DATA_HF_REPO = "open-reaction-database/ord-data"


def fetch_dataset(
    dataset_id: str,
    *,
    revision: str = "main",
    cache_dir: str | None = None,
    local_dir: str | None = None,
) -> str:
    """Downloads a dataset file from the ord-data Hugging Face mirror.

    Prefers the Parquet serialization and falls back to the legacy ``.pb.gz``
    binary format for datasets that have not been converted yet. The file is
    not parsed; dispatch on the returned suffix to choose a reader, e.g.
    ``parquet_dataset.read_dataset``/``parquet_dataset.DatasetView`` for
    ``.parquet`` and ``load_message`` for ``.pb.gz``.

    Downloads are cached and content-verified by ``huggingface_hub``: repeat
    calls reuse the cached file (re-validating against the remote for moving
    revisions like ``"main"``) instead of re-downloading. By default the cache
    lives under ``~/.cache/huggingface/hub``; set ``$HF_HOME`` (or pass
    ``cache_dir``) to relocate it, or set ``$HF_HUB_OFFLINE=1`` to reuse an
    already-warm cache without network access.

    Args:
        dataset_id: Dataset ID.
        revision: Branch, tag, or commit SHA to download (default: ``"main"``).
        cache_dir: Override the Hugging Face cache directory. Files are stored
            content-addressed and returned via a symlink into the cache.
        local_dir: If set, materialize the file as a real (non-symlink) copy
            under this directory instead of returning a cache symlink; useful
            for staging files into a project tree.

    Returns:
        Local filesystem path to the downloaded dataset file.

    Raises:
        RuntimeError: If no file exists for the dataset.
        ValueError: If the dataset ID is invalid.
    """
    if not validations.is_valid_dataset_id(dataset_id):
        raise ValueError(f"Invalid dataset ID: {dataset_id}")
    for suffix in (".parquet", ".pb.gz"):
        try:
            return hf_hub_download(
                repo_id=ORD_DATA_HF_REPO,
                repo_type="dataset",
                revision=revision,
                filename=id_filename(f"{dataset_id}{suffix}"),
                cache_dir=cache_dir,
                local_dir=local_dir,
            )
        except EntryNotFoundError:
            continue
    raise RuntimeError(
        f"Dataset {dataset_id} not found in {ORD_DATA_HF_REPO}@{revision}"
    )
