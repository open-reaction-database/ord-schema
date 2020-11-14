# Copyright 2020 Open Reaction Database Project Authors
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
"""Python API for the Open Reaction Database."""

import urllib.parse

import requests

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.interface import query
from ord_schema.proto import dataset_pb2

TARGET = 'https://client.open-reaction-database.org'
ORD_DATA_URL = 'https://github.com/Open-Reaction-Database/ord-data/raw/main/'


def fetch_dataset(dataset_id):
    """Fetches a single Dataset message.

    Datasets are first class objects in GitHub, so this method skips the
    client API and goes directly to the ord-data repository.

    Args:
        dataset_id: String dataset ID.

    Returns:
        Dataset message.

    Raises:
        RuntimeError: The dataset request failed.
        ValueError: The dataset ID is invalid.
    """
    if not validations.is_valid_dataset_id(dataset_id):
        raise ValueError(f'Invalid dataset ID: {dataset_id}')
    url = urllib.parse.urljoin(ORD_DATA_URL,
                               f'{message_helpers.id_filename(dataset_id)}.pb')
    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(
            f'Request {url} failed with status {response.status_code}')
    return dataset_pb2.Dataset.FromString(response.content)


class OrdClient:
    """Client for the Open Reaction Database."""

    def __init__(self, target=None):
        """Initializes the client.

        Args:
            target: Endpoint URL for client queries. Defaults to the public
                ORD client URL (defined in TARGET).
        """
        if not target:
            target = TARGET
        self._target = target

    def fetch_datasets(self, dataset_ids):
        """Fetches one or more Dataset messages.

        Args:
            dataset_ids: List of dataset IDs.

        Returns:
            List of Dataset messages.
        """
        # NOTE(kearnes): Return datasets in the same order as requested.
        datasets = {
            dataset_id: self.fetch_dataset(dataset_id)
            for dataset_id in set(dataset_ids)
        }
        return [datasets[dataset_id] for dataset_id in dataset_ids]

    def fetch_dataset(self, dataset_id):  # pylint: disable=no-self-use
        """Fetches a single Dataset message.

        Args:
            dataset_id: String dataset ID.

        Returns:
            Dataset message.
        """
        return fetch_dataset(dataset_id)

    def fetch_reactions(self, reaction_ids):
        """Fetches one or more Reaction messages.

        Args:
            reaction_ids: List of reaction IDs.

        Returns:
            List of Reaction messages.

        Raises:
            ValueError: A reaction ID is invalid.
        """
        for reaction_id in reaction_ids:
            if not validations.is_valid_reaction_id(reaction_id):
                raise ValueError(f'Invalid reaction ID: {reaction_id}')
        target = urllib.parse.urljoin(self._target, '/api/fetch_reactions')
        response = requests.post(target, json=list(set(reaction_ids)))
        dataset = dataset_pb2.Dataset.FromString(response.content)
        # NOTE(kearnes): Return reactions in the same order as requested.
        reactions = {
            reaction.reaction_id: reaction for reaction in dataset.reactions
        }
        return [reactions[reaction_id] for reaction_id in reaction_ids]

    def fetch_reaction(self, reaction_id):
        """Fetches a single Reaction message.

        Args:
            reaction_id: String reaction ID.

        Returns:
            Reaction message.
        """
        reactions = self.fetch_reactions([reaction_id])
        assert len(reactions) == 1
        return reactions[0]

    def query(self,
              reaction_ids=None,
              reaction_smarts=None,
              components=None,
              use_stereochemistry=None,
              similarity=None):
        """Executes a query against the Open Reaction Database.

        Args:
            reaction_ids: List of reaction IDs to fetch. This is provided for
                completeness, but fetch_reaction(s) should be preferred.
            reaction_smarts: Reaction SMARTS pattern.
            components: List of ComponentQuery instances.
            use_stereochemistry: Boolean whether to use stereochemistry when
                matching.
            similarity: Float similarity threshold for SIMILAR queries.

        Returns:
            A Dataset containing the matched reactions.
        """
        if components is not None:
            component_params = [
                component.get_params() for component in components
            ]
        else:
            component_params = None
        params = {
            'reaction_ids': ','.join(reaction_ids) if reaction_ids else None,
            'reaction_smarts': reaction_smarts,
            'component': component_params,
            'use_stereochemistry': use_stereochemistry,
            'similarity': similarity,
        }
        target = urllib.parse.urljoin(self._target, '/api/query')
        response = requests.get(target, params=params)
        return dataset_pb2.Dataset.FromString(response.content)


class ComponentQuery:
    """Client-side implementation of ReactionComponentPredicate."""

    _ALLOWED_SOURCES = list(
        query.ReactionComponentPredicate.SOURCE_TO_TABLE.keys())
    _MATCH_MODE = query.ReactionComponentPredicate.MatchMode

    def __init__(self, pattern, source, mode):
        """Initializes the query.

        Args:
            pattern: String SMILES or SMARTS pattern.
            source: String key into ReactionComponentPredicate.SOURCE_TO_TABLE.
            mode: String ReactionComponentPredicate.MatchMode value.
        """
        self._pattern = pattern
        if source not in self._ALLOWED_SOURCES:
            raise ValueError(f'source is not in {self._ALLOWED_SOURCES}')
        self._source = source
        try:
            self._MATCH_MODE.from_name(mode)
        except KeyError as error:
            raise ValueError(f'mode is not in {self._MATCH_MODE}') from error
        self._mode = mode

    def get_params(self):
        """Returns URL parameters for GET requests."""
        return ';'.join([self._pattern, self._source, self._mode])
