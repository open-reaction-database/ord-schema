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
"""Tests for ord_schema.scripts.list_dois."""

import os

from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.scripts import list_dois


class ValidateDatasetTest(absltest.TestCase):
    def test_simple(self):
        dataset = dataset_pb2.Dataset()
        dataset.dataset_id = "ord_dataset-1"
        dataset.reactions.add().provenance.doi = "foo/bar"
        dataset.reactions.add().provenance.doi = "doi:foo/bar"
        dataset.reactions.add().provenance.doi = "doi: foo/bar"
        dataset.reactions.add().provenance.doi = "DOI:foo/bar"
        dataset.reactions.add().provenance.doi = "DOI: \t foo/bar"
        tempdir = self.create_tempdir()
        message_helpers.write_message(dataset, os.path.join(tempdir, f"{dataset.dataset_id}.pb"))
        with flagsaver.flagsaver(input=os.path.join(tempdir, "*.pb")):
            list_dois.main(())

    def test_multiple_dois(self):
        dataset = dataset_pb2.Dataset()
        dataset.dataset_id = "ord_dataset-1"
        dataset.reactions.add().provenance.doi = "foo/bar"
        dataset.reactions.add().provenance.doi = "not/bar"
        tempdir = self.create_tempdir()
        message_helpers.write_message(dataset, os.path.join(tempdir, f"{dataset.dataset_id}.pb.gz"))
        with flagsaver.flagsaver(input=os.path.join(tempdir, "*.pb.gz")):
            list_dois.main(())


if __name__ == "__main__":
    absltest.main()
