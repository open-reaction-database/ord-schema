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
"""Tests for ord_schema.units."""

import urllib.error
from unittest import mock

import pytest
from rdkit import Chem

from ord_schema import resolvers
from ord_schema.proto import reaction_pb2


class TestNameResolvers:
    def test_resolve_names(self):
        def roundtrip_smi(smi):
            return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

        message = reaction_pb2.Reaction()
        message.inputs["test"].components.add().identifiers.add(type="NAME", value="aspirin")
        assert resolvers.resolve_names(message)
        resolved_smi = roundtrip_smi(message.inputs["test"].components[0].identifiers[1].value)
        assert resolved_smi == roundtrip_smi("CC(=O)Oc1ccccc1C(O)=O")
        assert (
            message.inputs["test"].components[0].identifiers[1].type
            == reaction_pb2.CompoundIdentifier.CompoundIdentifierType.SMILES
        )
        assert "NAME resolved" in message.inputs["test"].components[0].identifiers[1].details


class TestInputResolvers:
    def test_input_resolve(self):
        def roundtrip_smi(smi):
            return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

        string = "10 g of THF"
        reaction_input = resolvers.resolve_input(string)
        assert len(reaction_input.components) == 1
        assert reaction_input.components[0].amount.mass == reaction_pb2.Mass(value=10, units="GRAM")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(type="NAME", value="THF")
        assert reaction_input.components[0].identifiers[1].type == reaction_pb2.CompoundIdentifier.SMILES
        assert roundtrip_smi(reaction_input.components[0].identifiers[1].value) == roundtrip_smi("C1COCC1")

        string = "100 mL of 5.0uM sodium hydroxide in water"
        reaction_input = resolvers.resolve_input(string)
        assert len(reaction_input.components) == 2
        assert reaction_input.components[0].amount.moles == reaction_pb2.Moles(value=500, units="NANOMOLE")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(
            type="NAME", value="sodium hydroxide"
        )
        assert reaction_input.components[0].identifiers[1].type == reaction_pb2.CompoundIdentifier.SMILES
        assert roundtrip_smi(reaction_input.components[0].identifiers[1].value) == roundtrip_smi("[Na+].[OH-]")
        assert reaction_input.components[1].amount.volume == reaction_pb2.Volume(value=100, units="MILLILITER")
        assert reaction_input.components[1].amount.volume_includes_solutes
        assert reaction_input.components[1].identifiers[0] == reaction_pb2.CompoundIdentifier(
            type="NAME", value="water"
        )
        assert reaction_input.components[1].identifiers[1].type == reaction_pb2.CompoundIdentifier.SMILES
        assert roundtrip_smi(reaction_input.components[1].identifiers[1].value) == roundtrip_smi("O")

    @pytest.mark.parametrize(
        "string,expected",
        (
            ("100 g of 5.0uM sodium hydroxide in water", "amount of solution must be a volume"),
            ("100 L of 5 grapes in water", "String did not match template"),
        ),
    )
    def test_input_resolve_should_fail(self, string, expected):
        with pytest.raises((KeyError, ValueError), match=expected):
            resolvers.resolve_input(string)


def _http_error(code: int, msg: str) -> urllib.error.HTTPError:
    # urllib's HTTPError stub declares hdrs/fp as non-Optional, but the runtime
    # accepts None for both — fine for tests where we never read them.
    return urllib.error.HTTPError("", code, msg, hdrs=None, fp=None)  # ty: ignore[invalid-argument-type]


class TestPubChemRetry:
    """Tests for the tenacity-driven retry on PubChem 503 ServerBusy."""

    @staticmethod
    def _ok_response() -> mock.MagicMock:
        response = mock.MagicMock()
        response.read.return_value = b"CC(=O)Oc1ccccc1C(=O)O\n"
        response.__enter__.return_value = response
        return response

    def test_retries_then_succeeds(self, monkeypatch):
        # Two 503s then a successful response — the decorator should swallow both
        # 503s without surfacing them and return the eventual SMILES.
        responses = [
            _http_error(503, "PUGREST.ServerBusy"),
            _http_error(503, "PUGREST.ServerBusy"),
            self._ok_response(),
        ]
        urlopen = mock.MagicMock(side_effect=responses)
        monkeypatch.setattr("ord_schema.resolvers.urllib.request.urlopen", urlopen)
        # Drop tenacity's wait so the test runs in milliseconds, not seconds.
        # tenacity attaches a Retrying object to the wrapped function as .retry;
        # ty's stubs don't model this, hence the ignore.
        monkeypatch.setattr(resolvers._pubchem_resolve.retry, "wait", lambda *_, **__: 0)  # ty: ignore[unresolved-attribute]
        assert resolvers._pubchem_resolve("name", "aspirin") == "CC(=O)Oc1ccccc1C(=O)O"
        assert urlopen.call_count == 3

    def test_does_not_retry_404(self, monkeypatch):
        # 404 means "no such compound" — must not be retried.
        urlopen = mock.MagicMock(side_effect=_http_error(404, "Not Found"))
        monkeypatch.setattr("ord_schema.resolvers.urllib.request.urlopen", urlopen)
        with pytest.raises(urllib.error.HTTPError):
            resolvers._pubchem_resolve("name", "definitely-not-a-compound")
        assert urlopen.call_count == 1
