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

import logging
import os
import urllib.error
from unittest import mock

import pytest
from rdkit import Chem

from ord_schema import resolvers
from ord_schema.proto import reaction_pb2

# The live PubChem/OPSIN smoke tests below run on a single CI matrix entry (set by
# the workflow) to avoid tripping PubChem's per-IP throttling; locally they always
# run. They already retry 503 ServerBusy via tenacity in resolvers.py -- we do not
# stack additional reruns on top of that. If PubChem is still busy after those
# retries, skip rather than fail so a third-party outage never reddens the build.
_live_resolvers = pytest.mark.skipif(
    os.environ.get("ORD_LIVE_RESOLVERS", "true") != "true",
    reason="live resolver smoke tests run on a single CI matrix entry",
)


def _skip_if_pubchem_unavailable(caplog):
    """Skip the calling test when PubChem signalled it was overloaded.

    A 503 PUGREST.ServerBusy (or other 5xx) means PubChem is rate-limiting or
    down, not that resolution logic is broken; the resolver logs it at INFO. A
    wrong-but-resolved answer still fails via the test's own assertions.
    """
    if "ServerBusy" in caplog.text or "HTTP Error 5" in caplog.text:
        pytest.skip("PubChem unavailable (503 ServerBusy); skipping live resolver smoke test")


class TestNameResolvers:
    @_live_resolvers
    def test_resolve_names(self, caplog):
        def roundtrip_smi(smi):
            return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

        message = reaction_pb2.Reaction()
        message.inputs["test"].components.add().identifiers.add(type="NAME", value="aspirin")
        with caplog.at_level(logging.INFO, logger="ord_schema.resolvers"):
            modified = resolvers.resolve_names(message)
        _skip_if_pubchem_unavailable(caplog)
        assert modified
        resolved_smi = roundtrip_smi(message.inputs["test"].components[0].identifiers[1].value)
        assert resolved_smi == roundtrip_smi("CC(=O)Oc1ccccc1C(O)=O")
        assert (
            message.inputs["test"].components[0].identifiers[1].type
            == reaction_pb2.CompoundIdentifier.CompoundIdentifierType.SMILES
        )
        assert "NAME resolved" in message.inputs["test"].components[0].identifiers[1].details

    def test_resolve_names_mocked(self, monkeypatch):
        # Hermetic counterpart to test_resolve_names: covers the resolve_names
        # plumbing (adds a SMILES identifier with details) without the live call,
        # so the logic stays covered if the live smoke test is ever removed.
        monkeypatch.setattr(
            "ord_schema.resolvers.name_resolve",
            mock.MagicMock(return_value=("CC(=O)Oc1ccccc1C(=O)O", "PubChem API")),
        )
        message = reaction_pb2.Reaction()
        message.inputs["test"].components.add().identifiers.add(type="NAME", value="aspirin")
        assert resolvers.resolve_names(message)
        identifier = message.inputs["test"].components[0].identifiers[1]
        assert identifier.type == reaction_pb2.CompoundIdentifier.SMILES
        assert identifier.value == "CC(=O)Oc1ccccc1C(=O)O"
        assert identifier.details == "NAME resolved by the PubChem API"


class TestInputResolvers:
    @_live_resolvers
    def test_input_resolve(self, caplog):
        def roundtrip_smi(smi):
            return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

        string = "10 g of THF"
        with caplog.at_level(logging.INFO, logger="ord_schema.resolvers"):
            reaction_input = resolvers.resolve_input(string)
        _skip_if_pubchem_unavailable(caplog)
        assert len(reaction_input.components) == 1
        assert reaction_input.components[0].amount.mass == reaction_pb2.Mass(value=10, units="GRAM")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(type="NAME", value="THF")
        assert reaction_input.components[0].identifiers[1].type == reaction_pb2.CompoundIdentifier.SMILES
        assert roundtrip_smi(reaction_input.components[0].identifiers[1].value) == roundtrip_smi("C1COCC1")

        string = "100 mL of 5.0uM sodium hydroxide in water"
        with caplog.at_level(logging.INFO, logger="ord_schema.resolvers"):
            reaction_input = resolvers.resolve_input(string)
        _skip_if_pubchem_unavailable(caplog)
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

    def test_input_resolve_mocked(self, monkeypatch):
        # Hermetic counterpart to test_input_resolve: exercises the real string
        # parsing and amount handling while mocking name resolution, so the logic
        # stays covered if the live smoke test is ever removed.
        smiles_by_name = {"THF": "C1COCC1", "sodium hydroxide": "[Na+].[OH-]", "water": "O"}
        monkeypatch.setattr(
            "ord_schema.resolvers.name_resolve",
            lambda _value_type, value: (smiles_by_name[value], "PubChem API"),
        )

        reaction_input = resolvers.resolve_input("10 g of THF")
        assert len(reaction_input.components) == 1
        assert reaction_input.components[0].amount.mass == reaction_pb2.Mass(value=10, units="GRAM")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(type="NAME", value="THF")
        assert reaction_input.components[0].identifiers[1].type == reaction_pb2.CompoundIdentifier.SMILES
        assert reaction_input.components[0].identifiers[1].value == "C1COCC1"

        reaction_input = resolvers.resolve_input("100 mL of 5.0uM sodium hydroxide in water")
        assert len(reaction_input.components) == 2
        assert reaction_input.components[0].amount.moles == reaction_pb2.Moles(value=500, units="NANOMOLE")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(
            type="NAME", value="sodium hydroxide"
        )
        assert reaction_input.components[0].identifiers[1].value == "[Na+].[OH-]"
        assert reaction_input.components[1].amount.volume == reaction_pb2.Volume(value=100, units="MILLILITER")
        assert reaction_input.components[1].amount.volume_includes_solutes
        assert reaction_input.components[1].identifiers[0] == reaction_pb2.CompoundIdentifier(
            type="NAME", value="water"
        )
        assert reaction_input.components[1].identifiers[1].value == "O"

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


class TestCanonicalizeSmiles:
    def test_canonicalizes(self):
        # OCCO and C(CO)O are the same molecule written two ways; both must
        # produce the canonical form "OCCO".
        assert resolvers.canonicalize_smiles("OCCO") == "OCCO"
        assert resolvers.canonicalize_smiles("C(CO)O") == "OCCO"

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Could not parse SMILES"):
            resolvers.canonicalize_smiles("not a smiles")


class TestOpsinResolve:
    def test_opsin_resolve(self, monkeypatch):
        response = mock.MagicMock()
        response.read.return_value = b"OCCO\n"
        response.__enter__.return_value = response
        urlopen = mock.MagicMock(return_value=response)
        monkeypatch.setattr("ord_schema.resolvers.urllib.request.urlopen", urlopen)
        # value_type is ignored by OPSIN.
        assert resolvers._opsin_resolve("name", "ethane-1,2-diol") == "OCCO"
        assert urlopen.call_count == 1


class TestNameResolveFallback:
    def test_falls_back_to_next_resolver_on_http_error(self, monkeypatch):
        # First resolver raises an HTTPError; second resolver succeeds.
        first = mock.MagicMock(side_effect=_http_error(404, "Not Found"))
        second = mock.MagicMock(return_value="OCCO")
        monkeypatch.setattr("ord_schema.resolvers._NAME_RESOLVERS", {"first": first, "second": second})
        smiles, resolver_name = resolvers.name_resolve("name", "ethylene glycol")
        assert smiles == "OCCO"
        assert resolver_name == "second"
        first.assert_called_once_with("name", "ethylene glycol")
        second.assert_called_once_with("name", "ethylene glycol")

    def test_raises_when_no_resolver_succeeds(self, monkeypatch):
        only = mock.MagicMock(side_effect=_http_error(404, "Not Found"))
        monkeypatch.setattr("ord_schema.resolvers._NAME_RESOLVERS", {"only": only})
        with pytest.raises(ValueError, match="Could not resolve"):
            resolvers.name_resolve("name", "definitely-not-a-compound")


class TestResolveNamesSkipsStructuralIdentifiers:
    def test_skips_compound_with_existing_smiles(self, monkeypatch):
        # If a Compound already has a structural identifier, name resolution
        # must not be attempted.
        called = mock.MagicMock()
        monkeypatch.setattr("ord_schema.resolvers.name_resolve", called)
        message = reaction_pb2.Reaction()
        compound = message.inputs["test"].components.add()
        compound.identifiers.add(type="NAME", value="aspirin")
        compound.identifiers.add(type="SMILES", value="CC(=O)Oc1ccccc1C(=O)O")
        assert not resolvers.resolve_names(message)
        called.assert_not_called()

    def test_swallows_value_error(self, monkeypatch):
        # When name_resolve raises ValueError, resolve_names skips that NAME and
        # reports no modification.
        monkeypatch.setattr(
            "ord_schema.resolvers.name_resolve",
            mock.MagicMock(side_effect=ValueError("not found")),
        )
        message = reaction_pb2.Reaction()
        message.inputs["test"].components.add().identifiers.add(type="NAME", value="zzzzzz")
        assert not resolvers.resolve_names(message)


class TestResolveInputBadFormat:
    def test_missing_of_separator(self):
        with pytest.raises(ValueError, match="does not match template"):
            resolvers.resolve_input("just a description")
