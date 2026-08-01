# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the natural-language evaluation harness scoring logic.

These cover the offline parts only -- case loading and interpretation scoring -- so no
model, network, or database is required.
"""

from ord_schema.agent.nl_query import NLComponent, NLQuery
from ord_schema.agent.nl_query_eval import (
    CaseExpectation,
    ComponentExpectation,
    check_interpretation,
    load_cases,
)


def test_load_cases_returns_nonempty():
    cases = load_cases()
    assert cases
    assert all(case.question for case in cases)


def test_check_interpretation_passes_on_match():
    expect = CaseExpectation(
        components=[
            ComponentExpectation(identifier="benzene", target="INPUT", mode="EXACT")
        ],
        min_yield=70,
    )
    interpretation = NLQuery(
        components=[NLComponent(identifier="Benzene", target="INPUT", mode="EXACT")],
        min_yield=70,
    )
    # Matching is case- and whitespace-insensitive but otherwise exact.
    assert check_interpretation(expect, interpretation) == []


def test_check_interpretation_requires_exact_identifier():
    # A less-specific identifier must NOT match: "aminophenol" != "4-aminophenol".
    expect = CaseExpectation(
        components=[
            ComponentExpectation(
                identifier="4-aminophenol", target="INPUT", mode="EXACT"
            )
        ]
    )
    interpretation = NLQuery(
        components=[NLComponent(identifier="aminophenol", target="INPUT", mode="EXACT")]
    )
    mismatches = check_interpretation(expect, interpretation)
    assert any("missing component" in m for m in mismatches)
    assert any("unexpected component" in m for m in mismatches)


def test_check_interpretation_matches_equivalent_smarts():
    # Equivalent SMARTS that differ as strings should match once canonicalized.
    expect = CaseExpectation(
        components=[
            ComponentExpectation(identifier="cB(O)O", target="OUTPUT", mode="SMARTS")
        ]
    )
    interpretation = NLQuery(
        components=[NLComponent(identifier="[c]B(O)O", target="OUTPUT", mode="SMARTS")]
    )
    assert check_interpretation(expect, interpretation) == []


def test_check_interpretation_matches_equivalent_smiles():
    # Equivalent SMILES (different atom ordering) should match once canonicalized.
    expect = CaseExpectation(
        components=[
            ComponentExpectation(
                identifier="c1ccncc1", target="OUTPUT", mode="SUBSTRUCTURE"
            )
        ]
    )
    interpretation = NLQuery(
        components=[
            NLComponent(identifier="n1ccccc1", target="OUTPUT", mode="SUBSTRUCTURE")
        ]
    )
    assert check_interpretation(expect, interpretation) == []


def test_check_interpretation_flags_wrong_target():
    expect = CaseExpectation(
        components=[
            ComponentExpectation(identifier="ibuprofen", target="OUTPUT", mode="EXACT")
        ]
    )
    interpretation = NLQuery(
        components=[NLComponent(identifier="ibuprofen", target="INPUT", mode="EXACT")]
    )
    mismatches = check_interpretation(expect, interpretation)
    assert any("missing component" in m for m in mismatches)
    assert any("unexpected component" in m for m in mismatches)


def test_check_interpretation_flags_over_extracted_yield():
    expect = CaseExpectation(
        components=[
            ComponentExpectation(identifier="benzene", target="INPUT", mode="EXACT")
        ]
    )
    interpretation = NLQuery(
        components=[NLComponent(identifier="benzene", target="INPUT", mode="EXACT")],
        min_yield=70,
    )
    mismatches = check_interpretation(expect, interpretation)
    assert mismatches == ["unexpected min_yield=70.0"]


def test_check_interpretation_flags_wrong_yield_value():
    expect = CaseExpectation(min_yield=90)
    interpretation = NLQuery(min_yield=70)
    mismatches = check_interpretation(expect, interpretation)
    assert mismatches == ["min_yield: expected 90.0, got 70.0"]
