"""Cross-namespace safety helper — value_class + model-verification safety flag.

Every value returned by a precision_data_mcp tool MUST carry a ``value_class``
field that places it on the experimental-vs-theoretical spectrum, and a derived
``safe_for_model_verification`` boolean. When the value is unsafe for model
verification, a ``warning`` field is auto-attached explaining why.

The umbrella REFRESH_POLICY.md mandates this contract — see §1 guarantee 6
(added in this commit).

## Taxonomy

| value_class | safe_for_model_verification | Meaning |
|---|---|---|
| ``experimental`` | True | Single-source direct experimental measurement |
| ``experimental_world_average`` | True | PDG / equivalent world-average of multiple experiments |
| ``codata_adjusted`` | **False** | CODATA least-squares adjustment combines experiment + theory; do NOT use as experimental anchor for verifying QED / SM predictions |
| ``theoretical_lattice_qcd`` | **False** | Lattice-QCD calculation (e.g. FLAG averages, BMW result); use as theory comparator only |
| ``theoretical_sm_prediction`` | **False** | Perturbative or data-driven SM prediction (e.g. Theory Initiative a_mu) |
| ``theoretical_ritz_value`` | **False** | Computed from Rydberg formula + theoretical quantum-defect parameters (e.g. NIST ASD hydrogen levels) |
| ``definitional`` | True | A conventional definition (e.g. proton charge = +1e, electron spin = 1/2 hbar) — safe to use, not a measurement test target |

When a tool decorates a raw record with the safety contract:

* Reads ``value_class`` from the record (must be present).
* Sets ``safe_for_model_verification`` per the table above.
* If unsafe, attaches a ``warning`` field with the standard message for that class.

Consumers (verification documents, agents) should refuse to use a value with
``safe_for_model_verification == False`` as the experimental anchor when comparing
theoretical predictions; treat such values as theory comparators only.

Created under the umbrella safety-flag enforcement PR (post umbrella issue
[#92](https://github.com/temoTxt/PyPhysics/issues/92) closure).
"""

from __future__ import annotations

_SAFE_CLASSES = {"experimental", "experimental_world_average", "definitional"}

_UNSAFE_CLASSES = {
    "codata_adjusted": (
        "CODATA least-squares adjustment combines experimental measurements with "
        "theoretical inputs (one-loop QED corrections, SM expressions for derived "
        "quantities). Do NOT use as the experimental anchor when verifying QED or "
        "SM predictions — use the underlying primary experimental measurements directly."
    ),
    "theoretical_lattice_qcd": (
        "Lattice-QCD theoretical calculation. Use as a theory comparator only; never "
        "as an experimental anchor for model verification."
    ),
    "theoretical_sm_prediction": (
        "Standard Model prediction (perturbative or data-driven). Use as a theory "
        "comparator only; never as an experimental anchor for model verification."
    ),
    "theoretical_ritz_value": (
        "Ritz value derived from theoretical quantum-defect parameters (typical of "
        "NIST ASD hydrogen levels). Not an independent measurement; do NOT use for "
        "model verification."
    ),
    "theoretical_calculation": (
        "Theoretical calculation result. Use as a theory comparator only; never as "
        "an experimental anchor for model verification."
    ),
}

ALL_VALUE_CLASSES = _SAFE_CLASSES | set(_UNSAFE_CLASSES.keys())


def apply_safety_contract(record: dict) -> dict:
    """Decorate a value record with safe_for_model_verification + warning fields.

    Reads the ``value_class`` field from the record; raises ``ValueError`` if missing
    or unknown. Returns a new dict with the contract fields added (does not mutate).
    """
    if "value_class" not in record:
        raise ValueError(
            f"safety contract violation: record missing required value_class field. "
            f"Every precision_data_mcp value must declare its class per safety.py taxonomy. "
            f"Record keys present: {sorted(record.keys())}"
        )
    vc = record["value_class"]
    if vc not in ALL_VALUE_CLASSES:
        raise ValueError(
            f"safety contract violation: unknown value_class {vc!r}. "
            f"Valid classes: {sorted(ALL_VALUE_CLASSES)}"
        )

    out = dict(record)
    out["safe_for_model_verification"] = vc in _SAFE_CLASSES
    if vc in _UNSAFE_CLASSES:
        out["warning"] = _UNSAFE_CLASSES[vc]
    return out
