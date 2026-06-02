"""Curated, individually-cited Dirac-exploration target values.

These are the MHz / sub-MHz precision splittings the Bethe-Salpeter campaign
compares the dual-Dirac predictions against. They are deliberately NOT fetched
from the NIST ASD general tables: ASD does not carry them at this precision, and
for hydrogen the ASD level energies are theoretical/Ritz values (using them to
"corroborate" a theory prediction would be circular). They are therefore curated
here, each transcribed with its source.

Provenance: the values are transcribed verbatim from the measurement column of
``Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md`` (which
cites the primary precision-spectroscopy literature). ``reference`` points back
to that document rather than to invented author-year citations; the electron
g-factor is delegated to ``get_constant`` (CODATA) and is intentionally absent
here.

Positronium and muonium are curated here too: they are exotic atoms absent from
the NIST ASD, and the PDG machine-readable API (the ``pdg`` package) carries no
muonium/positronium entries either, so their hyperfine splittings are only
available from PDG review text / primary literature.

Per the umbrella safety-flag enforcement, every record carries ``value_class``
+ derived ``safe_for_model_verification`` per ``precision_data_mcp.safety``.
The Lamb shift here is ``codata_adjusted`` (CODATA-2018 combines experiment +
theory); the other five are direct experimental measurements transcribed from
their primary-source maser / spectroscopy papers.
"""

from precision_data_mcp.safety import apply_safety_contract

_DOC = "Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md"

DIRAC_TARGETS = [
    {
        "quantity": "hydrogen 2P_3/2 - 2P_1/2 fine-structure interval",
        "value": 10969.13,
        "uncertainty": 0.10,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; precision spectroscopy)",
        "reference": _DOC,
        "value_class": "experimental",
    },
    {
        "quantity": "hydrogen 2S_1/2 - 2P_1/2 Lamb shift",
        "value": 1057.845,
        "uncertainty": 0.009,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; precision spectroscopy)",
        "reference": _DOC,
        "value_class": "codata_adjusted",
    },
    {
        "quantity": "hydrogen 1S hyperfine splitting (21 cm)",
        "value": 1420.405751768,
        "uncertainty": 0.000000002,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; hydrogen maser)",
        "reference": _DOC,
        "value_class": "experimental",
    },
    {
        "quantity": "helium 3P_0 - 3P_1 fine-structure interval",
        "value": 29616.952,
        "uncertainty": None,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; precision spectroscopy)",
        "reference": _DOC,
        "value_class": "experimental",
    },
    {
        "quantity": "positronium 1S ortho-para hyperfine splitting",
        "value": 203389.0,
        "uncertainty": 2.0,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; PDG / primary literature)",
        "reference": _DOC,
        "value_class": "experimental",
    },
    {
        "quantity": "muonium 1S hyperfine splitting",
        "value": 4463.302776,
        "uncertainty": 0.000051,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; PDG / primary literature)",
        "reference": _DOC,
        "value_class": "experimental",
    },
]


def list_dirac_targets() -> list[dict]:
    """Return the curated Dirac-exploration target values, each safety-decorated."""
    return [apply_safety_contract(dict(t)) for t in DIRAC_TARGETS]
