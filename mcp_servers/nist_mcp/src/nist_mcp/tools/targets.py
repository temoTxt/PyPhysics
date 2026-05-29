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
"""

_DOC = "Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md"

DIRAC_TARGETS = [
    {
        "quantity": "hydrogen 2P_3/2 - 2P_1/2 fine-structure interval",
        "value": 10969.13,
        "uncertainty": 0.10,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; precision spectroscopy)",
        "reference": _DOC,
    },
    {
        "quantity": "hydrogen 2S_1/2 - 2P_1/2 Lamb shift",
        "value": 1057.845,
        "uncertainty": 0.009,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; precision spectroscopy)",
        "reference": _DOC,
    },
    {
        "quantity": "hydrogen 1S hyperfine splitting (21 cm)",
        "value": 1420.405751768,
        "uncertainty": 0.000000002,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; hydrogen maser)",
        "reference": _DOC,
    },
    {
        "quantity": "helium 3P_0 - 3P_1 fine-structure interval",
        "value": 29616.952,
        "uncertainty": None,
        "unit": "MHz",
        "source": "curated (Bethe-Salpeter cross-comparison; precision spectroscopy)",
        "reference": _DOC,
    },
]


def list_dirac_targets() -> list[dict]:
    """Return the curated Dirac-exploration target values (copies, in shared schema)."""
    return [dict(t) for t in DIRAC_TARGETS]
