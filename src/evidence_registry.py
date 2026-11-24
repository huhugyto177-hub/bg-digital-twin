"""Structured citation registry linking parameters to peer-reviewed literature.

The twin is calibrated to *directional* consensus across classical CBGT / PD /
dopaminergic pharmacology / QSP / LFP-beta sources. Triple-Control (Li & Jin)
is one entry among many — not the sole foundation.

Does NOT invent binding affinity, BBB permeability, or clinical efficacy for
novel chemical structures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

# Default registry (also written to configs/evidence_registry.yaml)
DEFAULT_EVIDENCE: dict[str, Any] = {
    "context_of_use": (
        "Mechanistic circuit-level translation of defined receptor/channel/pathway "
        "modulations and PK exposure curves to predicted basal ganglia activity, "
        "proxy LFP beta-band shifts, and motor-readiness readouts."
    ),
    "calibration_note": (
        "Calibrated to *directional* consensus from the literature listed below "
        "(pathway polarity, PD hypo-/hyperkinetic tendencies, L-DOPA restoration "
        "vs overdrive risk, D1 vs D2 engagement, STN–GPe beta associations). "
        "The twin does NOT invent binding affinity, BBB permeability, toxicity, "
        "or clinical efficacy for novel structures — those must be user-supplied "
        "assay/PK assumptions or left unknown."
    ),
    "non_claims": [
        "Does not predict de novo chemical binding affinity from SMILES/structures.",
        "Does not predict blood-brain barrier permeability.",
        "Does not predict peripheral ADMET/toxicity.",
        "Does not claim clinical efficacy without experimental validation.",
        "Does not invent EC50/Cmax/half-life for newly discovered molecules; user must supply mechanistic and PK assumptions.",
    ],
    "citations": [
        {
            "id": "albin_young_penney_1989",
            "parameter_keys": [
                "w_d1_gpi",
                "w_d2_gpe",
                "w_stn_gpi",
                "dopamine_tone",
                "pathway_balance",
            ],
            "citation": (
                "Albin RL, Young AB, Penney JB. The functional anatomy of basal "
                "ganglia disorders. Trends Neurosci. 1989;12(10):366–375."
            ),
            "assumption": (
                "Classical direct/indirect pathway polarity (D1→GPi inhibition; "
                "D2→GPe→STN→GPi) underpins the mean-field skeleton."
            ),
        },
        {
            "id": "delong_1990",
            "parameter_keys": ["w_d1_gpi", "w_stn_gpi", "w_gpe_stn", "dopamine_tone", "drive_gpi"],
            "citation": (
                "DeLong MR. Primate models of movement disorders of basal ganglia "
                "origin. Trends Neurosci. 1990;13(7):281–285."
            ),
            "assumption": (
                "Parkinsonian pathophysiology as altered rate balance across "
                "striatum–GPe–STN–GPi informs depleted-dopamine baseline."
            ),
        },
        {
            "id": "delong_wichmann_2007",
            "parameter_keys": ["w_d1_gpi", "w_stn_gpi", "w_gpe_stn", "dopamine_tone"],
            "citation": (
                "DeLong MR, Wichmann T. Circuits and circuit disorders of the basal "
                "ganglia. Arch Neurol. 2007;64(1):20–24."
            ),
            "assumption": "Updated classical direct/indirect/hyperdirect wiring used as mean-field skeleton.",
        },
        {
            "id": "wichmann_delong_2016",
            "parameter_keys": ["drive_stn", "w_stn_gpi", "w_gpe_stn", "beta_proxy", "gpi_snr"],
            "citation": (
                "Wichmann T, DeLong MR. Deep brain stimulation for movement "
                "disorders of basal ganglia origin: restoration of function or "
                "modulation of pathological activity? Neurotherapeutics. "
                "2016;13(2):264–283."
            ),
            "assumption": (
                "STN/GPi as key nodes whose pathological activity (including "
                "oscillatory patterns) tracks motor impairment."
            ),
        },
        {
            "id": "gerfen_surmeier_2011",
            "parameter_keys": ["d1_direct", "d2_indirect", "w_ctx_d1", "w_ctx_d2"],
            "citation": (
                "Gerfen CR, Surmeier DJ. Modulation of striatal projection systems "
                "by dopamine. Annu Rev Neurosci. 2011;34:441–466."
            ),
            "assumption": "D1 vs D2 MSN modulation differentially engages direct vs indirect pathways.",
        },
        {
            "id": "surmeier_2007_d1_d2",
            "parameter_keys": ["d1_direct", "d2_indirect", "dopamine_tone", "w_ctx_d1", "w_ctx_d2"],
            "citation": (
                "Surmeier DJ, Ding J, Day M, Wang Z, Shen W. D1 and D2 "
                "dopamine-receptor modulation of striatal glutamatergic signaling "
                "in striatal medium spiny neurons. Trends Neurosci. "
                "2007;30(5):228–235."
            ),
            "assumption": (
                "Opposite D1/D2 modulation of corticostriatal efficacy motivates "
                "separate d1_direct and d2_indirect gain channels."
            ),
        },
        {
            "id": "kravitz_2010_optogenetics",
            "parameter_keys": ["d1_direct", "d2_indirect", "motor_readiness", "pathway_balance"],
            "citation": (
                "Kravitz AV, et al. Regulation of parkinsonian motor behaviours by "
                "optogenetic control of basal ganglia circuitry. Nature. "
                "2010;466(7306):622–626."
            ),
            "assumption": (
                "Selective direct-pathway facilitation improves, and indirect "
                "facilitation worsens, Parkinsonian motor readouts — directional "
                "prior for pathway decomposition."
            ),
        },
        {
            "id": "nambu_2002_hyperdirect",
            "parameter_keys": ["w_ctx_stn", "stn_hyperdirect", "w_stn_gpi", "drive_stn"],
            "citation": (
                "Nambu A, Tokuno H, Takada M. Functional significance of the "
                "cortico–subthalamo–pallidal 'hyperdirect' pathway. Neurosci Res. "
                "2002;43(2):111–117."
            ),
            "assumption": "Hyperdirect cortex→STN→GPi route included as a distinct gain channel.",
        },
        {
            "id": "brown_2003_beta",
            "parameter_keys": ["beta_proxy", "drive_stn", "w_stn_gpi"],
            "citation": (
                "Brown P. Oscillatory nature of human basal ganglia activity: "
                "relationship to Parkinsonian pathophysiology. Mov Disord. "
                "2003;18(4):357–363."
            ),
            "assumption": (
                "Elevated beta-band activity associates with Parkinsonian circuit "
                "states; model uses a proxy, not a measured LFP spectrum."
            ),
        },
        {
            "id": "hammond_2007_beta_osc",
            "parameter_keys": ["beta_proxy", "w_gpe_stn", "w_stn_gpe", "drive_stn"],
            "citation": (
                "Hammond C, Bergman H, Brown P. Pathological synchronization in "
                "Parkinson's disease: networks, models and treatments. Trends "
                "Neurosci. 2007;30(7):357–364."
            ),
            "assumption": (
                "STN–GPe recurrent loops contribute to pathological beta "
                "synchronization; beta_proxy rises with STN/GPi drive imbalance."
            ),
        },
        {
            "id": "mallet_2008_stn_gpe_beta",
            "parameter_keys": ["beta_proxy", "w_gpe_stn", "w_stn_gpe", "gpe", "stn_hyperdirect"],
            "citation": (
                "Mallet N, et al. Parkinsonian beta oscillations in the external "
                "globus pallidus and their relationship with subthalamic nucleus "
                "activity. J Neurosci. 2008;28(52):14245–14258."
            ),
            "assumption": (
                "GPe–STN coupling strength influences Parkinsonian beta-related "
                "activity used as a directional prior for the LFP proxy."
            ),
        },
        {
            "id": "contin_levodopa_pk",
            "parameter_keys": ["dopamine_tone", "onset_h", "tmax_h", "half_life_h", "cmax_nM"],
            "citation": (
                "Contin M, Martinelli P. Pharmacokinetics of levodopa. J Neurol. "
                "2010;257(Suppl 2):S253–S261."
            ),
            "assumption": (
                "Empirical L-DOPA-like exposure envelopes (onset, tmax, short "
                "half-life) motivate default PK templates — not patient dosing advice."
            ),
        },
        {
            "id": "nutt_holford_1996_ldopa_pkpd",
            "parameter_keys": ["dopamine_tone", "motor_readiness", "overdrive_index", "ec50_nM"],
            "citation": (
                "Nutt JG, Holford NH. The response to levodopa in Parkinson's "
                "disease: imposing pharmacological law and order. Ann Neurol. "
                "1996;39(5):561–573."
            ),
            "assumption": (
                "Concentration–effect framing (EC50-like exposure–response) for "
                "motor benefit vs excess dopaminergic drive."
            ),
        },
        {
            "id": "obeso_2017_pd_essay",
            "parameter_keys": ["dopamine_tone", "motor_readiness", "overdrive_index"],
            "citation": (
                "Obeso JA, et al. Past, present, and future of Parkinson's disease: "
                "A special essay on the 200th Anniversary of the Shaking Palsy. "
                "Mov Disord. 2017;32(9):1264–1310."
            ),
            "assumption": (
                "Dopamine restoration can improve motor readiness but excessive "
                "drive may associate with dyskinesia risk proxies."
            ),
        },
        {
            "id": "calabresi_2010_dyskinesia",
            "parameter_keys": ["overdrive_index", "d1_direct", "dopamine_tone", "motor_readiness"],
            "citation": (
                "Calabresi P, Di Filippo M, Ghiglieri V, Tambasco N, Picconi B. "
                "Levodopa-induced dyskinesias in patients with Parkinson's "
                "disease: filling the bench-to-bedside gap. Lancet Neurol. "
                "2010;9(11):1106–1117."
            ),
            "assumption": (
                "Pulsatile / excessive dopaminergic stimulation links to "
                "dyskinesia risk; overdrive_index is a circuit proxy only."
            ),
        },
        {
            "id": "beste_qsp_pd_review",
            "parameter_keys": [
                "dopamine_tone",
                "d1_direct",
                "d2_indirect",
                "ec50_nM",
                "max_effect",
            ],
            "citation": (
                "Bakshi S, et al. Mathematical systems pharmacology and "
                "quantitative systems pharmacology approaches for Parkinson's "
                "disease (representative QSP framing). CPT Pharmacometrics Syst "
                "Pharmacol. / related QSP literature synthesizing DA–receptor–"
                "circuit response chains."
            ),
            "assumption": (
                "QSP-style separation of PK exposure → receptor occupancy → "
                "circuit effect motivates Hill/EC50 envelopes without claiming "
                "structure-based affinity prediction."
            ),
        },
        {
            "id": "li_jin_2023",
            "parameter_keys": [
                "d2_inverted_u_center",
                "d2_inverted_u_width",
                "pathway_balance",
                "w_ctx_d1",
                "w_ctx_d2",
            ],
            "citation": (
                "Li H, Jin X. Multiple pathways along the cortico-basal ganglia "
                "loops control distinct aspects of action. eLife. 2023. "
                "doi:10.7554/eLife.87644.3"
            ),
            "assumption": (
                "Triple-Control center–surround–context abstraction informs D2 "
                "inverted-U and pathway decomposition (one source among many)."
            ),
        },
    ],
}


def registry_path() -> Path:
    return Path(__file__).resolve().parents[1] / "configs" / "evidence_registry.yaml"


def ensure_registry_file() -> Path:
    path = registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(yaml.safe_dump(DEFAULT_EVIDENCE, sort_keys=False), encoding="utf-8")
    return path


def load_evidence_registry(path: Path | str | None = None) -> dict[str, Any]:
    p = Path(path) if path else ensure_registry_file()
    if not p.exists():
        return dict(DEFAULT_EVIDENCE)
    data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    # Merge defaults for missing top-level keys; prefer file citations when present
    out = dict(DEFAULT_EVIDENCE)
    out.update(data)
    if not out.get("citations"):
        out["citations"] = list(DEFAULT_EVIDENCE["citations"])
    return out


def evidence_table_rows(registry: dict[str, Any] | None = None) -> list[dict[str, str]]:
    reg = registry or load_evidence_registry()
    rows = []
    for item in reg.get("citations", []):
        rows.append(
            {
                "id": item.get("id", ""),
                "parameters": ", ".join(item.get("parameter_keys", [])),
                "citation": item.get("citation", ""),
                "assumption": item.get("assumption", ""),
            }
        )
    return rows


def references_markdown(registry: dict[str, Any] | None = None) -> str:
    """Full References / Literature Basis section for Decision Briefs and UI."""
    reg = registry or load_evidence_registry()
    lines = [
        "## References / Literature Basis",
        "",
        reg.get("calibration_note", DEFAULT_EVIDENCE["calibration_note"]),
        "",
        "| ID | Parameters | Citation | Assumption |",
        "|---|---|---|---|",
    ]
    for row in evidence_table_rows(reg):
        cite = row["citation"].replace("|", "\\|")
        assum = row["assumption"].replace("|", "\\|")
        lines.append(f"| {row['id']} | {row['parameters']} | {cite} | {assum} |")
    lines.append("")
    return "\n".join(lines)


def cou_and_limitations_markdown(registry: dict[str, Any] | None = None) -> str:
    reg = registry or load_evidence_registry()
    lines = [
        "## Context of Use & Limitations",
        "",
        f"**Context of Use:** {reg.get('context_of_use', DEFAULT_EVIDENCE['context_of_use'])}",
        "",
        f"**Calibration:** {reg.get('calibration_note', DEFAULT_EVIDENCE['calibration_note'])}",
        "",
        "**Non-claims / exclusions:**",
    ]
    for nc in reg.get("non_claims", DEFAULT_EVIDENCE["non_claims"]):
        lines.append(f"- {nc}")
    lines.append("")
    lines.append(
        "Every run is deterministic given its random seed; configuration JSON and "
        "P10–P90 uncertainty intervals should be retained with decision briefs."
    )
    return "\n".join(lines)
