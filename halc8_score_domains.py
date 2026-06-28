"""
halc8_score_domains.py
======================
HalC8 domain scoring pipeline — standalone, reproducible version.

Standalone implementation of the HalC8 v2.2 domain-scoring rules.
The script separates scoring from sequence retrieval: it accepts a CSV
with pre-curated sequences and produces scored CSV output plus a run
manifest.

QUICK START
-----------
  # Example domain-scoring command:
  python halc8_score_domains.py \\
      --input  input_domains.csv \\
      --output scored_domains.csv \\
      --input-type domain \\
      --ss-mode none

  # Run the bundled test suite:
  python halc8_score_domains.py \\
      --input-type domain \\
      --ss-mode none

INPUT CSV
---------
  Required columns:
    accession   unique identifier (WP_xxx.1 or any string)
    sequence    amino acid sequence (single-letter; gaps/stop codons stripped)
  Optional columns:
    category    free-text label (e.g. "novel_canonical", "9-cys_variant")
    notes       free-text; passed verbatim to output

INPUT TYPES
-----------
  domain        sequence IS the scored region (pre-extracted domain).
                No re-extraction. Recommended for pre-extracted HMM domains.
  full_protein  full precursor; extract_mature_region() applied before scoring
                (local alignment to HalC8 reference, UniProt P83716).

SS MODES
--------
  chou_fasman   local Chou-Fasman propensity prediction (Chou & Fasman 1978).
                Zero external dependencies. Default for published runs.
                Requires chou_fasman_ss.py in the same directory.
  jpred4        Jpred4 REST API (University of Dundee; Q3≈82%).
                Requires internet access; requires pip install requests.
  none          SS module skipped; 0 pts; ss_nterm_ok=True (conservative lb).
                Use when comparing scores without SS component.

OUTPUTS
-------
  <output>.csv                one row per sequence, all scoring columns
  <output>.manifest.json      run parameters + software versions (reviewer-proof)

SCORING (0–9.0 pts — frozen v2.2 rules)
---------------------------------------------------------------------------
  Hydrophobicity zone 1  (aa  1–30)  → 1.5 pts  (Kyte-Doolittle w=9)
  Hydrophobicity zone 2  (aa 31–50)  → 1.5 pts
  Hydrophobicity zone 3  (aa  >50)   → 1.0 pts
  Secondary structure                → 2.0 pts
  Conserved motif D/E/N/Q...G...C    → 3.0 pts  (4→3.0, 3→2.0, 2→1.0, 1→0.5)
  ────────────────────────────────────────────────────────────────
  TOTAL max                          → 9.0 pts

HARD GATES (exclude independently of score)
--------------------------------------------
  HG1  z1_pass=False         zone 1 not hydrophobic  → EXCLUDED
  HG2  ss_nterm_ok=False     N-terminal SS incompatible with HalC8 → EXCLUDED
  HG3  motif_has_anchor=False motif lacks D/E/N/Q anchor → EXCLUDED

VERDICTS (if all hard gates pass)
----------------------------------
  MOST LIKELY  score ≥ 6.5  AND  motif = 4/4
  PROBABLE     score ≥ 3.5  AND  motif ≥ 2/4
  MARGINAL     score ≥ 3.0  OR  (motif ≥ 3/4  AND  score ≥ 2.5)
  EXCLUDED     everything else

DEPENDENCIES
------------
  pip install pandas
  pip install biopython   # required only for --input-type full_protein
  pip install requests    # required only for --ss-mode jpred4

CITATION
--------
  This script (halc8_score_domains.py) is the canonical scorer used to
  run the frozen v2.2 scoring rules
  on pre-extracted HalC8-family candidate domains. The
  script writes per-domain scores and rule-based prioritization labels.
  
  for full performance metrics. The frozen v2.2 thresholds and hard-gate
  rules implemented here are the source of all published verdicts.
"""

import sys
import os
import re
import json
import time
import argparse
import platform
from datetime import datetime, timezone

# ── Third-party imports ──────────────────────────────────────────────────────
pd = None
try:
    import pandas as pd
except ImportError:
    pd = None

def _require_pandas():
    if pd is None:
        print("ERROR: pandas not installed. Run: pip install pandas", file=sys.stderr)
        sys.exit(1)

# BioPython — only needed for --input-type full_protein
try:
    from Bio.Align import PairwiseAligner
    _BIOPYTHON_OK = True
    import Bio as _bio_module
    _BIOPYTHON_VERSION = getattr(_bio_module, '__version__', 'unknown')
except ImportError:
    _BIOPYTHON_OK = False
    _BIOPYTHON_VERSION = 'not installed'

# requests — only needed for --ss-mode jpred4
try:
    import requests as _requests
    _REQUESTS_OK = True
except ImportError:
    _REQUESTS_OK = False


# =====================================================================
# VERSION / CONSTANTS
# =====================================================================

PIPELINE_VERSION = "halc8_score_domains_v1.0"
SCORE_VERSION    = "v2.2"   # frozen v2.2 scoring rules

# Reference mature HalC8 — UniProt P83716, Natrinema sp. AS7092, 74 aa
# Used only in --input-type full_protein for mature-region extraction
HALC8_MATURE = (
    "DIDITGCSACKYAAGQVCTIGCSAAGGFICGLLGITIPVAGLSLGFVEIVCTVADESYGC"
    "DAVAKEACNRAGLC"
)

# Kyte-Doolittle hydrophobicity scale
KD_SCALE = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

HYDRO_THRESH = {
    'z1_mean_min':  0.30,
    'z1_neg_max':   3,
    'z2_mean_min':  0.80,
    'z2_mean_max':  3.00,
    'z3_mean_min': -1.50,
    'z3_mean_max':  0.80,
}

SCORE_WEIGHTS = {
    'hydro_z1':  1.5,
    'hydro_z2':  1.5,
    'hydro_z3':  1.0,
    'structure': 2.0,
    'motif':     3.0,
}

VERDICT_THRESH = {
    'most_likely_min_score': 6.5,
    'probable_min_score':    3.5,
    'marginal_min_score':    3.0,
    'marginal_alt_motif':    3,
    'marginal_alt_score':    2.5,
}

MOTIF_POINTS = {4: 3.0, 3: 2.0, 2: 1.0, 1: 0.5, 0: 0.0}
VALID_AAS = set('ACDEFGHIKLMNPQRSTVWY')


# =====================================================================
# INPUT VALIDATION — fails loudly on errors, warns on soft issues
# =====================================================================

def validate_input(df, input_type):
    """
    Validate input DataFrame. Hard errors → sys.exit(1). Soft issues → stderr.
    Returns a cleaned copy of the DataFrame.
    """
    errors = []

    # ── Required columns ─────────────────────────────────────────────
    for col in ('accession', 'sequence'):
        if col not in df.columns:
            errors.append(f"Missing required column: '{col}'")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # ── Duplicate accessions ──────────────────────────────────────────
    dupes = df['accession'].duplicated()
    if dupes.any():
        dup_list = df.loc[dupes, 'accession'].tolist()
        print(
            f"WARNING: Duplicate accessions found (keeping first occurrence): "
            f"{dup_list}",
            file=sys.stderr,
        )
        df = df.drop_duplicates(subset='accession', keep='first').reset_index(drop=True)

    # ── Per-sequence checks ───────────────────────────────────────────
    cleaned_seqs = []
    for _, row in df.iterrows():
        acc = str(row['accession']).strip()
        raw = str(row['sequence']).strip().upper()

        # Remove gaps, stop codons, whitespace
        seq = re.sub(r'[-*.\s]', '', raw)

        # Check for non-standard amino acids
        invalid = set(seq) - VALID_AAS
        if invalid:
            print(
                f"WARNING: {acc}: non-standard characters stripped: "
                f"{sorted(invalid)!r}",
                file=sys.stderr,
            )
            seq = ''.join(c for c in seq if c in VALID_AAS)

        if len(seq) == 0:
            print(f"ERROR: {acc}: sequence is empty after cleaning.", file=sys.stderr)
            sys.exit(1)

        if input_type == 'domain':
            if len(seq) < 30:
                print(
                    f"WARNING: {acc}: domain is very short ({len(seq)} aa). "
                    f"Hydrophobicity zone 2/3 scoring may be unreliable.",
                    file=sys.stderr,
                )
            if len(seq) > 150:
                print(
                    f"WARNING: {acc}: domain is {len(seq)} aa — unusually long for "
                    f"--input-type domain (expected ≤120 aa). If this is a full "
                    f"precursor protein, use --input-type full_protein instead.",
                    file=sys.stderr,
                )

        cleaned_seqs.append(seq)

    df = df.copy()
    df['sequence'] = cleaned_seqs
    return df


# =====================================================================
# MATURE REGION EXTRACTION  (full_protein mode only)
# =====================================================================

def extract_mature_region(full_seq, reference=HALC8_MATURE, min_len=40):
    """
    Extract HalC8-homologous mature region from a full precursor protein
    via local pairwise alignment to the P83716 reference.
    Only called when --input-type full_protein.
    """
    if not _BIOPYTHON_OK:
        print(
            "ERROR: --input-type full_protein requires BioPython.\n"
            "       Run:  pip install biopython",
            file=sys.stderr,
        )
        sys.exit(1)

    if len(full_seq) <= len(reference) + 20:
        return full_seq

    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    try:
        alignments = aligner.align(full_seq, reference)
        best = alignments[0]
        coords = best.aligned[0]
        if coords:
            start  = coords[0][0]
            end    = coords[-1][1]
            region = full_seq[start:end]
            return region if len(region) >= min_len else full_seq
    except Exception:
        pass
    return full_seq


# =====================================================================
# SCORING MODULES
# =====================================================================

def _kd_profile(seq, window=9):
    """Kyte-Doolittle sliding-window hydrophobicity profile."""
    half   = window // 2
    scores = []
    for i in range(len(seq)):
        start = max(0, i - half)
        end   = min(len(seq), i + half + 1)
        vals  = [KD_SCALE.get(aa, 0.0) for aa in seq[start:end]]
        scores.append(sum(vals) / len(vals))
    return scores


def score_hydrophobicity(seq):
    """
    Evaluate three Kyte-Doolittle hydrophobicity zones.
    Returns (p_z1, p_z2, p_z3, z1_pass_bool, z1_mean, z2_mean, z3_mean, detail_str).
    """
    profile = _kd_profile(seq)
    n       = len(profile)

    z1      = profile[:min(30, n)]
    z1_mean = sum(z1) / len(z1) if z1 else 0.0
    z1_neg  = sum(1 for v in z1 if v < -0.3)

    z2      = profile[30:min(50, n)]
    z2_mean = sum(z2) / len(z2) if z2 else 0.0

    z3      = profile[50:]
    z3_mean = sum(z3) / len(z3) if z3 else 0.0

    p_z1 = SCORE_WEIGHTS['hydro_z1'] if (
        z1_mean >= HYDRO_THRESH['z1_mean_min'] and
        z1_neg  <= HYDRO_THRESH['z1_neg_max']
    ) else 0.0

    p_z2 = SCORE_WEIGHTS['hydro_z2'] if (
        HYDRO_THRESH['z2_mean_min'] <= z2_mean <= HYDRO_THRESH['z2_mean_max']
    ) else 0.0

    p_z3 = SCORE_WEIGHTS['hydro_z3'] if (
        HYDRO_THRESH['z3_mean_min'] <= z3_mean <= HYDRO_THRESH['z3_mean_max']
    ) else 0.0

    z1_pass = (p_z1 > 0.0)

    detail = (
        f"z1_mean={z1_mean:+.3f} ({z1_neg} neg peaks, thresh≥{HYDRO_THRESH['z1_mean_min']}); "
        f"z2_mean={z2_mean:+.3f}; "
        f"z3_mean={z3_mean:+.3f}"
    )
    return p_z1, p_z2, p_z3, z1_pass, z1_mean, z2_mean, z3_mean, detail


def score_conserved_motif(seq):
    """
    Search for the conserved D/E/N/Q...G...C motif in the first 30 aa.
    Scoring:
      +1  D/E/N/Q in positions 0–4  (s1 — N-terminal polar anchor)
      +1  G at positions 3–14       (s2 — invariant glycine)
      +1  C within 5 aa after G     (s3 — invariant cysteine)
      +1  A/V/L/I immediately before C  (s4 — small hydrophobic before Cys)
    Returns (motif_score_0_4, pts_float, pattern_str, g_pos_1indexed, s1_bool).
    """
    region     = seq[:min(30, len(seq))]
    best_score = 0
    best_pts   = 0.0
    best_patt  = "not found"
    best_gpos  = -1

    s1 = 1 if any(region[j] in 'DENQ' for j in range(min(5, len(region)))) else 0
    denq_char = next(
        (region[j] for j in range(min(5, len(region))) if region[j] in 'DENQ'), 'X'
    )

    for g_idx in range(3, min(15, len(region))):
        if region[g_idx] != 'G':
            continue

        c_idx = next(
            (j for j in range(g_idx + 1, min(g_idx + 6, len(region)))
             if region[j] == 'C'),
            None,
        )
        if c_idx is None:
            continue

        s4 = 1 if (c_idx > 0 and region[c_idx - 1] in 'AVLI') else 0
        avli_char = region[c_idx - 1] if s4 else 'X'
        score = s1 + 1 + 1 + s4   # s1 + s2(G) + s3(C) + s4
        gap   = 'X' * (c_idx - g_idx - 1)
        patt  = f"{denq_char}...G{gap}{avli_char}C"

        if score > best_score:
            best_score = score
            best_pts   = MOTIF_POINTS[score]
            best_patt  = patt
            best_gpos  = g_idx + 1

    return best_score, best_pts, best_patt, best_gpos, s1


def score_cysteines(seq):
    """Cysteine count and positions — informative variable, not scored."""
    positions = [i + 1 for i, aa in enumerate(seq) if aa == 'C']
    pos_str   = ','.join(str(p) for p in positions[:12])
    if len(positions) > 12:
        pos_str += '...'
    return len(positions), pos_str if pos_str else 'none'


def net_charge(seq):
    """Net charge at pH 7.0 — informative variable, not scored."""
    pos = seq.count('K') + seq.count('R') + seq.count('H') * 0.1
    neg = seq.count('D') + seq.count('E')
    return round(pos - neg, 1)


# ── Secondary structure ───────────────────────────────────────────────────────

def _evaluate_ss_pattern(pred_hec):
    """
    Given a per-residue H/E/C string, decide pts and ss_nterm_ok.
    Split: first 65% = N-terminal, last 35% = C-terminal.

    N-terminal criterion: (strand + coil) > helix
    C-terminal criterion: helix ≥ 25% AND helix > strand AND helix ≥ coil×0.6

    Returns (pts_float, verdict_str, detail_str, ss_nterm_ok_bool).
    """
    n     = len(pred_hec)
    if n < 20:
        return (
            0.0,
            "Sequence too short for SS evaluation",
            f"n={n} aa (< 20 aa minimum)",
            True,
        )

    split  = min(54, int(n * 0.65))
    nterm  = pred_hec[:split]
    cterm  = pred_hec[split:]

    n_h = nterm.count('H') / len(nterm)
    n_e = nterm.count('E') / len(nterm)
    n_c = nterm.count('C') / len(nterm)
    c_h = cterm.count('H') / len(cterm) if cterm else 0.0
    c_e = cterm.count('E') / len(cterm) if cterm else 0.0
    c_c = cterm.count('C') / len(cterm) if cterm else 0.0

    nterm_ok = (n_e + n_c) > n_h
    cterm_ok = (c_h >= 0.25 and c_h > c_e and c_h >= c_c * 0.6)

    detail = (
        f"split={split} aa | "
        f"Nterm: H={n_h:.2f} E={n_e:.2f} C={n_c:.2f} | "
        f"Cterm: H={c_h:.2f} E={c_e:.2f} C={c_c:.2f}"
    )

    if nterm_ok and cterm_ok:
        return SCORE_WEIGHTS['structure'], "Agrees", detail, True
    elif nterm_ok:
        return SCORE_WEIGHTS['structure'] * 0.5, "Agrees at N-term, no alpha-helix at C-term", detail, True
    elif cterm_ok:
        return SCORE_WEIGHTS['structure'] * 0.25, "Alpha-helix at C-term but N-term disagrees", detail, False
    else:
        return 0.0, "Disagrees at both ends", detail, False


def _predict_jpred4(seq, timeout_submit=30, timeout_poll=300):
    """
    Submit to Jpred4 REST API and return per-residue H/E/C string, or None.
    Requires pip install requests.
    """
    if not _REQUESTS_OK:
        return None

    SUBMIT_URL  = "https://www.compbio.dundee.ac.uk/jpred4/cgi-bin/rest/job"
    RESULTS_BASE = "https://www.compbio.dundee.ac.uk/jpred4/results"

    fasta = f">halc8_query\n{seq}\n"
    location = None

    for post_kwargs in [
        {'data': {'skipPDB': '1'}, 'files': {'seq': ('q.fasta', fasta.encode(), 'text/plain')}},
        {'data': {'format': 'single', 'seq': fasta}},
    ]:
        try:
            r = _requests.post(SUBMIT_URL, **post_kwargs,
                               timeout=timeout_submit, allow_redirects=False)
            if r.status_code in (201, 202):
                location = r.headers.get('Location', '')
                break
        except Exception:
            continue

    if not location:
        return None

    job_id  = location.rstrip('/').split('/')[-1]
    job_url = f"{SUBMIT_URL}/{job_id}"

    deadline = time.time() + timeout_poll
    while time.time() < deadline:
        time.sleep(15)
        try:
            poll = _requests.get(job_url, timeout=20)
        except Exception:
            continue

        if poll.status_code == 200:
            jnet_url = f"{RESULTS_BASE}/{job_id}/{job_id}.jnet"
            try:
                jnet = _requests.get(jnet_url, timeout=20)
                if jnet.status_code == 200:
                    pred = []
                    for line in jnet.text.splitlines():
                        if line.startswith('jnetpred:'):
                            raw = line.split(':', 1)[1]
                            for ch in raw:
                                if ch == 'H':
                                    pred.append('H')
                                elif ch == 'E':
                                    pred.append('E')
                                elif ch in ('-', 'C', ' '):
                                    pred.append('C')
                    return ''.join(pred) if pred else None
            except Exception:
                pass
            return None

        if poll.status_code == 410:   # job failed
            return None

    return None   # timeout


def score_secondary_structure(seq, ss_mode):
    """
    Score secondary structure.
    Returns (pts_float, verdict_str, detail_str, method_str, ss_nterm_ok_bool).
    """
    if ss_mode == 'none':
        detail = (
            "SS module skipped (--ss-mode none). "
            "0 pts assigned; ss_nterm_ok=True (conservative lower bound). "
            "HG2 gate not triggered."
        )
        return 0.0, "NOT_SCORED", detail, "none", True

    if ss_mode == 'chou_fasman':
        try:
            from chou_fasman_ss import predict_halc8_ss
        except ImportError:
            print(
                "ERROR: chou_fasman_ss.py not found.\n"
                "       Place chou_fasman_ss.py in the same directory as "
                "halc8_score_domains.py.",
                file=sys.stderr,
            )
            sys.exit(1)

        cf = predict_halc8_ss(seq)
        return (
            cf['pts_structure'],
            cf['secondary_structure'],
            cf['sec_detail'],
            'chou_fasman_1978',
            cf['ss_nterm_ok'],
        )

    if ss_mode == 'jpred4':
        if not _REQUESTS_OK:
            print(
                "ERROR: --ss-mode jpred4 requires the requests library.\n"
                "       Run:  pip install requests",
                file=sys.stderr,
            )
            sys.exit(1)

        print(f"    [Jpred4] submitting {seq[:12]}... (may take up to 5 min)")
        pred = _predict_jpred4(seq)
        if pred and len(pred) >= 20:
            pts, verdict, detail, ss_nterm_ok = _evaluate_ss_pattern(pred)
            detail = "[Jpred4] " + detail
            return pts, verdict, detail, "Jpred4", ss_nterm_ok
        else:
            detail = (
                "Jpred4 API unavailable or returned empty prediction. "
                "0 pts assigned; ss_nterm_ok=True (conservative)."
            )
            return 0.0, "NOT_PREDICTED (Jpred4 failed)", detail, "jpred4_failed", True

    raise ValueError(f"Unknown ss_mode: {ss_mode!r}")


# ── Verdict ───────────────────────────────────────────────────────────────────

def determine_verdict(total_score, motif_score_0_4,
                      z1_pass, ss_nterm_ok, motif_has_anchor):
    """
    Apply hard gates and continuous rules to assign a verdict.
    Implements the frozen v2.2 verdict logic used for all published results.
    """
    # Hard gates
    if not z1_pass:
        return "EXCLUDED"
    if not ss_nterm_ok:
        return "EXCLUDED"
    if not motif_has_anchor:
        return "EXCLUDED"

    # Continuous rules
    t = VERDICT_THRESH
    if total_score >= t['most_likely_min_score'] and motif_score_0_4 == 4:
        return "MOST LIKELY"
    if total_score >= t['probable_min_score'] and motif_score_0_4 >= 2:
        return "PROBABLE"
    if (total_score >= t['marginal_min_score'] or
            (motif_score_0_4 >= t['marginal_alt_motif'] and
             total_score >= t['marginal_alt_score'])):
        return "MARGINAL"
    return "EXCLUDED"


# =====================================================================
# MAIN SCORING FUNCTION
# =====================================================================

def score_one(accession, seq, input_type, ss_mode, category='', notes=''):
    """
    Score a single sequence. seq must already be cleaned (upper-case, no gaps).
    Returns a dict with all output columns.
    """
    # Step 1: extract mature region if needed
    scored_seq = extract_mature_region(seq) if input_type == 'full_protein' else seq

    # Step 2: hydrophobicity
    (p_z1, p_z2, p_z3,
     z1_pass, z1_mean, z2_mean, z3_mean,
     hydro_detail) = score_hydrophobicity(scored_seq)

    # Step 3: secondary structure
    (p_struct, sec_verdict, sec_detail,
     ss_method, ss_nterm_ok) = score_secondary_structure(scored_seq, ss_mode)

    # Step 4: conserved motif
    motif_score, p_motif, motif_patt, motif_gpos, motif_s1 = score_conserved_motif(scored_seq)

    # Step 5: informative variables
    cys_count, cys_pos = score_cysteines(scored_seq)
    charge              = net_charge(scored_seq)

    # Step 6: hard gate flags
    motif_has_anchor = (motif_score == 0) or (motif_s1 == 1)

    # Step 7: composite score and verdict
    total   = round(p_z1 + p_z2 + p_z3 + p_struct + p_motif, 2)
    verdict = determine_verdict(total, motif_score, z1_pass, ss_nterm_ok, motif_has_anchor)

    return {
        # Identity
        'accession':          accession,
        'sequence_length':    len(scored_seq),
        'input_type':         input_type,
        'category':           category,
        'notes':              notes,
        # Informative (not scored)
        'cysteine_count':     cys_count,
        'cysteine_positions': cys_pos,
        'net_charge':         charge,
        # Hydrophobicity
        'z1_mean':            round(z1_mean, 3),
        'z2_mean':            round(z2_mean, 3),
        'z3_mean':            round(z3_mean, 3),
        'pts_hydro_z1':       p_z1,
        'pts_hydro_z2':       p_z2,
        'pts_hydro_z3':       p_z3,
        'hydro_detail':       hydro_detail,
        # Secondary structure
        'secondary_structure': sec_verdict,
        'ss_method':          ss_method,
        'pts_structure':      p_struct,
        'sec_detail':         sec_detail,
        # Conserved motif
        'conserved_motif':    motif_patt,
        'motif_score_4':      motif_score,
        'pts_motif':          p_motif,
        # Hard gate flags
        'hg_z1_pass':         z1_pass,
        'hg_ss_nterm_ok':     ss_nterm_ok,
        'hg_motif_anchor':    motif_has_anchor,
        # Final score
        'total_score_9':      total,
        'verdict_v2':         verdict,
    }


# =====================================================================
# MANIFEST
# =====================================================================

def build_manifest(args, n_input, results_df, t_start, t_end):
    """Build a reviewer-proof run manifest dict."""
    manifest = {
        'pipeline_name':        PIPELINE_VERSION,
        'score_version':        SCORE_VERSION,
        'run_start_utc':        t_start,
        'run_end_utc':          t_end,
        'python_version':       sys.version.split()[0],
        'platform':             platform.platform(),
        'pandas_version':       pd.__version__,
        'biopython_version':    _BIOPYTHON_VERSION,
        'input_file':           os.path.abspath(args.input),
        'output_file':          os.path.abspath(args.output),
        'input_type':           args.input_type,
        'ss_mode':              args.ss_mode,
        'n_input_sequences':    n_input,
        'n_scored_sequences':   len(results_df),
        'verdict_summary':      results_df['verdict_v2'].value_counts().to_dict(),
        'thresholds': {
            'most_likely_min_score': VERDICT_THRESH['most_likely_min_score'],
            'probable_min_score':    VERDICT_THRESH['probable_min_score'],
            'marginal_min_score':    VERDICT_THRESH['marginal_min_score'],
            'z1_mean_min':           HYDRO_THRESH['z1_mean_min'],
            'z1_neg_max':            HYDRO_THRESH['z1_neg_max'],
            'z2_mean_min':           HYDRO_THRESH['z2_mean_min'],
            'z2_mean_max':           HYDRO_THRESH['z2_mean_max'],
            'z3_mean_min':           HYDRO_THRESH['z3_mean_min'],
            'z3_mean_max':           HYDRO_THRESH['z3_mean_max'],
        },
        'score_weights': SCORE_WEIGHTS,
    }
    return manifest


# =====================================================================
# MAIN
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        prog='halc8_score_domains.py',
        description='HalC8 domain scoring — standalone reproducible version.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # Example domain-scoring command:
  python halc8_score_domains.py \\
      --input  input_domains.csv \\
      --output scored_domains.csv \\
      --input-type domain --ss-mode none

  # Full precursor proteins with mature-region extraction:
  python halc8_score_domains.py \\
      --input  full_proteins.csv \\
      --output full_proteins_scores.csv \\
      --input-type full_protein --ss-mode none

""",
    )
    parser.add_argument(
        '--input', required=True, metavar='FILE',
        help='Input CSV. Required columns: accession, sequence. '
             'Optional: category, notes.',
    )
    parser.add_argument(
        '--output', required=True, metavar='FILE',
        help='Output CSV for scoring results.',
    )
    parser.add_argument(
        '--input-type', default='domain',
        choices=['domain', 'full_protein'],
        help='domain: score sequence as-is. '
             'full_protein: extract mature region first via local alignment. '
             'Default: domain.',
    )
    parser.add_argument(
        '--ss-mode', default='chou_fasman',
        choices=['chou_fasman', 'jpred4', 'none'],
        help='SS prediction method. '
             'chou_fasman: local, reproducible (default). '
             'jpred4: Jpred4 REST API (requires internet + pip install requests). '
             'none: 0 pts, ss_nterm_ok=True (conservative lower bound).',
    )
    parser.add_argument(
        '--manifest', default=None, metavar='FILE',
        help='Path for run manifest JSON. '
             'Default: <output_basename>.manifest.json',
    )

    args = parser.parse_args()

    _require_pandas()

    t_start = datetime.now(timezone.utc).isoformat()

    print("=" * 65)
    print(f"  {PIPELINE_VERSION}  |  score rules: {SCORE_VERSION}")
    print(f"  input-type: {args.input_type}  |  ss-mode: {args.ss_mode}")
    print(f"  started: {t_start}")
    print("=" * 65)

    # ── Load ──────────────────────────────────────────────────────────
    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"\nLoading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  {len(df)} rows, columns: {list(df.columns)}")

    # Ensure optional columns
    for col in ('category', 'notes'):
        if col not in df.columns:
            df[col] = ''

    # ── Validate ──────────────────────────────────────────────────────
    print("\nValidating input...")
    df = validate_input(df, args.input_type)
    n_input = len(df)
    print(f"  ✓ {n_input} sequences passed validation")

    # ── Score ─────────────────────────────────────────────────────────
    print(f"\nScoring {n_input} sequences...")
    results = []
    for _, row in df.iterrows():
        acc = str(row['accession']).strip()
        seq = str(row['sequence']).strip()
        cat = str(row.get('category', '')).strip()
        nts = str(row.get('notes', '')).strip()

        r = score_one(acc, seq, args.input_type, args.ss_mode, cat, nts)
        results.append(r)
        print(
            f"  {acc:<22}  len={r['sequence_length']:3d} aa  "
            f"Cys={r['cysteine_count']:2d}  "
            f"score={r['total_score_9']:4.1f}/9.0  "
            f"{r['verdict_v2']}"
        )

    results_df = pd.DataFrame(results)

    # Sort by verdict then score
    order = {'MOST LIKELY': 0, 'PROBABLE': 1, 'MARGINAL': 2, 'EXCLUDED': 3}
    results_df['_ord'] = results_df['verdict_v2'].map(
        lambda v: next((s for k, s in order.items() if k in str(v)), 4)
    )
    results_df = (
        results_df
        .sort_values(['_ord', 'total_score_9'], ascending=[True, False])
        .drop('_ord', axis=1)
        .reset_index(drop=True)
    )

    # ── Write outputs ─────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    results_df.to_csv(args.output, index=False)
    print(f"\n✓ Scores → {args.output}  ({len(results_df)} rows)")

    print("\nVerdict summary:")
    for v, c in results_df['verdict_v2'].value_counts().items():
        print(f"  {v}: {c}")

    t_end      = datetime.now(timezone.utc).isoformat()
    manifest   = build_manifest(args, n_input, results_df, t_start, t_end)
    mpath      = args.manifest or (
        re.sub(r'\.csv$', '', args.output) + '.manifest.json'
    )
    with open(mpath, 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(f"✓ Manifest → {mpath}")
    print("=" * 65)


if __name__ == '__main__':
    main()
