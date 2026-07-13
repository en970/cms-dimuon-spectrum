"""
Rediscovering known resonances in REAL CMS proton-proton collision data.

Source: the CMS 2012 "DoubleMuParked" dataset (Run2012B+C, 8 TeV), released
through the CERN Open Data Portal (https://opendata.cern.ch) as a pre-skimmed
NanoAOD-Outreach-Tool file -- the same file used in ROOT's official
`df102_NanoAODDimuonAnalysis` tutorial. This is REAL detector data recorded
by the CMS experiment, not simulation.

Method: select events with two or more reconstructed muons, take the first
two muon candidates as stored, require opposite electric charge, and compute
the invariant mass of the dimuon system via Lorentz-vector addition,

    M = sqrt( (E1+E2)^2 - |vec p1 + vec p2|^2 ),

using the `vector` library's Awkward-Array integration. Histogramming this
mass across millions of real collision events reproduces the famous "dimuon
spectrum": narrow peaks at the known masses of the J/psi, psi(2S), the
Upsilon triplet, and the Z boson -- rediscovering known particles directly
from real CERN detector data.

Data is streamed over HTTP directly from CERN's EOS Open Data storage; only
the derived (tiny) histogram and summary are saved, not the ~2.24 GB raw file.
"""
from __future__ import annotations
import os, json, time
import numpy as np
import uproot
import awkward as ak
import vector
vector.register_awkward()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

URL = ("https://opendata.cern.ch/eos/opendata/cms/derived-data/"
       "AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root:Events")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
OUT  = os.path.join(PROJ, "outputs")
os.makedirs(OUT, exist_ok=True)

STEP = 300_000               # events per streamed HTTP chunk (small -> resilient to drops)
CHECKPOINT_EVERY = 1          # write a checkpoint every N chunks (background runs get killed
                               # unpredictably -- checkpoint every chunk to minimise rework)
ENTRY_STOP = None            # None = process the full file; set an int to cap for testing
MAX_RETRIES = 6               # per-chunk retry budget (the remote HTTP stream is flaky)
RETRY_SLEEP_S = 8             # base backoff between retries
CKPT_PATH_TEMPLATE = None    # set in main()
N_BINS = 2000
MASS_LO, MASS_HI = 0.2, 200.0   # GeV, log-spaced (as in the classic dimuon plot)

RESONANCES = [
    ("eta",         0.548),
    ("rho/omega",   0.775),
    ("phi",         1.019),
    ("J/psi",       3.097),
    ("psi(2S)",     3.686),
    ("Upsilon(1S)", 9.460),
    ("Upsilon(2S)", 10.023),
    ("Upsilon(3S)", 10.355),
    ("Z",           91.19),
]


def process_chunk(tree, start, stop, edges, sample_masses):
    """Read one entry range and return (n_events, n_pairs, hist_counts). Raises on failure
    so the caller can retry."""
    chunk = tree.arrays(
        ["nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge"],
        entry_start=start, entry_stop=stop, library="ak",
    )
    sel = chunk[chunk.nMuon >= 2]
    mu = ak.zip({
        "pt": sel.Muon_pt, "eta": sel.Muon_eta, "phi": sel.Muon_phi,
        "mass": sel.Muon_mass,
    }, with_name="Momentum4D")
    q = sel.Muon_charge
    m1, m2 = mu[:, 0], mu[:, 1]
    q1, q2 = q[:, 0], q[:, 1]
    opp = q1 * q2 < 0
    mass = (m1[opp] + m2[opp]).mass
    mass = ak.to_numpy(mass)
    mass = mass[np.isfinite(mass) & (mass > 0)]
    if len(sample_masses) < 20000:
        need = 20000 - len(sample_masses)
        sample_masses.extend(mass[:need].tolist())
    h, _ = np.histogram(mass, bins=edges)
    return len(chunk), len(mass), h


def main():
    t0 = time.time()
    edges = np.logspace(np.log10(MASS_LO), np.log10(MASS_HI), N_BINS + 1)
    counts = np.zeros(N_BINS, dtype=np.int64)

    # Resume from a checkpoint if one exists (the remote HTTP stream is flaky; this
    # lets a re-run pick up where a previous attempt was interrupted).
    ckpt_path = os.path.join(OUT, "histogram_checkpoint.npz")
    start_at = 0
    n_pairs = 0
    sample_masses = []
    if os.path.exists(ckpt_path):
        d = np.load(ckpt_path)
        counts = d["counts"].copy()
        start_at = int(d["n_seen"])
        n_pairs = int(d["n_pairs"])
        print(f"Resuming from checkpoint: {start_at:,} events already processed", flush=True)

    tree = uproot.open(URL)
    n_total = tree.num_entries
    print(f"Opened remote file: {n_total:,} total events", flush=True)

    stop_at = ENTRY_STOP if ENTRY_STOP is not None else n_total
    n_seen = start_at
    n_chunks = 0

    start = start_at
    while start < stop_at:
        stop = min(start + STEP, stop_at)
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                n_ev, n_pr, h = process_chunk(tree, start, stop, edges, sample_masses)
                break
            except Exception as exc:  # noqa: BLE001 -- the remote stream fails in varied ways
                wait = RETRY_SLEEP_S * attempt
                print(f"  [retry {attempt}/{MAX_RETRIES}] chunk [{start:,}:{stop:,}) failed: "
                      f"{type(exc).__name__}: {exc} -- retrying in {wait}s", flush=True)
                if attempt == MAX_RETRIES:
                    raise
                time.sleep(wait)
                tree = uproot.open(URL)  # fresh connection for the retry

        n_seen += n_ev
        n_pairs += n_pr
        counts += h
        n_chunks += 1
        start = stop

        elapsed = time.time() - t0
        rate = (n_seen - start_at) / elapsed if elapsed > 0 else 0
        print(f"  chunk {n_chunks}: {n_seen:,}/{n_total:,} events "
              f"({100*n_seen/n_total:.1f}%), {n_pairs:,} opp-charge dimuon pairs, "
              f"{elapsed:.0f}s elapsed, {rate:,.0f} evt/s", flush=True)

        if n_chunks % CHECKPOINT_EVERY == 0:
            np.savez(ckpt_path, edges=edges, counts=counts,
                     n_seen=n_seen, n_pairs=n_pairs, n_total=n_total)
            with open(os.path.join(OUT, "progress.json"), "w") as fh:
                json.dump({"n_seen": n_seen, "n_total": n_total,
                          "n_pairs": n_pairs, "elapsed_s": round(elapsed, 1)}, fh, indent=2)

    np.savez(os.path.join(OUT, "histogram_final.npz"),
             edges=edges, counts=counts, n_seen=n_seen, n_pairs=n_pairs, n_total=n_total)
    with open(os.path.join(OUT, "sample_dimuon_masses.json"), "w") as fh:
        json.dump(sample_masses[:2000], fh)

    # ---- plot ----
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.step(centers, np.where(counts > 0, counts, np.nan), where="mid", color="C0", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    ax.set_ylabel("events / bin")
    ax.set_title(f"Dimuon invariant-mass spectrum -- real CMS 2012 collision data\n"
                 f"({n_pairs:,} opposite-charge pairs from {n_seen:,} events)", fontsize=11)
    ymax = counts.max()
    for name, m in RESONANCES:
        if MASS_LO < m < MASS_HI:
            ax.axvline(m, color="C3", alpha=0.3, lw=1)
            ax.text(m, ymax, name, rotation=90, fontsize=8, va="top", ha="right", color="C3")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "dimuon_spectrum.png"), dpi=140)
    plt.close(fig)

    summary = {
        "source": "CMS 2012 DoubleMuParked (Run2012B+C, 8 TeV), CERN Open Data Portal",
        "url": URL,
        "n_events_total_in_file": int(n_total),
        "n_events_processed": int(n_seen),
        "coverage_fraction": round(n_seen / n_total, 4),
        "n_opposite_charge_dimuon_pairs": int(n_pairs),
        "runtime_seconds": round(time.time() - t0, 1),
        "mass_range_GeV": [MASS_LO, MASS_HI],
        "n_bins": N_BINS,
        "resonances_expected_GeV": RESONANCES,
    }
    with open(os.path.join(OUT, "run_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\n=== SUMMARY ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
