# Rediscovering Known Resonances in Real CMS Collision Data

**Live report:** https://en970.github.io/cms-dimuon-spectrum/
**PDF (LaTeX, 6 pp):** [`reports/cms-dimuon-spectrum-report.pdf`](reports/cms-dimuon-spectrum-report.pdf)

**One line:** stream real 2012 CMS proton-proton collision data from the CERN Open Data Portal,
compute the dimuon invariant-mass spectrum across all 61.5 million events, and recover eight
known resonances (rho/omega, phi, J/psi, psi(2S), the Upsilon triplet, and the Z boson) to
better than 1% mass accuracy against the Particle Data Group values.

- **Data type:** Real detector data (not simulation) · **Compute:** CPU, HTTP streaming, no GPU
  needed · **Runtime:** ~93 minutes for the full 61.5M-event file · **Status:** complete, full
  dataset processed (100% coverage)

## Why this project

This project uses **the actual CERN Open Data Portal** and **real LHC collision data recorded by
CMS** -- the same data, format, and analysis pattern (`uproot`/`awkward`/`vector`, the
Scikit-HEP stack) that DESY's CMS group works with directly (DESY is one of CMS's largest member
institutes). It is a direct, independent reproduction of the philosophy behind ROOT's official
`df102_NanoAODDimuonAnalysis` tutorial, done with the Pythonic HEP stack.

## Data

Source: CMS **2012 DoubleMuParked** dataset (Run2012B+C, 8 TeV proton-proton collisions),
released as a pre-skimmed NanoAOD-Outreach-Tool ROOT file on CERN's public EOS storage:

```
https://opendata.cern.ch/eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root
```

- **61,540,413 events**, each with per-event arrays of muon `pt`, `eta`, `phi`, `mass`, `charge`.
- 2.24 GB, streamed directly over HTTP with `uproot` in ~300,000-event chunks -- **the raw file
  is never downloaded or persisted**; only the derived histogram (tens of KB) and a small event
  sample are kept.
- Public domain (CC0), citable via the CERN Open Data Portal.

## Method

For each event with two or more reconstructed muons, take the first two muon candidates as
stored, require opposite electric charge, and compute the invariant mass of the dimuon system by
Lorentz-vector addition (`vector` library, `Momentum4D` behavior on Awkward Arrays):

```
M = sqrt( (E1+E2)^2 - |p1_vec + p2_vec|^2 )
```

Masses are histogrammed on a log scale (2000 bins, 0.2-200 GeV) across the full dataset. This is
a direct reproduction of a classic HEP analysis pattern: known-mass resonances that decay to
mu+mu- (or that are misidentified combinatorial background near well-known masses) produce sharp
peaks wherever real muon pairs from that decay accumulate.

### Engineering notes (the remote HTTP stream is flaky)
The naive approach (`TTree.iterate()` over the whole file in one call) crashed with a
`ClientPayloadError` before finishing even the first chunk. The working version:
- reads **manual, small (300k-event) entry ranges** via `TTree.arrays(entry_start, entry_stop)`,
- **retries** a failed chunk up to 6 times with backoff and a fresh connection,
- **checkpoints** the running histogram after every chunk (`outputs/histogram_checkpoint.npz` +
  `progress.json`), and **resumes** from the checkpoint on restart.

In practice the background process was killed by the environment's execution-time limits several
times over the course of processing; each restart picked up automatically from the last
checkpoint with no data loss, until 100% of the file was processed.

## Results

| Quantity | Value |
|---|---|
| Events processed | 61,540,413 (100% of the file) |
| Opposite-charge dimuon pairs | 39,047,770 |
| Total wall-clock processing time | ~93 minutes (across several resumed runs) |

### Peak validation (measured vs. Particle Data Group)

Each peak's mass was located by a local search around the literature value followed by parabolic
interpolation between the peak bin and its two neighbours:

| Resonance | PDG mass (GeV) | Measured (GeV) | Residual | Residual (%) |
|---|---:|---:|---:|---:|
| rho/omega | 0.775 | 0.7813 | +0.0063 | +0.81% |
| phi | 1.019 | 1.0166 | -0.0024 | -0.24% |
| J/psi | 3.097 | 3.0938 | -0.0032 | -0.10% |
| psi(2S) | 3.686 | 3.6810 | -0.0050 | -0.13% |
| Upsilon(1S) | 9.460 | 9.4519 | -0.0081 | -0.09% |
| Upsilon(2S) | 10.023 | 10.0144 | -0.0086 | -0.09% |
| Upsilon(3S) | 10.355 | 10.3370 | -0.0180 | -0.17% |
| Z | 91.190 | 90.9217 | -0.2683 | -0.29% |

All eight resonances are recovered to **better than 1% mass accuracy** with no calibration
corrections applied -- a strong quantitative confirmation that the spectrum is genuine physics,
not an artefact.

See `outputs/dimuon_spectrum.png` for the full spectrum and `outputs/run_summary.json` /
`outputs/peak_validation.json` for the raw numbers.

## Caveats

- **No muon momentum-scale calibration.** Real CMS analyses apply per-event momentum corrections
  (Rochester corrections) to remove small detector-alignment biases; this is why the Z peak (most
  sensitive to high-pT scale effects) shows the largest residual (-0.29%). The uncorrected
  agreement is still sub-percent.
- **"First two muons" is a simplification.** A production analysis selects the highest-quality
  (often highest-pT) opposite-charge pair and applies muon identification/isolation quality cuts;
  this reproduction uses the muons as ordered in the file.
- **Educational/derived data.** This is the NanoAOD Outreach Tool's simplified format, not the
  full research-grade AOD/MiniAOD used in a real CMS publication.

## Reproduce

```bash
python3 src/dimuon_spectrum.py     # streams the full file; resumable via outputs/histogram_checkpoint.npz
python3 src/analyze_peaks.py       # quantitative peak-position validation against PDG values
```

Requires `uproot`, `awkward`, `vector`, `numpy`, `matplotlib` (all in the Scikit-HEP stack).

## Repository layout

```
index.html        formal technical report (open locally or via GitHub Pages)
assets/           figure used by the report
src/              analysis code
outputs/          result figures, histograms, and summary JSON
reports/          combined LaTeX technical report (PDF)
```

## References
- CERN Open Data Portal: https://opendata.cern.ch/
- CMS 2012 DoubleMuParked dataset record: https://opendata.cern.ch/record/17
- ROOT `df102_NanoAODDimuonAnalysis` tutorial (the reference analysis this reproduces): https://root.cern/doc/master/df102__NanoAODDimuonAnalysis_8py.html
- Particle Data Group (PDG), Review of Particle Physics: https://pdg.lbl.gov/
- uproot: https://github.com/scikit-hep/uproot5 · awkward: https://github.com/scikit-hep/awkward · vector: https://github.com/scikit-hep/vector

## Licence

MIT (see `LICENSE`). The analysed data are the property of the CMS Collaboration and CERN, used
under the CERN Open Data CC0 terms.
