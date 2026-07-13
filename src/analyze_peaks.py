"""
Quantitative validation: locate each resonance peak in the measured histogram
(local search + parabolic refinement around the known literature mass) and
compare to the PDG value. This turns the qualitative "the peaks are visible"
claim into a quantitative one.
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "outputs")

d = np.load(os.path.join(OUT, "histogram_final.npz"))
edges, counts = d["edges"], d["counts"].astype(np.float64)
centers = 0.5 * (edges[:-1] + edges[1:])

# (name, PDG mass GeV, search half-window in log10(mass))
RESONANCES = [
    ("rho/omega",   0.775, 0.03),
    ("phi",         1.019, 0.02),
    ("J/psi",       3.097, 0.01),
    ("psi(2S)",     3.686, 0.01),
    ("Upsilon(1S)", 9.460, 0.01),
    ("Upsilon(2S)", 10.023, 0.01),
    ("Upsilon(3S)", 10.355, 0.01),
    ("Z",           91.19, 0.03),
]

results = []
for name, m_pdg, hw in RESONANCES:
    lo, hi = 10 ** (np.log10(m_pdg) - hw), 10 ** (np.log10(m_pdg) + hw)
    mask = (centers > lo) & (centers < hi)
    idx_local = np.where(mask)[0]
    if len(idx_local) < 3:
        continue
    sub_c, sub_y = centers[idx_local], counts[idx_local]
    i_peak = np.argmax(sub_y)
    # parabolic refinement using the peak bin and its two neighbours
    if 0 < i_peak < len(sub_y) - 1:
        y0, y1, y2 = sub_y[i_peak - 1], sub_y[i_peak], sub_y[i_peak + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        delta = np.clip(delta, -1, 1)
        x0, x1, x2 = sub_c[i_peak - 1], sub_c[i_peak], sub_c[i_peak + 1]
        m_meas = x1 + delta * (x2 - x0) / 2
    else:
        m_meas = sub_c[i_peak]
    resid = m_meas - m_pdg
    resid_pct = 100 * resid / m_pdg
    results.append({
        "name": name, "pdg_mass_GeV": m_pdg, "measured_mass_GeV": round(float(m_meas), 5),
        "residual_GeV": round(float(resid), 5), "residual_pct": round(float(resid_pct), 3),
        "peak_bin_count": int(sub_y[i_peak]),
    })
    print(f"{name:14s}  PDG={m_pdg:8.4f} GeV   measured={m_meas:8.4f} GeV   "
          f"residual={resid:+.4f} GeV ({resid_pct:+.2f}%)")

with open(os.path.join(OUT, "peak_validation.json"), "w") as f:
    json.dump(results, f, indent=2)
