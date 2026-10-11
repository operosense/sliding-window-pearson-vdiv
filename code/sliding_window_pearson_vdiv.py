# -*- coding: utf-8 -*-
"""
sliding_window_pearson_vdiv.py
==============================
Sliding-window Pearson correlation for operando battery data: detect the
divergence voltage(s) V_div at which the correlation between capacity Q(V) and a
mechanical signal (strain from a fibre/gauge sensor, or load-cell force) breaks down.

METHOD
    For each voltage V_i the Pearson correlation r_p between Q and the mechanical
    signal is computed over a centred voltage window of width dV, stepped one sample
    at a time. Ideal elastic coupling gives r_p ~= 1; a chemo-mechanical regime change
    makes r_p fall away from unity, often inverting toward -1. The r_p ~= 1 baseline is
    detected automatically (see `autodetect_baseline`). Contiguous troughs where
    r_p < tau = mu_base - K_SIGMA * sigma_base (persisting >= N_CONSEC windows and
    plunging to r_p <= DEEP_RP) are flagged, and V_div is reported at the r_p MINIMUM
    of each trough -- one V_div per inversion. A Fisher z-transform on non-overlapping
    bins provides a 95% significance flag for each V_div.

INPUT  (plain text; format configurable in USER SETTINGS below)
    One row per (voltage, mechanical, capacity) sample, already interpolated onto a
    common voltage grid. Defaults: TAB-separated, one header row, columns in the order
    Voltage, Mechanical, Capacity. Edit DELIM / HAS_HEADER / COL_* for other layouts.

OUTPUT (written to the ../results folder of this repository)
    <name>_rp.txt     voltage vs r_p trace
    <name>_Vdiv.txt   table: V_div, r_p_min, trough range, n_windows, 95% significance
    <name>_rp.png     r_p(V) plot with baseline, threshold and V_div markers

USAGE
    # process every matching file in ../data (INPUT_GLOB):
    python sliding_window_pearson_vdiv.py
    # or name specific files:
    python sliding_window_pearson_vdiv.py data_charge.txt data_discharge.txt
    In an IDE (e.g. Spyder) just press Run; it processes matching files in ../data.
    The functions below can also be imported and reused as a library.

PARAMETERS
    Defaults in USER SETTINGS reproduce the reference publication (see CITATION).
    Adjust them for other cells, sensors or sampling.

Requires : Python 3, numpy, matplotlib (matplotlib only for the optional plot).
License  : MIT (see LICENSE).
Author   : Quyen Do, Norwegian University of Science and Technology (NTNU)
Citation : Q. Do, H. Q. Nguyen, D. Rettenwander, S. T. Boles, "Decoupling the
           chemo-mechanical response during battery cycling using sliding-window
           Pearson correlation", J. Phys. Energy (2026).
           https://doi.org/10.1088/2515-7655/aeb06e
"""
import sys, os, glob
import numpy as np
# matplotlib is only needed for the plot. If it fails to import (e.g. a NumPy 1.x/2.x
# ABI mismatch in a given install), the numeric results and .txt tables are still
# written and only the .png plot is skipped.
try:
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception as _e:
    HAVE_MPL = False
    print("[note] matplotlib could not be imported (%s: %s)." % (type(_e).__name__, _e))
    print("       Results and .txt files are still written; the .png plot is skipped.")
    print('       To enable plots, fix the environment, e.g.:  pip install "numpy<2"')

# ============================== USER SETTINGS ==============================
# --- input file format ---
DELIM       = "\t"     # column delimiter: "\t", ",", or None for any whitespace
HAS_HEADER  = True     # True if the first line is a header row to skip
COL_V       = 0        # column index of cell voltage
COL_MECH    = 1        # column index of the mechanical signal (strain or force)
COL_CAP     = 2        # column index of capacity
INPUT_GLOB  = "*Good correlation.txt"  # files in DATA_DIR processed when none are named
                       # on the command line (this tool's own _rp.txt / _Vdiv.txt
                       # outputs are skipped automatically)
# --- analysis parameters (defaults follow the reference paper) ---
DV              = 0.02  # V   sliding-window width
MIN_PTS         = 10    #     minimum samples per window (windows with fewer are skipped)
BASELINE_N      = 6     #     baseline = first N consecutive windows with r_p >= BASELINE_RP_MIN
BASELINE_RP_MIN = 0.95  #     r_p threshold defining the near-unity baseline plateau
K_SIGMA         = 5     #     trough threshold tau = mu_base - K_SIGMA * sigma_base
N_CONSEC        = 5     #     a trough must persist at least this many consecutive windows
DEEP_RP         = -0.5  #     keep a trough as a V_div only if its minimum r_p plunges at
                        #     least this deep (toward anti-correlation, r_p -> -1);
                        #     raise toward 0 to also capture shallow, non-inverting dips
Z_CRIT          = 1.96  #     Fisher-z critical value (1.96 = two-sided 95%, 1.645 = one-sided)
# ==========================================================================

# Folder holding this script. Falls back to the current working directory when
# __file__ is undefined (e.g. running a console cell).
try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    HERE = os.getcwd()
DATA_DIR = os.path.join(HERE, os.pardir, "data")      # input files
OUT_DIR  = os.path.join(HERE, os.pardir, "results")   # generated outputs


def load_raw(path):
    """Read (voltage, mechanical, capacity) columns from a text file per USER SETTINGS."""
    V = []; M = []; Q = []
    need = max(COL_V, COL_MECH, COL_CAP)
    with open(path, encoding="utf-8", errors="replace") as f:
        if HAS_HEADER:
            next(f, None)
        for ln in f:
            ln = ln.replace("\r", "").rstrip("\n")
            if not ln.strip():
                continue
            c = ln.split(DELIM) if DELIM else ln.split()
            if len(c) <= need:
                continue
            try:
                V.append(float(c[COL_V])); M.append(float(c[COL_MECH])); Q.append(float(c[COL_CAP]))
            except ValueError:
                pass          # skip non-numeric / malformed rows
    return np.array(V), np.array(M), np.array(Q)


def sliding_rp(V, Q, M, dv=DV, min_pts=MIN_PTS):
    """Centred, one-sample-step Pearson r_p between capacity Q and mechanical signal M.
    Window membership is by |V - V_i| <= dv/2 (order-independent). r_p is evaluated in
    the record's native order, so the first samples are the start of the half-cycle."""
    order = np.argsort(V, kind="mergesort")      # sorted view, used only for fast lookup
    Vs, Qs, Ms = V[order], Q[order], M[order]
    half = dv / 2.0
    lo = np.searchsorted(Vs, V - half, side="left")
    hi = np.searchsorted(Vs, V + half, side="right")
    rp = np.full(len(V), np.nan)
    for i in range(len(V)):
        a, b = lo[i], hi[i]
        if b - a >= min_pts:
            q, m = Qs[a:b], Ms[a:b]
            if q.std() > 0 and m.std() > 0:
                rp[i] = np.corrcoef(q, m)[0, 1]
    return rp


def autodetect_baseline(rp, n=BASELINE_N, rp_min=BASELINE_RP_MIN):
    """Locate the baseline automatically: the first run of `n` consecutive windows with
    r_p >= `rp_min`. Returns (start_idx, end_idx, mu_base, sigma_base)."""
    ok = np.isfinite(rp) & (rp >= rp_min)
    for i in range(len(rp) - n + 1):
        if ok[i:i + n].all():
            seg = rp[i:i + n]
            return i, i + n, seg.mean(), seg.std()
    raise RuntimeError("No baseline plateau (>= %d consecutive windows with r_p >= %.2f)"
                       % (n, rp_min))


def detect_vdiv(V, rp):
    """V_div = the voltage at each r_p MINIMUM (anti-correlation trough, r_p -> -1).
    A trough is a contiguous excursion where r_p < tau (= mu_base - K_SIGMA*sigma_base)
    that persists >= N_CONSEC windows and whose minimum reaches r_p <= DEEP_RP. Within
    each trough, V_div is the voltage of the lowest r_p. Several troughs -> several V_div."""
    i0, i1, mu, sd = autodetect_baseline(rp)
    tau = mu - K_SIGMA * sd
    below = np.isfinite(rp) & (rp < tau)
    vdivs = []                                   # one dict per trough
    i, n = 0, len(rp)
    while i < n:
        if below[i]:
            j = i
            while j < n and below[j]:
                j += 1
            if (j - i) >= N_CONSEC:
                kmin = i + int(np.argmin(rp[i:j]))
                if rp[kmin] <= DEEP_RP:          # keep only genuine plunges toward -1
                    vdivs.append(dict(vdiv=V[kmin], rp_min=rp[kmin],
                                      v_lo=V[i], v_hi=V[j - 1], n=j - i))
            i = j
        else:
            i += 1
    return dict(i0=i0, i1=i1, mu=mu, sigma=sd, tau=tau, baseV=(V[i0], V[i1 - 1]), vdivs=vdivs)


def fisher_z_bins(V, Q, M, n=BASELINE_N, rp_min=BASELINE_RP_MIN):
    """Non-overlapping dV bins (statistically independent) for the Fisher-z stage.
    Bins are ordered from the start of the half-cycle (ascending V for charge,
    descending V for discharge). The baseline is the first `n` consecutive bins with
    r_p >= `rp_min`, and z_base = arctanh(mean r_p of those bins) (SI Section S7).
    Returns (edges_lo, edges_hi, r, npts, z_base); z_base is None if no baseline."""
    vmin, vmax = V.min(), V.max()
    edges = np.arange(vmin, vmax + DV, DV)
    lo, hi = edges[:-1], edges[1:]
    if V[0] > V[-1]:                                # discharge: start at high voltage
        lo, hi = lo[::-1], hi[::-1]
    r = np.full(len(lo), np.nan); npts = np.zeros(len(lo), dtype=int)
    for k in range(len(lo)):
        m = (V >= lo[k]) & (V < hi[k])
        npts[k] = m.sum()
        if npts[k] >= 4 and Q[m].std() > 0 and M[m].std() > 0:
            r[k] = np.corrcoef(Q[m], M[m])[0, 1]
    ok = np.isfinite(r) & (r >= rp_min)
    z_base = None
    for k in range(len(r) - n + 1):
        if ok[k:k + n].all():
            z_base = np.arctanh(np.clip(r[k:k + n].mean(), -0.999999, 0.999999))
            break
    return lo, hi, r, npts, z_base


def fisher_z_significance(bins, vdiv):
    """Significance of a V_div at the non-overlapping dV bin containing it:
    r_p,lower = tanh(z_base - Z_CRIT * SE), SE = 1/sqrt(n - 3); the departure is
    significant when the bin's r_p falls below r_p,lower (SI Section S7)."""
    lo, hi, r, npts, z_base = bins
    if vdiv is None or z_base is None:
        return None
    for k in range(len(lo)):
        if lo[k] <= vdiv < hi[k] or (k == np.argmax(hi) and vdiv == hi[k]):
            if not np.isfinite(r[k]) or npts[k] < 4:
                return None
            se = 1.0 / np.sqrt(npts[k] - 3)
            lb = np.tanh(z_base - Z_CRIT * se)
            return dict(bin=(lo[k], hi[k]), n=int(npts[k]), r=r[k], lower95=lb,
                        significant=bool(r[k] < lb))
    return None


def run(path):
    """Process one input file: compute r_p(V), detect V_div, and write the outputs."""
    name = os.path.splitext(os.path.basename(path))[0]
    V, M, Q = load_raw(path)
    rp = sliding_rp(V, Q, M)
    res = detect_vdiv(V, rp)
    vdivs = res["vdivs"]

    # --- write the r_p(V) trace ---
    os.makedirs(OUT_DIR, exist_ok=True)
    out_txt = os.path.join(OUT_DIR, name + "_rp.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("Voltage\tr_p\n")
        for v, r in zip(V, rp):
            f.write(f"{v:.4f}\t{'' if not np.isfinite(r) else f'{r:.6f}'}\n")

    # --- write the V_div table ---
    bins = fisher_z_bins(V, Q, M)
    for d in vdivs:
        d["sig"] = fisher_z_significance(bins, d["vdiv"])
    out_vdiv = os.path.join(OUT_DIR, name + "_Vdiv.txt")
    with open(out_vdiv, "w", encoding="utf-8") as f:
        f.write("V_div\trp_min\ttrough_V_lo\ttrough_V_hi\tn_windows\t"
                "bin_V_lo\tbin_V_hi\tbin_rp\tbin_rp_lower95\tsignificant_95\n")
        for d in vdivs:
            sig = d["sig"]
            if sig is None:
                cols = "\t\t\t\t"
            else:
                cols = (f"{sig['bin'][0]:.4f}\t{sig['bin'][1]:.4f}\t{sig['r']:.4f}\t"
                        f"{sig['lower95']:.4f}\t{'yes' if sig['significant'] else 'no'}")
            f.write(f"{d['vdiv']:.4f}\t{d['rp_min']:.4f}\t{d['v_lo']:.4f}\t"
                    f"{d['v_hi']:.4f}\t{d['n']}\t{cols}\n")

    # --- plot (only if matplotlib imported OK) ---
    out_png = "(plot skipped: matplotlib unavailable)"
    if HAVE_MPL:
        fig, ax = plt.subplots(figsize=(6.6, 3.9))
        good = np.isfinite(rp)
        ax.plot(V[good], rp[good], ".", ms=2.0, color="#0b3d64", label=r"$r_p(V)$")
        ax.axhline(res["tau"], ls="--", lw=1.0, color="#b00020",
                   label=fr"$\tau=\mu_{{base}}-{K_SIGMA}\sigma_{{base}}$ = {res['tau']:.3f}")
        ax.axvspan(res["baseV"][0], res["baseV"][1], color="#2a7a2a", alpha=0.12,
                   label=f"auto baseline {res['baseV'][0]:.2f}-{res['baseV'][1]:.2f} V")
        for k, d in enumerate(vdivs):
            ax.axvline(d["vdiv"], color="#c07d00", lw=1.5)
            ax.plot(d["vdiv"], d["rp_min"], "v", color="#c07d00", ms=7,
                    label=(r"$V_{div}$ (r$_p$ minima)" if k == 0 else None))
            ax.annotate(f"{d['vdiv']:.3f} V", (d["vdiv"], d["rp_min"]),
                        textcoords="offset points", xytext=(3, 10), fontsize=7, color="#7a4e00")
        ax.set_xlabel("Voltage / V"); ax.set_ylabel(r"$r_p$"); ax.set_ylim(-1.18, 1.18)
        ax.legend(fontsize=7, loc="center left"); ax.set_title(name, fontsize=9)
        out_png = os.path.join(OUT_DIR, name + "_rp.png")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")   # left open so it shows in an IDE

    print(f"\n=== {os.path.basename(path)} ===")
    print(f" rows={len(V)}  V {V[0]:.3f}->{V[-1]:.3f}  valid r_p windows={np.isfinite(rp).sum()}")
    print(f" auto baseline {res['baseV'][0]:.3f}-{res['baseV'][1]:.3f} V  "
          f"mu={res['mu']:.4f} sigma={res['sigma']:.4f}  tau={res['tau']:.4f}")
    print(f" V_div (r_p minima, one per trough, >= {N_CONSEC} windows, r_p<=DEEP_RP):")
    for d in vdivs:
        print(f"    V_div = {d['vdiv']:.3f} V   r_p,min = {d['rp_min']:+.3f}   "
              f"trough {d['v_lo']:.3f}-{d['v_hi']:.3f} V ({d['n']} windows)  "
              f"significant_95: "
              f"{'n/a' if d['sig'] is None else ('yes' if d['sig']['significant'] else 'no')}")
    if not vdivs:
        print("    (none)")
    print(f" wrote: {os.path.basename(out_txt)}, {os.path.basename(out_vdiv)}, "
          f"{os.path.basename(out_png)}")
    return res


def find_inputs():
    """Default input list: files matching INPUT_GLOB in DATA_DIR (or cwd),
    excluding this tool's own generated outputs."""
    hits = sorted(glob.glob(os.path.join(DATA_DIR, INPUT_GLOB)))
    if not hits:
        hits = sorted(glob.glob(os.path.join(os.getcwd(), INPUT_GLOB)))
    return [p for p in hits if not p.endswith(("_rp.txt", "_Vdiv.txt"))]


def main(files=None):
    if not files:
        files = find_inputs()
    if not files:
        print("No input files found (pattern %r)." % INPUT_GLOB)
        print("Looked in:", os.path.abspath(DATA_DIR))
        return
    print("Found %d input file(s):" % len(files))
    for p in files:
        print("   -", os.path.basename(p))
    for p in files:
        if not os.path.isabs(p) and not os.path.exists(p):
            p = os.path.join(DATA_DIR, p)
        run(p)
    if HAVE_MPL:
        plt.show()


# Runs on "Run file" in an IDE, or `python this_script.py [files...]` from a shell.
# Importing this module (e.g. `from sliding_window_pearson_vdiv import sliding_rp`)
# does not trigger a run.
if __name__ == "__main__":
    main(sys.argv[1:])
