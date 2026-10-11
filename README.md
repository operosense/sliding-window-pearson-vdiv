# Sliding-window Pearson correlation for operando battery chemo-mechanics

Data and code supporting the paper:

> Q. Do, H. Q. Nguyen, D. Rettenwander, S. T. Boles,
> **"Decoupling the chemo-mechanical response during battery cycling using sliding-window Pearson correlation"**,
> *Journal of Physics: Energy* (2026). https://doi.org/10.1088/2515-7655/aeb06e

The method computes a voltage-resolved Pearson correlation, r_p(V), between operando
capacity Q(V) and a mechanical signal (fibre-Bragg-grating strain ε(V) or load-cell
force F(V)). When the two signals stay elastically coupled, r_p ≈ 1. A chemo-mechanical
regime change pulls r_p away from 1, often all the way to −1. The **divergence voltage
V_div** is the voltage of each r_p minimum. Its statistical significance is tested with a
Fisher z-transform on non-overlapping voltage bins.

The full workflow (the "Operational V_div criterion") is shown as a flowchart in
Figure 1b of the [paper](https://doi.org/10.1088/2515-7655/aeb06e). It is written out
step by step in [`docs/Algorithm1_Vdiv_pseudocode.md`](docs/Algorithm1_Vdiv_pseudocode.md).

## Repository contents

```
.
├── code/
│   └── sliding_window_pearson_vdiv.py   # pipeline: (V, ε, Q) → r_p(V) → V_div + significance
├── data/
│   ├── Ch2 Good correlation.txt         # Case 1, 2nd charge:    Voltage, Strain, Capacity
│   ├── DCh2 Good correlation.txt        # Case 1, 2nd discharge: Voltage, Strain, Capacity
│   ├── correlation_results Ch2 Good.txt # Case 1, 2nd charge: non-overlapping r_p(V) shown in Figure 2k
│   └── README.md                        # column definitions
├── docs/
│   └── Algorithm1_Vdiv_pseudocode.md    # Algorithm 1 (plain text + LaTeX)
├── results/                             # example outputs of the pipeline on the Case 1 data
├── requirements.txt
├── CITATION.cff
├── LICENSE                              # MIT (code)
└── LICENSE-DATA                         # CC BY 4.0 (data)
```

## Quick start

Requires Python ≥ 3.8 with `numpy` and `matplotlib`.

```bash
pip install -r requirements.txt
cd code
python sliding_window_pearson_vdiv.py        # processes data/*Good correlation.txt → results/
```

To analyse your own data, pass the file(s) explicitly:

```bash
python sliding_window_pearson_vdiv.py path/to/my_charge.txt path/to/my_discharge.txt
```

Input files are plain text with one header row and tab-separated columns
`Voltage, Mechanical (strain or force), Capacity`, already synchronised onto a common
voltage grid. You can change the delimiter, header and column order under
`USER SETTINGS` at the top of the script. The functions (`sliding_rp`,
`detect_vdiv`, `fisher_z_bins`, `fisher_z_significance`) can also be imported and used
as a library.

## Outputs

For each input `<name>.txt` the pipeline writes to `results/`:

| file | content |
|------|---------|
| `<name>_rp.txt`   | voltage vs. sliding-window r_p |
| `<name>_Vdiv.txt` | one row per V_div: r_p minimum, trough range, number of windows, Fisher-z bin, bin r_p, 95 % lower bound and significance flag |
| `<name>_rp.png`   | r_p(V) with baseline, threshold τ and V_div markers |

## Analysis parameters

These defaults reproduce the paper (Section 3.1 and SI Section S7):

| parameter | value | meaning |
|-----------|-------|---------|
| `DV` | 0.02 V | sliding-window width (centred, one-sample step) |
| `MIN_PTS` | 10 | minimum samples per window |
| `BASELINE_N`, `BASELINE_RP_MIN` | 6, 0.95 | baseline = first 6 consecutive windows with r_p ≥ 0.95 |
| `K_SIGMA` | 5 | threshold τ = μ_base − 5 σ_base |
| `N_CONSEC` | 5 | a trough must last at least 5 consecutive windows |
| `DEEP_RP` | −0.5 | the trough minimum must reach r_p ≤ −0.5 to count as a V_div |
| `Z_CRIT` | 1.96 | Fisher-z 95 % lower bound: r_p,lower = tanh(z_base − 1.96 / √(n − 3)) |

## Example results (Case 1, NMC811 all-solid-state cell, second cycle)

Running `sliding_window_pearson_vdiv.py` on the files in `data/` gives:

| half-cycle | V_div (r_p minima) | Fisher-z 95 % |
|------------|--------------------|---------------|
| charge     | 3.903 V, 3.982 V, 4.184 V | all significant |
| discharge  | 4.104 V, 3.969 V, 3.828 V, 2.726 V | all significant |

The three charge inversions correspond to the ≈3.9 V, ≈4.0 V and H2→H3 (≈4.16 V) features
discussed in the paper. The last charge trough and the last discharge trough are flat
r_p ≈ −1 plateaus that run to the cut-off voltage. For these, the "minimum" is
degenerate and its position along the plateau is not physically meaningful, so read them
together with the trough range written in `<name>_Vdiv.txt`.

## Citation

If you use this code or data, please cite the paper above. See also `CITATION.cff`.

```bibtex
@article{Do2026SlidingPearson,
  author  = {Do, Quyen and Nguyen, Hung Quoc and Rettenwander, Daniel and Boles, Steven T.},
  title   = {Decoupling the chemo-mechanical response during battery cycling using sliding-window Pearson correlation},
  journal = {Journal of Physics: Energy},
  year    = {2026},
  doi     = {10.1088/2515-7655/aeb06e}
}
```

## License

- Code (`code/`): MIT License, see [`LICENSE`](LICENSE).
- Data, results and documentation (`data/`, `results/`, `docs/`):
  Creative Commons Attribution 4.0 International (CC BY 4.0), see [`LICENSE-DATA`](LICENSE-DATA).

## Contact

Quyen Do, Department of Energy and Process Engineering, Norwegian University of Science
and Technology (NTNU), Trondheim, Norway.
Corresponding author: Steven T. Boles (steven.boles@ntnu.no).
