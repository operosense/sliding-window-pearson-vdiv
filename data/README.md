# Data

Operando data from **Case 1** of the paper: a LiNi0.8Mn0.1Co0.1O2 (NMC811)
all-solid-state cell with an embedded fibre-Bragg-grating (FBG) strain sensor, second
cycle. The strain and electrochemical channels were synchronised by their shared
acquisition timestamps and interpolated onto a common voltage grid. Each row is one
synchronised (V, ε, Q) sample.

| file | half-cycle | columns (tab-separated, 1 header row) | rows |
|------|------------|----------------------------------------|------|
| `Ch2 Good correlation.txt`  | 2nd charge    | `Voltage` (V), `Strain` (FBG mechanical channel), `Capacity` | 14 375 |
| `DCh2 Good correlation.txt` | 2nd discharge | `Voltage` (V), `Strain` (FBG mechanical channel), `Capacity` | 13 713 |

Rows are in acquisition order, so the first row is the start of the half-cycle
(ascending voltage for charge, descending for discharge).

`correlation_results Ch2 Good.txt` is the r_p(V) trace for the 2nd charge computed with
non-overlapping 0.02 V windows, as plotted in Figure 2k. The first two columns are
`Voltage` and `Correlation`. The other columns hold the voltage/strain record and its
first derivative as exported from the plotting software.
`--` marks windows with too few points.

License: CC BY 4.0 (see `../LICENSE-DATA`).
