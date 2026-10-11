# Algorithm 1 — Operational V_div detection

Pseudocode for the Figure 1b flowchart, implemented in
[`code/sliding_window_pearson_vdiv.py`](../code/sliding_window_pearson_vdiv.py).

Parameters (as used in the paper and its Supporting Information):

- window width **ΔV = 0.02 V**
- minimum **≥ 10 points per window**
- persistence **N = 5 consecutive windows**
- threshold **k = 5 standard deviations** below the baseline mean
- depth cutoff **DEEP_RP = −0.5**: a trough counts as a V_div only if r_p falls at least this far (toward −1)
- significance: **95 % Fisher-z lower bound** on non-overlapping ΔV bins (z_crit = 1.96)
- reported resolution **±0.01 V**

**V_div is the voltage of each r_p minimum** (the anti-correlation point, where r_p → −1).
A trough is a contiguous excursion below τ that lasts at least N windows and reaches
r_p ≤ DEEP_RP. V_div is the voltage of the lowest r_p in that trough. A half-cycle can
contain several inversions, so several V_div values can be returned.

---

## Plain-text version

```
Algorithm 1  Operational V_div detection

Input : synchronized Q(V) and eps(V) (or F(V)) on a common voltage grid
        (shared acquisition timestamps)
        window width           dV       = 0.02 V
        persistence count       N        = 5
        threshold multiple      k        = 5
        depth cutoff            DEEP_RP  = -0.5
        baseline plateau test:  n_b      = 6 windows with r_p >= 0.95 (auto-detected)
Output: list of V_div (one per r_p-minimum inversion), each with a significance flag

 1  interpolate Q and eps onto a common voltage grid
 2  for i = 1 .. M do                          # M windows, one-sample step, centered
 3      W_i <- samples with voltage in [V_i - dV/2, V_i + dV/2]
 4      if |W_i| >= 10 then
 5          r_p[i] <- Pearson( Q[W_i], eps[W_i] )
 6  # ---- autonomous baseline: no chemistry-specific voltage window ----
 7  i0 <- smallest i such that r_p[i .. i+n_b-1] are ALL >= 0.95
 8  (mu_base, sigma_base) <- mean, std of r_p[i0 .. i0+n_b-1]
 9  tau <- mu_base - k * sigma_base
10  # ---- significance stage: non-overlapping dV bins ----
11  r_bin[b], n_bin[b] <- Pearson r and point count of each non-overlapping dV bin
12  z_base <- arctanh( mean r_bin over the first n_b consecutive bins with r_bin >= 0.95 )
13  # ---- find troughs and take V_div at each r_p MINIMUM ----
14  Vdiv_list <- [ ]
15  i <- 1
16  while i <= M do
17      if r_p[i] < tau then
18          j <- end of the contiguous run where r_p < tau          # the trough
19          if (j - i) >= N then                                    # lasts >= N windows
20              m <- argmin r_p over [i, j)                         # deepest point
21              if r_p[m] <= DEEP_RP then
22                  V_div <- V_m                                    # voltage of the r_p minimum
23                  b     <- non-overlapping bin containing V_m
24                  SE    <- 1 / sqrt(n_bin[b] - 3)
25                  LB    <- tanh( z_base - z_crit * SE )           # 95 % lower bound
26                  flag  <- "confirmed" if r_bin[b] < LB else "qualitative"
27                  append (V_div +/- 0.01 V, r_p[m], flag) to Vdiv_list
28          i <- j
29      else
30          i <- i + 1
31  return Vdiv_list
```

---

## LaTeX version (algorithm2e)

```latex
\begin{algorithm}[t]
\caption{Operational $V_{div}$ detection}
\KwIn{synchronized $Q(V)$ and $\varepsilon(V)$ (or $F(V)$) on a common voltage grid;
      $\Delta V = 0.02$~V; $N=5$; $k=5$; depth cutoff $r_\mathrm{deep}=-0.5$;
      baseline = first 6 windows with $r_p\ge0.95$}
\KwOut{list of $V_{div}$ (one per $r_p$-minimum inversion), each with a significance flag}
interpolate $Q,\varepsilon$ onto a common voltage grid\;
\ForEach{voltage $v$ (one-sample step, centered window)}{
  $W \leftarrow$ samples with voltage in $[v-\tfrac{\Delta V}{2},\,v+\tfrac{\Delta V}{2}]$\;
  \If{$|W|\ge 10$}{ $r_p(v) \leftarrow \mathrm{Pearson}(Q_W,\varepsilon_W)$\; }
}
$(\mu_\mathrm{base},\sigma_\mathrm{base}) \leftarrow$ mean, std of $r_p$ over the auto baseline\;
$\tau \leftarrow \mu_\mathrm{base} - k\,\sigma_\mathrm{base}$\;
$z_\mathrm{base} \leftarrow \operatorname{arctanh}(\bar r_\mathrm{base})$ over the baseline non-overlapping bins\;
$\mathcal{V}\leftarrow[\,]$\;
\ForEach{maximal run $[a,b)$ with $r_p<\tau$}{
  \If{$(b-a)\ge N$}{
    $m \leftarrow \arg\min_{a\le i<b} r_p(i)$\;
    \If{$r_p(m)\le r_\mathrm{deep}$}{
       $V_{div}\leftarrow V_{m}$\tcp*{voltage of the $r_p$ minimum}
       $LB \leftarrow \tanh\!\big(z_\mathrm{base} - z_\mathrm{crit}/\sqrt{n_{b^\star}-3}\big)$\tcp*{bin $b^\star \ni V_m$}
       flag $\leftarrow$ ``confirmed'' \textbf{if} $r_{b^\star}<LB$ \textbf{else} ``qualitative''\;
       append $(V_{div}\pm0.01\text{ V},\,r_p(m),\,\text{flag})$ to $\mathcal{V}$\;
    }
  }
}
\Return $\mathcal{V}$\;
\end{algorithm}
```
