# Threshold Estimation — Method and Mathematical Background

How the error-correction **threshold** $p_{\mathrm{th}}$ is extracted from the Monte-Carlo
logical-error-rate data, the **finite-size-scaling (FSS)** hypothesis behind it, the statistical
model of the data, the fitting and uncertainty procedure, and the assumptions and limitations.

Files covered:
- [threshold.py](./threshold.py) — `SamplePoint`, `FitResult`, `fss`, `estimate_threshold`
- [simulation.py](./simulation.py) — `estimate_logical_error_rate`, `wilson_interval`
- [main.py](./main.py) — sweep, threshold call, `result_<eta>.png` / `collapse_<eta>.png`

This document is the companion to the project [README.md](../README.md) "Threshold Estimation"
section and to [circuit_builder/README.md](./circuit_builder/README.md) (how the circuits are built).
Symbols ($p$, $p_{\mathrm{th}}$, $d$, $\nu$, $p_L$) match those used there and in the code.

---

## 1. The threshold as a phase transition

Fix the noise model and the bias $\eta$, and sweep the physical error rate $p$ and the code
distance $d$. The decoder either succeeds or fails on each shot, giving a **logical error rate**
$p_L(p, d)$. The defining property of a good code family is a sharp **threshold** $p_{\mathrm{th}}$:

$$
\lim_{d \to \infty} p_L(p, d) =
\begin{cases}
0 & p < p_{\mathrm{th}} \quad (\text{increasing } d \text{ helps}) \\
\text{const} > 0 & p > p_{\mathrm{th}} \quad (\text{increasing } d \text{ hurts})
\end{cases}
$$

So below threshold, making the code bigger drives the logical error rate down exponentially; above
threshold, it drives it up. The crossover is a **continuous (second-order) phase transition**:
$p_{\mathrm{th}}$ plays the role of the critical point, $p - p_{\mathrm{th}}$ the role of the
reduced temperature, and $d$ the role of the system size. Right at $p = p_{\mathrm{th}}$, $p_L$ is
(asymptotically) independent of $d$ — which is why the raw $p_L$-vs-$p$ curves for different $d$
all **cross at one point**, the first visual estimate of the threshold (the `result_<eta>.png`
figure). The XZZX-vs-CSS comparison is exactly a comparison of where these crossings sit, and how
they move with the bias $\eta$.

---

## 2. The finite-size-scaling hypothesis

Near a continuous transition the only relevant length is the **correlation length**

$$
\xi(p) \sim |p - p_{\mathrm{th}}|^{-\nu}
$$

which diverges at the critical point with a **critical exponent** $\nu$. Universality says that a
finite system of size $d$ does not depend on $p$ and $d$ separately but only on the dimensionless
ratio $d / \xi \sim (p - p_{\mathrm{th}})\, d^{1/\nu}$. Applied to the logical error rate this is
the **finite-size-scaling ansatz**

$$
p_L(p, d) = F\big( \underbrace{(p - p_{\mathrm{th}}) d^{1/\nu}}_{x} \big)
$$

with a single **universal scaling function** $F$ shared by all distances. Two consequences drive
the whole method:

- **Collapse.** If we plot $p_L$ against the rescaled variable $x = (p - p_{\mathrm{th}}) d^{1/\nu}$
  using the correct $(p_{\mathrm{th}}, \nu)$, the points from *every* distance fall on the *same*
  curve $F(x)$. A good collapse is the qualitative evidence that the threshold is real (this is the
  `collapse_<eta>.png` figure).
- **Crossing.** At $x = 0$, i.e. $p = p_{\mathrm{th}}$, $p_L = F(0)$ is independent of $d$ — the
  crossing point of §1.

$\nu$ is a by-product of the fit, not the target, but it is worth reporting: for the 2D surface code
under independent noise it sits around $\nu \approx 1.5$ (random-bond Ising / percolation
universality), so a fitted $\nu$ wildly outside $\sim 1$ – $1.5$ is a warning that the data or the fit
window are poor.

---

## 3. Polynomial approximation of the scaling function (and why the fit must be windowed)

The universal function $F$ has no closed form, so following the standard practice (Wang–Harrington–
Preskill, Harrington's thesis) we approximate it **near the threshold** by a low-order Taylor
expansion in $x$ — here a quadratic:

$$
F(x) \approx a + b\,x + c\,x^2 .
$$

This is implemented as `fss(X, p_th, nu, a, b, c)` in [threshold.py](./threshold.py), with five free
parameters $(p_{\mathrm{th}}, \nu, a, b, c)$.

The Taylor truncation is **only valid for small $|x|$**, i.e. close to the crossing. Far below
threshold $p_L$ decays roughly exponentially in $d$ and far above it saturates — neither is a
quadratic in $x$. Including those far-from-threshold points would force the parabola to bend to fit
data it cannot describe, pulling the apparent $p_{\mathrm{th}}$ and inflating $\chi^2$. **This is the
reason `estimate_threshold` windows the fit** (§5): the quadratic is a near-threshold model and must
be fed near-threshold data.

> Caveat noted in code (`# TODO: take the error bias into account`): the same scaling form is used
> for every bias $\eta$. A more refined analysis could let the scaling function or $\nu$ depend on
> the bias; the present quadratic collapse is bias-agnostic and is used only to *locate*
> $p_{\mathrm{th}}$, which is adequate for the CSS-vs-XZZX comparison.

---

## 4. The statistical model of each data point

Each `SamplePoint` is a Monte-Carlo estimate. For a configuration $(p, d)$ we take $N$ shots
([simulation.py](./simulation.py), `estimate_logical_error_rate`) and count $k$ logical failures.
The failures are independent Bernoulli trials, so

$$
k \sim \mathrm{Binomial}(N, p_L), \qquad \hat{p}_L = \frac{k}{N}.
$$

Sampling is **adaptive**: shots are drawn in batches until either $k$ reaches `target_errors` or $N$
reaches `max_shots`. Stopping on a fixed error count keeps the *relative* statistical uncertainty
roughly constant (≈ $1/\sqrt{k}$) across points spanning orders of magnitude in $p_L$.

### Wilson score interval

The naive Bernoulli error $\sqrt{\hat{p}_L(1-\hat{p}_L)/N}$ collapses to zero when $k = 0$ and is
badly miscalibrated for small $p_L$ — exactly the deep-sub-threshold regime that dominates this
sweep. We instead use the **Wilson score interval** (`wilson_interval`), which inverts the binomial
score test. For $z$ standard deviations and $\hat{p} = k/N$:

$$
\text{center} = \frac{\hat{p} + \dfrac{z^2}{2N}}{1 + \dfrac{z^2}{N}}, \qquad
\text{half-width} = \frac{z}{1 + \dfrac{z^2}{N}} \sqrt{ \frac{\hat{p}(1-\hat{p})}{N} + \frac{z^2}{4N^2} } .
$$

It stays inside $[0,1]$, gives a sensible non-degenerate bound even at $k = 0$, and is well-behaved
as $p_L \to 0$. We take the **half-width at $z = 1$** as the 1σ uncertainty $\sigma$ stored on each
`SamplePoint`; the *same* 1σ interval is drawn as the plot error bars in [main.py](./main.py), so the
figure and the fit weights use one consistent convention. Points with $k = 0$ (no failures observed)
carry no usable $\sigma$ for a two-sided fit and are excluded from the FSS fit (§5); in the
`result_<eta>.png` plot they are instead shown as Wilson upper-bound arrows.

The relation between $z$ and confidence: $z = 1 \Rightarrow$ ≈ 68 % (1σ), $z = 1.96 \Rightarrow$
95 %. We standardize on 1σ everywhere so that "error bar" means the same thing in the data, the fit,
and the threshold.

---

## 5. Fitting procedure and uncertainty on the threshold

`estimate_threshold(sample_points, window_frac, n_boot, seed)` in [threshold.py](./threshold.py)
performs the following.

**(0) Drop zero-error points.** Points with $k = 0$ are removed (§4): the collapse ansatz is invalid
that far below threshold and their uncertainty is one-sided.

**(1) Locate the threshold (pass 1).** A coarse weighted least-squares fit of `fss` over all
remaining points gives a first estimate $p_{\mathrm{th}}^{(0)}$. The objective is the usual
chi-square,

$$
\chi^2 = \sum_i \left( \frac{\hat{p}_{L,i} - F(x_i)}{\sigma_i} \right)^2,
$$

minimized by `scipy.optimize.curve_fit`, with each residual weighted by the point's 1σ Wilson
$\sigma_i$.

**(2) Window and refit (pass 2).** Keep only points within
$|p - p_{\mathrm{th}}^{(0)}| \le w$, with half-width $w = \text{window\_frac}\cdot(p_{\max}-p_{\min})$,
and refit. Restricting to the near-threshold window is what makes the quadratic model legitimate
(§3). The refit passes **`absolute_sigma=True`**, which tells `curve_fit` to treat the $\sigma_i$ as
*real* uncertainties so that the returned covariance `pcov` is a genuine statistical covariance —
**not** rescaled by the reduced $\chi^2$, as it would be by default. The covariance-based error bar
is then

$$
\delta p_{\mathrm{th}}^{\text{cov}} = \sqrt{(\mathrm{pcov})_{00}}, \qquad
\delta \nu = \sqrt{(\mathrm{pcov})_{11}} .
$$

If the window holds too few points for the five parameters, a warning is emitted and the fit falls
back to all measured points; the window count is reported as `n_points` so this condition is visible.

**(3) Goodness of fit.** The reduced chi-square over the windowed points,

$$
\chi^2_{\mathrm{red}} = \frac{1}{N_{\text{win}} - 5} \sum_{i \in \text{win}}
\left( \frac{\hat{p}_{L,i} - F(x_i)}{\sigma_i} \right)^2 ,
$$

is reported. $\chi^2_{\mathrm{red}} \approx 1$ means the quadratic FSS model is consistent with the
Monte-Carlo scatter; $\gg 1$ signals model misfit (window too wide, distances too small, or genuine
scaling violations) and $\ll 1$ over-fitting or over-estimated error bars.

**(4) Bootstrap confidence interval.** The covariance error assumes the linearized model is exact at
the optimum, which is optimistic for a non-linear fit on noisy data. We therefore also run a
**parametric bootstrap**: for `n_boot` replicas we resample each windowed point's failure count from
its Binomial law, $k_i^\ast \sim \mathrm{Binomial}(N_i, \hat{p}_{L,i})$, refit, and collect the
distribution of $p_{\mathrm{th}}^\ast$. Its standard deviation gives
$\delta p_{\mathrm{th}}^{\text{boot}}$ and its 16/84 percentiles give a confidence interval. The
bootstrap RNG is seeded (`np.random.default_rng(seed)`) so the interval is **reproducible**.

The reported `p_th_err` is the **more conservative** (larger) of the covariance and bootstrap
estimates. All of this is returned in a `FitResult`:

```text
FitResult(p_th, p_th_err, p_th_err_cov, p_th_err_boot, p_th_ci, nu, nu_err, chi2_red, n_points, window, popt, pcov)
```

[main.py](./main.py) uses it to label the threshold line as $p_{\mathrm{th}} \pm \delta$ with a 1σ
band (`axvspan`), and to title each panel of the collapse figure with $p_{\mathrm{th}}$, $\nu$, and
$\chi^2_{\mathrm{red}}$.

---

## 6. Reading the data-collapse figure

`collapse_<eta>.png` (built by `draw_collapse` in [main.py](./main.py)) is the primary validation
plot. Each measured point is rescaled to $x = (p - p_{\mathrm{th}})d^{1/\nu}$ and plotted against
$p_L$, with the fitted parabola $a + bx + cx^2$ overlaid:

- **Good threshold:** points from all distances $d$ land on the single overlaid curve (a clean
  collapse), and the curves cross near $x = 0$.
- **Poor fit:** distances peel apart into separate branches, or the parabola misses the cloud — read
  together with a large $\chi^2_{\mathrm{red}}$, this means the window is too wide, the distances are
  too small, or the bias breaks the simple scaling form.

The collapse plot and $\chi^2_{\mathrm{red}}$ together are what let a reader judge the quoted
$p_{\mathrm{th}} \pm \delta$ rather than taking it on faith.

---

## 7. Assumptions and limitations

- **X-memory threshold only.** The experiment holds and reads out a single logical-X observable
  (see [circuit_builder/README.md](./circuit_builder/README.md) §0). The quoted $p_{\mathrm{th}}$ is
  the X-memory threshold under Z-biased noise — the regime where XZZX is expected to win — not a
  full (X *and* Z) code threshold.
- **Decoder is MWPM.** Thresholds are decoder-dependent; MWPM under correlated / high-bias noise is
  sub-optimal, so the values are a lower bound for what a correlated-matching or BP+OSD decoder would
  give. See the project [README.md](../README.md) "Decoder choice and its limitations".
- **Small distances.** The sweep uses $d \in \{3, 5, 7\}$. Small distances carry large finite-size
  corrections beyond the leading FSS form, so $p_{\mathrm{th}}$ has a real systematic on top of the
  reported statistical error. Adding $d = 9, 11$ (the diagram code already builds them) and sampling
  more densely near the crossing are the most effective ways to shrink both the statistical and the
  systematic error — watch `n_points` and $\chi^2_{\mathrm{red}}$ to know when the window is starved.
- **Quadratic, bias-agnostic collapse.** Only the leading quadratic term of $F$ is kept and the same
  form is used for every $\eta$ (the `# TODO` in `fss`). This is sufficient to locate the threshold
  for the CSS-vs-XZZX comparison but is not a high-precision determination of $p_{\mathrm{th}}$ or of
  the critical exponent $\nu$.
- **Per-model $p$-axis.** Thresholds from different noise models live on different $p$-axes and must
  not be compared at the same numeric $p$ (see [circuit_builder/README.md](./circuit_builder/README.md)
  §5, "On the meaning of `p` across models"). The comparison that *is* valid is CSS vs XZZX within one
  model and one $\eta$.