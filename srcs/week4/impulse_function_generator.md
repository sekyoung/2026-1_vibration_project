# `impulse_function_generator.m` — Discrete Tuned (1-cosine) Gust

Week 4. This replaces the Week-3 mathematical Dirac impulse (`impulse_response.m`, which used MATLAB's `impulse()`) with the **certification "impulsive" gust**: the discrete tuned 1-cosine gust defined by **FAA 14 CFR §25.341(a)** and **EASA CS-25.341(a)**. The model input is gust velocity `w_g` [m/s], so the generated profile feeds the Week-3 state-space model (`wing_2dof.mat`) directly.

---

## 1. Why a discrete tuned gust, not a Dirac impulse

A Dirac delta is a mathematical idealization — infinite amplitude, zero duration, infinite bandwidth. It is convenient for extracting the impulse-response function but it is not a load case any airworthiness authority recognizes. The regulatory analog of an "impulsive" atmospheric disturbance is the **discrete gust**, and both the FAA and EASA prescribe a specific **1-cosine** velocity profile for it:

$$U(s) = \frac{U_{ds}}{2}\left[\,1 - \cos\!\left(\frac{\pi s}{H}\right)\right], \qquad 0 \le s \le 2H$$

where `s` is the distance the aircraft has penetrated into the gust, `H` is the **gust gradient** (the distance over which the gust rises to its peak), and `U_ds` is the **design gust velocity**. Converting distance to time along the flight path with `s = U·t` gives the time signal the script produces:

$$w_g(t) = \frac{U_{ds}}{2}\left[\,1 - \cos\!\left(\frac{\pi U t}{H}\right)\right], \qquad 0 \le t \le \frac{2H}{U}, \quad \text{else } 0.$$

The design gust velocity scales with the sixth root of the gradient:

$$U_{ds} = U_{ref}\, F_g \left(\frac{H}{107}\right)^{1/6} \quad [\text{H in metres}] \qquad\left(= U_{ref}\, F_g \left(\tfrac{H}{350}\right)^{1/6} \text{ for H in ft}\right)$$

with `U_ref` the reference gust velocity (17.07 m/s ≈ 56 ft/s EAS at sea level, between V_B and V_C) and `F_g` the flight-profile alleviation factor.

**The criterion is not a single gust.** §25.341(a)(3) requires that *a sufficient number of gust gradients be investigated to find the critical response for each load quantity*. That "tuning" sweep is the heart of the regulation, so the generator produces a **family** of gusts over a range of `H`, runs each through the wing model, and reports the gradient that maximizes the peak plunge `h` and pitch `θ`.

---

## 2. Why the *scaled* gradient band (default), with literal CS-25 as a toggle

The Week-3 model is a **bench/section-scale** wing, not a full-scale transport:

| quantity | model value | CS-25 design world |
|---|---|---|
| flight speed `U` | 12 m/s | ~230–250 m/s cruise |
| chord | 0.30 m | several metres |
| structural modes | ≈ 5.5 Hz (plunge), ≈ 9.8 Hz (pitch) | ~1–10 Hz but at full scale |
| CS-25 gradient range | — | 9.1–107 m (30–350 ft) |

Applying the CS-25 gradient range literally at `U = 12 m/s`: even the **shortest** 30-ft gust (H = 9.1 m) takes `2H/U ≈ 1.5 s` to traverse, against a 5.5 Hz (0.18 s period) structural mode. Every CS-25 gradient therefore lands in the **quasi-static** regime — the wing follows the gust without any dynamic amplification, and the tuning sweep finds nothing. The standard's numbers are calibrated to full-scale geometry and speed.

The aerodynamically meaningful parameter is the **gust frequency relative to the structural modes**, not the absolute gradient. A 1-cosine gust of total length `2H` has a dominant frequency content around

$$f_{gust} \approx \frac{U}{2H}.$$

So the regulation's "tune `H` to find the critical response" maps cleanly onto "sweep the gust so its frequency passes through the wing's resonances" — the same resonance argument used for the engine excitation in Week 2, now applied to the gust.

**Default `standard = 'scaled'`** keeps the §25.341 *form and tuning methodology verbatim* but chooses the gradient band so `f_gust` brackets the structural modes:

$$f_{gust} \in [\,0.5\,f_{n,\min},\; 2\,f_{n,\max}\,] \;\Rightarrow\; H = \frac{U}{2 f_{gust}} \in \big[\,\tfrac{U}{4 f_{n,\max}},\; \tfrac{U}{f_{n,\min}}\,\big].$$

For this model that is roughly `H ≈ 0.3–2.2 m` — physically meaningful and dynamically active.

**`standard = 'cs25'`** restores the literal 9.1–107 m range and `U_ref = 17.07 m/s` for a by-the-book demonstration of the regulation (expect a quasi-static result at this scale). Both paths use the identical 1-cosine shape, `U_ds` formula, and critical-gradient search — only the swept band and `U_ref` differ.

For the scaled mode, `U_ref` defaults to `0.10·U` as a documented, physically modest gust intensity. Because the wing model is **linear**, the *location* of the critical gradient is independent of `U_ref` (amplitude scales linearly), so this choice affects only the absolute magnitude of the reported response, not which `H` is critical.

---

## 3. Parameters (name–value overrides)

All defaults live in `default_config()`; override any as `('name', value)` pairs.

| field | default | meaning |
|---|---|---|
| `standard` | `'scaled'` | `'scaled'` (mode-bracketing band) or `'cs25'` (literal 9.1–107 m, U_ref = 17.07) |
| `model_file` | `''` | path to `wing_2dof.mat`; empty → auto-locate beside script, then `../week3/` |
| `fs` | `2000` | output sample rate [Hz] |
| `nH` | `25` | number of gust gradients in the sweep (log-spaced) |
| `settle` | `6` | quiet time [s] appended after the longest gust so transients decay |
| `Fg` | `1.0` | flight-profile alleviation factor (1.0 = conservative, sea-level limit) |
| `Uref` | `[]` | reference gust velocity [m/s EAS]; empty → set by `standard` |
| `band` | `[0.5 2.0]` | scaled mode only: `f_gust` limits as multiples of `[min fn, max fn]` |
| `plot` | `true` | draw the 3-panel summary figure |

## 4. Output `gust` struct

| field | description |
|---|---|
| `t` | common time vector [s] |
| `wg` | `[numel(t) × nH]` gust-velocity family `w_g(t)` [m/s], one column per gradient |
| `H` | swept gust gradients [m] |
| `Uds` | design gust velocity per gradient [m/s] |
| `f_gust` | equivalent gust frequency `U/(2H)` per gradient [Hz] |
| `resp.h_peak`, `resp.th_peak` | peak \|h\| [mm] and \|θ\| [deg] per gradient |
| `crit` | critical gradient and peak for each output (`H_h`, `H_th`, `idx_h`, `idx_th`, `h_max`, `th_max`) |
| `cfg` | the resolved configuration used |

The figure shows: (1) the 1-cosine gust family with the critical gust highlighted, (2) the **tuning curve** of peak response vs `H` with critical gradients marked, and (3) the wing's `h`/`θ` time response to the critical gust.

## 5. Usage

```matlab
gust = impulse_function_generator;                       % scaled, mode-bracketing (default)
gust = impulse_function_generator('standard','cs25');    % literal CS-25 range
gust = impulse_function_generator('band',[0.25 4], 'nH',60);   % wider, finer sweep
```

Requires the Control System Toolbox (`ss`, `lsim`) — same dependency as the Week-3 scripts. Place `wing_2dof.mat` beside this file or leave it in `week3/` (auto-located).

---

## References

1. **FAA — 14 CFR §25.341(a)**, *Gust and turbulence loads — Discrete Gust Design Criteria.* Defines the 1-cosine shape `U = (U_ds/2)[1 − cos(πs/H)]` for `0 ≤ s ≤ 2H`; the `U_ds = U_ref·F_g·(H/350)^{1/6}` formula; the 30–350 ft gradient range; `U_ref = 56 ft/s` EAS at sea level (V_B–V_C); and the requirement to investigate a sufficient number of gradients to find the critical response.
2. **EASA — CS-25.341(a)** and **AMC 25.341**, *Easy Access Rules for Large Aeroplanes (CS-25).* European equivalent; states the gradient range as 9.1–107 m, `U_ref = 17.07 m/s` EAS at sea level, and the note that gust gradients beyond 350 ft be considered if 12.5× the mean geometric chord exceeds 350 ft.
3. **FAA — AC 25.341-1**, *Dynamic Gust Loads.* Guidance material: typical 1-cosine design gust velocity profiles, definitions of `U_ref`, `F_g`, `H`, and the discrete gust response procedure.
4. **FAA — 14 CFR §23.341 / EASA CS-23** (normal-category). The small-aircraft counterpart, with lower reference velocities and no `F_g`; the more natural regulatory home for a model of this scale, noted here as the alternative to CS-25.
5. **DOT/FAA/AR-99/62**, *Studies of Time-Phased Vertical and Lateral Gusts: Development of Multiaxis One-Minus-Cosine Gust Model* (Oct 1999). Origin of the multiaxis 1-cosine formulation referenced by AC 25.341-1.
