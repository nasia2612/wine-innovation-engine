# Wine Innovation Engine

*Computational Food Science portfolio project*

A three-phase project modeling and optimizing white-wine fermentation: an ODE-based
kinetic "digital twin" (Phase 1), a 1D spatial heat-gradient model layered on top of it
(Phase 2), and a Design-of-Experiments / Response Surface Methodology (RSM) optimization
of the process (Phase 3).

---

## Phase 1 — Fermentation Digital Twin

### The model

The simulator solves a system of 5 coupled ODEs:

```
dX/dt  = (μ - k_d) · X                # viable biomass
dXt/dt = μ · X                        # total biomass
dS/dt  = -(β·X / Y_es) - MNT·X        # sugar consumption
dN/dt  = -μ·X / Y_xn                  # nitrogen depletion
dE/dt  = β·X + MNT·X                  # ethanol production
```

### Key findings

- Ethanol turns out to be almost entirely non-growth-associated.
- Sugar is not the limiting nutrient in grape juice , it's 150–300 g/L, hugely
  abundant. What actually runs out first is nitrogen.
- Coleman (2007) uses an nth-order polynomial `kd_prime(T)` that rises sharply above
  25°C, capturing thermal death at high temperatures.
- We also studied stuck fermentation, caused by low YAN (<80 mg N/L ,this can also
  cause sluggish, very slow fermentation). Stuck fermentation can also be caused by
  high temperature → high `kd`.

### Conclusions

- **Baseline behavior (14°C, S0=220 g/L, N0=300 mg/L):** the model reaches 65.8 g/L
  final ethanol (≈8.3% v/v ABV) with 0 g/L residual sugar in ≈212 h , a clean,
  complete fermentation under adequate-nitrogen conditions.
- **Validation against Coleman Dataset 2** (stuck-fermentation conditions: S0=265 g/L,
  N0=80 mg/L, T=15°C): the model predicts a max viable biomass of 3.59 g/L against the
  paper's reported 3.5 g/L — a ≈2.6% deviation — supporting the model's kinetic
  parameters under nitrogen-limited conditions.
- Ethanol being non-growth-associated and nitrogen (not sugar) being the limiting
  nutrient are the two structural facts that carry through into every later phase:
  they're why Phase 3's DoE varies N0 but not X0 (inoculum size).

---

## Phase 2 — CFD-Lite Spatial Heat Gradient

### The model

Boulton (1979) showed that wine fermentors can be modeled as well-mixed heat
balances — an assumption valid for small white wine tanks. Miller & Block (2020)
demonstrated that large-scale fermentations violate this assumption, exhibiting
measurable temperature gradients. Coleman et al. (2007) established that these
gradients directly alter local Monod kinetics, with `kd` rising sharply above 25°C.
This module implements a 1D discretization of the Boulton heat equation coupled to
the Coleman kinetic model, serving as a computationally lightweight alternative to
full CFD (Miller 2019a).

We have a gradient because of the size of the tank, not because of a cap (this is
white wine, so there's no cap). Real tanks are not well-mixed:

- Fermentation is exothermic — it generates heat.
- That heat is generated wherever yeast is eating sugar, which is mostly in the
  upper and middle liquid.
- The cooling jacket sits at the base and walls, pulling heat out from below.
- The top of the tank is open air — effectively insulated.

As a result, the liquid develops a vertical temperature gradient: the top is warmer
than the bottom.

We model a jacketed cylindrical white-wine fermentor, discretized into N horizontal
layers ("CFD-lite"). Convection happens inside the liquid wine: yeast activity
creates heat, the hot wine in the middle becomes lighter and floats to the top; when
it hits the cold steel wall it cools down, gets heavier, and sinks to the bottom.
This creates a continuous mixing loop called a convection current.

Molecular thermal diffusivity α = 1.33×10⁻⁷ m²/s (derived from Bogard et al. 2020,
Table 1) gives a conduction timescale τ = H²/α ≈ 10,000 h, approximately 40× longer
than fermentation duration. An effective diffusivity α_eff ≈ 1.33×10⁻⁴ m²/s is
therefore used to approximate buoyancy-driven natural convection, following the
Nusselt-number scaling implicit in Bogard et al. 2020.

At 14°C the simulated fermentation is slow (~500–700 h) and retains modest residual
sugar, consistent with Coleman et al. (2007), whose measured 15°C fermentation
(their Fig. 2a) similarly remained above 70 g/L sugar after 600 h. Cold white-wine
fermentations are inherently sluggish. μ_max was set to 0.05 /h following Coleman's
reported value near 11–15°C (their Fig. 3a).

### Conclusions

- **Grid independence:** peak temperature gradient converges as layer count
  increases (N=10, 20, 40, 80 → 6.768, 6.769, 6.799, 6.789 °C respectively) — the
  N=40 vs. N=80 difference is negligible, confirming the result is grid-independent
  rather than a discretization artifact.
- **Interpolation choice is not a sensitive assumption:** using linear vs. cubic
  interpolation for `dS/dt` gives peak gradients of 6.7851°C vs. 6.7822°C , close
  enough that the faster linear interpolation is justified.
- The model supports the qualitative claim that ignoring spatial gradients is a real
  simplification in Phase 1/3: a well-mixed, scalar-temperature assumption
  understates how warm the top layer of a large tank can get relative to the
  jacketed base.

---

## Phase 3 — RSM & Multi-Response Optimization

### Design of Experiments

What we want:
- Final ethanol → want higher.
- Residual sugar → want low (we want it dry, with no stuck fermentation).
- Fermentation time → want low (economic/operational goal: less time in the tank
  means more production cycles and less contamination risk).

These goals conflict: warmer finishes faster and drier but risks overheating; more
initial sugar means more ethanol but a slower, stickier fermentation.

We already know the screening factors (T, S0, N0), so no screening experiment is
needed. The three factors were selected a priori based on the mechanistic structure
of the Boulton ODE model (Arrhenius temperature dependence, dual-substrate Monod
kinetics), not through empirical screening , so the DoE here serves to map the
response surface around conditions we already know have an effect, not to discover
which ones do.

The DoE layer varies temperature, initial sugar (S₀), and initial assimilable
nitrogen (N₀), but not initial yeast concentration (inoculum size, X₀). In this
model the yeast population is not static: biomass grows logistically toward a
carrying capacity set by the limiting nutrient (nitrogen, secondarily sugar), not by
the starting cell count ,a small inoculum reaches the same population ceiling as a
large one, just after a longer lag. Final ethanol is therefore governed by nutrient
availability, not inoculum size. This is consistent with Thuy et al. (2023, Food
Sci. Technol), whose three-factor Box–Behnken optimization of fruit-wine
fermentation found dry yeast concentration (0.15–0.25 g/L) was the least influential
of three factors on ethanol, with no further gain above ~0.2 g/L.

There is no explicit enzyme state, because grape must consists of directly
fermentable hexoses (glucose, fructose) and needs no saccharification. Enzymatic
capacity is absorbed into βmax and the viable-biomass term X_V , biomass is the
model's proxy for enzyme concentration.

A three-factor Box–Behnken design (pyDOE2) is used over T, S0, N0. Factor levels are
coded (-1, 0, +1) with symmetric spacing so the quadratic terms are not distorted by
unequal step sizes. All levels are kept inside the validated domain of Coleman et
al. (2007), but the temperature range is deliberately kept away from the tested
extremes (11°C, 35°C), since all four Coleman validation datasets resulted in stuck
or incomplete fermentations.



### Multi-response optimization (desirability)

Each fitted response is mapped to an individual desirability d(Yᵢ) ∈ [0,1] using
one-sided Derringer–Suich functions (`desirability` in R):

- ethanol → `dMax`
- residual sugar → `dMin`
- fermentation time → `dMin`


The overall desirability is the geometric mean:

```
D = (d_ethanol · d_sugar · d_time)^(1/3)
```
#in the future
Because the Phase 2 CFD-lite model shows a real spatial temperature gradient
(Dirichlet at tank base, Neumann at top), the scalar-T assumption is validated
**after** the optimization, not inside the DoE loop. We follow these steps:

1. Take the optimum setpoint T* from the desirability optimization.
2. Run the Phase 2 spatial model **once** with T* as the boundary condition.
3. Extract the resulting temperature gradient across tank height (ΔT, base vs. top)
   over the fermentation window.

Removing the maintenance term from ethanol production had a larger-than-expected
effect on stuck/complete outcomes (5/15 → 3/15 stuck points), due to an indirect
feedback: lower ethanol accumulation reduces the ethanol-driven death rate (`kd`),
allowing viable biomass to persist longer and consume more sugar before dying —
illustrating the coupled nature of the ODE system beyond the directly edited term.

### Conclusions

- **Optimizer result** (`results/tables/optimal_conditions.csv`): the
  desirability-weighted optimum is **T = 15°C, S0 = 253.64 g/L, N0 = 321.18 mg/L**,
  with predicted time 230.75 h, predicted peak_rate 1.7557 g/L/h, predicted
  conversion 1.0089, and overall desirability D = 0.8386.
- The coded temperature factor sits exactly at its lower bound (x1 = −1), matching
  the earlier observation that the model "wanted to go further with the
  temperature but can't" , the optimizer is boundary-clamped, not sitting at an
  interior optimum, which is consistent with neither `peak_rate` nor `time` having
  an in-range stationary point (both are saddle/ridge surfaces, per the
  eigenanalyses above).
- Two of the 15 Box–Behnken design points (13%) failed to reach dryness within the
  simulated window, and both were at the low-nitrogen level (N0 = 140 mg/L) — this
  is the same nitrogen-limited stuck-fermentation mechanism identified in Phase 1,
  now showing up quantitatively in the DoE responses.
- At a low fermentation temperature, yeasts are less sensitive to the toxic effects
  of ethanol concentration. Cell growth and fermentation rate slow down, but cell
  viability improves (Shener et al., 2007); wines fermented at 16°C were reported as
  fresher and lighter in taste. Wines with a higher initial sugar content were more
  voluminous with a fruity character, while those with lower initial sugar were
  light, airy, with a strong freshness and fruity character. On the other hand,
  alcoholic fermentation at a higher temperature (>15°C) possibly reduces thiol
  aromas due to increased ester content; lower fermentation temperatures (10–15°C)
  can improve the aromatic profile of wines (Georgiev et al., 2024). The optimizer's
  low-temperature recommendation is therefore consistent with both the process
  objectives (dryness, avoiding stuck fermentation) and the sensory literature.

---

## Limitations

- **Phase 2 heat transport:** the 1D model captures conductive heat transport only.
- **Deterministic pure error:** the Box–Behnken center point is replicated 3 times,
  but the simulator is fully deterministic (no stochastic noise), so those replicate
  runs are numerically identical. This makes the "pure error" and lack-of-fit tests
  in the RSM output above near-degenerate (F-values on the order of 10²⁹–10³⁰) ,
  they should not be read as evidence of a genuine model-fit problem, since there is
  essentially no real replicate variance to compare against.
- **Ethanol not in the reported optimum:** an `ethanol` RSM model is fit in
  `rsm_fit.Rmd`, but `desirability_optimize.R` only builds desirability functions for
  `time`, `peak_rate`, and `conversion` (as a constraint). Ethanol yield is not
  actually part of the multi-response optimum reported above.
- **Extrapolation past the physical ceiling:** the optimizer's predicted conversion
  (1.0089) exceeds the physical maximum of 1.0 (100%) ,a symptom of fitting a
  second-order polynomial near the boundary of the design space.
-


