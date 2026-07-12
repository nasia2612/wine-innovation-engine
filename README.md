
Computational Food Science portfolio project
# Wine Innovation Engine



-----Phase 1 — Fermentation digital twin-----

The model

The simulator solves a system of 5 coupled ODEs:

dX/dt  = (μ - k_d) · X                         # viable biomass
dXt/dt = μ · X                                  # total biomass  
dS/dt  = -(β·X / Y_es) - MNT·X                 # sugar consumption
dN/dt  = -μ·X / Y_xn                           # nitrogen depletion
dE/dt  = β·X + MNT·X                           # ethanol production


# ethanol turns out to be almost entirely non-growth-associated 
# sugar is not the limiting nutrient in grape juice. Sugar is 150–300 g/L — hugely abundant. What actually runs out first is nitrogen.

# Coleman 2007 uses an nth-order polynomial kd_prime(T) that rises sharply above 25°C, capturing thermal death at high temperatures.

We studied  also the stuck fermentation that is caused by low YAN <80 mg N/L (it might also cause sluggish fermentation,very slow fermantation).Stuck fermentation might also be caused by high temperature ->high kd

For phase 2 :
Boulton (1979) showed that wine fermentors can be modeled as well-mixed heat balances  an assumption valid for small white wine tanks. 
Miller & Block (2020) demonstrated that large-scale fermentations violate this assumption, exhibiting measurable temperature gradients. 
Coleman et al. (2007) established that these gradients directly alter local Monod kinetics, with kd rising sharply above 25°C. This module implements a 1D discretization of the Boulton heat equation coupled to the Coleman kinetic model, serving as a computationally lightweight alternative to full CFD (Miller 2019a).




#we have gradient because of the size of the tank not because of the cap (because we have white wine we dont have cap)

real tanks are not well-mixed:
Fermentation is exothermic — it generates heat
That heat is generated wherever yeast is eating sugar, which is mostly in the upper and middle liquid
The cooling jacket sits at the base and walls — pulling heat out from below
The top of the tank is open air — effectively insulated



#: the liquid develops a vertical temperature gradient. The top is warmer than the bottom.


lets say we use a Jacketed cylindrical white wine fermentor (we discretized it into N horizontal layers for CFD-lite)

Convection happens inside the liquid wine. During fermentation, yeast creates heat. The hot wine in the middle becomes lighter and floats to the top. When it hits the cold steel wall, it cools down, gets heavier, and sinks to the bottom. This creates a continuous 'mixing loop' called a convection current.

#Limitation: The 1D model captures conductive heat transport only. 

Molecular thermal diffusivity α = 1.33×10⁻⁷ m²/s (derived from  Bogard et al. 2020, Table 1) gives a conduction timescale 
τ = H²/α ≈ 10,000 h , approximately 40× longer than fermentation  duration. An effective diffusivity α_eff ≈ 1.33×10⁻⁴ m²/s is 
therefore used to approximate buoyancy-driven natural convection,  following the Nusselt-number scaling implicit in Bogard et al. 2020. 


At 14°C the simulated fermentation is slow (~500–700 h) and retains
modest residual sugar, consistent with Coleman et al. (2007), whose
measured 15°C fermentation (their Fig. 2a) similarly remained above
70 g/L sugar after 600 h. Cold white-wine fermentations are inherently
sluggish. μ_max was set to 0.05 /h following Coleman's reported value
near 11–15°C (their Fig. 3a).



#PHASE 3 RSM -OPTIMIZATION

What we want :
-final ethanol-> want higher 
-residual sugar->want low (we want it dry and to complete the fermentation)
-fermentation time -> want low (Economic/operational goal: less time in the tank = more production cycle, less risk of contamination)

We have conflicts Warmer finishes faster and drier but can overheat; more sugar means more ethanol but slower, stickier fermentation

we  know the screening factors(T,S0,N0) so we do not need a screening experiment.
The three factors were selected a priori based on the mechanistic structure of the Boulton ODE model (Arrhenius temperature dependence, dual-substrate Monod kinetics), not through empirical screening SO the DoE here serves to map the response surface around conditions that we already know have an effect, not to find out which ones do.

The DoE layer varies temperature, initial sugar (S₀), and initial assimilable nitrogen (N₀), but not initial yeast concentration (inoculum size, X₀).
. In this model the yeast population is not static: biomass grows logistically toward a carrying capacity set by the limiting nutrient (nitrogen, secondarily sugar), not by the starting cell count ,a small inoculum reaches the same population ceiling as a large one, just after a longer lag. Final ethanol is therefore governed by nutrient availability, not inoculum size.
 in a three-factor Box–Behnken optimisation of fruit-wine fermentation (Thuy et al., 2023, Food Sci. Technol), dry yeast concentration (0.15–0.25 g/L) was the least influential of the three factors on ethanol, with no further gain above ~0.2 g/L.
 there is no explicit enzyme state because grape must consists of directly fermentable hexoses (glucose, fructose) and needs no saccharification.
 So, enzymatic capacity is absorbed into βmax and the viable-biomass term X_V —> biomass is the model's proxy for enzyme concentration.

A three-factor Box–Behnken design (pyDOE2) is used over T, S0, N0. Factor levels
are coded (-1, 0, +1) with symmetric spacing so the quadratic terms are not
distorted by unequal step sizes.

Also  All levels are kept inside the validated domain of the Coleman et al. (2007) but The temperature range is deliberately kept away from thetested extremes (11 °C, 35 °C), since all four Coleman validation datasets resulted in stuck or incomplete fermentations.


Each model is evaluated by ANOVA; non-significant terms are dropped (model
reduction) so the retained surface reflects only statistically supported effects.
Model adequacy is judged on both R² and residual behaviour.
### Multi-response optimization (desirability)
Each fitted response is mapped to an individual desirability d(Yᵢ) ∈ [0,1] using
one-sided Derringer–Suich functions (`desirability` in R):

- i will add importance  weight to sugar because we care mostly for dryness/no stuck fermentation

- ethanol → `dMax`
- residual sugar → `dMin`
- fermentation time → `dMin`

The overall desirability is the geometric mean

    D = (d_ethanol · d_sugar · d_time)^(1/3)
Because the Phase 2 CFD-lite model shows a real spatial temperature gradient
(Dirichlet at tank base, Neumann at top), the scalar-T assumption is validated
**after** the optimization, not inside the DoE loop->so we follow these steps :
1. Take the optimum setpoint T* from the desirability optimization.
2. Run the Phase 2 spatial model **once** with T* as boundary condition.
3. Extract the resulting temperature gradient across tank height (ΔT, base vs top)
   over the fermentation window.