
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


