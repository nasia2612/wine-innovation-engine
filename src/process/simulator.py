# because we know the rate of change we solve the system of equations using solve_ivp

# state vector :the numbers that we want to track during the simulation
STATE_NAMES = ["X", "Xt", "S", "N", "E"]
#   X  = viable yeast biomass (g/L)
#   Xt = total yeast biomass  (g/L)
#   S  = sugar (g/L)                  ← glucose + fructose
#   N  = assimilable nitrogen (mg/L)  ← if it finishes → stuck
#   E  = ethanol (g/L)                ← product + inhibitor


# Parameters (white wine preset, literature values) ─────
params = {
    # Growth
    "mu_max": 0.10,  # max specific growth rate (1/h) — λευκό ~14°C
    "KN": 50.0,  # half-saturation for nitrogen (mg/L) — Monod KS
    "kE": 0.035,  # ethanol inhibition constant — Coleman et al. fixed
    "kd": 0.001,  # death rate (1/h) — viability decay
    # Stoichiometry
    "Yxn": 0.004,  # nitrogen consumed per biomass (mg N / g X)
    "Yes": 0.47,  # ethanol yield on sugar (g E / g S) — Gay-Lussac ~0.51
    # Sugar uptake
    "Ks": 2.0,  # half-saturation for sugar (g/L)
    "vmax_s": 2.0,  # max sugar uptake rate (g S / g X / h)
    # Temperature
    "T_opt": 14.0,  # optimal temp for white wine (°C)
    "T_std": 4.0,  # Gaussian width (°C)
    # Initial conditions (white wine)
    "X0": 0.25,  # initial biomass (g/L)
    "Xt0": 0.25,
    "S0": 220.0,  # initial sugar (g/L) — ~13% ABV potential
    "N0": 200.0,  # initial YAN (mg/L) — adequate nitrogen
    "E0": 0.0,  # no ethanol at start
}
