# Phase 2: CFD-lite 1D heat transport in fermentation tank
# - 1D heat equation (method of lines) along tank height z
#   - Dirichlet BC at base (z=0): T = T_coolant (cooling jacket)
#   - Neumann BC at top (z=L):  dT/dz = 0 (adiabatic, no heat loss)
#   - Source term: Q_ferm = DeltaH_ferm * |dS/dt| (exothermic fermentation)

import numpy as np

# Sources:
#  Bogard et al. 2020, Foods 9:865 — tank geometry, wine thermal props
#   Nelson & Boulton 2024, Fermentation


params = {
    "H": 2.2,  # m, tank height (Bogard 2020)
    "N": 20,  # layers — start here; dz = H/(N-1) ≈ 0.116 m
    "alpha": 1.9e-5,  # m²/s, effective (convective) (we make it bigger than pure diffusion to account for convective mixing)
    "rho": 1080.0,  # kg/m³, fermenting must (Nelson & Boulton),density of the fermenting liquid
    "Cp": 3800.0,  # J/kg/K, approximated  ,specific heat capacity of the fermenting liquid
    "T_coolant": 12.0,  # °C, jacket/cellar (Bogard 2020)
    "T_initial": 12.0,  # °C
    "DELTA_H_FERM": 561.6,  # J/g sugar (101.2 kJ/mol ÷ 180.16),Enthalpy of Fermentation ,the heat released per gram of sugar fermented, calculated from the enthalpy change of the fermentation reaction.
}


# derived quantities
params["dz"] = params["H"] / (params["N"] - 1)  # N-1 gaps
params["alpha_h"] = (
    params["alpha"] * 3600
)  # convert it to hours, because the ODE solver works in hours, and the heat equation is in seconds. So we multiply by 3600 to convert from seconds to hours. This way, the units are consistent when we use the heat equation in our simulation.


# we are gonna use the method of lines to solve the heat equation. The method of lines is a technique for solving partial differential equations (PDEs) by discretizing the spatial domain into a finite number of points (or layers) and then solving the resulting system of ordinary differential equations (ODEs) in time. In our case, we are going to discretize the height of the fermentation tank into N layers, and then we will solve the heat equation for each layer as a function of time.
# Converts the heat PDE  ∂T/∂t = α·∂²T/∂z²  into N coupled ODEs


def heat_rhs(t, T, params):

    N = params["N"]
    dz = params["dz"]
    alpha = params["alpha_h"]

    # dTdt-> ndarray of length N, representing the rate of change of temperature at each layer

    dTdt = np.zeros(N)

    # base (Dirichlet BC): T[0] = T_coolant
    dTdt[0] = 0.0
    # finite difference approximation for the interior layers (1 to N-2)
    for i in range(
        1, N - 1
    ):  # the loop covers the interior layers of the tank, from layer 1 to layer N-2 (0-indexed). The first and last layers are handled separately due to boundary conditions.
        d2T_dz2 = (T[i + 1] - 2 * T[i] + T[i - 1]) / dz**2  # the central difference
        dTdt[i] = alpha * d2T_dz2

    # we invent a ghost node above node 19 -> T[N] = T[N-1] (Neumann BC: dT/dz=0)
    d2T_top = 2 * (T[N - 2] - T[N - 1]) / dz**2
    dTdt[N - 1] = alpha * d2T_top

    return dTdt
