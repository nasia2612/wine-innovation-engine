import numpy as np
from scipy.integrate import solve_ivp


# Phase 2: CFD-lite 1D heat transport in fermentation tank
# - 1D heat equation (method of lines) along tank height z
#   - Dirichlet BC at base (z=0): T = T_coolant (cooling jacket)
#   - Neumann BC at top (z=L):  dT/dz = 0 (adiabatic, no heat loss)
#   - Source term: Q_ferm = DeltaH_ferm * |dS/dt| (exothermic fermentation)


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

    # we invent a ghost node above node 19 -> T[N] = T[N-2] (it makes a mirror,Neumann BC: dT/dz=0)
    d2T_top = 2 * (T[N - 2] - T[N - 1]) / dz**2
    dTdt[N - 1] = alpha * d2T_top

    return dTdt


# no fermentation heat source term, just pure diffusion, to test the heat equation solver, a linear gradient (cold base,warm top) is expected to flatten over time, and the temperature should converge to a uniform value (the average of the initial temperatures) as time progresses.
#  The test checks if the final temperature profile is close to this expected uniform value, within a specified tolerance.
# without source term, the heat equation reduces to a pure diffusion problem, and the temperature should eventually become uniform throughout the tank. The test checks if the final temperature profile is close to this expected uniform value, within a specified tolerance.
def run_diffusion_test(params, T_top=24, t_end=100, n_save=200):

    N = params["N"]

    # Initial condition : Linear gradient from base to top
    T0 = np.linspace(
        params["T_coolant"], T_top, N
    )  # linear gradient from T_coolant at the base to T_top at the top
    # base (node 0) =T_coolant=12
    # top  (node N-1)=T_top=24

    # enforce dirichelt BC at the base (node 0)
    T0[0] = params["T_coolant"]

    # time pointa to record the solution
    t_eval = np.linspace(0, t_end, n_save)

    # solve the heat equation using solve_ivp
    sol = solve_ivp(
        fun=heat_rhs,
        t_span=(0, t_end),
        y0=T0,
        args=(params,),
        method="BDF",  # stiff solver(stiff we mean by when fast and slow processes coexist), because the heat equation can be stiff, especially for small dz and large alpha. The BDF method is suitable for stiff problems and can handle the rapid changes in temperature that may occur in the early stages of the simulation.
        t_eval=t_eval,
        dense_output=True,
    )

    return sol


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    sol = run_diffusion_test(params)

    # plot the temperature profile at a few times
    z = np.linspace(0, params["H"], params["N"])
    for k in [0, 10, 40, 100, 199]:  # indices into the saved snapshots
        plt.plot(sol.y[:, k], z, label=f"t = {sol.t[k]:.0f} h")

    plt.xlabel("Temperature (°C)")
    plt.ylabel("Height z (m)")
    plt.title("Stage 1: gradient relaxation (no source term)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("results/figures/phase2_diffusion_test.png", dpi=150)
    plt.show()
