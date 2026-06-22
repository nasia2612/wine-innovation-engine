# because we know the rate of change we solve the system of equations using solve_ivp
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# state vector :the numbers that we want to track during the simulation
STATE_NAMES = ["X", "Xt", "S", "N", "E"]
#   X  = viable yeast biomass (g/L)
#   Xt = total yeast biomass  (g/L)
#   S  = sugar (g/L)                  ← glucose + fructose
#   N  = assimilable nitrogen (mg/L)  ← if it finishes → stuck
#   E  = ethanol (g/L)                ← product + inhibitor

#'The ethanol yield is not 0.51 but ~0.47-0.49 because carbon is shared between ethanol, biomass and glycerol

# Parameters (white wine preset, literature values) ─────
params = {
    # Growth
    "mu_max": 0.10,  # max specific growth rate (1/h) — λευκό ~14°C
    "KN": 50.0,  # half-saturation for nitrogen (mg/L) — Monod KS
    "kE": 0.035,  # ethanol inhibition constant — Coleman et al.
    "kd": 0.001,  # death rate (1/h) — viability decay
    # Stoichiometry
    "Yxn": 83.3,  # nitrogen consumed per biomass (mg N / g X)
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
    "N0": 300.0,  # initial YAN (mg/L) — adequate nitrogen
    "E0": 0.0,  # no ethanol at start
    "T": 14.0,  # fermentation temperature (°C)
}


# rhs() takes the current state and parameters and returns the 5 derivatives (dX/dt, dXt/dt, dS/dt, dN/dt, dE/dt) as a list or array.


def rhs(t, y, p):  # t as time,y the state vector,p as the parameters
    X, Xt, S, N, E = y  # unpack the state vector

    # modulating factors : nitrogen limitation, ethanol inhibition, temperature effect
    # AS F_N increases, growth rate increases, as F_N decreases, growth rate decreases
    # nitrogen limitation (monod type)    , N=nitrogen concentration, KN=half-saturation constant for nitrogen
    f_N = N / (
        p["KN"] + N
    )  # as N increases, f_N approaches 1, as N decreases, f_N approaches 0

    # 2. Ethanol inhibition — όσο ανεβαίνει το E, πέφτει η δράση
    f_E = np.exp(-p["kE"] * E)

    # 3. Temperature — Gaussian γύρω από το βέλτιστο
    f_T = np.exp(-(((p["T"] - p["T_opt"]) / p["T_std"]) ** 2))  # ✅

    #   # RATES

    # specific growth rate  -how fast the biomass grows per unit of existing biomass
    mu = p["mu_max"] * f_N * f_E * f_T  # modulated by the three factors

    # sugar uptake rate per biomass
    # wwe model the sugar uptake as a Monod-type function of sugar concentration (S) modulated by ethanol and temperature, but not directly by nitrogen, because nitrogen affects growth more than uptake.

    v_s = (
        p["vmax_s"] * (S / (p["Ks"] + S)) * f_E * f_T
    )  # modulated by ethanol and temperature

    # dervatives of the state variables
    dX = mu * X - p["kd"] * X
    dXt = mu * X
    dS = -v_s * X
    dN = -p["Yxn"] * mu * X
    dE = p["Yes"] * v_s * X

    return [dX, dXt, dS, dN, dE]


# so i will solve the system of equations using solve_ivp, which is a numerical integrator for ordinary differential equations (ODEs). This will allow us to simulate the fermentation process over time and observe how the state variables change.


def run_simulation(p=params, t_end=300):
    # t=0 to t_end hours

    y0 = [p["X0"], p["Xt0"], p["S0"], p["N0"], p["E0"]]  # initial state vector

    sol = solve_ivp(
        fun=rhs,
        t_span=(0, t_end),
        y0=y0,
        args=(p,),
        method="RK45",  # Runge-Kutta
        dense_output=True,  #  λείες καμπύλες
        max_step=1.0,  # step max 1 hour
    )

    return sol


def plot_simulation(sol):
    """Σχεδιάζει τις 5 καμπύλες."""
    t = sol.t
    X, Xt, S, N, E = sol.y

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("White Wine Fermentation — Digital Twin", fontsize=14)

    axes[0, 0].plot(t, X, color="green")
    axes[0, 0].set_title("Viable Biomass X (g/L)")
    axes[0, 1].plot(t, S, color="orange")
    axes[0, 1].set_title("Sugar S (g/L)")
    axes[0, 2].plot(t, E, color="purple")
    axes[0, 2].set_title("Ethanol E (g/L)")
    axes[1, 0].plot(t, N, color="blue")
    axes[1, 0].set_title("Nitrogen N (mg/L)")
    axes[1, 1].plot(t, Xt, color="gray")
    axes[1, 1].set_title("Total Biomass Xt (g/L)")

    for ax in axes.flat:
        ax.set_xlabel("Time (h)")
        ax.grid(True, alpha=0.3)

    axes[1, 2].axis("off")  # κενό panel
    plt.tight_layout()
    plt.savefig("results/figures/fermentation_simulation.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    sol = run_simulation()
    plot_simulation(sol)
