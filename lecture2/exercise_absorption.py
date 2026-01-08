# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import numpy as np
import scipy.sparse as sparse
import matplotlib.pyplot as plt

from mesohops.util.bath_corr_functions import bcf_convert_dl_to_exp
from mesohops.trajectory.dyadic_spectra import DyadicSpectra as DHOPS
from mesohops.trajectory.dyadic_spectra import (
    prepare_spectroscopy_input_dict,
    prepare_chromophore_input_dict,
    prepare_convergence_parameter_dict,
)
from mesohops.util.spectroscopy_analysis import zeropadded_fft
from mesohops.timing.helper_functions.loperator_generation import (
    generate_spectroscopy_loperators,
)


# ---------------------------------------------------------------------
# Parameter Control Center
# ---------------------------------------------------------------------
nsite = 10                                  # Number of sites in system


# ---------------------------------------------------------------------
# Defining the System
# ---------------------------------------------------------------------
V = -100                               # Inter-pigment coupling [cm^-1]
E_A = 200                                       # Site energy A [cm^-1]
E_B = -200                                      # Site energy B [cm^-1]


def generate_alt_energy_spectroscopy_hamiltonian(nsites, V, E_A, E_B):
    """
    Generates a spectroscopy Hamiltonian (with global ground at index 0) for a linear
    chain with nearest-neighbor coupling V and alternating site energies E_A, E_B
    across the excited-state manifold.

    Parameters
    ----------
    1. nsites: int
               Number of chromophore sites in the system.
    2. V: float
          Inter-pigment coupling [cm^-1].
    3. E_A: float
            On-site energy for sublattice A [cm^-1].
    4. E_B: float
            On-site energy for sublattice B [cm^-1].

    Returns
    -------
    1. H2_sys_hamiltonian: sparse.coo_matrix
                           The system Hamiltonian (shape (nsites+1, nsites+1)).
    """
    # Off-diagonal couplings in the excited-state block (indices 1..nsites)
    list_row = list(range(2, nsites + 1)) + list(range(1, nsites))
    list_col = list(range(1, nsites)) + list(range(2, nsites + 1))
    list_val = [V] * (2 * (nsites - 1))

    # Start from coupling-only Hamiltonian
    H = sparse.coo_matrix((list_val, (list_row, list_col)), shape=(nsites + 1, nsites + 1))

    # Add alternating on-site energies on the excited-state diagonal
    diag_rows = []
    diag_cols = []
    diag_vals = []
    for site in range(1, nsites + 1):
        diag_rows.append(site)
        diag_cols.append(site)
        diag_vals.append(E_A if (site % 2 == 1) else E_B)
    H_diag = sparse.coo_matrix((diag_vals, (diag_rows, diag_cols)), shape=(nsites + 1, nsites + 1))

    return (H + H_diag).tocoo()


H2_sys_hamiltonian = generate_alt_energy_spectroscopy_hamiltonian(nsite, V, E_A, E_B)

# Set number of chromophores to be equal to number of sites
nchrom = nsite

# ---------------------------------------------------------------------
# Defining the Bath
# ---------------------------------------------------------------------
e_lambda = 100     # Reorganization energy [cm^-1]
gamma = 500       # Reorganization timescale [cm^-1]
temp = 300        # Temperature [K]

# Set GW and L-operators
list_modes = bcf_convert_dl_to_exp(e_lambda, gamma, temp)
list_loperators = generate_spectroscopy_loperators(nsite)

# ---------------------------------------------------------------------
# Convergence Parameters
# ---------------------------------------------------------------------
ntraj = 500                         # Number of trajectories in ensemble
kmax = 3                                   # Hierarchy truncation depth
dt = 2                                      # Simulation time step [fs]
da = 0.01                            # Auxiliary adaptivity error bound
ds = 0.001                                # State adaptivity error bound
upstep = 10                           # Steps between adaptivity checks



# ---------------------------------------------------------------------
# Monte Carlo Site-Sampling Parameters (for absorption)
# ---------------------------------------------------------------------
site_sampling = True                         # If True, randomly select excitation site per trajectory
site_seed_base = 5000                        # Base seed for site RNG (independent of dynamics RNG)
site_weights = None                          # Optional probability weights over sites (length = nsite); None -> uniform

# ---------------------------------------------------------------------
# Prepare Spectroscopy Dictionaries
# ---------------------------------------------------------------------
# Time dictionary for absorption (uses t_1 only)
abs_time_dict = {"t_1": 500}

# Field and site definitions
E_vec = np.array([0.0, 0.0, 1.0])            # Defined in xyz coordinates
field_dict_abs = {"E_1": E_vec}

# Chromophore input dict: transition dipoles aligned along E_vec for all sites
M2_mu_ge = np.tile(E_vec, (nchrom, 1))

bath_dict = {
    "list_lop": list_loperators,
    "list_modes": list_modes,
}
chromophore_input = prepare_chromophore_input_dict(
    M2_mu_ge, H2_sys_hamiltonian, bath_dict=bath_dict
)

# Convergence parameters dict
convergence_dict = prepare_convergence_parameter_dict(
    dt, kmax, delta_a=da, delta_s=ds, set_update_step=upstep
)

# ---------------------------------------------------------------------
# Run Ensemble Calculation and Compute Spectrum
# ---------------------------------------------------------------------
print("Running ensemble absorption calculation (linear chain, alternating site energies)...")

# Collect trajectory response functions
list_R_abs = []
list_sampled_sites = []
base_seed = 1000
for i in range(ntraj):
    dyn_seed = base_seed + i

    # Monte Carlo sample excitation site (1..nsite) for this trajectory
    if site_sampling:
        rng_site = np.random.default_rng(site_seed_base + i)
        if site_weights is not None:
            probs = np.asarray(site_weights, dtype=float)
            if probs.shape[0] != nsite:
                raise ValueError("site_weights must have length nsite")
            probs = probs / probs.sum()
            ket_site_i = int(rng_site.choice(np.arange(1, nsite + 1), p=probs))
        else:
            ket_site_i = int(rng_site.integers(1, nsite + 1))
    else:
        ket_site_i = 1  # default fixed site

    list_sampled_sites.append(ket_site_i)

    # Build spectroscopy input with sampled site
    site_dict_abs_i = {"list_ket_sites": [ket_site_i]}
    spec_input_abs_i = prepare_spectroscopy_input_dict(
        "ABSORPTION", abs_time_dict, field_dict_abs, site_dict_abs_i
    )

    # Run trajectory with deterministic dynamics seed
    dyad_abs = DHOPS(spec_input_abs_i, chromophore_input, convergence_dict, dyn_seed)
    R_t_abs = np.asarray(dyad_abs.calculate_spectrum())
    list_R_abs.append(R_t_abs)

# Diagnostics: show counts of sampled sites
if site_sampling:
    _u, _c = np.unique(list_sampled_sites, return_counts=True)
    print("Sampled excitation site counts:", dict(zip(_u.tolist(), _c.tolist())))

# Ensemble average (trim to shortest length if needed)
min_len_abs = min(len(arr) for arr in list_R_abs)
if any(len(arr) != min_len_abs for arr in list_R_abs):
    print(f"WARNING: Absorption trajectory lengths differ; trimming all to {min_len_abs} points before averaging.")
R_abs_avg = np.mean(np.vstack([arr[:min_len_abs] for arr in list_R_abs]), axis=0)

# ---------------------------------------------------------------------
# FFT to spectrum using zero-padded FFT and plot
# ---------------------------------------------------------------------
S_abs, w_abs = zeropadded_fft(R_abs_avg, dt)

plt.figure(figsize=(8, 5))
plt.plot(w_abs, S_abs / np.max(np.abs(S_abs)))
plt.xlabel("Frequency (cm$^{-1}$)")
plt.ylabel("Intensity")
plt.xlim(-1000, 1000)
plt.title(f'Ensemble-Averaged Absorption Spectrum ({ntraj} Trajectories)')
plt.tight_layout()
plt.show()
