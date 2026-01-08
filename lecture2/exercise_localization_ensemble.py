# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import numpy as np
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from mesohops.trajectory.exp_noise import bcf_exp
from mesohops.util.bath_corr_functions import bcf_convert_dl_to_exp
from mesohops.trajectory.hops_trajectory import HopsTrajectory

# ---------------------------------------------------------------------
# Parameter Control Center
# ---------------------------------------------------------------------
V = 0                                  # Inter-pigment coupling [cm^-1]
N_TRAJECTORIES = 300                  # Number of Trajectories Averaged

# ---------------------------------------------------------------------
# Defining the System
# ---------------------------------------------------------------------
N_SITES = 4                       # Number of Sites in the Linear Chain

def generate_1_particle_hamiltonian(nsite, V):
    """
    Generates the Hamiltonian for one particle Holstein calculations for a linear chain
    with nsite sites and inter-site coupling V.

    Parameters
    ----------
    1. nsite: int
               Number of sites in the system.
    2. V: float
           Inter-site coupling strength [units: cm^-1].

    Returns
    -------
    1. H2_sys_hamiltonian: sparse.coo_matrix
                           The Hamiltonian in COO format.
    """
    list_row = list(range(1, nsite))
    list_row += list(range(0, nsite-1))
    list_col = list(range(0, nsite-1))
    list_col += list(range(1, nsite))
    list_val = [V]*2*(nsite-1)

    # Create the Hamiltonian in COO Format
    H2_sys_hamiltonian = sparse.coo_matrix((list_val, (list_row, list_col)),
                                           shape=(nsite, nsite))
    return H2_sys_hamiltonian


H2_sys_hamiltonian = generate_1_particle_hamiltonian(N_SITES, V)

# ---------------------------------------------------------------------
# Defining the Bath
# ---------------------------------------------------------------------
E_LAMBDA = 100                           # Reorganization energy [cm^-1]
GAMMA = 500                           # Reorganization timescale [cm^-1]
TEMPERATURE = 300                                      # Temperature [K]


def generate_holstein_1_particle_loperators(nsite):
    """
    Generates the L-operators for one particle Holstein calculations for a linear chain
    with nsite sites.

    Parameters
    ----------
    1. nsite: int
               Number of sites in the system.

    Returns
    -------
    1. list_loperators: list(sparse.coo_matrix)
                         List of L-operators in COO format.
    """
    list_loperators = [sparse.coo_matrix(([1], ([site_ind], [site_ind])),
                       shape=(nsite, nsite)) for site_ind in range(nsite)]
    return list_loperators


list_loperators = generate_holstein_1_particle_loperators(N_SITES)
list_dl_modes = bcf_convert_dl_to_exp(E_LAMBDA, GAMMA, TEMPERATURE)


def prepare_gw_sysbath(list_lop, list_modes):
    """
    A helper function that builds the lists taken in by the HopsTrajectory object as
    the parameters of the correlation functions.

    Parameters
    ----------
    1. list_lop: list(sparse matrix)
                 A list of site-projection operators for the system-bath interaction.
    2. list_modes: list(complex)
                   A list of complex exponential modes in (g, w) form (see below).
                   Assumed here to be the same for all baths.

    Returns
    -------
    1. list_gw_noise: list(tuple)
                      A list of (g, w) modes that make up the bath correlation
                      function for the noise in the form g*np.exp(-w*t/hbar) 
                      [units: cm^-2, cm^-1].
    2. list_l_noise: list(sparse matrix)
                      A list of site-projection operators for the noise, matched to
                      the modes in list_gw_noise.

    """
    list_gw_noise = []
    list_l_noise = []

    (g_0, w_0) = list_modes

    for L2_ind_bath in list_lop:
        list_gw_noise.append([g_0, w_0])
        list_l_noise.append(L2_ind_bath)

    return list_gw_noise, list_l_noise


# ---------------------------------------------------------------------
# Convergence Parameters
# ---------------------------------------------------------------------
K_MAX = 2                                  # Hierarchy truncation depth
DT = 2                                      # Simulation time step [fs]
TAU = 0.5 if DT % 1 == 0 else DT / 2             # Noise time step [fs]


# ---------------------------------------------------------------------
# Building Dictionaries for HOPS
# ---------------------------------------------------------------------
def build_HOPS_dictionaries(list_lop, list_modes, T_MAX, seed):
    """
    Constructs the dictionaries that will define the HOPS trajectory object.

    Parameters
    ----------
    1. list_lop: list(sparse matrix)
                  A list of site-projection operators for the system-bath interaction.
    2. list_modes: list(complex)
                    A list of complex exponential modes in (g, w) form. Assumed here to
                    be the same for all baths.
    3. T_MAX: float
               Final propagation time [fs]. Used to size the noise length.
    4. seed: int
               Seed used to initialize the stochastic noise realization.

    Returns
    -------
    1. sys_param: dictionary
                   The parameters of the system, bath, and overall simulation
    2. noise_param: dictionary
                     The parameters that define how noise is handled
    3. hierarchy_param: dictionary
                        The parameters that define how the hierarchy is handled
    4. eom_param: dictionary
                   The parameters that define the equation-of-motion
    """
    # Define Bath Parameters
    list_gw_noise, list_l_noise = prepare_gw_sysbath(list_lop, list_modes)

    # Define System Parameters
    sys_param = {'HAMILTONIAN': H2_sys_hamiltonian,
                 'GW_SYSBATH': list_gw_noise,
                 'L_HIER': list_l_noise,
                 'L_NOISE1': list_l_noise,
                 'ALPHA_NOISE1': bcf_exp,
                 'PARAM_NOISE1': list_gw_noise,
                 }

    # Define Noise Parameters
    noise_param = {'SEED': seed,
                   'MODEL': 'FFT_FILTER',
                   'TLEN': T_MAX + 1000,  # Units: fs
                   'TAU': TAU,            # Units: fs
                   }

    # Define Hierarchy Parameters
    hierarchy_param = {'MAXHIER': K_MAX}

    # Define Equation-of-Motion Parameters
    eom_param = {'EQUATION_OF_MOTION': "NORMALIZED NONLINEAR"}

    # Return all Dictionaries
    return sys_param, noise_param, hierarchy_param, eom_param


# ---------------------------------------------------------------------
# Running the Trajectory
# ---------------------------------------------------------------------
def bootstrap_mean_ci(data: np.ndarray, axis: int = 0, n_boot: int = 1000,
                      ci: float = 95.0, rng=None):
    """
    Compute a bootstrap percentile confidence interval for the mean along `axis`.

    Parameters
    ----------
    1. data: np.ndarray
             Input array.

    2. axis: int, default 0
             Axis along which to bootstrap-resample the data.

    3. n_boot: int, default 1000
               Number of bootstrap resamples.

    4. ci: float, default 95.0
           Confidence level as a percentage (e.g., 95 for 95% CI).

    5. rng: np.random.Generator | None
            Random generator. If None, uses NumPy default Generator with a fixed seed
            for reproducibility.

    Returns
    -------
    1. lower_bound : np.ndarray
                     Lower bound estimate of mean population dynamics.

    2. upper_bound : np.ndarray
                     Upper bound estimate of mean population dynamics.
    """
    # Assign default RNG if necessary
    if rng is None:
        rng = np.random.default_rng(12345)

    data = np.asarray(data)
    n = data.shape[axis]

    # Move the bootstrapped axis to 0 for convenience
    data_swapped = np.moveaxis(data, axis, 0)  # shape: (n, ...rest)

    # Collect bootstrap means with shape (n_boot, ...rest)
    means = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = data_swapped[idx, ...]
        means.append(sample.mean(axis=0))
    means = np.stack(means, axis=0)

    # Define confidence interval
    alpha = 100.0 - ci
    lower = np.percentile(means, alpha/2.0, axis=0)
    upper = np.percentile(means, 100.0 - alpha/2.0, axis=0)
    return lower, upper


def run_single_seed(seed_val: int, T_MAX: float):
    sys_param, noise_param, hierarchy_param, eom_param = build_HOPS_dictionaries(
        list_loperators, list_dl_modes, T_MAX=T_MAX, seed=seed_val)

    trajectory = HopsTrajectory(
        sys_param,
        noise_param=noise_param,
        hierarchy_param=hierarchy_param,
        eom_param=eom_param,
    )

    # Begin Fully Delocalized
    psi_0 = np.ones(N_SITES, dtype=np.complex128)/2
    trajectory.initialize(psi_0)

    trajectory.propagate(T_MAX, DT)

    # Retrieve Time Axis and Psi Trajectory from Storage
    t_axis = np.asarray(trajectory.storage['t_axis'])
    psi_traj = np.asarray(trajectory.storage['psi_traj'])  # shape [nt, N_SITES]
    pops = np.abs(psi_traj)**2
    return t_axis, pops


def main():
    # Time of Simulation
    T_MAX = 300.0  # fs

    # Seeds to Iterate Over
    seeds = list(range(N_TRAJECTORIES))

    list_pops = []
    t_axis_ref = None

    for s in seeds:
        t_axis, pops = run_single_seed(int(s), T_MAX)
        if t_axis_ref is None:
            t_axis_ref = t_axis
        else:
            # Basic sanity: make sure the same length and values
            if len(t_axis) != len(t_axis_ref) or not np.allclose(t_axis, t_axis_ref):
                raise RuntimeError("Time axes differ across seeds; "
                                   "align before averaging.")
        list_pops.append(pops)

    pops_arr = np.stack(list_pops, axis=0)
    pops_mean = pops_arr.mean(axis=0)
    # Bootstrap 95% confidence intervals for the mean across seeds
    pops_ci_low, pops_ci_high = bootstrap_mean_ci(pops_arr, axis=0, n_boot=1000, ci=95.0)

    # Compute Participation Ratio (PR): PR(t) = 1 / sum_i p_i(t)^2
    pr_all = 1.0 / np.sum(pops_arr**2, axis=2)  # shape: [nseed, nt]
    pr_mean = pr_all.mean(axis=0)               # shape: [nt]

    # Combined Plot: Ensemble Averaged Population Dynamics (left) and PR (right)
    fig, (ax_pop, ax_ipr) = plt.subplots(1, 2, figsize=(14, 5),
                                         sharex=True)

    # Left: Ensemble Averaged Population Dynamics
    colors = ['#A72608', '#F49D37', '#68A357', '#0197F6']
    for i_site in range(N_SITES):
        ax_pop.plot(t_axis_ref, pops_mean[:, i_site], color=colors[i_site], lw=2.0,
                    label=f'Site {i_site+1}')
        # Shaded error band (95% bootstrap CI) around the mean population
        y_low = pops_ci_low[:, i_site]
        y_high = pops_ci_high[:, i_site]
        ax_pop.fill_between(
            t_axis_ref,
            y_low,
            y_high,
            color=colors[i_site],
            alpha=0.25,
            linewidth=0,
        )
    ax_pop.set_xlabel('Time [fs]')
    ax_pop.set_ylabel('Population')
    ax_pop.set_title(f'Ensemble-Averaged Population Dynamics ({pops_arr.shape[0]} '
                     f'Trajectories)')
    ax_pop.set_xlim(t_axis_ref.min(), t_axis_ref.max())
    ax_pop.set_ylim(0.0, 1.05)
    ax_pop.spines['top'].set_visible(False)
    ax_pop.spines['right'].set_visible(False)
    ax_pop.spines['bottom'].set_linewidth(1.5)
    ax_pop.spines['left'].set_linewidth(1.5)
    ax_pop.tick_params(axis='both', which='major', labelsize=10, width=1.5, length=4.0)
    # Legend 1: site colors (keep this when adding a second legend)
    leg1 = ax_pop.legend(ncol=2)
    ax_pop.add_artist(leg1)

    # Legend 2: meaning of line and shaded band (mean vs 95% bootstrap CI)
    proxy_mean = Line2D([0], [0], color='k', lw=2.0, label='Mean')
    proxy_ci = Patch(facecolor='0.5', alpha=0.25, edgecolor='none', label='95% Confidence Interval')
    ax_pop.legend(handles=[proxy_mean, proxy_ci], loc='upper left',
                  frameon=False)

    # Right: Ensemble Averaged PR
    ax_ipr.plot(t_axis_ref, pr_mean, color='k', lw=2.0)
    ax_ipr.set_xlabel('Time [fs]')
    ax_ipr.set_ylabel('PR')
    ax_ipr.set_title(f'Ensemble-Averaged Participation Ratio ('
                     f'{pops_arr.shape[0]} Trajectories)')
    ax_ipr.set_xlim(t_axis_ref.min(), t_axis_ref.max())
    ax_ipr.set_ylim(1.0, float(N_SITES) + 0.05)
    ax_ipr.spines['top'].set_visible(False)
    ax_ipr.spines['right'].set_visible(False)
    ax_ipr.spines['bottom'].set_linewidth(1.5)
    ax_ipr.spines['left'].set_linewidth(1.5)
    ax_ipr.tick_params(axis='both', which='major', labelsize=10, width=1.5, length=4.0)

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
