
# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import time

import numpy as np
import scipy.sparse as sparse
import matplotlib.pyplot as plt

from mesohops.trajectory.exp_noise import bcf_exp
from mesohops.util.bath_corr_functions import bcf_convert_dl_to_exp
from mesohops.trajectory.hops_trajectory import HopsTrajectory


# ---------------------------------------------------------------------
# Parameter Control Center
# ---------------------------------------------------------------------
nsites = 5                         # Number of sites in system
dist_init = 2                      # Distance of initial site to quenching site
tmax = 250.0                       # Time length of simulation [fs]


# ---------------------------------------------------------------------
# Defining the System
# ---------------------------------------------------------------------
V = 100                            # Inter-pigment coupling [cm^-1]
dE = -1000.0                       # Difference in energy level of quenching site [cm^-1]
init_site = nsites//2 + dist_init  # Initial site index (begins at 0)
quenching_site = nsites // 2       # Index of quenching site

if init_site >= nsites:
    print(f"Distance too large, the initial site {init_site} is past the end of the chain {nsites}.")
if dist_init+nsites//2 < 0:
    print(f"Distance too small, the initial site {init_site} is less than 0.")


def generate_1_particle_hamiltonian_w_quenching(nsites, V, dE):
    """
    Generates the Hamiltonian for a 1-particle system with quenching on one site to a
    lower energy state.

    Parameters
    ----------
    1. nsites: int
               Number of sites in the system.

    2. V: float
          Inter-pigment coupling [cm^-1].

    3. dE: float
           Difference in energy level of quenching site [cm^-1].

    Returns
    -------
    1. H2_sys_hamiltonian: np.array
                           The system Hamiltonian matrix.
    """

    # Explicitly define the Hamiltonian
    H2_sys_hamiltonian = (np.diag([0] + [V] * (nsites - 1), k=1) +
                          np.diag([0] + [V] * (nsites - 1), k=-1))
    H2_sys_hamiltonian[0, 0] = dE
    return H2_sys_hamiltonian

H2_sys_hamiltonian = generate_1_particle_hamiltonian_w_quenching(nsites, V, dE)

# Define initial state
psi_0 = np.zeros(nsites+1)
psi_0[init_site+1] = 1.0


# ---------------------------------------------------------------------
# Defining the Bath
# ---------------------------------------------------------------------
e_lambda_holstein = 50          # Holstein Mode Reorganization energy [cm^-1]
gamma_holstein = 500            # Holstein Mode Reorganization timescale [cm^-1]
e_lambda_quench = 200           # Quenching Coupling Mode Reorganization energy [cm^-1]
gamma_quench = 2000           # Quenching Coupling Mode Reorganization timescale [cm^-1]
temp = 300                  # Temperature [K]

# Holstein Mode GW
gw_holstein = bcf_convert_dl_to_exp(e_lambda_holstein, gamma_holstein, temp, 0)

# Quenching site coupling mode
gw_quench = bcf_convert_dl_to_exp(e_lambda_quench, gamma_quench, temp, 0)


def generate_holstein_1_particle_loperators_w_quench(nsite, quenching_site):
    """
    Generates the L-operators for one particle Holstein calculations for a linear chain
    with nsite sites.

    Parameters
    ----------
    1. nsite: int
              Number of sites in the system.

    2. quenching_site: int
                       Index of the quenching site in the system.

    Returns
    -------
    1. list_loperators: list(sparse.coo_matrix)
                        List of L-operators in COO format.
    """

    list_loperators = []
    # Iterate across Holstein modes
    for i in range(0, nsite+1):
        list_loperators.append(sparse.coo_matrix(([1.0], ([i], [i])),
                                                 (nsite+1, nsite+1)))

    # Add coupling mode for population decay from quenching site to lower energy site
    list_loperators.append(sparse.coo_matrix(([1.0, 1.0], ([quenching_site+1, 0],
                                                           [0, quenching_site+1])),
                                             (nsite+1, nsite+1)))

    return list_loperators

# Index lists of matched L-operators and associated exponential coefficients, g and w
list_lop = generate_holstein_1_particle_loperators_w_quench(nsites, quenching_site)
gw_sysbath = [(gw_holstein)] * (nsites+1) + [(gw_quench)]


# ---------------------------------------------------------------------
# Convergence Parameters
# ---------------------------------------------------------------------
ntraj = 100                     # Number of trajectories to run
kmax = 2                        # Hierarchy truncation depth
dt = 1.0                        # Simulation time step [fs]
tau = 0.25                      # Noise time step [fs]
da = 0.01                       # Auxiliary adaptivity error bound
ds = 0.01                       # State adaptivity error bound
upstep = 10                     # Steps between adaptivity checks

# ---------------------------------------------------------------------
# Initializing the Trajectories
# ---------------------------------------------------------------------
def build_HOPS_dictionaries(list_lop, T_MAX, TAU, seed=None):
    """
    Constructs the dictionaries that will define the HOPS trajectory object.

    Parameters
    ----------
    1. list_lop: list(sparse matrix)
                 A list of site-projection operators for the system-bath interaction.
    2. T_MAX: float
              The maximum time of the simulation [fs].
    3. TAU: float
            The noise time step [fs].
    4. seed: int
             The integer seed that makes the calculation reproducible.

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
    # Get the description of the bath and system-bath interaction terms
    list_gw_noise, list_l_noise = gw_sysbath, list_lop

    # Define sys_param
    sys_param = {'HAMILTONIAN': H2_sys_hamiltonian,
                 'GW_SYSBATH': list_gw_noise,
                 'L_HIER': list_l_noise,
                 'L_NOISE1': list_l_noise,
                 'ALPHA_NOISE1': bcf_exp,
                 'PARAM_NOISE1': list_gw_noise,
                 }

    # Define noise parameters.
    noise_param = {'SEED': seed,
                   'MODEL': 'FFT_FILTER',
                   'TLEN': T_MAX + 1000,  # Units: fs
                   'TAU': TAU,  # Units: fs
                   }

    # Define hierarchy parameters
    hierarchy_param = {'MAXHIER': kmax}

    # Define equation-of-motion parameters.
    eom_param = {'EQUATION_OF_MOTION': "NORMALIZED NONLINEAR"}

    # Return all dictionaries in the proper order.
    return sys_param, noise_param, hierarchy_param, eom_param


def main():
    # ---------------------------------------------------------------------
    # Running the Trajectories
    # ---------------------------------------------------------------------
    total_propagation_time = 0
    list_traj = []
    for traj_idx in range(ntraj):
        sys_param, noise_param, hierarchy_param, eom_param = (
            build_HOPS_dictionaries(list_lop, tmax, tau, seed=traj_idx))
        hops_traj = HopsTrajectory(system_param=sys_param,
                                   eom_param=eom_param,
                                   noise_param=noise_param,
                                   hierarchy_param=hierarchy_param)
        hops_traj.make_adaptive(da, ds, upstep)
        hops_traj.initialize(psi_0)
        hops_traj.propagate(tmax, dt)
        traj_result = hops_traj.storage['psi_traj']
        traj_time = hops_traj.storage.metadata['LIST_PROPAGATION_TIME'][0]
        list_traj.append(traj_result)
        total_propagation_time += traj_time

    print('Total Propagation Time:', total_propagation_time)


    # ---------------------------------------------------------------------
    # Ensemble Analysis
    # ---------------------------------------------------------------------
    # Convert trajectories from wavefunction to population
    list_pop = np.abs(list_traj) ** 2

    # Average over trajectories
    mean_pop = np.mean(list_pop, axis=0)

    # Collect all states to the left or right of the quenching site
    left_pop = np.sum(mean_pop[:, 1:quenching_site+1], axis=1)
    right_pop = np.sum(mean_pop[:, quenching_site+2:], axis=1)

    # Collect the mean population of the quenching site
    quench_site_pop = mean_pop[:, quenching_site+1]

    # Collect total chain population
    total_pop = left_pop + right_pop + quench_site_pop

    # Collect the mean population of the low energy state
    quenched_pop = mean_pop[:, 0]

    # ---------------------------------------------------------------------
    # Plotting
    # ---------------------------------------------------------------------
    # Plot the populations of the states to the left of the quenching site, the right of
    # the quenching site, and the quenching site
    t_ax = np.arange(0, tmax+dt, dt)

    fig, (ax_pop, ax_quench) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    colors = ['#A72608', '#F49D37', '#68A357', '#0197F6']

    # Plot the populations on the chain
    ax_pop.plot(t_ax, left_pop, label="Population Left of Quenching Site", color=colors[0])
    ax_pop.plot(t_ax, right_pop, label="Population Right of Quenching Site", color=colors[1])
    ax_pop.plot(t_ax, quench_site_pop, label="Quenching Site Population", color=colors[2])

    ax_pop.set_title(f'Chain Population Dynamics from {t_ax.min():.1f} fs to '
                     f'{t_ax.max():.1f} fs')
    ax_pop.set_xlim(t_ax.min(), t_ax.max())
    ax_pop.set_ylim(0.0, 1.05)
    ax_pop.set_xlabel('Time (fs)')
    ax_pop.set_ylabel('Population')
    ax_pop.spines['top'].set_visible(False)
    ax_pop.spines['right'].set_visible(False)
    ax_pop.spines['bottom'].set_linewidth(1.5)
    ax_pop.spines['left'].set_linewidth(1.5)
    ax_pop.tick_params(axis='both', which='major', labelsize=10, width=1.5, length=4.0)
    ax_pop.legend(ncol=1, frameon=False)

    # Plot the quenched population
    ax_quench.plot(t_ax, total_pop, label="Total Chain Population", color=colors[3])
    ax_quench.plot(t_ax, quenched_pop, label="Quenched Population", color=colors[0])

    ax_quench.set_title(f'Quenching Dynamics from {t_ax.min():.1f} fs to {t_ax.max():.1f} '
                        f'fs')
    ax_quench.set_xlim(t_ax.min(), t_ax.max())
    ax_quench.set_ylim(0.0, 1.05)
    ax_quench.set_xlabel('Time (fs)')
    ax_quench.set_ylabel('Population')
    ax_quench.spines['top'].set_visible(False)
    ax_quench.spines['right'].set_visible(False)
    ax_quench.spines['bottom'].set_linewidth(1.5)
    ax_quench.spines['left'].set_linewidth(1.5)
    ax_quench.tick_params(axis='both', which='major', labelsize=10, width=1.5, length=4.0)
    ax_quench.legend(ncol=1, frameon=False)

    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
