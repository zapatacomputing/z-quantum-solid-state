
"""Module to generate 1D FHM Hamiltonians at half-filling."""

import numpy as np
from zquantum.solid_state.fermi_hubbard import get_fermi_hubbard_hamiltonian
from openfermion.transforms import jordan_wigner


def get_1d_fhm_hamiltonian(fhm_model_specs_list, **kwargs):
    """Generates 1D FHM Hamiltonians from a list of dictionaries, with each
    dictionary specifying a problem/Hamiltonian.

    Notes:
        Assumes t=1 (kinetic energy) and assumes half-filling (i.e. chemical potential = U/2).

    Args:
        fhm_model_specs (list[Dict]): List of pecifications for the 1D FHM problem. 
                                        Each dictionary should contain
                                        keys 'n_sites', 'n_sites'
                                        where:
                                        U is the potential energy
                                        n_sites is the number of lattice sites
    Returns:
        List of zquantum.core.qubitoperator.QubitOperator object describing the 
        1D FHM Hamiltonians at half-filling
    """
    hamiltonians = []

    for i in range(len(fhm_model_specs_list)):

        hamiltonian = get_fermi_hubbard_hamiltonian(x_dimension=fhm_model_specs_list[i]['n_sites'],
                                                    y_dimension=1,
                                                    tunneling=1,
                                                    coulomb=fhm_model_specs_list[i]['U'],
                                                    chemical_potential=fhm_model_specs_list[i]['U']/2.,
                                                    magnetic_field=0.0,
                                                    periodic=False, 
                                                    spinless=False,
                                                    particle_hole_symmetry=False)
        hamiltonian = jordan_wigner(hamiltonian)
        hamiltonians.append(hamiltonian)
    return hamiltonians


# def get_fhm_hamiltonians_over_U_range(fhm_model_specs, n_instances=10, **kwargs):
#     """Generates 1D FHM Hamiltonians with potential energy ranging from U_low 
#     to U_high with `n_instances` number of points.

#     Notes:
#         Assumes t=1 (kinetic energy) and assumes half-filling (i.e. chemical potential = U/2).

#     Args:
#         fhm_model_specs (dict): Specifications of the 1D FHM problem. It should contain
#                                 keys 'U_range', 't', 'chem_potential', 'x_dim', and 'y_dim'
#                                 where:
#                                 U_range is a list of potential energies, U
#                                 t is the kinetic energy
#                                 chem_potential is the chemical potential (U/2 for half-filling)
#                                 x_dim is the width of the grid
#                                 y_dim is the height of the grid
#         n_instances (int): Number of points to choose over range of potential energies

#     Returns:
#         List of zquantum.core.qubitoperator.QubitOperator object describing the 
#         1D FHM Hamiltonians at half-filling
#     """
#     assert len(fhm_model_specs['U_range']) == n_instances

#     hamiltonians = []

#     for i in range(n_instances):

#         hamiltonian = get_fermi_hubbard_hamiltonian(x_dimension=fhm_model_specs['x_dim'],
#                                                     y_dimension=fhm_model_specs['y_dim'],
#                                                     tunneling=fhm_model_specs['t'],
#                                                     coulomb=fhm_model_specs['U_range'][i],
#                                                     chemical_potential=fhm_model_specs['chem_potential'],
#                                                     magnetic_field=0.0,
#                                                     periodic=False, 
#                                                     spinless=False,
#                                                     particle_hole_symmetry = False)
#         hamiltonians.append(hamiltonian)

#     return hamiltonians
