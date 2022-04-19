################################################################################
# © Copyright 2020-2022 Zapata Computing Inc.
################################################################################
"""Module to generate 1D FHM Hamiltonians at half-filling."""

import numpy as np
from zquantum.core.openfermion.transforms import jordan_wigner
from zquantum.solid_state.fermi_hubbard import get_fermi_hubbard_hamiltonian


def get_1d_fhm_hamiltonian(fhm_model_specs_list, **kwargs):
    """Generates 1D FHM Hamiltonians from a list of dictionaries, with each
    dictionary specifying a problem/Hamiltonian.

    Notes:
        Assumes t=1 (kinetic energy) and assumes half-filling
        (i.e. chemical potential = U/2).

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

        hamiltonian = get_fermi_hubbard_hamiltonian(
            x_dimension=fhm_model_specs_list[i]["n_sites"],
            y_dimension=1,
            tunneling=1,
            coulomb=fhm_model_specs_list[i]["U"],
            chemical_potential=fhm_model_specs_list[i]["U"] / 2.0,
            magnetic_field=0.0,
            periodic=False,
            spinless=False,
            particle_hole_symmetry=False,
        )
        hamiltonian = jordan_wigner(hamiltonian)
        hamiltonians.append(hamiltonian)
    return hamiltonians
