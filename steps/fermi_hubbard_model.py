
from zquantum.solid_state.fermi_hubbard import (
    calculate_exact_density_of_energy_for_2_D_fermi_hubbard,
    get_fermi_hubbard_hamiltonian
)
from zquantum.core.utils import ValueEstimate, save_value_estimate
from qeopenfermion import save_interaction_operator


def calculate_and_save_exact_density_of_energy_for_2_D_fermi_hubbard(tunneling_energy, 
                                                                    coulomb_interaction_energy, 
                                                                    x_dimension=None,
                                                                    y_dimension=1, 
                                                                    magnetic_field=0):
    """
    Calculates and saves the exact density of energy for 1D Fermi-Hubbard model of finite length.
    It works only for a half-filling case.
    
    Args:
        tunneling_energy (float): Tunneling energy
        coulomb_interaction_energy (float): Coulomb interaction energy.
        x_dimension (int): x dimension of the FH model
        y_dimension (int): y dimension of the FH model
        magnetic_field (float): strength of the magnetic field
    Returns:
        (float): energy density
    """
    energy_density = calculate_exact_density_of_energy_for_2_D_fermi_hubbard(tunneling_energy, 
                                                                             coulomb_interaction_energy, 
                                                                             x_dimension,
                                                                             y_dimension, 
                                                                             magnetic_field)

    val_estimate = ValueEstimate(energy_density)
    save_value_estimate(val_estimate, 'value_estimate.json')


def build_and_save_fermi_hubbard_hamiltonian(x_dimension, y_dimension, tunneling, coulomb,
                                             chemical_potential=0.0, magnetic_field=0.0, 
                                             periodic=True, spinless=False, particle_hole_symmetry=False):
    """
        Generates and saves the Hamiltonian corresponding to the Fermi-Hubbard Model for 
        a number of interacting fermions on a rectangular lattice of dimensions
        x_dimension times y_dimension.
        
        The Hamiltonian has the form:
        
        H = -t \sum_{<i,j>} \sum_{\sigma} (a^{\dagger}_{i,\sigma} a_{j,\sigma}
            + a^{\dagger}_{j,\sigma} a_{i,\sigma}) 
            + U \sum_{i} a^{\dagger}_{i,up} a_{i,up}a^{\dagger}_{i,down} a_{i,down}
            + \mu \sum_{i} \sum_{\sigma} a^{\dagger}_{i,\sigma} a_{i,\sigma}
            - h \sum_{i} (a^{\dagger}_{i,up} a_{i,up} - a^{\dagger}_{i,down} a_{i,down})
            where \sigma is the spin (up or down)
            t is the tunneling amplitude
            U is the Coulomb potential
            \mu is the chemical potential
            h is the magnetic field
        
        Args:
            x_dimension (int): The width of the grid.
            y_dimension (int): The height of the grid.
            tunneling (float): The tunneling amplitude :math:`t`.
            coulomb (float): The attractive local interaction strength :math:`U`.
            chemical_potential (float, optional): The chemical potential
                :math:`\mu` at each site. Default value is 0.
            magnetic_field (float, optional): The magnetic field :math:`h`
                at each site. Default value is 0. Ignored for the spinless case.
            periodic (bool, optional): If True, add periodic boundary conditions,
                in both directions. Default is True.
            spinless (bool, optional): If True, return a spinless Fermi-Hubbard
                model. Default is False.
            particle_hole_symmetry (bool, optional): If False, the repulsion
                term corresponds to:
                \sum_{i} a^{\dagger}_{i,up} a_{i,up}a^{\dagger}_{i,down} a_{i,down}
                If true, it corresponds to:
                \sum_{i} (a^{\dagger}_{i,up} a_{i,up} - 1/2) (a^{\dagger}_{i,down}a_{i,down} - 1/2)
    """
    hamiltonian = get_fermi_hubbard_hamiltonian(x_dimension, y_dimension, tunneling, coulomb,
                                                chemical_potential, magnetic_field, periodic, 
                                                spinless, particle_hole_symmetry)

    save_interaction_operator(hamiltonian, "hamiltonian.json")
