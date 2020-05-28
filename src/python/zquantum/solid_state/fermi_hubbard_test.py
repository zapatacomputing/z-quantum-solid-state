import unittest
from .fermi_hubbard import (calculate_exact_density_of_energy_for_2_D_fermi_hubbard, calculate_exact_density_of_energy_for_infinite_1D_fermi_hubbard, 
                            get_fermi_hubbard_hamiltonian, get_fermi_hubbard_ordering)

from openfermion import FermionOperator, get_interaction_operator

class TestFermiHubbard(unittest.TestCase):

    def test_calculate_exact_density_of_energy_for_2_D_fermi_hubbard(self):
        # Given
        tunneling_energy = 1
        coulomb_interaction_energy = 8
        x_dimension = 4
        y_dimension = 2
        magnetic_field = 0
        correct_energy_density = -0.3782403507113399

        # When
        energy_density = calculate_exact_density_of_energy_for_2_D_fermi_hubbard(tunneling_energy, 
                coulomb_interaction_energy, 
                x_dimension, 
                y_dimension, 
                magnetic_field)

        # Then
        self.assertAlmostEqual(energy_density, correct_energy_density, 5)

        # Given
        x_dimension = None
        y_dimension = 1
        magnetic_field = 0
        correct_energy_density = -0.32753053437955576

        # When
        energy_density = calculate_exact_density_of_energy_for_2_D_fermi_hubbard(tunneling_energy, 
                coulomb_interaction_energy, 
                x_dimension, 
                y_dimension, 
                magnetic_field)

        # Then
        self.assertAlmostEqual(energy_density, correct_energy_density, 5)

        # Given
        x_dimension = 2
        y_dimension = 2
        magnetic_field = 0
        correct_energy_density = -0.33005873956798837

        # When
        energy_density = calculate_exact_density_of_energy_for_2_D_fermi_hubbard(tunneling_energy, 
                coulomb_interaction_energy, 
                x_dimension, 
                y_dimension, 
                magnetic_field)

        # Then
        self.assertAlmostEqual(energy_density, correct_energy_density, 5)


    def test_calculate_exact_density_of_energy_for_infinite_1D_fermi_hubbard(self):
        # Given
        tunneling_energy = 1
        coulomb_interaction_energy = 8
        correct_energy_density = -0.32753053437955576

        # When
        energy_density = calculate_exact_density_of_energy_for_infinite_1D_fermi_hubbard(tunneling_energy, coulomb_interaction_energy)

        # Then
        self.assertAlmostEqual(energy_density, correct_energy_density, 5)

    def test_get_fermi_hubbard_hamiltonian(self):
        # Given
        fh_2by2_spinless = FermionOperator("""
        -0.5 [0^ 0] +
        4.0 [0^ 0 1^ 1] +
        4.0 [0^ 0 2^ 2] +
        -1.0 [0^ 1] +
        -1.0 [0^ 2] +
        -1.0 [1^ 0] +
        -0.5 [1^ 1] +
        4.0 [1^ 1 3^ 3] +
        -1.0 [1^ 3] +
        -1.0 [2^ 0] +
        -0.5 [2^ 2] +
        4.0 [2^ 2 3^ 3] +
        -1.0 [2^ 3] +
        -1.0 [3^ 1] +
        -1.0 [3^ 2] +
        -0.5 [3^ 3]
        """)

        interaction_op = get_interaction_operator(fh_2by2_spinless)

        # When
        test_interaction_op = get_fermi_hubbard_hamiltonian(2, 2, 1.0, 4.0, chemical_potential=0.5, spinless=True)

        # Then
        self.assertEqual(interaction_op, test_interaction_op)

        # Given
        fh_2by2_spinful = FermionOperator("""
        -0.8 [0^ 0] +
        4.0 [0^ 0 1^ 1] +
        -1.0 [0^ 2] +
        -1.0 [0^ 4] +
        -0.2 [1^ 1] +
        -1.0 [1^ 3] +
        -1.0 [1^ 5] +
        -1.0 [2^ 0] +
        -0.8 [2^ 2] +
        4.0 [2^ 2 3^ 3] +
        -1.0 [2^ 6] +
        -1.0 [3^ 1] +
        -0.2 [3^ 3] +
        -1.0 [3^ 7] +
        -1.0 [4^ 0] +
        -0.8 [4^ 4] +
        4.0 [4^ 4 5^ 5] +
        -1.0 [4^ 6] +
        -1.0 [5^ 1] +
        -0.2 [5^ 5] +
        -1.0 [5^ 7] +
        -1.0 [6^ 2] +
        -1.0 [6^ 4] +
        -0.8 [6^ 6] +
        4.0 [6^ 6 7^ 7] +
        -1.0 [7^ 3] +
        -1.0 [7^ 5] +
        -0.2 [7^ 7]
        """)
        
        interaction_op = get_interaction_operator(fh_2by2_spinful)

        # When
        test_interaction_op = get_fermi_hubbard_hamiltonian(2, 2, 1.0, 4.0,
            chemical_potential=0.5,
            magnetic_field=0.3,
            spinless=False)
        # Then
        self.assertEqual(interaction_op, test_interaction_op)

        # Given
        fh_2by2_spinful_aperiodic = FermionOperator("""
        -0.8 [0^ 0] +
        4.0 [0^ 0 1^ 1] +
        -1.0 [0^ 2] +
        -1.0 [0^ 4] +
        -0.2 [1^ 1] +
        -1.0 [1^ 3] +
        -1.0 [1^ 5] +
        -1.0 [2^ 0] +
        -0.8 [2^ 2] +
        4.0 [2^ 2 3^ 3] +
        -1.0 [2^ 6] +
        -1.0 [3^ 1] +
        -0.2 [3^ 3] +
        -1.0 [3^ 7] +
        -1.0 [4^ 0] +
        -0.8 [4^ 4] +
        4.0 [4^ 4 5^ 5] +
        -1.0 [4^ 6] +
        -1.0 [5^ 1] +
        -0.2 [5^ 5] +
        -1.0 [5^ 7] +
        -1.0 [6^ 2] +
        -1.0 [6^ 4] +
        -0.8 [6^ 6] +
        4.0 [6^ 6 7^ 7] +
        -1.0 [7^ 3] +
        -1.0 [7^ 5] +
        -0.2 [7^ 7]
        """)

        interaction_op = get_interaction_operator(fh_2by2_spinful_aperiodic)
        
        # When
        test_interaction_op = get_fermi_hubbard_hamiltonian(2, 2, 1.0, 4.0,
            chemical_potential=0.5,
            magnetic_field=0.3,
            spinless=False,
            periodic=False)

        # Then
        self.assertEqual(interaction_op, test_interaction_op)

        # Given
        fh_2by2_spinful_phs = FermionOperator("""
        4.0 [] +
        -2.8 [0^ 0] +
        4.0 [0^ 0 1^ 1] +
        -1.0 [0^ 2] +
        -1.0 [0^ 4] +
        -2.2 [1^ 1] +
        -1.0 [1^ 3] +
        -1.0 [1^ 5] +
        -1.0 [2^ 0] +
        -2.8 [2^ 2] +
        4.0 [2^ 2 3^ 3] +
        -1.0 [2^ 6] +
        -1.0 [3^ 1] +
        -2.2 [3^ 3] +
        -1.0 [3^ 7] +
        -1.0 [4^ 0] +
        -2.8 [4^ 4] +
        4.0 [4^ 4 5^ 5] +
        -1.0 [4^ 6] +
        -1.0 [5^ 1] +
        -2.2 [5^ 5] +
        -1.0 [5^ 7] +
        -1.0 [6^ 2] +
        -1.0 [6^ 4] +
        -2.8 [6^ 6] +
        4.0 [6^ 6 7^ 7] +
        -1.0 [7^ 3] +
        -1.0 [7^ 5] +
        -2.2 [7^ 7]
        """)

        interaction_op = get_interaction_operator(fh_2by2_spinful_phs)

        # When
        test_interaction_op = get_fermi_hubbard_hamiltonian(2, 2, 1.0, 4.0,
            chemical_potential=0.5,
            magnetic_field=0.3,
            spinless=False,
            particle_hole_symmetry=True)

        # Then
        self.assertEqual(interaction_op, test_interaction_op)
    
    def test_get_fermi_hubbard_ordering(self):
        # Given/When
        horizontal_ordering = get_fermi_hubbard_ordering(2, 2, 2, 4, ordering_type='sycamore-horizontal')
        vertical_ordering = get_fermi_hubbard_ordering(2, 2, 2, 4, ordering_type='sycamore-vertical')
        # Then
        self.assertEqual([0, 2, 1, 3, 5, 7, 4, 6], horizontal_ordering)
        self.assertEqual([0, 2, 4, 6, 5, 7, 1, 3], vertical_ordering)

        # Given/When
        horizontal_ordering = get_fermi_hubbard_ordering(2, 2, 2, 4, ordering_type='sycamore-horizontal',
                                                         initial_x=1, initial_y=1)
        vertical_ordering = get_fermi_hubbard_ordering(2, 2, 2, 4, ordering_type='sycamore-vertical',
                                                         initial_x=1, initial_y=1)
        # Then
        self.assertEqual([7, 10, 8, 11, 14, 17, 13, 16], horizontal_ordering)
        self.assertEqual([7, 10, 13, 16, 14, 17, 8, 11], vertical_ordering)

        # Given/When
        interleaved_ordering = get_fermi_hubbard_ordering(3, 3, 4, 5, ordering_type='sycamore-interleaved')
        # Then
        self.assertEqual(interleaved_ordering, [0, 1, 2, 3, 10, 11, 8, 9, 16, 17, 18, 19, 14, 15, 12, 13, 4, 5])

        # Given/When
        interleaved_ordering = get_fermi_hubbard_ordering(3, 4, 7, 4, ordering_type='sycamore-interleaved')
        # Then
        self.assertEqual(interleaved_ordering, [0, 1, 2, 3, 4, 5, 18, 19, 16, 17, 14, 15, 21, 22, 23, 24, 25, 26, 11, 12, 9, 10, 7, 8])

