from openfermion import get_interaction_operator
from openfermion import get_sparse_operator, get_ground_state
from openfermion.hamiltonians import fermi_hubbard
import scipy.integrate as integrate
import scipy
import numpy as np


def compute_energy_density(energy, x_dimension, y_dimension, chemical_potential):
    """Computes energy density from raw energy."""
    return (energy / (x_dimension * y_dimension)) + chemical_potential


def compute_energy_from_density(density, x_dimension, y_dimension, chemical_potential):
    """Computes energy from energy density."""
    return (x_dimension * y_dimension) * (density - chemical_potential)


def calculate_ground_state_for_2_D_fermi_hubbard(
    tunneling_energy,
    coulomb_interaction_energy,
    x_dimension=None,
    y_dimension=1,
    magnetic_field=0,
):
    """
    Calculates the exact density of energy for 1D Fermi-Hubbard model of finite length.
    It works only for a half-filling case.

    For the general case
        Args:
            tunneling_energy (float): Tunneling energy
            coulomb_interaction_energy (float): Coulomb interaction energy.
            x_dimension (int): x dimension of the FH model
            y_dimension (int): y dimension of the FH model
            magnetic_field (float): strength of the magnetic field
        Returns:
            (float): energy density
    """
    if x_dimension is None and y_dimension == 1:
        print("WARNING!")
        print(
            "Value of the magnetic field is not taken into account in calculation for infite 1D Fermi Hubbard model."
        )
        return calculate_exact_density_of_energy_for_infinite_1D_fermi_hubbard(
            tunneling_energy, coulomb_interaction_energy
        )

    chemical_potential = coulomb_interaction_energy / 2
    energy_data = []

    hubbard_model = fermi_hubbard(
        x_dimension,
        y_dimension,
        tunneling_energy,
        coulomb_interaction_energy,
        chemical_potential,
        magnetic_field,
        False,
        False,
    )

    # Get scipy.sparse.csc representation.
    sparse_operator = get_sparse_operator(hubbard_model)

    return get_ground_state(sparse_operator)[1]


def calculate_exact_density_of_energy_for_2_D_fermi_hubbard(
    tunneling_energy,
    coulomb_interaction_energy,
    x_dimension=None,
    y_dimension=1,
    magnetic_field=0,
):
    """
    Calculates the exact density of energy for 1D Fermi-Hubbard model of finite length.
    It works only for a half-filling case.

    For the general case
        Args:
            tunneling_energy (float): Tunneling energy
            coulomb_interaction_energy (float): Coulomb interaction energy.
            x_dimension (int): x dimension of the FH model
            y_dimension (int): y dimension of the FH model
            magnetic_field (float): strength of the magnetic field
        Returns:
            (float): energy density
    """
    if x_dimension is None and y_dimension == 1:
        print("WARNING!")
        print(
            "Value of the magnetic field is not taken into account in calculation for infite 1D Fermi Hubbard model."
        )
        return calculate_exact_density_of_energy_for_infinite_1D_fermi_hubbard(
            tunneling_energy, coulomb_interaction_energy
        )

    chemical_potential = coulomb_interaction_energy / 2
    energy_data = []

    hubbard_model = fermi_hubbard(
        x_dimension,
        y_dimension,
        tunneling_energy,
        coulomb_interaction_energy,
        chemical_potential,
        magnetic_field,
        False,
        False,
    )

    # Get scipy.sparse.csc representation.
    sparse_operator = get_sparse_operator(hubbard_model)

    # It works only for the case with one electron per site.
    energy_density = (
        get_ground_state(sparse_operator)[0] / (x_dimension * y_dimension)
        + chemical_potential
    )

    return energy_density


def calculate_exact_density_of_energy_for_infinite_1D_fermi_hubbard(
    tunneling_energy, coulomb_interaction_energy
):
    """
    Calculates the exact density of energy for an infinite 1D Fermi-Hubbard model.
        Args:
            tunneling_energy (float): Tunneling energy
            coulomb_interaction_energy (float): Coulomb interaction energy.
        Returns:
            (float): energy density
    """

    def integrand(x, t, u):
        return (
            -4
            * t
            * scipy.special.jv(0, x)
            * scipy.special.jv(1, x)
            / (x * (1 + np.exp(x * u / (2 * t))))
        )

    return integrate.quad(
        integrand, 0, np.inf, args=(tunneling_energy, coulomb_interaction_energy)
    )[0]


def get_fermi_hubbard_hamiltonian(
    x_dimension,
    y_dimension,
    tunneling,
    coulomb,
    chemical_potential=0.0,
    magnetic_field=0.0,
    periodic=True,
    spinless=False,
    particle_hole_symmetry=False,
):
    """
    Generates the Hamiltonian corresponding to the Fermi-Hubbard Model for
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
    Returns:

    """
    hamiltonian = fermi_hubbard(
        x_dimension,
        y_dimension,
        tunneling,
        coulomb,
        chemical_potential=chemical_potential,
        magnetic_field=magnetic_field,
        periodic=periodic,
        spinless=spinless,
        particle_hole_symmetry=particle_hole_symmetry,
    )

    return get_interaction_operator(hamiltonian)


def get_fermi_hubbard_ordering(
    x_dimension,
    y_dimension,
    qpu_x_dimension,
    qpu_y_dimension,
    ordering_type="sycamore-horizontal",
    initial_x=0,
    initial_y=0,
):
    """
    Generate an ordering of the qubits (mapping of spin-orbitals to qubits)
    for the Fermi-Hubbard model with dimensions x_dimension by y_dimension
    on a qubit array of size qpu_x_dimension by qpu_y_dimension.

    The ordering of the qubits in the qubit array is organized by rows from
    left to right. For example an array of 2 x 4 has qubit numbers:

    0 1
    2 3
    4 5
    6 7

    The set of qubits above could represent a Fermi-Hubbard grid of size
    2 x 2 by a horizontal mapping:

    S0 - S1
    |    |
    S3 - S2

    S0: 0, 2
    S1: 1, 3
    S2: 5, 7
    S3: 4, 6

    Ordering: [0, 2, 1, 3, 5, 7, 4, 6] (Horizontal)

    The same grid of size 2 x 2 with initial_x=1 and initial_y=1 will be
    mapped based on the qubit mapping of a 3 x 3 Fermi-Hubbard grid,
    taking the subset corresponding to S4, S3, S8 and S7, in that order:

    S0 - S1 - S2
    |    |     |
    S5 - S4 - S3
    |    |     |
    S6 - S7 - S8

    S4: 7, 10
    S3: 8, 11
    S8: 14, 17
    S7: 13, 16

    Ordering: [7, 10, 8, 11, 14, 17, 13, 16] (Horizontal)

    Args:
        x_dimension (int): The width of the Fermi-Hubbard grid.
        y_dimension (int): The height of the Fermi-Hubbard grid.
        qpu_x_dimension (int): maximum length of the array of qubits in the x direction
        qpu_y_dimension (int): maximum length of the array of qubits in the y direction
        ordering_type (str): ordering type
        initial_x (int): Starting in the x dimension point of the Fermi-Hubbard grid.
        initial_y (int): Starting in the y dimension point of the Fermi-Hubbard grid.

    Return:
        ordering (list)
    """
    n_qubits = 2 * x_dimension * y_dimension

    ordering = []
    # Horizontal with oblique sites
    if ordering_type == "sycamore-horizontal":
        # compute how many rows of sites will be needed.
        # In the horizontal-oblique case, the use of the x dimension
        # is maximized
        n_rows = (x_dimension * y_dimension) // qpu_x_dimension
        # The residue is the number of sites appearing in the last
        # row, which will not be full
        residue = (x_dimension * y_dimension) % qpu_x_dimension
        if residue > 0:
            n_rows += 1
        # generates dimers row by row
        for counter in range(initial_y, n_rows + initial_y):
            # the shift is the number of qubits that have been already utilized
            shift = 2 * (counter) * (qpu_x_dimension + initial_x)
            # limit_x is the number of sites that will get assigned in this passing
            if residue > 0 and counter == n_rows - 1:
                limit_x = residue
            else:
                limit_x = qpu_x_dimension
            # even rows are numbered from right to left
            if (counter - initial_y) % 2 == 0:
                for n in range(initial_x, limit_x + initial_x):
                    # the first index in the tuple is the upper qubit
                    # the second one is the bottom qubit.
                    # Each row of sites requires two rows of qubits in the
                    # oblique mapping.
                    ordering.extend(
                        [shift + n, n + (qpu_x_dimension + initial_x) + shift]
                    )
            # odd rows are numbered from left to right
            elif (counter - initial_y) % 2 == 1:
                for n in reversed(range(initial_x, limit_x + initial_x)):
                    ordering.extend(
                        [shift + n, n + (qpu_x_dimension + initial_x) + shift]
                    )
        # Computes the final size of the grid of qubits required to map all the
        # desired sites.
        required_x = qpu_x_dimension + initial_x
        required_y = 2 * (n_rows + initial_y)

    # Vertical with oblique sites
    elif ordering_type == "sycamore-vertical":
        # Compute the number of columns of sites. For this mapping, we saturate the
        # y-dimension of the qpu
        n_columns = (x_dimension * y_dimension) // (qpu_y_dimension // 2)
        # The residue is the number of sites appearing in the last
        # column, which will not be full
        residue = (x_dimension * y_dimension) % (qpu_y_dimension // 2)
        if residue > 0:
            n_columns += 1

        for counter in range(initial_x, n_columns + initial_x):
            # limit_y is the number of sites that will get assigned in this passing
            if residue > 0 and counter == n_columns - 1:
                limit_y = residue
            else:
                limit_y = qpu_y_dimension // 2
            # even columns are numbered from top to bottom
            if (counter - initial_x) % 2 == 0:
                for n in range(initial_y, limit_y + initial_y):
                    ordering.extend(
                        [
                            (2 * n) * (qpu_x_dimension + initial_x) + counter,
                            (2 * n + 1) * (qpu_x_dimension + initial_x) + counter,
                        ]
                    )
            # odd columns are number from bottom to top
            elif (counter - initial_x) % 2 == 1:
                for n in reversed(range(initial_y, limit_y + initial_y)):
                    ordering.extend(
                        [
                            (2 * n) * (qpu_x_dimension + initial_x) + counter,
                            (2 * n + 1) * (qpu_x_dimension + initial_x) + counter,
                        ]
                    )

        required_x = n_columns + initial_x
        required_y = qpu_y_dimension + 2 * initial_y

    # Horizontal sites with interleaved ordering
    elif ordering_type == "sycamore-interleaved":

        # the number of sites that can be fitted into the x dimension of the qpu
        half_effective_x = qpu_x_dimension // 2
        x_residue = qpu_x_dimension % 2
        # the number of sites that can be fitted in the y dimension of the qpu
        effective_y = (2 * x_dimension * y_dimension) // (2 * half_effective_x)
        y_residue = (2 * x_dimension * y_dimension) % (2 * half_effective_x)

        if y_residue > 0:
            effective_y += 1

        # number of qubits that would be needed for the requested fermi-hubbard dimension
        required_x = 2 * (half_effective_x + initial_x)
        required_y = effective_y + 2 * initial_y

        # In this ordering, the sites are assigned in two moments:
        # First, we follow a horizontal mapping skipping one row of sites each pass
        # until it reaches the maximum dimension of the qpu. Then, during the second moment,
        # tt goes in a reversed horizontal through the rows of sites that were skipped in the first moment.
        counter = -1
        # First moment
        for y in range(initial_y, effective_y + initial_y, 2):
            counter += 1
            if y_residue > 0 and counter == effective_y - 1:
                limit_x = y_residue // 2
            else:
                limit_x = half_effective_x
            # from right to left
            if counter % 2 == 0:
                x_range = range(initial_x, limit_x + initial_x)
            # from left to right
            else:
                x_range = reversed(range(initial_x, limit_x + initial_x))
            for x in x_range:
                ordering.extend(
                    [y * qpu_x_dimension + 2 * x, y * qpu_x_dimension + 2 * x + 1]
                )
        # second moment
        for y in reversed(range(initial_y + 1, effective_y + initial_y, 2)):
            counter += 1
            if y_residue > 0 and counter == effective_y - 1:
                limit_x = y_residue // 2
            else:
                limit_x = half_effective_x
            # from right to left
            if counter % 2 == 0:
                x_range = range(initial_x, limit_x + initial_x)
            # from left to right
            else:
                x_range = reversed(range(initial_x, limit_x + initial_x))
            for x in x_range:
                ordering.extend(
                    [y * qpu_x_dimension + 2 * x, y * qpu_x_dimension + 2 * x + 1]
                )

    # check if the required size of the qpu is larger than the available size
    if required_x * required_y > (qpu_x_dimension + initial_x) * (
        qpu_y_dimension + 2 * initial_y
    ):
        raise Warning(
            "The number of required qubits is larger than the available number of qubits using the {0} mapping".format(
                ordering_type
            )
        )

    return ordering
