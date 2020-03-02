"""Python interface to C++ HMC code in `src`
"""
from typing import Optional
import numpy as np

from qubo_ising_hmc._qubo_ising_hmc import Ising as _Ising  # pylint: disable=E0611


class Ising(_Ising):  # pylint: disable=R0903
    """Class wraps HMC solver for spin-glass Hamiltonian such that

    H[psi] = psi@J@psi + h@psi + offset

    and the partition function corresponds to

    Z = sum(exp(- beta H[psi]), psi)


    Attributes:
        J: Spin coupling matrix
        h: External magentic field
        offset: Offset of Hamiltonian
        beta: Inverse termperature
        md_steps: Number of molecular dynamics steps
        ergodicity_jumps: Parameter which helps with erogdicity. Defaults to -100

    Note:
        Multiplies with minus one before inserting in CPP class.

    Todo:
        * Test
        * Think about matrix symmerty (Tom uses symmertrization)
    """
