# distutils: language = c++

import numpy as np
cimport numpy as cnp
from libcpp.vector cimport vector


cdef extern from "quboIsingHMC.h":
    cdef cppclass ising:
        ising( # initialize with parameters to avoid reading files
            const int,  # total number of sites (= L^dim)
            const vector[vector[double]],  # connectivity matrix (actually stores inverse of connectivity matrix)
            const vector[double], # external h-field
            const double,  # this is an overall shift to the Hamiltonian
            const double, # mass term to regulate the connectivity matrix
            const double, # self explanatory (incorporates beta)
            const int,
            const int
        ) except +


cdef class Ising:
    cdef ising* _cobj


    def __cinit__(
        self,
        const vector[vector[double]] K_in,  # connectivity matrix (actually stores inverse of connectivity matrix)
        const vector[double] h_in, # external h-field
        const double mathcalE_in,  # this is an overall shift to the Hamiltonian
        const double beta_in, # self explanatory (incorporates beta)
        const int MDsteps_in = 10,
        const int ergJumps_in = -100
    ):
        Lambda_in = len(h_in)

        assert len(K_in) == Lambda_in
        assert np.sum(np.abs(np.transpose(K_in) - K_in)) < 1.e-12

        eigs = np.linalg.eigvalsh(K_in)

        C_in = max(0.0, eigs.min()) + 0.1

        self._cobj = new ising(
            Lambda_in,
            K_in,
            h_in,
            mathcalE_in,
            C_in,
            beta_in,
            MDsteps_in,
            ergJumps_in
        )
