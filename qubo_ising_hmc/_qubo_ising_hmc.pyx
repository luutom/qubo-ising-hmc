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
            const double,
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
        const double mass_in,
        const int MDsteps_in,
        const int ergJumps_in
    ):
        Lambda_in = len(h_in)

        assert len(K_in) == Lambda_in

        C_in = 6

        self._cobj = new ising(
            Lambda_in,
            K_in,
            h_in,
            mathcalE_in,
            C_in,
            beta_in,
            mass_in,
            MDsteps_in,
            ergJumps_in
        )
