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
        double *k; # used in various places (most important for version II)
        double kappa; # constant for version II (mean value of k[i])
        double beta;
        double C;
        double mathcalE;
        vector[vector[double]] K_mat();
        vector[double] h;
        int Lambda;
        vector[double]  energy;
        vector[double]  acceptance;
        vector[vector[double]] configs;
        void thermalize(const size_t, const int, const int);
        void run_hmc(const size_t, const size_t);


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


    @property
    def kappa(self) -> float:
        """The kappa parameter
        """
        return self._cobj.kappa

    @property
    def k(self) -> np.ndarray:
        """The kappa parameter
        """
        return np.array(self._cobj.k[0])

    @property
    def c(self) -> float:
        """
        """
        return self._cobj.C

    @property
    def epsilon(self) ->float:
        """
        """
        return self._cobj.mathcalE

    @property
    def Lambda(self) ->int:
        """
        """
        return self._cobj.Lambda

    @property
    def K(self) -> np.ndarray:
        """
        """
        return self._cobj.K_mat()

    @property
    def h(self) -> np.ndarray:
        """
        """
        return np.array(self._cobj.h[0])


    def thermalize(self,const size_t numOfTherm, const int numberOfMDSteps, const int ergJumpFrequency=-100):
        """
        """
        self._cobj.thermalize(numOfTherm,numberOfMDSteps, ergJumpFrequency)

    def run_hmc(self, const size_t numOfTrajs, const size_t saveFrequency=10):
        """
        """
        self._cobj.run_hmc(numOfTrajs, saveFrequency)

    @property
    def configs(self) -> np.ndarray:
        """
        """
        return np.array(self._cobj.configs)

    @property
    def acceptance(self) -> np.ndarray:
        """
        """
        return np.array(self._cobj.acceptance)


    @property
    def energy(self) -> np.ndarray:
        """
        """
        return np.array(self._cobj.energy)
