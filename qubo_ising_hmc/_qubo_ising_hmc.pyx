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
        const vector[vector[double]] J,  # connectivity matrix (actually stores inverse of connectivity matrix)
        const vector[double] h, # external h-field
        const double offset,  # this is an overall shift to the Hamiltonian
        const double beta, # self explanatory (incorporates beta)
        const int md_steps = 10,
        const int ergodicity_jumps = -100
    ):
        Lambda_in = len(h)

        self._J = J
        self._h = h
        self._offset = offset
        self._beta = beta
        self._md_steps = md_steps
        self._ergodicity_jumps = ergodicity_jumps


        assert len(J) == Lambda_in
        assert np.abs(np.diag(J)).sum() < 1.e-12

        if not np.sum(np.abs(np.transpose(J) - J)) < 1.e-12:
            print("Symmetrisizing J. Ask Tom about the factor of 1/2")
            J_sym = np.transpose(J) + J
        else:
            J_sym = J

        # Use the minus sign definition
        J_sym *= -1
        h_sym = np.array(h) * (-1)

        eigs = np.linalg.eigvalsh(J_sym)

        # Slightly shift C_in to avoid zero
        C_in = - min(0.0, eigs.min()) + 0.1

        print("J =", J_sym)
        print("h_sym =", h_sym)
        print("offset =", offset)
        print("C_in =", C_in)
        print("beta =", beta)

        self._cobj = new ising(
            Lambda_in,
            J_sym,
            h_sym,
            offset,
            C_in,
            beta,
            md_steps,
            ergodicity_jumps
        )

    @property
    def J(self):
        return self._J

    @property
    def h(self):
        return self._h

    @property
    def offset(self):
        return self._offset

    @property
    def beta(self):
        return self._beta

    @property
    def md_steps(self):
        return self._md_steps

    @property
    def ergodicity_jumps(self):
        return self._ergodicity_jumps

    @property
    def kappa(self) -> float:
        """The kappa parameter
        """
        return self._cobj.kappa

    @property
    def k(self) -> np.ndarray:
        """The kappa parameter
        """
        return np.array([self._cobj.k[i] for i in range(self.Lambda)])

    @property
    def Lambda(self) ->int:
        """
        """
        return self._cobj.Lambda

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
