#include <cmath>
#include <complex>
#include <cstdlib>
#include <string>
#include <time.h>
#include <random>
#include <iostream>

#ifndef QUBO_ISING_HMC
#define	QUBO_ISING_HMC

class ising {
public:
  ising();
  ising(const ising& orig);
  virtual ~ising();

  std::default_random_engine generator;
  std::discrete_distribution<int> discrete {1,1};
  std::uniform_real_distribution<double> uniform; //(0.0,1.0);
  std::normal_distribution<double> normal; //(0.0,1.0);  /* for some reason, I can't initialize in the class.  But default has mean=0 and stddev=1.0 */

  int dim; // dimension of system
  int L; // number of sites in in one direction (this code will assume that the volume is symmetric in all directions)
  int Lambda;  // total number of sites (= L^dim)
  double Hshift;  // this is an overall shift to the Hamiltonian
  double beta;  // self explanatory (incorporates beta)
  double sqrtBeta; // sqrt(beta)
  double **K = NULL;  // connectivity matrix (actually stores inverse of connectivity matrix)
  double mass; // mass term to regulate the connectivity matrix
  double *phi;  // hubbard auxilliary field
  double *phiNew; // new proposed field
  double *ppsi; // used only in version II
  double *ppsi2;
  double *ppsi3;
  double *h; // external h-field
  double *p;  // conjugate field to phi
  double *pdot; // time derivative (MD) of p
  double actS;  // action S[phi]
  double artH;  // artificial Hamiltonian: p^2/2 + S[phi]
  double *k; // used in various places (most important for version II)
  double kappa; // constant for version II (mean value of k[i])

  std::string version; // stores vesion of HMC calculation
  
  double acceptP[100]; // array that keeps tally of accept/reject. . .
  
  // other functions
  int initialize(int L, int dim ,double beta, double mass, int MDsteps, int ergJumps);  // this routine allocates all arrays and sets up the connectivity matrix
  int reset(double beta, int MDsteps, int ergJumps); // resets certain parameters

  int ergJumpFreq; // frequency in which to do ergJumps (i.e. phi --> -phi)
  int nMD; // number of MD steps per trajectory
  double epsilon; // = 1/nMD
  double calcH();// calculates artificial Hamiltonian
  double calcH(double *p, double *phi);// calculates artificial Hamiltonian (overloaded)
  double calcS(); // calculates effective potential S[phi]
  double calcS(double *phi); // calculates effective potential S[phi] (overloaded)
  int calcPdot(); // calculates pdot
  int calcPdot(double *phi); // calculates pdot (overloaded)
  int leapfrog();
  int hmcTraj(int traj); // perform one HMC trajectory (w/ accept/reject)
  int ergJump(); // do phi ---> -phi jump

  double calcSi(int i); // calcualte polarization at site i
  double calcM();   // <  m  >  (unsigned)
  double M1;
  double calcE(); // average energy
  double E1;

  //private:

  // arrays, functions, and variables needed for sparse matrix format, including inversion and multiplication . . .
  int nmax;
  double thresh;
  double *var = NULL;
  int *rld = NULL;
  int setUpSparseMatrix(double **m, int n);
  void sprsin(double **a, int n, double thresh, int nmax,double sa[], int ija[]);
  int getNmax(double **a, int n);
  void sprsax(double sa[], int ija[], double x[], double b[],int n);
  void sprstx(double sa[], int ija[], double x[], double b[],int n);
  void asolve(int n, double *b, double *x, int useless);
  void atimes(int n, double *x, double *r, int itrnsp);
  double snrm(int n, double *sx, int itol);
  void linbcg(int n, double b[], double x[], int itol, double tol,int itmax, int *iter, double *err);
  
  // this is a recursive routine that calculates the index for spin at site {x_0, x_1, . . ., x_(dim-1)}
  int getIndexFromCoordinates(int *coordinates, int dim);  
  int getCoordinatesFromIndex(int *coordinates, int Index); // gets d-dimensional coordinates give index

  // set up the K matrix
  int setupK();
  int setupJasonsK();

  // sampling and simple statistical routines
  int sampleGaussian(double *p);
  double mean(double *p, int dim);
  double variance(double *e, int dim);
};

#endif	/* QUBO_ISING_HMC */

