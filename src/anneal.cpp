#include <cstdlib>
#include <iostream>
#include <complex>
#include <cstdlib>
#include "quboIsingHMC.h"

using namespace std;

/*
 * 
 */

int main(int argc, char** argv)
{
  ising ising;
  int numberOfMDSteps;
  double beta;
  double mass;
  int numOfTrajs,saveFrequency,traj,ergJumpFrequency;
  int index;
  int numOfTherm;

  if(argc !=2){
    std::cout << "# Missing command line argument!" << std::endl;
    std::cout << "# Usage:: ./anneal KH_inputFile" << std::endl;
    exit(1);
  }

  // the first thing we do is load in the connectivity matrix and external field.  We also symmetrize the connectivity matrix in the process
  ising.readKandH(argv[1]);

  // since Jason and Christopher define the hamiltonian using a different sign convention, we need to multiply the K and h terms with an overall minus sign
  for(int i=0; i < ising.Lambda;++i) {
    for(int j=0;j < ising.Lambda; ++j) ising.K[i][j] *= -1;
    ising.h[i] *= -1;
  }

  // Now this parameter is IMPORTANT!  This is the overall mass term that must be added to the symmetrized connectivity matrix that makes the HS
  // transformation stable.  In principle, one needs to know what this parameters is before running the code, and this parameter depends on
  // eigenvalues of the symmetrized connectivity matrix.  We should probably allow the user to load this value in through the input file.  Maybe later. . .
  mass = 6.9; //   This is the coefficient "C" in my notes. Must be greater than 6.76978649856507 for simple toy problem.

  // set number of MD steps.  Trajectory length is always set to 1
  numberOfMDSteps = 10;  // this gives a pretty high acceptance rate
  ergJumpFrequency = -100; // please ignore this for now.  I did not discuss this in the notes.  Setting negative means no erg jumps.

  // set inverse "temperature"
  beta = 1.0;

  // initialize remaining parameters before running. . .
  ising.initialize(beta,mass,numberOfMDSteps,ergJumpFrequency);
  
  // some run parameters. . .
  numOfTherm = 500000;
  numOfTrajs = 500000;
  saveFrequency=10;

  double *m = new double [numOfTrajs/saveFrequency];
  double *e = new double [numOfTrajs/saveFrequency];
  double *psi = new double [numOfTrajs/saveFrequency];
  double *accP = new double [numOfTrajs/saveFrequency];
  index = 0;

  double betaStart,betaEnd,deltaBeta;

  // here we start the annealing process. 
  betaStart=.2;  // we start at some high temperature
  betaEnd=beta;  // and this is our ending temperature
  deltaBeta=(betaEnd-betaStart)/numOfTherm;  // and we change in these small increments
  
  for(traj=0;traj<=numOfTherm;traj++) {
    ising.reset(betaStart+traj * deltaBeta, numberOfMDSteps,ergJumpFrequency);  // this call resets the temperature, num of MD steps, and ergJump frequency
    ising.hmcThermTraj(traj);  // this does one hmc thermal trajectory  (i.e. it always accepts)
  }
  ising.reset(beta, numberOfMDSteps, ergJumpFrequency);
  for(traj=0;traj<=numOfTherm;traj++) {
    ising.hmcTraj(traj);  // now run with accept/reject
  }

  // ok, now start taking statistics. . .
  for(traj=0;traj<=numOfTrajs;traj++){
    ising.hmcTraj(traj);
    if(traj%saveFrequency==0) {
      m[index] = ising.calcM();  // average magnetization per site
      e[index] = ising.calcE();  // extensive energy
      psi[index] = ising.mean(ising.psi,ising.Lambda);  // this is the mean value of the HS aux. field
      accP[index] = ising.mean(ising.acceptP,100);
      std::cout << traj << " " << ising.mean(ising.acceptP,100) << " " << psi[index] << " " << m[index] << " " << e[index];
      for (int j = 0; j<ising.Lambda; j++) std::cout << " " << ising.calcSi(j);  // this loops through each spin site
      std::cout << std::endl;
      index += 1;
    }
  }

  // now do error analysis on observables. . .
  
  double M,E;

  // NOTE:  I didn't code up any bootstrap routines.  YOU need to do this!
  
  M = ising.mean(m,numOfTrajs/saveFrequency);
  E = ising.mean(e,numOfTrajs/saveFrequency);
  std::cout << "# beta = " << beta << std::endl;
  std::cout << "#----------------" << std::endl;
  std::cout << "# M = " << M << std::endl;
  std::cout << "# E = " << E << std::endl;
  std::cout << "# accept. rate = " << ising.mean(accP,numOfTrajs/saveFrequency) << std::endl;
  
  return 0;
}
