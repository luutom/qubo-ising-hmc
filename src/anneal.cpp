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
  int L,dim,numberOfMDSteps;
  double beta;
  double mass;
  int numOfTrajs,saveFrequency,traj,ergJumpFrequency;
  int index;
  int numOfTherm;

  L=7;  // How many sites per side.  Total number of sites is L^d, where d is dimension of the system
  dim=2;  // this is the dimension of the problem
  mass = 6.9; //   This is the coefficient "C" in my notes. Must be greater than 6.76978649856507 for Jason's problem

  // set number of MD steps.  Trajectory length is always set to 1
  numberOfMDSteps = 4;
  ergJumpFrequency = 100; // please ignore this for now.  I did not discuss this in the notes

  // set inverse "temperature"
  beta = 1.5;

  // initialize
  ising.initialize(L,dim,beta,mass,numberOfMDSteps,ergJumpFrequency);
  
  // some run parameters. . .
  numOfTherm = 300000;
  numOfTrajs = 500000;
  saveFrequency=10;

  double *m = new double [numOfTrajs/saveFrequency];
  double *e = new double [numOfTrajs/saveFrequency];
  double *phi = new double [numOfTrajs/saveFrequency];
  double *accP = new double [numOfTrajs/saveFrequency];
  index = 0;

  double betaStart,betaEnd,deltaBeta;

  betaStart=.2;
  betaEnd=2.0;
  deltaBeta=(betaEnd-betaStart)/numOfTherm;
  
  for(traj=0;traj<=numOfTherm;traj++) {
    ising.reset(betaStart+traj * deltaBeta, numberOfMDSteps,-1);
    ising.hmcTraj(traj);
  }
  ising.reset(beta, numberOfMDSteps, ergJumpFrequency);
  for(traj=0;traj<=numOfTherm;traj++) {
    ising.hmcTraj(traj);
  }

  
  for(traj=0;traj<=numOfTrajs;traj++){
    ising.hmcTraj(traj);
    if(traj%saveFrequency==0) {
      m[index] = ising.calcM();
      e[index] = ising.calcE()/beta;
      phi[index] = ising.mean(ising.phi,ising.Lambda);
      accP[index] = ising.mean(ising.acceptP,100);
      std::cout << traj << " " << ising.mean(ising.acceptP,100) << " " << phi[index] << " " << m[index] << " " << e[index];
      for (int j = 0; j<L; j++) std::cout << " " << L*ising.calcSi(j);
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
