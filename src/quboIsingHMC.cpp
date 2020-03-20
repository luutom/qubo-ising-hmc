#include "quboIsingHMC.h"
#include <cmath>
#include <iostream>
#include <vector>

ising::ising() {
}

ising::ising(
  const int Lambda_in,  // total number of sites (= L^dim)
  const std::vector<std::vector<double> > K_in,  // connectivity matrix (actually stores inverse of connectivity matrix)
  const std::vector<double> h_in, // external h-field
  const double mathcalE_in,  // this is an overall shift to the Hamiltonian
  const double C_in, // mass term to regulate the connectivity matrix
  const double beta, // self explanatory (incorporates beta)
  const int MDsteps,
  const int ergJumps
) : Lambda{Lambda_in}, mathcalE{mathcalE_in} {
    K = new double* [Lambda];   // set up temp array to store connectivity matrix
    h = new double [Lambda];
    for(int i=0; i < Lambda;++i) {
      K[i] = new double [Lambda];
      for(int j=0;j < Lambda; ++j) K[i][j] = K_in[i][j];
      h[i] = h_in[i];
    }
    initialize(beta, C_in, MDsteps, ergJumps);
}

ising::ising(const ising& orig) {
}

ising::~ising() {
}

int ising::reset(double Beta, int MDsteps, int ergJumps){
  // resets certain parameters
  beta = Beta;
  sqrtBeta = sqrt(beta);

  nMD = MDsteps;
  epsilon = 1./nMD;
  ergJumpFreq = ergJumps;

  vu = 1.0; // adiabatic parameter for external field

  return 0;
}

int ising::reset(double Beta, int MDsteps, int ergJumps, double lambda){
  // resets certain parameters
  beta = Beta;
  sqrtBeta = sqrt(beta);

  nMD = MDsteps;
  epsilon = 1./nMD;
  ergJumpFreq = ergJumps;

  vu = lambda; // adiabatic parameter for external field

  return 0;
}

std::vector<std::vector<double> > const ising::K_mat(){
    std::vector<std::vector<double> > kk(Lambda, std::vector<double>(Lambda, 0));
    std::cout<<"Not implemented yet"<<std::endl;
    return kk;
}

int ising::initialize(double Beta, double mass, int MDsteps, int ergJumps) {
  /*
    Sets up arrays and so on
   */

  beta = Beta;
  sqrtBeta = sqrt(beta);
  C = mass;

  ergJumpFreq = ergJumps;
  nMD = MDsteps;
  epsilon = 1./nMD;

  totNumErgJumps=0;
  ergJumpTallies[0]=0;
  ergJumpTallies[1]=0;
  
  // allocate arrays. . .
  psi = new double [Lambda]; sampleGaussian(psi);
  psiNew = new double [Lambda];
  p = new double [Lambda];   sampleGaussian(p);
  pdot = new double [Lambda];
  for(int i=0;i<100;i++) acceptP[i]=0.0;

  // now add the C term to the connectivity matrix
  for(int i=0;i < Lambda;++i) K[i][i] += C;

  thresh = .000001; // set threshold for 0 matrix elements
  setUpSparseMatrix(K, Lambda);  // ok, this was assuming that K is sparse.  In the future, this most likely will not be needed anymore.

  // after putting K in sparse format, no longer need K itself (again, this part will be removed in the future)
  for (int i = 0; i < Lambda; i++){
    delete [] K[i];
  }
  delete [] K;

  /*
    Calculate constants here. . . only calculate if H != 0
  */
  double *b = new double [Lambda];
  k = new double [Lambda];

  for (int i=0;i<Lambda;i++) {
    b[i]=h[i];
    k[i]=0.0;
  }

  int itcount;
  double err;
  std::cout << "# Solve K.k=h . . ." << std::endl;
  linbcg(Lambda,b,k,3,1.e-9,10000,&itcount,&err);
  std::cout << "# . . . success!" << std::endl;

  kappa = 0.0;
  for (int i=0;i<Lambda;i++) kappa += k[i];
  kappa /= Lambda;
  std::cout << "# Kappa = " << kappa << std::endl;

  varphi = new double [Lambda]; // this stores K * psi
  varphiNew = new double [Lambda]; // proposed K*psiNew
  varphi2 = new double [Lambda]; // used for pdot calculation
  varphi3 = new double [Lambda]; // also used for pdot calculation

  return 0;
}

double ising::calcSi(int i){
  // calculates polarization on site i

  return psi[i]/sqrtBeta-vu*k[i];
}

double ising::calcM(){
  // calculates magnetization

  M1 = mean(psi,Lambda)/sqrtBeta - vu*kappa;
  return M1;

}

double ising::calcE(){
  // calculates energy

  E1 =  mathcalE + Lambda*C/2.;
  for (int i=0;i< Lambda; i++)
    E1 +=  vu*h[i]*vu*k[i]/2. - vu*h[i]*psi[i]/2./sqrtBeta - varphi[i]*tanh(sqrtBeta*varphi[i])/2.0/sqrtBeta;
  return E1;

}

int ising::ergJump(){
  // three possible erg jumps:  
  //        psi --> -psi
  //        psi --> -psi + 2*sqrtBeta*k
  // recall that k = invK*h
  int ergJ;

  ergJ = discrete(generator); // choose the cases randomly (but with equal weigth)

  switch(ergJ) {
  case 0:
    for (int i = 0; i < Lambda; i++)
      psiNew[i] = -psi[i];
    break;
  case 1:
    for (int i = 0; i < Lambda; i++)
      psiNew[i] = -psi[i] + 2*sqrtBeta*vu*k[i];
    break;
  }
    
  return ergJ;
}

int ising::hmcTraj(int traj){
  // perform one HMC trajectory (w/ accept/reject)
  double Hstart,Hend;
  double Sstart;
  int whichErgJump = 0;

  // first sample momentum with gaussian
  sampleGaussian(p);
  Hstart = calcH();
  Sstart = actS;

  if(traj%ergJumpFreq == 0 && ergJumpFreq > 0){
    // do ergodicity jump
    whichErgJump = ergJump();
    totNumErgJumps += 1;
  } else {
    // now do leapfrog integration
    leapfrog();
  }

  Hend = calcH(p,psiNew);
  if(uniform(generator) <= exp(-(Hend-Hstart))) {
    // accept!
    for(int i=0;i<Lambda;i++) {
      psi[i] = psiNew[i];
      varphi[i] = varphiNew[i];
    }
    acceptP[traj%100]=1.0;
    if(traj%ergJumpFreq == 0 && ergJumpFreq > 0) ergJumpTallies[whichErgJump] += 1;
  } else {
    actS = Sstart; // reset action to old value
    acceptP[traj%100]=0.0;
  }

  return 0;
}

int ising::hmcThermTraj(int traj){
  // perform one HMC trajectory and always accept! (thermalization)

  // first sample momentum with gaussian
  sampleGaussian(p);
  calcH();

  if(traj%ergJumpFreq == 0 && ergJumpFreq > 0){
    // do ergodicity jump
    ergJump();
  } else {
    // now do leapfrog integration
    leapfrog();
  }

  calcH(p,psiNew);  // always accept during thermalization
  for(int i=0;i<Lambda;i++) {
    psi[i] = psiNew[i];
    varphi[i] = varphiNew[i];
  }
  acceptP[traj%100]=1.0;

  return 0;
}

int ising::leapfrog() {

  for(int i=0;i<Lambda;i++)
    psiNew[i] = psi[i];

  // first half step
  calcPdot(psiNew);
  for (int i=0;i<Lambda;i++)
    p[i] += (epsilon/2.0)*pdot[i];

  for (int n=0;n<nMD;n++) {
    for (int i = 0; i<Lambda;i++)
      psiNew[i] += epsilon * p[i];
    calcPdot(psiNew);
    for (int i = 0;i<Lambda;i++)
      p[i] += epsilon * pdot[i];
  }

  // correct for overshoot. . .
  calcPdot(psiNew);
  for(int i=0;i<Lambda;i++)
    p[i] -= (epsilon/2.0)*pdot[i];

  return 0;
}

void ising::thermalize(
		       const double beta,
		       const size_t numOfTherm,
		       const int numberOfMDSteps,
		       const int ergJumpFrequency
		       ){
  // here we just thermalize, accepting everything 

  reset(beta, numberOfMDSteps,ergJumpFrequency);  // this call resets the temperature, num of MD steps, and ergJump frequency
  for(size_t traj=0;traj<=numOfTherm;traj++) {
    hmcThermTraj(traj);  // this does one hmc thermal trajectory  (i.e. it always accepts)
  }
}

void ising::anneal(
		   const double initBeta,
		   const double finalBeta,
		   const size_t numOfTherm,
		   const int numberOfMDSteps,
		   const int ergJumpFrequency
){
    // here we start the annealing process.
    const double betaStart= initBeta;  // we start at some high temperature
    const double betaEnd= finalBeta;  // and this is our ending temperature
    const double deltaBeta=(betaEnd-betaStart)/numOfTherm;  // and we change in these small increments

    for(size_t traj=0;traj<=numOfTherm;traj++) {
      reset(betaStart+traj * deltaBeta, numberOfMDSteps,ergJumpFrequency);  // this call resets the temperature, num of MD steps, and ergJump frequency
      hmcThermTraj(traj);  // this does one hmc thermal trajectory  (i.e. it always accepts)
    }
    //    reset(finalBeta, numberOfMDSteps, ergJumpFrequency,lambda);
    //    for(size_t traj=0;traj<=numOfTherm;traj++) {
    //      hmcTraj(traj);  // now run with accept/reject
    //    }
}

void ising::run_hmc(
		    const double beta,
		    const size_t numOfTrajs,
		    const int numberOfMDSteps,
		    const int ergJumpFrequency,
		    const size_t saveFrequency
		    ){
    // ok, now start taking statistics.
    std::vector<double> config(Lambda, 0);

    energy.clear();
    acceptance.clear();
    configs.clear();

    reset(beta, numberOfMDSteps, ergJumpFrequency);

    for(size_t traj=0;traj<=numOfTrajs;traj++){
      hmcTraj(traj);
      if(traj%saveFrequency==0) {
        energy.push_back(calcE());  // extensive energy
        acceptance.push_back(mean(acceptP,100));
        for(int i=0;i<Lambda;i++){
            config[i] = psi[i];
        }
        configs.push_back(config);
      }
    }
}

void ising::annealNoH(
		   const double initBeta,
		   const double finalBeta,
		   const size_t numOfTherm,
		   const int numberOfMDSteps,
		   const int ergJumpFrequency
){
    // here we start the annealing process but set H=0
    const double betaStart= initBeta;  // we start at some high temperature
    const double betaEnd= finalBeta;  // and this is our ending temperature
    const double deltaBeta=(betaEnd-betaStart)/numOfTherm;  // and we change in these small increments

    for(size_t traj=0;traj<=numOfTherm;traj++) {
      reset(betaStart+traj * deltaBeta, numberOfMDSteps,ergJumpFrequency,0.0);  // this call resets the temperature, num of MD steps, and ergJump frequency
      hmcTraj(traj);  // this does one hmc trajectory 
    }
    //    reset(finalBeta, numberOfMDSteps, ergJumpFrequency,lambda);
    //    for(size_t traj=0;traj<=numOfTherm;traj++) {
    //      hmcTraj(traj);  // now run with accept/reject
    //    }
}

void ising::turnOnH(
		    const double beta,
		    const size_t numOfTherm,
		    const int numberOfMDSteps,
		    const int ergJumpFrequency
		    ){
  // here we adiabatically turn on H
  const double vuStart= 0.0;  // we start at vu * H = 0
  const double vuEnd= 1.0;  // and we end at vu * H = H
  const double deltaVu=(vuEnd-vuStart)/numOfTherm;  // and we change in these small increments

  for(size_t traj=0;traj<=numOfTherm;traj++) {
    reset(beta, numberOfMDSteps,ergJumpFrequency,vuStart+traj * deltaVu);  // this call resets the temperature, num of MD steps, and ergJump frequency
    hmcTraj(traj);  // this does one hmc trajectory
  }
}

void ising::annealAndTurnOnH(
		   const double initBeta,
		   const double finalBeta,
		   const size_t numOfTherm,
		   const int numberOfMDSteps,
		   const int ergJumpFrequency
){
    // here we start the annealing process and adiabatically turn on H as well
    const double betaStart= initBeta;  // we start at some high temperature
    const double betaEnd= finalBeta;  // and this is our ending temperature
    const double deltaBeta=(betaEnd-betaStart)/numOfTherm;  // and we change in these small increments
    const double vuStart= 0.0;  // we start at vu * H = 0
  const double vuEnd= 1.0;  // and we end at vu * H = H
  const double deltaVu=(vuEnd-vuStart)/numOfTherm;  // and we change in these small increments

    for(size_t traj=0;traj<=numOfTherm;traj++) {
      reset(betaStart+traj * deltaBeta, numberOfMDSteps,ergJumpFrequency,vuStart+traj * deltaVu);  // this call resets the temperature, num of MD steps, and ergJump frequency
      hmcTraj(traj);  // this does one hmc trajectory 
    }
    //    reset(finalBeta, numberOfMDSteps, ergJumpFrequency,lambda);
    //    for(size_t traj=0;traj<=numOfTherm;traj++) {
    //      hmcTraj(traj);  // now run with accept/reject
    //    }
}

int ising::calcPdot(){

  sprsax(var,rld,psi,varphi2,Lambda);  // this gives varphi2[i] = K[i][j]*psi[j]
  for(int i=0; i < Lambda; i++)   // this makes the vector varphi3[i] = sqrtBeta*tanh(sqrtBeta*K[i][j]*psi[j])
    varphi3[i]=sqrtBeta*tanh(sqrtBeta*varphi2[i]);
  sprsax(var,rld,varphi3,pdot,Lambda); // this gives pdot[i] = sqrtBeta*K[i][j]*tanh(spqrtJ*varphi2[j])
  for(int i=0; i < Lambda; i++)
    pdot[i] += -varphi2[i] + vu*h[i]*sqrtBeta;

  return 0;
}

int ising::calcPdot(double *psi){

  sprsax(var,rld,psi,varphi2,Lambda);  // this gives varphi2[i] = K[i][j]*psi[j]
  for(int i=0; i < Lambda; i++)   // this makes the vector varphi3[i] = sqrtBeta*tanh(sqrtBeta*K[i][j]*psi[j])
    varphi3[i]=sqrtBeta*tanh(sqrtBeta*varphi2[i]);
  sprsax(var,rld,varphi3,pdot,Lambda); // this gives pdot[i] = sqrtBeta*K[i][j]*tanh(sqrtBeta*varphi2[j])
  for(int i=0; i < Lambda; i++)     // now add the remaining terms
    pdot[i] += -varphi2[i] + vu*h[i]*sqrtBeta;

  return 0;
}

double ising::calcH(double *p, double *psi){
  // calculates artificial hamiltonian

  artH = 0.0;
  for(int i=0; i < Lambda; i++)
    artH += p[i]*p[i]/2.0;

  artH += calcS(psi);
  return artH;
}

double ising::calcH(){
  // calculates artificial hamiltonian

  artH = 0.0;
  for(int i=0; i < Lambda; i++)
    artH += p[i]*p[i]/2.0;

  artH += calcS();
  return artH;
}

double ising::calcS(double *psi){
  // calculates effective potential S[psi]

  actS = 0.0;

  sprsax(var,rld,psi,varphiNew,Lambda);  // this gives varphi[i] = K[i][j]*psi[j]
  for(int i=0;i < Lambda; i++){
    actS += psi[i]*varphiNew[i]/2.0;
    actS -= vu*h[i]*psi[i]*sqrtBeta;
    actS -= log(2.*cosh(sqrtBeta*varphiNew[i]));
  }

  return actS;
}

double ising::calcS(){
  // calculates effective potential S[psi]

  actS = 0.0;

  sprsax(var,rld,psi,varphi,Lambda);  // this gives varphi[i] = K[i][j]*psi[j]
  for(int i=0;i < Lambda; i++){
    actS += psi[i]*varphi[i]/2.0;
    actS -= vu*h[i]*psi[i]*sqrtBeta;
    actS -= log(2.*cosh(sqrtBeta*varphi[i]));
  }

  return actS;
}

int ising::readKandH(std::string Kfile) {
  using namespace std;
  ifstream inputFile;
  string value;
  int i,j;

  inputFile.open(Kfile.c_str());
  if(inputFile.is_open()){
    getline(inputFile,value);
    sscanf(value.c_str(), "%d", &Lambda);
    K = new double* [Lambda];   // set up temp array to store connectivity matrix
    for (i = 0; i < Lambda; i++){
      K[i] = new double [Lambda];
      for (j=0; j < Lambda; j++)
	K[i][j] = 0.0;  // set all entries to be zero initially
    }
    for (i=0;i<Lambda;++i) {
      for (j=0;j<Lambda;++j) {
	inputFile >> value; K[i][j] = atof(value.c_str());
      }
      getline(inputFile,value);
    }
    h = new double [Lambda];
    for (j=0;j<Lambda;++j) {
      inputFile >> value; h[j] = atof(value.c_str());
    }
    getline(inputFile,value);
    getline(inputFile,value);
    sscanf(value.c_str(), "%lf", &mathcalE);
  } else {
    cout << "# Could not load connectivity matrix!"<< endl;
    return 1;
  }
  // now symmetrize the K matrix
  for(i=0;i<Lambda;++i) {
    for(j=i;j<Lambda;++j) {
      K[i][j] += K[j][i];
      K[j][i]=K[i][j];
    }
  }

  inputFile.close();

  return 0;
}

int ising::sampleGaussian(double *p) {

  for (int i=0;i<Lambda; i++) p[i]=normal(generator);

  return 0;
}

double ising::mean(double *m, int dim){
  double answer=0.0;

  for (int i=0;i<dim;i++)
    answer += m[i];

  return answer/dim;
}

double ising::variance(double *e, int dim){
  double em;
  double e2=0.0;

  em=mean(e,dim);

  for(int i=0;i<dim;i++)
    e2 += e[i]*e[i];
  e2 /= dim;

  return e2-em*em;
}

int ising::setUpSparseMatrix(double **m, int n)
{
  //  int nmax;

  nmax = getNmax(m, n);
  var  = new double [nmax];
  rld  = new int [nmax];
  sprsin(m, n, thresh, nmax, var, rld);

  return 0;
}

void ising::sprsin(double **a, int n, double thresh, int nmax,
		   double sa[], int ija[])
{
  for (int j=1;j<=n;j++) sa[j-1]=a[j-1][j-1];
  ija[1-1]=n+2;
  int k=n+1;
  for (int i=1;i<=n;i++) {
    for (int j=1;j<=n;j++) {
      if (fabs(a[i-1][j-1]) >= thresh && i != j) {
	if (++k > nmax) std::cout << "sprsin: nmax too small" << std::endl;
	sa[k-1]=a[i-1][j-1];
	ija[k-1]=j;
      }
    }
    ija[i+1-1]=k+1;
  }
}

int ising::getNmax(double **a, int n){
  int i,j;
  unsigned long k;

  k=n+1;
  for (i=1;i<=n;i++) {
    for (j=1;j<=n;j++) {
      if (std::abs(a[i-1][j-1]) >= thresh && i != j) {
	++k;
      }
    }
  }
  return k;
}

void ising::sprsax(double sa[], int ija[], double x[], double b[],
		   int n)
{
  int i,k;

  if (ija[1 - 1] != n+2) std::cout << "sprsax: mismatched vector and matrix" << std::endl;
  for(i=1;i<=n;i++) {
    b[i - 1]=sa[i - 1]*x[i - 1];
    for(k=ija[i - 1];k<=ija[i+1 - 1]-1;k++)
      b[i - 1] += sa[k - 1]*x[ija[k - 1] - 1];

  }
}

void ising::sprstx(double sa[], int ija[], double x[], double b[],
		   int n)
{
  int i,j,k;

  if(ija[1 - 1] != n+2) std::cout << "sprstx: mismatched vector and matrix" << std::endl;
  for(i=1;i<=n;i++) b[i - 1]=sa[i - 1]*x[i - 1];
  for(i=1;i<=n;i++) {
    for(k=ija[i - 1];k<=ija[i+1 - 1]-1;k++) {
      j=ija[k - 1];
      b[j - 1] += sa[k - 1]*x[i - 1];
    }
  }
}

void ising::asolve(int n, double *b, double *x, int useless){

  int i;

  for(i = 0; i < n; i++) x[i]=(var[i] != 0.0 ? b[i]/var[i] : b[i]);
}

void ising::atimes(int n, double *x, double *r, int itrnsp) {

  if(itrnsp) sprstx(var,rld,x,r,n);
  else sprsax(var,rld,x,r,n);
}

double ising::snrm(int n, double *sx, int itol){
  int i, isamax;
  double ans;

  if(itol <=3) {
    ans = 0.0;
    for (i=0;i<n;i++) ans += sx[i]*sx[i];
    return sqrt(ans);
  } else {
    isamax=0;
    for(i=0;i<n;i++){
      if(fabs(sx[i]) > fabs(sx[isamax])) isamax = i;
    }
    return fabs(sx[isamax]);
  }
}

#define EPS 1.0e-14

void ising::linbcg(int n, double b[], double x[], int itol, double tol,
		   int itmax, int *iter, double *err)
{
  int j;
  double bnrm,dxnrm,xnrm,zm1nrm,znrm;
  double akden,bknum,bk,bkden,ak;
  double *p, *pp, *r, *rr, *z, *zz;

  p  = new double [n];
  pp = new double [n];
  r  = new double [n];
  rr = new double [n];
  z  = new double [n];
  zz = new double [n];

  bkden = bnrm = 1.0;

  *iter=0;
  atimes(n,x,r,0);
  for (j=1;j<=n;j++) {
    r[j - 1]=b[j - 1]-r[j - 1];
    rr[j - 1]=r[j - 1];
  }
  znrm=1.0;
  if (itol == 1) {
    bnrm=snrm(n,b,itol);
  } else if (itol == 2) {
    asolve(n,b,z,0);
    bnrm=snrm(n,z,itol);
  }
  else if (itol == 3 || itol == 4) {
    asolve(n,b,z,0);
    bnrm=snrm(n,z,itol);
    asolve(n,r,z,0);
    znrm=snrm(n,z,itol);
  }
  asolve(n,r,z,0);
  while (*iter <= itmax) {
    ++(*iter);
    zm1nrm=znrm;
    asolve(n,rr,zz,1);
    for (bknum=0.0,j=1;j<=n;j++)
      {bknum += z[j - 1]*rr[j - 1];
      }
    if (*iter == 1) {
      for (j=1;j<=n;j++) {
	p[j - 1]=z[j - 1];
	pp[j - 1]=zz[j - 1];
      }
    }
    else {
      if(bkden==bknum) {
	bk = 1.0;
      } else {
	bk=bknum/bkden;
      }
      for (j=1;j<=n;j++) {
	p[j - 1]=bk*p[j - 1]+z[j - 1];
	pp[j - 1]=bk*pp[j - 1]+zz[j - 1];
      }
    }
    bkden=bknum;
    atimes(n,p,z,0);
    for (akden=0.0,j=1;j<=n;j++) {
      akden += z[j - 1]*pp[j - 1];
    }
    if(akden==bknum) {
      ak = 1.0;
    } else {
      ak=bknum/akden;
    }
    atimes(n,pp,zz,1);
    for (j=1;j<=n;j++) {
      x[j - 1] += ak*p[j - 1];
      r[j - 1] -= ak*z[j - 1];
      rr[j - 1] -= ak*zz[j - 1];
    }
    asolve(n,r,z,0);
    if (itol == 1 || itol == 2) {
      znrm=1.0;
      *err=snrm(n,r,itol)/bnrm;
    } else if (itol == 3 || itol == 4) {
      znrm=snrm(n,z,itol);
      if (fabs(zm1nrm-znrm) > EPS*znrm) {
	dxnrm=std::abs(ak)*snrm(n,p,itol);
	*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
      } else {
	*err=znrm/bnrm;
	continue;
      }
      xnrm=snrm(n,x,itol);
      if (*err <= 0.5*xnrm) *err /= xnrm;
      else {
	*err=znrm/bnrm;
	continue;
      }
    }
    if (*err <= tol) {
      std::cout << "# bicg :: iter = "<< *iter << " : err = "<< *err << std::endl;
      break;
    }
  }

  delete []  p;
  delete [] pp;
  delete []  r;
  delete [] rr;
  delete []  z;
  delete [] zz;
}

#undef EPS
