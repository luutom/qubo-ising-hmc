#include "quboIsingHMC.h"
#include <cmath>
#include <iostream>
#include "lapacke.h"

ising::ising() {
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

  return 0;
}

  
int ising::initialize(int nx, int d ,double Beta, double C, int MDsteps, int ergJumps) {
  /*
    Sets up arrays and so on
   */

  beta = Beta;
  sqrtBeta = sqrt(beta);
  dim = d;
  L = nx;
  mass = C;

  ergJumpFreq = ergJumps;
  nMD = MDsteps;
  epsilon = 1./nMD;

  Lambda = nx; // pow(L,d);  this is specific for Jason's problem

  // allocate arrays. . .
  phi = new double [Lambda]; sampleGaussian(phi);
  phiNew = new double [Lambda];
  p = new double [Lambda];   sampleGaussian(p);
  pdot = new double [Lambda];
  for(int i=0;i<100;i++) acceptP[i]=0.0;
  

  std::cout << "# Seting up Jason's connectivity . . .";
  K = new double* [Lambda];   // in the future I want to get rid of this step and just put the matrix directly into sparse format. . .
  for (int i = 0; i < Lambda; i++){
    K[i] = new double [Lambda];
    for (int j=0; j < Lambda; j++)
      K[i][j] = 0.0;  // set all entries to be zero initially
  }  
  setupJasonsK(); // sets up the connectivity matrix with variable mass on the diagonal. . .
  thresh = .000001; // set threshold for 0 matrix elements
  setUpSparseMatrix(K, Lambda);
  std::cout << " success!" << std::endl;
  // after putting K in sparse format, no longer need K itself (again, this part will be removed in the future)
  for (int i = 0; i < Lambda; i++){
    delete [] K[i]; 
  }
  delete [] K;
  /*
    Calculate constants here. . . only calculate if H != 0
  */

  double *b = new double [Lambda];
  h = new double [Lambda];
  k = new double [Lambda];
  // -2.5, -3.5, -2.5,  1. ,  2. ,  4. ,  1.
  b[0]= 2.5; b[1]= 3.5; b[2]= 2.5; b[3]=-1.0; b[4]=-2.0; b[5]=-4.0; b[6]=-1.0;  // again, this is specific to Jason's problem
  Hshift = 5.5;  // again, specific to Jason's problem
  for (int i=0;i<Lambda;i++) h[i]=b[i];
  int itcount;
  double err;
  std::cout << "# Solve K.x=(-2.5, -3.5, -2.5,  1. ,  2. ,  4. ,  1.)" << std::endl;
  linbcg(Lambda,b,k,3,1.e-9,10000,&itcount,&err);
  kappa = 0.0;
  // for (int i=0;i<Lambda;i++) std::cout << k[i]  << std::endl;
  // exit(0);
  for (int i=0;i<Lambda;i++) kappa += k[i];
  kappa /= Lambda;
  std::cout << "# Kappa = " << kappa << std::endl;

  // after constants are calculated, don't need inverse K anymore, but only K.  So re-calculate it
  ppsi = new double [Lambda]; // this stores K * phi
  ppsi2 = new double [Lambda]; // used for pdot calculation
  ppsi3 = new double [Lambda]; // also used for pdot calculation
    
  return 0;
}

double ising::calcSi(int i){
  // calculates polarization on site i
  
  return phi[i]/Lambda/sqrtBeta-k[i]/Lambda;
}

double ising::calcM(){
  // calculates magnetization

  M1 = mean(phi,Lambda)/sqrtBeta - kappa;
  return M1;

}

double ising::calcE(){
  // calculates energy

  E1 =  mass*beta/2.+Hshift*beta/Lambda;
  for (int i=0;i< Lambda; i++)
    E1 += -sqrtBeta*h[i]*phi[i]/2./Lambda - sqrtBeta*ppsi[i]*tanh(sqrtBeta*ppsi[i])/2.0/Lambda + beta*h[i]*k[i]/2./Lambda;
  return E1;
    
}

int ising::ergJump(){

  for (int i = 0; i < Lambda; i++)
    phiNew[i] = -phi[i]; //(2*discrete(generator)-1)*phi[i];

  return 0;
}

int ising::hmcTraj(int traj){
  // perform one HMC trajectory (w/ accept/reject)
  double Hstart,Hend;
  double Sstart;
  
  // first sample momentum with gaussian
  sampleGaussian(p);
  Hstart = calcH();
  Sstart = actS;

  if(traj%ergJumpFreq == 0 && ergJumpFreq > 0){
    // do ergodicity jump
    ergJump();
  } else {
    // now do leapfrog integration
    leapfrog();
  }

  Hend = calcH(p,phiNew);
  if(uniform(generator) <= exp(-(Hend-Hstart))) {
    // accept!
    for(int i=0;i<Lambda;i++) phi[i]= phiNew[i];
    acceptP[traj%100]=1.0;
  } else {
    actS = Sstart; // reset action to old value
    acceptP[traj%100]=0.0;
  }
  
  return 0;
}

int ising::leapfrog() {

  for(int i=0;i<Lambda;i++)
    phiNew[i] = phi[i];

  // first half step
  calcPdot(phiNew);
  for (int i=0;i<Lambda;i++)
    p[i] += (epsilon/2.0)*pdot[i];

  for (int n=0;n<nMD;n++) {
    for (int i = 0; i<Lambda;i++)
      phiNew[i] += epsilon * p[i];
    calcPdot(phiNew);
    for (int i = 0;i<Lambda;i++)
      p[i] += epsilon * pdot[i];
  }

  // correct for overshoot. . .
  calcPdot(phiNew);
  for(int i=0;i<Lambda;i++)
    p[i] -= (epsilon/2.0)*pdot[i];
  
  return 0;
}

int ising::calcPdot(){

  sprsax(var,rld,phi,ppsi2,Lambda);  // this gives ppsi2[i] = K[i][j]*phi[j]
  for(int i=0; i < Lambda; i++)   // this makes the vector ppsi3[i] = sqrtBeta*tanh(sqrtBeta*K[i][j]*phi[j])
    ppsi3[i]=sqrtBeta*tanh(sqrtBeta*ppsi2[i]);
  sprsax(var,rld,ppsi3,pdot,Lambda); // this gives pdot[i] = sqrtBeta*K[i][j]*tanh(spqrtJ*ppsi2[j])
  for(int i=0; i < Lambda; i++)
    pdot[i] += -ppsi2[i] + h[i]*sqrtBeta;

  return 0;
}

int ising::calcPdot(double *phi){

  sprsax(var,rld,phi,ppsi2,Lambda);  // this gives ppsi2[i] = K[i][j]*phi[j]
  for(int i=0; i < Lambda; i++)   // this makes the vector ppsi3[i] = sqrtBeta*tanh(sqrtBeta*K[i][j]*phi[j])
    ppsi3[i]=sqrtBeta*tanh(sqrtBeta*ppsi2[i]);
  sprsax(var,rld,ppsi3,pdot,Lambda); // this gives pdot[i] = sqrtBeta*K[i][j]*tanh(sqrtBeta*ppsi2[j])
  for(int i=0; i < Lambda; i++)     // now add the remaining terms
    pdot[i] += -ppsi2[i] + h[i]*sqrtBeta;

  return 0;
}

double ising::calcH(double *p, double *phi){
  // calculates artificial hamiltonian
  
  artH = 0.0;
  for(int i=0; i < Lambda; i++)
    artH += p[i]*p[i]/2.0;

  artH += calcS(phi);
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

double ising::calcS(double *phi){
  // calculates effective potential S[phi]

  actS = 0.0;

  sprsax(var,rld,phi,ppsi,Lambda);  // this gives ppsi[i] = K[i][j]*phi[j]
  for(int i=0;i < Lambda; i++){
    actS += phi[i]*ppsi[i]/2.0;
    actS -= h[i]*phi[i]*sqrtBeta;
    actS -= log(2.*cosh(sqrtBeta*ppsi[i]));
  }
  
  return actS;
}

double ising::calcS(){
  // calculates effective potential S[phi]

  actS = 0.0;

  sprsax(var,rld,phi,ppsi,Lambda);  // this gives ppsi[i] = K[i][j]*phi[j]
  for(int i=0;i < Lambda; i++){
    actS += phi[i]*ppsi[i]/2.0;
    actS -= h[i]*phi[i]*sqrtBeta;
    actS -= log(2.*cosh(sqrtBeta*ppsi[i]));
  }
  
  return actS;
}

int ising::setupJasonsK() {

  // K =
  //   {
  //    {0., 2., 1, -1, -1, -2., 0.},
  //    {2., 0., 2., -1, -1, -2., -1},
  //    {1, 2., 0., 0., -1, -2., -1},
  //    {-1, -1, 0., 0., 0., 0., 0.},
  //    {-1, -1, -1, 0., 0., 2., 0.},
  //    {-2., -2., -2., 0., 2., 0., 0.},
  //    {0., -1, -1, 0., 0., 0., 0.}
  //   };

  K[0][1] = -2.0; K[1][0] = K[0][1];
  K[0][2] = -1.0; K[2][0] = K[0][2];
  K[0][3] = 1.0; K[3][0] = K[0][3];
  K[0][4] = 1.0; K[4][0] = K[0][4];
  K[0][5] = 2.0; K[5][0] = K[0][5];

  K[1][2] = -2.0; K[2][1] = K[1][2];
  K[1][3] = 1.0; K[3][1] = K[1][3];
  K[1][4] = 1.0; K[4][1] = K[1][4];
  K[1][5] = 2.0; K[5][1] = K[1][5];
  K[1][6] = 1.0; K[6][1] = K[1][6];

  K[2][4] = 1.0; K[4][2] = K[2][4];
  K[2][5] = 2.0; K[5][2] = K[2][5];
  K[2][6] = 1.0; K[6][2] = K[2][6];

  K[4][5] = -2.0; K[5][4] = K[4][5];

  for (int i = 0; i< Lambda; i++) K[i][i] += mass;
  
  return 0;
}

int ising::setupK() {
  int *coordinates;
  int *nn;
  int nnIndex;
  
  coordinates = new int [dim];
  nn = new int [dim];

  for (int i=0; i< Lambda; i++) { // loop thru all spin sites
    getCoordinatesFromIndex(coordinates,i);  // extract coordinate of this spin
    for(int j=0; j< dim; j++) { // now loop thru nearest neighbors of spin at this coordinate
      for( int k = 0; k < dim; k++) nn[k]=coordinates[k];
      nn[j] = (nn[j]+1)%L;  // takes into account periodic boundaries
      nnIndex = getIndexFromCoordinates(nn,dim);
      K[nnIndex][i] += 1.0;
      K[i][nnIndex] += 1.0;
    };
    K[i][i] += mass;  // put mass term on diagonal. . .
  }

  delete [] coordinates;  // garbage collection. . .
  delete [] nn;

  return 0;
}

int ising::getIndexFromCoordinates(int *coordinates, int d){

  if(d > 1) {
    return L*getIndexFromCoordinates(coordinates, d-1)+coordinates[d-1];
  } else {
    return coordinates[d-1];
  }
}

int ising::getCoordinatesFromIndex(int *coordinates, int Index){
  int n;

  n=Index;
  for (int i = 0; i < dim; i++) coordinates[i]=0;

  for (int i = dim-1; i >= 0; i--) {
    coordinates[i]=n%L;
    n -= coordinates[i];
    n /= L;
  }
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

  if (ija[1 - 1] != n+2) std::cout << "sprsax: mismatched vector and matrix" << std::endl;
  for(int i=1;i<=n;i++) {
    b[i - 1]=sa[i - 1]*x[i - 1];
    for(int k=ija[i - 1];k<=ija[i+1 - 1]-1;k++)
      b[i - 1] += sa[k - 1]*x[ija[k - 1] - 1];

  }
}

void ising::sprstx(double sa[], int ija[], double x[], double b[],
		   int n)
{

  if(ija[1 - 1] != n+2) std::cout << "sprstx: mismatched vector and matrix" << std::endl;
  for(int i=1;i<=n;i++) b[i - 1]=sa[i - 1]*x[i - 1];
  for(int i=1;i<=n;i++) {
    for(int k=ija[i - 1];k<=ija[i+1 - 1]-1;k++) {
      int j=ija[k - 1];
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
