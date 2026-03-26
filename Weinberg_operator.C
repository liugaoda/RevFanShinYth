#include <stdio.h>
#include <string> 
#include <iostream>  
#include <complex>
#include <math.h>
#include <fstream>
#include <cmath>
//#include <complex.h>

using namespace std;

typedef long double ldouble;
typedef complex<ldouble> cx_ldouble;
typedef vector<cx_ldouble> cx_lvec;
typedef vector<cx_lvec> cx_lmat;
typedef vector<ldouble> lvec;

typedef vector<int> in_vec;
typedef vector<in_vec> in_mat;


long double gammarecursion(long double n, long double s, long double z)
{
	long double res;
	if(n > 20.0L)
	{
		res = z;
	}
	else
	{
		res = z + (n-s) / ( 1.0L + n / gammarecursion(n+1.0L, s, z) );
	}
	return res;
}




ldouble gamma(int z, ldouble x){
	ldouble res = 0.L;
	if (z == 1){
	res = exp(-x);
	} else if (z == 2){
	res = (x + 1.L) * exp(-x);
	} else if (z == 3){
	res = (pow(x,2) + 2.L * x + 2.L) * exp(-x);
	} else if (z == 4){
	res = (pow(x,3) + 3.L*pow(x,2) + 6.L * x + 6.L) * exp(-x);
	} else if (z == 5){
	res = (pow(x,4) + 4.L*pow(x,3) + 12.L*pow(x,2) + 24.L * x + 24.L) * exp(-x);
	} else if (z == 6){
	res = (120.L + 120.L * x + 60.L * pow(x, 2) + 20.L * pow(x, 3) + 5.L * pow(x, 4) + pow(x, 5)) * exp(-x);
	} else if (z == 7){
	res = (720.L + 720.L * x + 360.L * pow(x, 2) + 120.L * pow(x, 3) + 30.L * pow(x, 4) + 6.L * pow(x, 5) + pow(x, 6)) * exp(-x);
	} else if (z == 8){
	res = (5040.L + 5040.L * x + 2520.L * pow(x, 2) + 840.L * pow(x, 3) + 210.L * pow(x, 4) + 42.L * pow(x, 5) + 7.L * pow(x, 6) + pow(x, 7)) * exp(-x);
	} else if (z == 9){
	res = (40320.L + 40320.L * x + 20160.L * pow(x, 2) + 6720.L * pow(x, 3) + 1680.L * pow(x, 4) + 336.L * pow(x, 5) + 56.L * pow(x, 6) + 8.L * pow(x, 7) + pow(x, 8)) * exp(-x);
	} else if (z == 10){
	res = (362880.L + 362880.L * x + 181440.L * pow(x, 2) + 60480.L * pow(x, 3) + 15120.L * pow(x, 4) + 3024.L * pow(x, 5) + 504.L * pow(x, 6) + 72.L * pow(x, 7) + 9.L * pow(x, 8) + pow(x, 9)) * exp(-x);
	} else if (z == 11){
	res = (3628800.L + 3628800.L * x + 1814400.L * pow(x, 2) + 604800.L * pow(x, 3) + 151200.L * pow(x, 4) + 30240.L * pow(x, 5) + 5040.L * pow(x, 6) + 720.L * pow(x, 7) + 90.L * pow(x, 8) + 10.L * pow(x, 9) + pow(x, 10)) * exp(-x);
	} else if (z == 0){
	res = 1.0L*exp(-x) / gammarecursion(1.0L, 0.0L, x);
	}
	return res;
}


complex<long double>* Matrixsq(complex<long double> *M, complex<long double> *res)
{
  
    //first row of matrixsq
    res[0] = M[0] * conj(M[0]) + M[3] * conj(M[3]) + M[6] * conj(M[6]);
    res[1] = M[1] * conj(M[0]) + M[4] * conj(M[3]) + M[7] * conj(M[6]);
    res[2] = M[2] * conj(M[0]) + M[5] * conj(M[3]) + M[8] * conj(M[6]);

    //second row of matrixsq
    res[3] = M[0] * conj(M[1]) + M[3] * conj(M[4]) + M[6] * conj(M[7]);
    res[4] = M[1] * conj(M[1]) + M[4] * conj(M[4]) + M[7] * conj(M[7]);
    res[5] = M[2] * conj(M[1]) + M[5] * conj(M[4]) + M[8] * conj(M[7]);

    //third row of matrixsq
    res[6] = M[0] * conj(M[2]) + M[3] * conj(M[5]) + M[6] * conj(M[8]);
    res[7] = M[1] * conj(M[2]) + M[4] * conj(M[5]) + M[7] * conj(M[8]);
    res[8] = M[2] * conj(M[2]) + M[5] * conj(M[5]) + M[8] * conj(M[8]);
  
    return res;
  }

complex<long double>* invmatrix(complex<long double> *A, complex<long double> *res)
{
    long double det = abs(A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6])); 
    //first row of res
    res[0] = (A[4] * A[8] - A[5] * A[7])/det;
    res[1] = (A[2] * A[7] - A[1] * A[8])/det;
    res[2] = (A[1] * A[5] - A[2] * A[4])/det;
    
    //second row of res
    res[3] = (A[5] * A[6] - A[3] * A[8])/det;
    res[4] = (A[0] * A[8] - A[2] * A[6])/det;
    res[5] = (A[2] * A[3] - A[0] * A[5])/det;
    
    //third row of res
    res[6] = (A[3] * A[7] - A[4] * A[6])/det;
    res[7] = (A[1] * A[6] - A[0] * A[7])/det;
    res[8] = (A[0] * A[4] - A[1] * A[3])/det;
  
  return res;
}


complex<long double>* multimatrix(complex<long double> *A, complex<long double> *B, complex<long double> *res)
{
    //first row of res
    res[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    res[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    res[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
    
    //second row of res
    res[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    res[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    res[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
    
    //third row of res
    res[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    res[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    res[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
  
  return res;
}

struct Coef
{
     long double h1, h4, h6, a, b, c, d, alpha, beta, lambda1, lambda2, lambda3, S, cer;
     complex<long double> h2, h3, h5, c1, c2, c3, c11, c21, c31, v11, v12, v13, v21, v22, v23, v31, v32, v33;
};

struct Coef extractcoef(complex<long double> M[9], long double cert)
{
     struct Coef res;
     
     res.h1 = real(M[0]);
     res.h4 = real(M[4]);
     res.h6 = real(M[8]);
     res.h2 = M[1];
     res.h3 = M[2];
     res.h5 = M[5];
     
     res.a = 1.;
     res.b = -(res.h1 + res.h4 + res.h6);

     res.c = res.h1 * res.h4 + res.h1 * res.h6
             +res.h4 * res.h6 - real((res.h2 * conj(res.h2)))
             - real(res.h3 * conj(res.h3)) - real(res.h5 * conj(res.h5));

     res.d = -res.h1 * res.h4 * res.h6
             + real(res.h2 * conj(res.h2)) * res.h6
             + real(res.h3 * conj(res.h3)) * res.h4
             - real(res.h2 * conj(res.h3) * res.h5 + conj(res.h2) * res.h3 * conj(res.h5)) 
             + res.h1 * real(res.h5 * conj(res.h5));
     
     res.alpha = -pow(res.b, 3) / (27 * pow(res.a, 3))
                 -res.d / (2 * res.a) + res.b * res.c / (6 * pow(res.a, 2));
     res.beta = res.c / (3 * res.a) - pow(res.b, 2) / (9 * pow(res.a, 2));
     
     if (pow(res.alpha, 2) + pow(res.beta, 3) < 0.)
     {
          res.S = cos(acos(res.alpha / pow(-res.beta, 1.5)) / 3);

	  if (cert > 0.)
	  {
          res.lambda1 = -res.b / (3 * res.a) - pow(-res.beta, 0.5) * (res.S + pow(3, 0.5) * pow(1 - pow(res.S, 2), 0.5));
          res.lambda2 = -res.b / (3 * res.a) - pow(-res.beta, 0.5) * (res.S - pow(3, 0.5) * pow(1 - pow(res.S, 2), 0.5));
          res.lambda3 = -res.b / (3 * res.a) + 2 * pow(-res.beta, 0.5) * res.S;
	  }
	  else
	    {
             res.lambda3 = -res.b / (3 * res.a) - pow(-res.beta, 0.5) * (res.S + pow(3, 0.5) * pow(1 - pow(res.S, 2), 0.5));
             res.lambda1 = -res.b / (3 * res.a) - pow(-res.beta, 0.5) * (res.S - pow(3, 0.5) * pow(1 - pow(res.S, 2), 0.5));
             res.lambda2 = -res.b / (3 * res.a) + 2 * pow(-res.beta, 0.5) * res.S;  
	    }  
	  
      if (res.lambda1>0. && res.lambda2>0. && res.lambda3>0.)
      {
          //*************************************************************************
          res.c1 = pow(pow(abs((res.lambda1-res.h6) * conj(res.h2) + res.h5 * conj(res.h3)), 2)
		   + pow(abs((res.lambda1 - res.h4) * conj(res.h3) + conj(res.h2) * conj(res.h5)),2)
		   + pow((res.lambda1 - res.h4) * (res.lambda1 - res.h6) - pow(abs(res.h5), 2), 2), 0.5);
      
          res.c2 = pow(pow(abs((res.lambda2 - res.h6) * res.h2 + res.h3 * conj(res.h5)), 2)
		   + pow(abs((res.lambda2 - res.h1) * conj(res.h5) + res.h2 * conj(res.h3)),2)
		   + pow((res.lambda2 - res.h1) * (res.lambda2 - res.h6) - pow(abs(res.h3), 2), 2), 0.5);
	  
          res.c3 = pow(pow(abs((res.lambda3 - res.h4) * res.h3 + res.h2 * res.h5), 2)
		   + pow(abs((res.lambda3 - res.h1) * res.h5 + res.h3 * conj(res.h2)), 2)
		   + pow((res.lambda3 - res.h1) * (res.lambda3 - res.h4) - pow(abs(res.h2), 2), 2), 0.5);
      
          res.c11 = real(res.c1);
          res.c21 = real(res.c2);
          res.c31 = real(res.c3);
      
          //*************************************************************************

          res.v11 = ((res.lambda1 - res.h4) * (res.lambda1 - res.h6) - pow(abs(res.h5), 2)) / res.c11;
          res.v12 = ((res.lambda1 - res.h6) * conj(res.h2) + res.h5 * conj(res.h3)) / res.c11;
          res.v13 = ((res.lambda1 - res.h4) * conj(res.h3) + conj(res.h2) * conj(res.h5)) / res.c11;

          res.v21 = ((res.lambda2-res.h6) * res.h2 + res.h3 * conj(res.h5)) / res.c21;
          res.v22 = ((res.lambda2-res.h1) * (res.lambda2 - res.h6) - pow(abs(res.h3), 2)) / res.c21;
          res.v23 = ((res.lambda2-res.h1) * conj(res.h5) + res.h2 * conj(res.h3)) / res.c21;

          res.v31 = ((res.lambda3-res.h4) * res.h3 + res.h2 * res.h5) / res.c31;
          res.v32 = ((res.lambda3-res.h1) * res.h5 + res.h3 * conj(res.h2)) / res.c31;
          res.v33 = ((res.lambda3-res.h1) * (res.lambda3 - res.h4) - pow(abs(res.h2),2)) / res.c31;
	  
          res.cer = 1.;
      }
      else
	  res.cer = 0.;
     }
  else
      res.cer = 0.;
     
  return res;
}


long double *calculate_chisq(const double *par1, long double *reschisq)
{

   long double par[12] = {par1[0], par1[1], par1[2], par1[3], par1[4], par1[5], par1[6], par1[7], par1[8], par1[9], par1[10], par1[11]};

  
  
    long double Pi = M_PI;
  
  const complex<long double> i(0.,1.);
  const complex<long double> I(0.,1.);
  complex<long double> q,q31,q32,q51,q52,q53,q54 ;
  q = exp(2.L * Pi * i * complex<long double>(par[0],par[1]));
  q31 = exp(2.L * Pi * i * complex<long double>(par[0],par[1])/3.L);
  q32 = exp(4.L * Pi * i * complex<long double>(par[0],par[1])/3.L);
  q51 = exp(2.L * Pi * i * complex<long double>(par[0],par[1])/5.L);
  q52 = exp(4.L * Pi * i * complex<long double>(par[0],par[1])/5.L);
  q53 = exp(6.L * Pi * i * complex<long double>(par[0],par[1])/5.L);
  q54 = exp(8.L * Pi * i * complex<long double>(par[0],par[1])/5.L);
  ldouble y= par[1];

//--------------------------------------//
//             定义模形式                //
//--------------------------------------//

complex<long double> Y0411,Y0431,Y0432,Y0433;
complex<long double> Y03211,Y03212,Y03231,Y03232;
complex<long double> Y0211,Y0231,Y0232,Y0233;
complex<long double> Y01211,Y01212,Y01231,Y01232;
complex<long double> Y0011,Y0031,Y0032,Y0033;
complex<long double> Y1211,Y1212,Y1231,Y1232;
complex<long double> Y211,Y231,Y232,Y233;
complex<long double> Y3211,Y3212,Y3231,Y3232;
complex<long double> Y4111,Y4121,Y431,Y432,Y433;
complex<long double> Y5211,Y5212,Y5221,Y5222,Y5231,Y5232;
complex<long double> Y611,Y6311,Y6312,Y6313,Y6321,Y6322,Y6323;


//--------------------------------------//
//              定义参数                 //
//--------------------------------------//


complex<long double> alpha, beta,gammainput, g1, g2,g3, Lamda, v;

  alpha = complex<long double>(1.L, 0.);
  beta = complex<long double>(par[2], 0.);
  gammainput = complex<long double>(par[3], 0.);
  g1 = complex<long double>(1.L, 0.);
 g2 = complex<long double>(par[4], 0.);






 const ldouble eulergamma = 0.577215664901533;

 const ldouble zeta3 = 1.2020569031595942;
 const ldouble zeta4=1.08232323371114;
const ldouble zeta5 = 1.03692775514337;
const ldouble zeta6 = 1.01734306198445;
const ldouble m4r3c1=0.017247857248905454;
const ldouble m4r3pc1=-0.017247857248906224;
const ldouble m4r5c1=-0.007992880519988832;
const ldouble cwm2r1c1=-0.29076134702187606;
const ldouble m2r3c1=-0.10769300935905161;
const ldouble m2r3pc1=0.10769300935905152;
const ldouble m2r5c1=0.05591564365805286;
const ldouble m0r3c1=0.856271381317112;
 const ldouble m0r3pc1=-0.8562713813171127;
 const ldouble m0r5c1=-0.6403749984084703;
 const ldouble cw1r6II1=-8.85224051831748;
 const ldouble cw1r6II2=2.11363312370849;

const ldouble ploygamma9and1o5 = 3.54375005874724*pow(10,12.L);
const ldouble ploygamma9and2o5 = 3.46070596386813*pow(10,9.L);
const ldouble ploygamma9and3o5 = 6.00170445943624*pow(10,7.L);
const ldouble ploygamma9and4o5 = 3.38061259232518*pow(10,6.L);

 
const cx_ldouble dirichletL522 = cx_ldouble(0.958716, 0.145566);
const cx_ldouble dirichletL542 = cx_ldouble(0.958716, -0.145566);
const cx_ldouble dirichletL526 = cx_ldouble(0.999776, 0.0142581 );
const cx_ldouble dirichletL546 = cx_ldouble(0.999776, -0.0142581);
const cx_ldouble dirichletL5210 = cx_ldouble(0.999999,0.00095963 );
const cx_ldouble dirichletL5410 = cx_ldouble(0.999999,-0.00095963);



  
// ------------------------------------------------------------------//
//                         Modular forms                             //
// ------------------------------------------------------------------//

// Weight: w=-5



 




// Weight: w=-4

Y0411=(189.L*q)/(16.L*pow(Pi, 5)) + (6237.L*pow(q, 2))/(512.L*pow(Pi, 5)) 
+ (427.L*pow(q, 3))/(36.L*pow(Pi, 5)) + (199773.L*pow(q, 
4))/(16384.L*pow(Pi, 5)) + (295407.L*pow(q, 5))/(25000.L*pow(Pi, 5)) 
+ (4697.L*pow(q, 6))/(384.L*pow(Pi, 5)) + (56727.L*pow(q, 
7))/(4802.L*pow(Pi, 5)) + (6392925.L*pow(q, 8))/(524288.L*pow(Pi, 5)) 
+ (415051.L*pow(q, 9))/(34992.L*pow(Pi, 5)) + (9748431.L*pow(q, 
10))/(800000.L*pow(Pi, 5)) + (7609707.L*pow(q, 11))/(644204.L*pow(Pi, 
5)) + (451339.L*pow(q, 12))/(36864.L*pow(Pi, 5)) + (35087283.L*pow(q, 
13))/(2970344.L*pow(Pi, 5)) + (1871991.L*pow(q, 
14))/(153664.L*pow(Pi, 5)) + (222467.L*pow(q, 15))/(18750.L*pow(Pi, 
5)) + (204573789.L*pow(q, 16))/(16777216.L*pow(Pi, 5)) + 
(134176581.L*pow(q, 17))/(11358856.L*pow(Pi, 5)) + (4565561.L*pow(q, 
18))/(373248.L*pow(Pi, 5)) + (116995725.L*pow(q, 
19))/(9904396.L*pow(Pi, 5)) + (312245199.L*pow(q, 
20))/(25600000.L*pow(Pi, 5)) + pow(y, 5) / (5.L) + (63.L*gamma(5.L, 
4.L*Pi*y))/(128.L*pow(Pi, 5)*q) + (2079.L*gamma(5.L, 
8.L*Pi*y))/(4096.L*pow(Pi, 5)*pow(q, 2)) + (427.L*gamma(5.L, 
12.L*Pi*y))/(864.L*pow(Pi, 5)*pow(q, 3)) + (66591.L*gamma(5.L, 
16.L*Pi*y))/(131072.L*pow(Pi, 5)*pow(q, 4)) + (98469.L*gamma(5.L, 
20.L*Pi*y))/(200000.L*pow(Pi, 5)*pow(q, 5)) + (4697.L*gamma(5.L, 
24.L*Pi*y))/(9216.L*pow(Pi, 5)*pow(q, 6)) + (18909.L*gamma(5.L, 
28.L*Pi*y))/(38416.L*pow(Pi, 5)*pow(q, 7)) + (2130975.L*gamma(5.L, 
32.L*Pi*y))/(4194304.L*pow(Pi, 5)*pow(q, 8)) + (415051.L*gamma(5.L, 
36.L*Pi*y))/(839808.L*pow(Pi, 5)*pow(q, 9)) + (3249477.L*gamma(5.L, 
40.L*Pi*y))/(6400000.L*pow(Pi, 5)*pow(q, 10)) + (2536569.L*gamma(5.L, 
44.L*Pi*y))/(5153632.L*pow(Pi, 5)*pow(q, 11)) + (451339.L*gamma(5.L, 
48.L*Pi*y))/(884736.L*pow(Pi, 5)*pow(q, 12)) + (11695761.L*gamma(5.L, 
52.L*Pi*y))/(23762752.L*pow(Pi, 5)*pow(q, 13)) + (623997.L*gamma(5.L, 
56.L*Pi*y))/(1229312.L*pow(Pi, 5)*pow(q, 14)) + (222467.L*gamma(5.L, 
60.L*Pi*y))/(450000.L*pow(Pi, 5)*pow(q, 15)) + (68191263.L*gamma(5.L, 
64.L*Pi*y))/(134217728.L*pow(Pi, 5)*pow(q, 16)) + 
(44725527.L*gamma(5.L, 68.L*Pi*y))/(90870848.L*pow(Pi, 5)*pow(q, 17)) 
+ (4565561.L*gamma(5.L, 72.L*Pi*y))/(8957952.L*pow(Pi, 5)*pow(q, 18)) 
+ (38998575.L*gamma(5.L, 76.L*Pi*y))/(79235168.L*pow(Pi, 5)*pow(q, 
19)) + (104081733.L*gamma(5.L, 80.L*Pi*y))/(204800000.L*pow(Pi, 
5)*pow(q, 20)) + (189.L*zeta5)/(16.L*pow(Pi, 5));
 

 Y0431=(-1647.L*q)/(416.L*pow(Pi, 5)) - (54351.L*pow(q, 2))/(13312.L*pow(Pi, 
5)) - (14641.L*pow(q, 3))/(3744.L*pow(Pi, 5)) - (1740879.L*pow(q, 
4))/(425984.L*pow(Pi, 5)) - (2574261.L*pow(q, 5))/(650000.L*pow(Pi, 
5)) - (161051.L*pow(q, 6))/(39936.L*pow(Pi, 5)) - (3460347.L*pow(q, 
7))/(873964.L*pow(Pi, 5)) - (55709775.L*pow(q, 
8))/(13631488.L*pow(Pi, 5)) - (3557581.L*pow(q, 9))/(909792.L*pow(Pi, 
5)) - (84950613.L*pow(q, 10))/(20800000.L*pow(Pi, 5)) - 
(66313161.L*pow(q, 11))/(16749304.L*pow(Pi, 5)) - (15475537.L*pow(q, 
12))/(3833856.L*pow(Pi, 5)) - (305760609.L*pow(q, 
13))/(77228944.L*pow(Pi, 5)) - (114191451.L*pow(q, 
14))/(27966848.L*pow(Pi, 5)) - (7627961.L*pow(q, 
15))/(1950000.L*pow(Pi, 5)) - (1782714447.L*pow(q, 
16))/(436207616.L*pow(Pi, 5)) - (1169253063.L*pow(q, 
17))/(295330256.L*pow(Pi, 5)) - (39133391.L*pow(q, 
18))/(9704448.L*pow(Pi, 5)) - (1019534175.L*pow(q, 
19))/(257514296.L*pow(Pi, 5)) - (2720993877.L*pow(q, 
20))/(665600000.L*pow(Pi, 5)) + pow(y, 5) / (5.L) - (549.L*gamma(5.L, 
4.L*Pi*y))/(3328.L*pow(Pi, 5)*q) - (18117.L*gamma(5.L, 
8.L*Pi*y))/(106496.L*pow(Pi, 5)*pow(q, 2)) - (14641.L*gamma(5.L, 
12.L*Pi*y))/(89856.L*pow(Pi, 5)*pow(q, 3)) - (580293.L*gamma(5.L, 
16.L*Pi*y))/(3407872.L*pow(Pi, 5)*pow(q, 4)) - (858087.L*gamma(5.L, 
20.L*Pi*y))/(5200000.L*pow(Pi, 5)*pow(q, 5)) - (161051.L*gamma(5.L, 
24.L*Pi*y))/(958464.L*pow(Pi, 5)*pow(q, 6)) - (1153449.L*gamma(5.L, 
28.L*Pi*y))/(6991712.L*pow(Pi, 5)*pow(q, 7)) - (18569925.L*gamma(5.L, 
32.L*Pi*y))/(109051904.L*pow(Pi, 5)*pow(q, 8)) - 
(3557581.L*gamma(5.L, 36.L*Pi*y))/(21835008.L*pow(Pi, 5)*pow(q, 9)) - 
(28316871.L*gamma(5.L, 40.L*Pi*y))/(166400000.L*pow(Pi, 5)*pow(q, 
10)) - (22104387.L*gamma(5.L, 44.L*Pi*y))/(133994432.L*pow(Pi, 
5)*pow(q, 11)) - (15475537.L*gamma(5.L, 
48.L*Pi*y))/(92012544.L*pow(Pi, 5)*pow(q, 12)) - 
(101920203.L*gamma(5.L, 52.L*Pi*y))/(617831552.L*pow(Pi, 5)*pow(q, 
13)) - (38063817.L*gamma(5.L, 56.L*Pi*y))/(223734784.L*pow(Pi, 
5)*pow(q, 14)) - (7627961.L*gamma(5.L, 
60.L*Pi*y))/(46800000.L*pow(Pi, 5)*pow(q, 15)) - 
(594238149.L*gamma(5.L, 64.L*Pi*y))/(3489660928.L*pow(Pi, 5)*pow(q, 
16)) - (389751021.L*gamma(5.L, 68.L*Pi*y))/(2362642048.L*pow(Pi, 
5)*pow(q, 17)) - (39133391.L*gamma(5.L, 
72.L*Pi*y))/(232906752.L*pow(Pi, 5)*pow(q, 18)) - 
(339844725.L*gamma(5.L, 76.L*Pi*y))/(2060114368.L*pow(Pi, 5)*pow(q, 
19)) - (906997959.L*gamma(5.L, 80.L*Pi*y))/(5324800000.L*pow(Pi, 
5)*pow(q, 20)) - (405.L*zeta5)/(104.L*pow(Pi, 5));

 Y0432=q31*(6561.L/(832.L*pow(Pi, 5)) + (6934977.L*q)/(851968.L*pow(Pi, 5)) + 
(13784661.L*pow(q, 2))/(1747928.L*pow(Pi, 5)) + (338409819.L*pow(q, 
3))/(41600000.L*pow(Pi, 5)) + (1218029967.L*pow(q, 
4))/(154457888.L*pow(Pi, 5)) + (7101632961.L*pow(q, 
5))/(872415232.L*pow(Pi, 5)) + (4061423025.L*pow(q, 
6))/(515028592.L*pow(Pi, 5)) + (792496629.L*pow(q, 
7))/(97450496.L*pow(Pi, 5)) + (64092775311.L*pow(q, 
8))/(8125000000.L*pow(Pi, 5)) + (2081483811.L*pow(q, 
9))/(255696896.L*pow(Pi, 5)) + (5869870821.L*pow(q, 
10))/(744357926.L*pow(Pi, 5)) + (153708857577.L*pow(q, 
11))/(18901136384.L*pow(Pi, 5)) + (227482854219.L*pow(q, 
12))/(28847086112.L*pow(Pi, 5)) + (13874802579.L*pow(q, 
13))/(1703936000.L*pow(Pi, 5)) + (241130600271.L*pow(q, 
14))/(30577756144.L*pow(Pi, 5)) + (174194018559.L*pow(q, 
15))/(21420149504.L*pow(Pi, 5)) + (1853430385977.L*pow(q, 
16))/(235019407168.L*pow(Pi, 5)) + (1287457675119.L*pow(q, 
17))/(158164877312.L*pow(Pi, 5)) + (412890743709.L*pow(q, 
18))/(52341575000.L*pow(Pi, 5)) + (2220465309975.L*pow(q, 
19))/(273044415488.L*pow(Pi, 5)) + (72171.L*gamma(5.L, (8.L*Pi*y) / 
(3.L)))/(212992.L*pow(Pi, 5)*q) + (3418281.L*gamma(5.L, (20.L*Pi*y) / 
(3.L)))/(10400000.L*pow(Pi, 5)*pow(q, 2)) + (73975275.L*gamma(5.L, 
(32.L*Pi*y) / (3.L)))/(218103808.L*pow(Pi, 5)*pow(q, 3)) + 
(88055181.L*gamma(5.L, (44.L*Pi*y) / (3.L)))/(267988864.L*pow(Pi, 
5)*pow(q, 4)) + (151631271.L*gamma(5.L, (56.L*Pi*y) / 
(3.L)))/(447469568.L*pow(Pi, 5)*pow(q, 5)) + (1552614723.L*gamma(5.L, 
(68.L*Pi*y) / (3.L)))/(4725284096.L*pow(Pi, 5)*pow(q, 6)) + 
(3613123017.L*gamma(5.L, (80.L*Pi*y) / (3.L)))/(10649600000.L*pow(Pi, 
5)*pow(q, 7)) + (1759535541.L*gamma(5.L, (92.L*Pi*y) / 
(3.L)))/(5355037376.L*pow(Pi, 5)*pow(q, 8)) + 
(13398329637.L*gamma(5.L, (104.L*Pi*y) / 
(3.L)))/(39541219328.L*pow(Pi, 5)*pow(q, 9)) + 
(22428942525.L*gamma(5.L, (116.L*Pi*y) / 
(3.L)))/(68261103872.L*pow(Pi, 5)*pow(q, 10)) + 
(75750753771.L*gamma(5.L, (128.L*Pi*y) / 
(3.L)))/(223338299392.L*pow(Pi, 5)*pow(q, 11)) + 
(7181808381.L*gamma(5.L, (140.L*Pi*y) / 
(3.L)))/(21849100000.L*pow(Pi, 5)*pow(q, 12)) + 
(44675653275.L*gamma(5.L, (152.L*Pi*y) / 
(3.L)))/(131847319552.L*pow(Pi, 5)*pow(q, 13)) + 
(126688756887.L*gamma(5.L, (164.L*Pi*y) / 
(3.L)))/(385569436928.L*pow(Pi, 5)*pow(q, 14)) + 
(93074326317.L*gamma(5.L, (176.L*Pi*y) / 
(3.L)))/(274420596736.L*pow(Pi, 5)*pow(q, 15)) + 
(31348595781.L*gamma(5.L, (188.L*Pi*y) / 
(3.L)))/(95407522912.L*pow(Pi, 5)*pow(q, 16)) + 
(705020528421.L*gamma(5.L, (200.L*Pi*y) / 
(3.L)))/(2080000000000.L*pow(Pi, 5)*pow(q, 17)) + 
(457296772689.L*gamma(5.L, (212.L*Pi*y) / 
(3.L)))/(1391754600704.L*pow(Pi, 5)*pow(q, 18)) + 
(155422052775.L*gamma(5.L, (224.L*Pi*y) / 
(3.L)))/(458208837632.L*pow(Pi, 5)*pow(q, 19)) + 
(390884861025.L*gamma(5.L, (236.L*Pi*y) / 
(3.L)))/(1189634033536.L*pow(Pi, 5)*pow(q, 20)));
 

Y0433=q32*(216513.L/(26624.L*pow(Pi, 5)) + (10254843.L*q)/(1300000.L*pow(Pi, 
5)) + (221925825.L*pow(q, 2))/(27262976.L*pow(Pi, 5)) + 
(264165543.L*pow(q, 3))/(33498608.L*pow(Pi, 5)) + (454893813.L*pow(q, 
4))/(55933696.L*pow(Pi, 5)) + (4657844169.L*pow(q, 
5))/(590660512.L*pow(Pi, 5)) + (10839369051.L*pow(q, 
6))/(1331200000.L*pow(Pi, 5)) + (5278606623.L*pow(q, 
7))/(669379672.L*pow(Pi, 5)) + (40194988911.L*pow(q, 
8))/(4942652416.L*pow(Pi, 5)) + (67286827575.L*pow(q, 
9))/(8532637984.L*pow(Pi, 5)) + (227252261313.L*pow(q, 
10))/(27917287424.L*pow(Pi, 5)) + (21545425143.L*pow(q, 
11))/(2731137500.L*pow(Pi, 5)) + (134026959825.L*pow(q, 
12))/(16480914944.L*pow(Pi, 5)) + (380066270661.L*pow(q, 
13))/(48196179616.L*pow(Pi, 5)) + (279222978951.L*pow(q, 
14))/(34302574592.L*pow(Pi, 5)) + (94045787343.L*pow(q, 
15))/(11925940364.L*pow(Pi, 5)) + (2115061585263.L*pow(q, 
16))/(260000000000.L*pow(Pi, 5)) + (1371890318067.L*pow(q, 
17))/(173969325088.L*pow(Pi, 5)) + (466266158325.L*pow(q, 
18))/(57276104704.L*pow(Pi, 5)) + (1172654583075.L*pow(q, 
19))/(148704254192.L*pow(Pi, 5)) + (2187.L*gamma(5.L, (4.L*Pi*y) / 
(3.L)))/(6656.L*pow(Pi, 5)*q) + (2311659.L*gamma(5.L, (16.L*Pi*y) / 
(3.L)))/(6815744.L*pow(Pi, 5)*pow(q, 2)) + (4594887.L*gamma(5.L, 
(28.L*Pi*y) / (3.L)))/(13983424.L*pow(Pi, 5)*pow(q, 3)) + 
(112803273.L*gamma(5.L, (40.L*Pi*y) / (3.L)))/(332800000.L*pow(Pi, 
5)*pow(q, 4)) + (406009989.L*gamma(5.L, (52.L*Pi*y) / 
(3.L)))/(1235663104.L*pow(Pi, 5)*pow(q, 5)) + 
(2367210987.L*gamma(5.L, (64.L*Pi*y) / (3.L)))/(6979321856.L*pow(Pi, 
5)*pow(q, 6)) + (1353807675.L*gamma(5.L, (76.L*Pi*y) / 
(3.L)))/(4120228736.L*pow(Pi, 5)*pow(q, 7)) + (264165543.L*gamma(5.L, 
(88.L*Pi*y) / (3.L)))/(779603968.L*pow(Pi, 5)*pow(q, 8)) + 
(21364258437.L*gamma(5.L, (100.L*Pi*y) / 
(3.L)))/(65000000000.L*pow(Pi, 5)*pow(q, 9)) + 
(693827937.L*gamma(5.L, (112.L*Pi*y) / (3.L)))/(2045575168.L*pow(Pi, 
5)*pow(q, 10)) + (1956623607.L*gamma(5.L, (124.L*Pi*y) / 
(3.L)))/(5954863408.L*pow(Pi, 5)*pow(q, 11)) + 
(51236285859.L*gamma(5.L, (136.L*Pi*y) / 
(3.L)))/(151209091072.L*pow(Pi, 5)*pow(q, 12)) + 
(75827618073.L*gamma(5.L, (148.L*Pi*y) / 
(3.L)))/(230776688896.L*pow(Pi, 5)*pow(q, 13)) + 
(4624934193.L*gamma(5.L, (160.L*Pi*y) / 
(3.L)))/(13631488000.L*pow(Pi, 5)*pow(q, 14)) + 
(80376866757.L*gamma(5.L, (172.L*Pi*y) / 
(3.L)))/(244622049152.L*pow(Pi, 5)*pow(q, 15)) + 
(58064672853.L*gamma(5.L, (184.L*Pi*y) / 
(3.L)))/(171361196032.L*pow(Pi, 5)*pow(q, 16)) + 
(617810128659.L*gamma(5.L, (196.L*Pi*y) / 
(3.L)))/(1880155257344.L*pow(Pi, 5)*pow(q, 17)) + 
(429152558373.L*gamma(5.L, (208.L*Pi*y) / 
(3.L)))/(1265319018496.L*pow(Pi, 5)*pow(q, 18)) + 
(137630247903.L*gamma(5.L, (220.L*Pi*y) / 
(3.L)))/(418732600000.L*pow(Pi, 5)*pow(q, 19)) + 
(740155103325.L*gamma(5.L, (232.L*Pi*y) / 
(3.L)))/(2184355323904.L*pow(Pi, 5)*pow(q, 20)));







//weight:w=-3



Y03211=q31*(-729.L/(64.L*sqrt(2.L)*pow(Pi, 4)) - 
(175689.L*q)/(16384.L*sqrt(2.L)*pow(Pi, 4)) - (875529.L*pow(q, 
2))/(76832.L*sqrt(2.L)*pow(Pi, 4)) - (85293.L*pow(q, 
3))/(8000.L*sqrt(2.L)*pow(Pi, 4)) - (10410849.L*pow(q, 
4))/(913952.L*sqrt(2.L)*pow(Pi, 4)) - (44965449.L*pow(q, 
5))/(4194304.L*sqrt(2.L)*pow(Pi, 4)) - (47502369.L*pow(q, 
6))/(4170272.L*sqrt(2.L)*pow(Pi, 4)) - (10005525.L*pow(q, 
7))/(937024.L*sqrt(2.L)*pow(Pi, 4)) - (284310729.L*pow(q, 
8))/(25000000.L*sqrt(2.L)*pow(Pi, 4)) - (211002489.L*pow(q, 
9))/(19668992.L*sqrt(2.L)*pow(Pi, 4)) - (336623769.L*pow(q, 
10))/(29552672.L*sqrt(2.L)*pow(Pi, 4)) - (14270175.L*pow(q, 
11))/(1336336.L*sqrt(2.L)*pow(Pi, 4)) - (683132049.L*pow(q, 
12))/(59973152.L*sqrt(2.L)*pow(Pi, 4)) - (21920301.L*pow(q, 
13))/(2048000.L*sqrt(2.L)*pow(Pi, 4)) - (1246153329.L*pow(q, 
14))/(109401632.L*sqrt(2.L)*pow(Pi, 4)) - (95626575.L*pow(q, 
15))/(8954912.L*sqrt(2.L)*pow(Pi, 4)) - (4204290987.L*pow(q, 
16))/(368947264.L*sqrt(2.L)*pow(Pi, 4)) - (2509014609.L*pow(q, 
17))/(233971712.L*sqrt(2.L)*pow(Pi, 4)) - 
(10405746.L*sqrt(2.L)*pow(q, 18))/(1830125.L*pow(Pi, 4)) - 
(483381675.L*pow(q, 19))/(45265984.L*sqrt(2.L)*pow(Pi, 4)) - 
(3645.L*gamma(4.L, (8.L*Pi*y) / (3.L)))/(2048.L*sqrt(2.L)*pow(Pi, 
4)*q) - (9477.L*gamma(4.L, (20.L*Pi*y) / 
(3.L)))/(5000.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 2)) - 
(936765.L*gamma(4.L, (32.L*Pi*y) / 
(3.L)))/(524288.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 3)) - 
(222345.L*gamma(4.L, (44.L*Pi*y) / 
(3.L)))/(117128.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 4)) - 
(4377645.L*gamma(4.L, (56.L*Pi*y) / 
(3.L)))/(2458624.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 5)) - 
(317115.L*gamma(4.L, (68.L*Pi*y) / 
(3.L)))/(167042.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 6)) - 
(2283957.L*gamma(4.L, (80.L*Pi*y) / 
(3.L)))/(1280000.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 7)) - 
(2125035.L*gamma(4.L, (92.L*Pi*y) / 
(3.L)))/(1119364.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 8)) - 
(52054245.L*gamma(4.L, (104.L*Pi*y) / 
(3.L)))/(29246464.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 9)) - 
(10741815.L*gamma(4.L, (116.L*Pi*y) / 
(3.L)))/(5658248.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 10)) - 
(239815485.L*gamma(4.L, (128.L*Pi*y) / 
(3.L)))/(134217728.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 11)) - 
(11381877.L*gamma(4.L, (140.L*Pi*y) / 
(3.L)))/(6002500.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 12)) - 
(237511845.L*gamma(4.L, (152.L*Pi*y) / 
(3.L)))/(133448704.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 13)) - 
(21458115.L*gamma(4.L, (164.L*Pi*y) / 
(3.L)))/(11303044.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 14)) - 
(53585145.L*gamma(4.L, (176.L*Pi*y) / 
(3.L)))/(29984768.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 15)) - 
(18527535.L*gamma(4.L, (188.L*Pi*y) / 
(3.L)))/(9759362.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 16)) - 
(284310729.L*gamma(4.L, (200.L*Pi*y) / 
(3.L)))/(160000000.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 17)) - 
(119836665.L*gamma(4.L, (212.L*Pi*y) / 
(3.L)))/(63123848.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 18)) - 
(1125054765.L*gamma(4.L, (224.L*Pi*y) / 
(3.L)))/(629407744.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 19)) - 
(184032405.L*gamma(4.L, (236.L*Pi*y) / 
(3.L)))/(96938888.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 20)));

Y03212=(369.L*q)/(64.L*pow(Pi, 4)) + (675.L*pow(q, 2))/(128.L*pow(Pi, 4)) + 
(3281.L*pow(q, 3))/(576.L*pow(Pi, 4)) + (88929.L*pow(q, 
4))/(16384.L*pow(Pi, 4)) + (702.L*pow(q, 5))/(125.L*pow(Pi, 4)) + 
(1025.L*pow(q, 6))/(192.L*pow(Pi, 4)) + (443169.L*pow(q, 
7))/(76832.L*pow(Pi, 4)) + (173475.L*pow(q, 8))/(32768.L*pow(Pi, 4)) 
+ (265721.L*pow(q, 9))/(46656.L*pow(Pi, 4)) + (43173.L*pow(q, 
10))/(8000.L*pow(Pi, 4)) + (82350.L*pow(q, 11))/(14641.L*pow(Pi, 4)) 
+ (790721.L*pow(q, 12))/(147456.L*pow(Pi, 4)) + (5269689.L*pow(q, 
13))/(913952.L*pow(Pi, 4)) + (810675.L*pow(q, 14))/(153664.L*pow(Pi, 
4)) + (2132.L*pow(q, 15))/(375.L*pow(Pi, 4)) + (22760289.L*pow(q, 
16))/(4194304.L*pow(Pi, 4)) + (469800.L*pow(q, 17))/(83521.L*pow(Pi, 
4)) + (166075.L*pow(q, 18))/(31104.L*pow(Pi, 4)) + (24044409.L*pow(q, 
19))/(4170272.L*pow(Pi, 4)) + (84591.L*pow(q, 20))/(16000.L*pow(Pi, 
4)) + pow(y, 4) / (4.L) + (729.L*0.9400256809)/(128.L*pow(Pi, 4)) + (15.L*gamma(4.L, 4.L*Pi*y))/(16.L*pow(Pi, 
4)*q) + (1845.L*gamma(4.L, 8.L*Pi*y))/(2048.L*pow(Pi, 4)*pow(q, 2)) + 
(205.L*gamma(4.L, 12.L*Pi*y))/(216.L*pow(Pi, 4)*pow(q, 3)) + 
(3615.L*gamma(4.L, 16.L*Pi*y))/(4096.L*pow(Pi, 4)*pow(q, 4)) + 
(4797.L*gamma(4.L, 20.L*Pi*y))/(5000.L*pow(Pi, 4)*pow(q, 5)) + 
(16405.L*gamma(4.L, 24.L*Pi*y))/(18432.L*pow(Pi, 4)*pow(q, 6)) + 
(18015.L*gamma(4.L, 28.L*Pi*y))/(19208.L*pow(Pi, 4)*pow(q, 7)) + 
(474165.L*gamma(4.L, 32.L*Pi*y))/(524288.L*pow(Pi, 4)*pow(q, 8)) + 
(33215.L*gamma(4.L, 36.L*Pi*y))/(34992.L*pow(Pi, 4)*pow(q, 9)) + 
(351.L*gamma(4.L, 40.L*Pi*y))/(400.L*pow(Pi, 4)*pow(q, 10)) + 
(112545.L*gamma(4.L, 44.L*Pi*y))/(117128.L*pow(Pi, 4)*pow(q, 11)) + 
(49405.L*gamma(4.L, 48.L*Pi*y))/(55296.L*pow(Pi, 4)*pow(q, 12)) + 
(214215.L*gamma(4.L, 52.L*Pi*y))/(228488.L*pow(Pi, 4)*pow(q, 13)) + 
(2215845.L*gamma(4.L, 56.L*Pi*y))/(2458624.L*pow(Pi, 4)*pow(q, 14)) + 
(42653.L*gamma(4.L, 60.L*Pi*y))/(45000.L*pow(Pi, 4)*pow(q, 15)) + 
(925215.L*gamma(4.L, 64.L*Pi*y))/(1048576.L*pow(Pi, 4)*pow(q, 16)) + 
(160515.L*gamma(4.L, 68.L*Pi*y))/(167042.L*pow(Pi, 4)*pow(q, 17)) + 
(1328605.L*gamma(4.L, 72.L*Pi*y))/(1492992.L*pow(Pi, 4)*pow(q, 18)) + 
(977415.L*gamma(4.L, 76.L*Pi*y))/(1042568.L*pow(Pi, 4)*pow(q, 19)) + 
(1156077.L*gamma(4.L, 80.L*Pi*y))/(1280000.L*pow(Pi, 4)*pow(q, 20));

Y03231=(-45.L*q)/(8.L*pow(Pi, 4)) - (5535.L*pow(q, 2))/(1024.L*pow(Pi, 4)) - 
(205.L*pow(q, 3))/(36.L*pow(Pi, 4)) - (10845.L*pow(q, 
4))/(2048.L*pow(Pi, 4)) - (14391.L*pow(q, 5))/(2500.L*pow(Pi, 4)) - 
(16405.L*pow(q, 6))/(3072.L*pow(Pi, 4)) - (54045.L*pow(q, 
7))/(9604.L*pow(Pi, 4)) - (1422495.L*pow(q, 8))/(262144.L*pow(Pi, 4)) 
- (33215.L*pow(q, 9))/(5832.L*pow(Pi, 4)) - (1053.L*pow(q, 
10))/(200.L*pow(Pi, 4)) - (337635.L*pow(q, 11))/(58564.L*pow(Pi, 4)) 
- (49405.L*pow(q, 12))/(9216.L*pow(Pi, 4)) - (642645.L*pow(q, 
13))/(114244.L*pow(Pi, 4)) - (6647535.L*pow(q, 
14))/(1229312.L*pow(Pi, 4)) - (42653.L*pow(q, 15))/(7500.L*pow(Pi, 
4)) - (2775645.L*pow(q, 16))/(524288.L*pow(Pi, 4)) - (481545.L*pow(q, 
17))/(83521.L*pow(Pi, 4)) - (1328605.L*pow(q, 18))/(248832.L*pow(Pi, 
4)) - (2932245.L*pow(q, 19))/(521284.L*pow(Pi, 4)) - 
(3468231.L*pow(q, 20))/(640000.L*pow(Pi, 4)) + pow(y, 4) / (4.L) - 
(729.L*0.9400256809)/(128.L*pow(Pi, 4)) - 
(123.L*gamma(4.L, 4.L*Pi*y))/(128.L*pow(Pi, 4)*q) - (225.L*gamma(4.L, 
8.L*Pi*y))/(256.L*pow(Pi, 4)*pow(q, 2)) - (3281.L*gamma(4.L, 
12.L*Pi*y))/(3456.L*pow(Pi, 4)*pow(q, 3)) - (29643.L*gamma(4.L, 
16.L*Pi*y))/(32768.L*pow(Pi, 4)*pow(q, 4)) - (117.L*gamma(4.L, 
20.L*Pi*y))/(125.L*pow(Pi, 4)*pow(q, 5)) - (1025.L*gamma(4.L, 
24.L*Pi*y))/(1152.L*pow(Pi, 4)*pow(q, 6)) - (147723.L*gamma(4.L, 
28.L*Pi*y))/(153664.L*pow(Pi, 4)*pow(q, 7)) - (57825.L*gamma(4.L, 
32.L*Pi*y))/(65536.L*pow(Pi, 4)*pow(q, 8)) - (265721.L*gamma(4.L, 
36.L*Pi*y))/(279936.L*pow(Pi, 4)*pow(q, 9)) - (14391.L*gamma(4.L, 
40.L*Pi*y))/(16000.L*pow(Pi, 4)*pow(q, 10)) - (13725.L*gamma(4.L, 
44.L*Pi*y))/(14641.L*pow(Pi, 4)*pow(q, 11)) - (790721.L*gamma(4.L, 
48.L*Pi*y))/(884736.L*pow(Pi, 4)*pow(q, 12)) - (1756563.L*gamma(4.L, 
52.L*Pi*y))/(1827904.L*pow(Pi, 4)*pow(q, 13)) - (270225.L*gamma(4.L, 
56.L*Pi*y))/(307328.L*pow(Pi, 4)*pow(q, 14)) - (1066.L*gamma(4.L, 
60.L*Pi*y))/(1125.L*pow(Pi, 4)*pow(q, 15)) - (7586763.L*gamma(4.L, 
64.L*Pi*y))/(8388608.L*pow(Pi, 4)*pow(q, 16)) - (78300.L*gamma(4.L, 
68.L*Pi*y))/(83521.L*pow(Pi, 4)*pow(q, 17)) - (166075.L*gamma(4.L, 
72.L*Pi*y))/(186624.L*pow(Pi, 4)*pow(q, 18)) - (8014803.L*gamma(4.L, 
76.L*Pi*y))/(8340544.L*pow(Pi, 4)*pow(q, 19)) - (28197.L*gamma(4.L, 
80.L*Pi*y))/(32000.L*pow(Pi, 4)*pow(q, 20));

Y03232=q32*(-10935.L/(1024.L*sqrt(2.L)*pow(Pi, 4)) - 
(28431.L*q)/(2500.L*sqrt(2.L)*pow(Pi, 4)) - (2810295.L*pow(q, 
2))/(262144.L*sqrt(2.L)*pow(Pi, 4)) - (667035.L*pow(q, 
3))/(58564.L*sqrt(2.L)*pow(Pi, 4)) - (13132935.L*pow(q, 
4))/(1229312.L*sqrt(2.L)*pow(Pi, 4)) - (951345.L*pow(q, 
5))/(83521.L*sqrt(2.L)*pow(Pi, 4)) - (6851871.L*pow(q, 
6))/(640000.L*sqrt(2.L)*pow(Pi, 4)) - (6375105.L*pow(q, 
7))/(559682.L*sqrt(2.L)*pow(Pi, 4)) - (156162735.L*pow(q, 
8))/(14623232.L*sqrt(2.L)*pow(Pi, 4)) - (32225445.L*pow(q, 
9))/(2829124.L*sqrt(2.L)*pow(Pi, 4)) - (719446455.L*pow(q, 
10))/(67108864.L*sqrt(2.L)*pow(Pi, 4)) - (34145631.L*pow(q, 
11))/(3001250.L*sqrt(2.L)*pow(Pi, 4)) - (712535535.L*pow(q, 
12))/(66724352.L*sqrt(2.L)*pow(Pi, 4)) - (64374345.L*pow(q, 
13))/(5651522.L*sqrt(2.L)*pow(Pi, 4)) - (160755435.L*pow(q, 
14))/(14992384.L*sqrt(2.L)*pow(Pi, 4)) - (55582605.L*pow(q, 
15))/(4879681.L*sqrt(2.L)*pow(Pi, 4)) - (852932187.L*pow(q, 
16))/(80000000.L*sqrt(2.L)*pow(Pi, 4)) - (359509995.L*pow(q, 
17))/(31561924.L*sqrt(2.L)*pow(Pi, 4)) - (3375164295.L*pow(q, 
18))/(314703872.L*sqrt(2.L)*pow(Pi, 4)) - (552097215.L*pow(q, 
19))/(48469444.L*sqrt(2.L)*pow(Pi, 4)) - (243.L*gamma(4.L, (4.L*Pi*y) 
/ (3.L)))/(128.L*sqrt(2.L)*pow(Pi, 4)*q) - (58563.L*gamma(4.L, 
(16.L*Pi*y) / (3.L)))/(32768.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 2)) - 
(291843.L*gamma(4.L, (28.L*Pi*y) / 
(3.L)))/(153664.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 3)) - 
(28431.L*gamma(4.L, (40.L*Pi*y) / (3.L)))/(16000.L*sqrt(2.L)*pow(Pi, 
4)*pow(q, 4)) - (3470283.L*gamma(4.L, (52.L*Pi*y) / 
(3.L)))/(1827904.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 5)) - 
(14988483.L*gamma(4.L, (64.L*Pi*y) / 
(3.L)))/(8388608.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 6)) - 
(15834123.L*gamma(4.L, (76.L*Pi*y) / 
(3.L)))/(8340544.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 7)) - 
(3335175.L*gamma(4.L, (88.L*Pi*y) / 
(3.L)))/(1874048.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 8)) - 
(94770243.L*gamma(4.L, (100.L*Pi*y) / 
(3.L)))/(50000000.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 9)) - 
(70334163.L*gamma(4.L, (112.L*Pi*y) / 
(3.L)))/(39337984.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 10)) - 
(112207923.L*gamma(4.L, (124.L*Pi*y) / 
(3.L)))/(59105344.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 11)) - 
(4756725.L*gamma(4.L, (136.L*Pi*y) / 
(3.L)))/(2672672.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 12)) - 
(227710683.L*gamma(4.L, (148.L*Pi*y) / 
(3.L)))/(119946304.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 13)) - 
(7306767.L*gamma(4.L, (160.L*Pi*y) / 
(3.L)))/(4096000.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 14)) - 
(415384443.L*gamma(4.L, (172.L*Pi*y) / 
(3.L)))/(218803264.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 15)) - 
(31875525.L*gamma(4.L, (184.L*Pi*y) / 
(3.L)))/(17909824.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 16)) - 
(1401430329.L*gamma(4.L, (196.L*Pi*y) / 
(3.L)))/(737894528.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 17)) - 
(836338203.L*gamma(4.L, (208.L*Pi*y) / 
(3.L)))/(467943424.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 18)) - 
(1734291.L*sqrt(2.L)*gamma(4.L, (220.L*Pi*y) / 
(3.L)))/(1830125.L*pow(Pi, 4)*pow(q, 19)) - (161127225.L*gamma(4.L, 
(232.L*Pi*y) / (3.L)))/(90531968.L*sqrt(2.L)*pow(Pi, 4)*pow(q, 20)));










// Weight: w=-2


Y0211=(-15.L*q)/(2.L*pow(Pi, 3)) - (135.L*pow(q, 2))/(16.L*pow(Pi, 3)) - 
(70.L*pow(q, 3))/(9.L*pow(Pi, 3)) - (1095.L*pow(q, 4))/(128.L*pow(Pi, 
3)) - (189.L*pow(q, 5))/(25.L*pow(Pi, 3)) - (35.L*pow(q, 
6))/(4.L*pow(Pi, 3)) - (2580.L*pow(q, 7))/(343.L*pow(Pi, 3)) - 
(8775.L*pow(q, 8))/(1024.L*pow(Pi, 3)) - (3785.L*pow(q, 
9))/(486.L*pow(Pi, 3)) - (1701.L*pow(q, 10))/(200.L*pow(Pi, 3)) - 
(9990.L*pow(q, 11))/(1331.L*pow(Pi, 3)) - (2555.L*pow(q, 
12))/(288.L*pow(Pi, 3)) - (16485.L*pow(q, 13))/(2197.L*pow(Pi, 3)) - 
(5805.L*pow(q, 14))/(686.L*pow(Pi, 3)) - (196.L*pow(q, 
15))/(25.L*pow(Pi, 3)) - (70215.L*pow(q, 16))/(8192.L*pow(Pi, 3)) - 
(36855.L*pow(q, 17))/(4913.L*pow(Pi, 3)) - (3785.L*pow(q, 
18))/(432.L*pow(Pi, 3)) - (51450.L*pow(q, 19))/(6859.L*pow(Pi, 3)) - 
(13797.L*pow(q, 20))/(1600.L*pow(Pi, 3)) + pow(y, 3) / (3.L) - 
(15.L*gamma(3.L, 4.L*Pi*y))/(4.L*pow(Pi, 3)*q) - (135.L*gamma(3.L, 
8.L*Pi*y))/(32.L*pow(Pi, 3)*pow(q, 2)) - (35.L*gamma(3.L, 
12.L*Pi*y))/(9.L*pow(Pi, 3)*pow(q, 3)) - (1095.L*gamma(3.L, 
16.L*Pi*y))/(256.L*pow(Pi, 3)*pow(q, 4)) - (189.L*gamma(3.L, 
20.L*Pi*y))/(50.L*pow(Pi, 3)*pow(q, 5)) - (35.L*gamma(3.L, 
24.L*Pi*y))/(8.L*pow(Pi, 3)*pow(q, 6)) - (1290.L*gamma(3.L, 
28.L*Pi*y))/(343.L*pow(Pi, 3)*pow(q, 7)) - (8775.L*gamma(3.L, 
32.L*Pi*y))/(2048.L*pow(Pi, 3)*pow(q, 8)) - (3785.L*gamma(3.L, 
36.L*Pi*y))/(972.L*pow(Pi, 3)*pow(q, 9)) - (1701.L*gamma(3.L, 
40.L*Pi*y))/(400.L*pow(Pi, 3)*pow(q, 10)) - (4995.L*gamma(3.L, 
44.L*Pi*y))/(1331.L*pow(Pi, 3)*pow(q, 11)) - (2555.L*gamma(3.L, 
48.L*Pi*y))/(576.L*pow(Pi, 3)*pow(q, 12)) - (16485.L*gamma(3.L, 
52.L*Pi*y))/(4394.L*pow(Pi, 3)*pow(q, 13)) - (5805.L*gamma(3.L, 
56.L*Pi*y))/(1372.L*pow(Pi, 3)*pow(q, 14)) - (98.L*gamma(3.L, 
60.L*Pi*y))/(25.L*pow(Pi, 3)*pow(q, 15)) - (70215.L*gamma(3.L, 
64.L*Pi*y))/(16384.L*pow(Pi, 3)*pow(q, 16)) - (36855.L*gamma(3.L, 
68.L*Pi*y))/(9826.L*pow(Pi, 3)*pow(q, 17)) - (3785.L*gamma(3.L, 
72.L*Pi*y))/(864.L*pow(Pi, 3)*pow(q, 18)) - (25725.L*gamma(3.L, 
76.L*Pi*y))/(6859.L*pow(Pi, 3)*pow(q, 19)) - (13797.L*gamma(3.L, 
80.L*Pi*y))/(3200.L*pow(Pi, 3)*pow(q, 20)) - 
(15.L*zeta3)/(2.L*pow(Pi, 3));

Y0231=(21.L*q)/(8.L*pow(Pi, 3)) + (189.L*pow(q, 2))/(64.L*pow(Pi, 3)) + 
(169.L*pow(q, 3))/(72.L*pow(Pi, 3)) + (1533.L*pow(q, 
4))/(512.L*pow(Pi, 3)) + (1323.L*pow(q, 5))/(500.L*pow(Pi, 3)) + 
(169.L*pow(q, 6))/(64.L*pow(Pi, 3)) + (129.L*pow(q, 7))/(49.L*pow(Pi, 
3)) + (12285.L*pow(q, 8))/(4096.L*pow(Pi, 3)) + (4543.L*pow(q, 
9))/(1944.L*pow(Pi, 3)) + (11907.L*pow(q, 10))/(4000.L*pow(Pi, 3)) + 
(6993.L*pow(q, 11))/(2662.L*pow(Pi, 3)) + (12337.L*pow(q, 
12))/(4608.L*pow(Pi, 3)) + (23079.L*pow(q, 13))/(8788.L*pow(Pi, 3)) + 
(1161.L*pow(q, 14))/(392.L*pow(Pi, 3)) + (1183.L*pow(q, 
15))/(500.L*pow(Pi, 3)) + (98301.L*pow(q, 16))/(32768.L*pow(Pi, 3)) + 
(51597.L*pow(q, 17))/(19652.L*pow(Pi, 3)) + (4543.L*pow(q, 
18))/(1728.L*pow(Pi, 3)) + (36015.L*pow(q, 19))/(13718.L*pow(Pi, 3)) 
+ (96579.L*pow(q, 20))/(32000.L*pow(Pi, 3)) + pow(y, 3) / (3.L) + 
(21.L*gamma(3.L, 4.L*Pi*y))/(16.L*pow(Pi, 3)*q) + (189.L*gamma(3.L, 
8.L*Pi*y))/(128.L*pow(Pi, 3)*pow(q, 2)) + (169.L*gamma(3.L, 
12.L*Pi*y))/(144.L*pow(Pi, 3)*pow(q, 3)) + (1533.L*gamma(3.L, 
16.L*Pi*y))/(1024.L*pow(Pi, 3)*pow(q, 4)) + (1323.L*gamma(3.L, 
20.L*Pi*y))/(1000.L*pow(Pi, 3)*pow(q, 5)) + (169.L*gamma(3.L, 
24.L*Pi*y))/(128.L*pow(Pi, 3)*pow(q, 6)) + (129.L*gamma(3.L, 
28.L*Pi*y))/(98.L*pow(Pi, 3)*pow(q, 7)) + (12285.L*gamma(3.L, 
32.L*Pi*y))/(8192.L*pow(Pi, 3)*pow(q, 8)) + (4543.L*gamma(3.L, 
36.L*Pi*y))/(3888.L*pow(Pi, 3)*pow(q, 9)) + (11907.L*gamma(3.L, 
40.L*Pi*y))/(8000.L*pow(Pi, 3)*pow(q, 10)) + (6993.L*gamma(3.L, 
44.L*Pi*y))/(5324.L*pow(Pi, 3)*pow(q, 11)) + (12337.L*gamma(3.L, 
48.L*Pi*y))/(9216.L*pow(Pi, 3)*pow(q, 12)) + (23079.L*gamma(3.L, 
52.L*Pi*y))/(17576.L*pow(Pi, 3)*pow(q, 13)) + (1161.L*gamma(3.L, 
56.L*Pi*y))/(784.L*pow(Pi, 3)*pow(q, 14)) + (1183.L*gamma(3.L, 
60.L*Pi*y))/(1000.L*pow(Pi, 3)*pow(q, 15)) + (98301.L*gamma(3.L, 
64.L*Pi*y))/(65536.L*pow(Pi, 3)*pow(q, 16)) + (51597.L*gamma(3.L, 
68.L*Pi*y))/(39304.L*pow(Pi, 3)*pow(q, 17)) + (4543.L*gamma(3.L, 
72.L*Pi*y))/(3456.L*pow(Pi, 3)*pow(q, 18)) + (36015.L*gamma(3.L, 
76.L*Pi*y))/(27436.L*pow(Pi, 3)*pow(q, 19)) + (96579.L*gamma(3.L, 
80.L*Pi*y))/(64000.L*pow(Pi, 3)*pow(q, 20)) + 
(9.L*zeta3)/(4.L*pow(Pi, 3));

Y0232=q31*(-81.L/(16.L*pow(Pi, 3)) - (5913.L*q)/(1024.L*pow(Pi, 3)) - 
(3483.L*pow(q, 2))/(686.L*pow(Pi, 3)) - (45927.L*pow(q, 
3))/(8000.L*pow(Pi, 3)) - (89019.L*pow(q, 4))/(17576.L*pow(Pi, 3)) - 
(379161.L*pow(q, 5))/(65536.L*pow(Pi, 3)) - (138915.L*pow(q, 
6))/(27436.L*pow(Pi, 3)) - (242757.L*pow(q, 7))/(42592.L*pow(Pi, 3)) 
- (1275831.L*pow(q, 8))/(250000.L*pow(Pi, 3)) - (254259.L*pow(q, 
9))/(43904.L*pow(Pi, 3)) - (150822.L*pow(q, 10))/(29791.L*pow(Pi, 3)) 
- (1791153.L*pow(q, 11))/(314432.L*pow(Pi, 3)) - (2051487.L*pow(q, 
12))/(405224.L*pow(Pi, 3)) - (597051.L*pow(q, 13))/(102400.L*pow(Pi, 
3)) - (1610037.L*pow(q, 14))/(318028.L*pow(Pi, 3)) - 
(1108809.L*pow(q, 15))/(194672.L*pow(Pi, 3)) - (9557433.L*pow(q, 
16))/(1882384.L*pow(Pi, 3)) - (6498387.L*pow(q, 
17))/(1124864.L*pow(Pi, 3)) - (1699299.L*pow(q, 
18))/(332750.L*pow(Pi, 3)) - (8890155.L*pow(q, 
19))/(1560896.L*pow(Pi, 3)) - (729.L*gamma(3.L, (8.L*Pi*y) / 
(3.L)))/(256.L*pow(Pi, 3)*q) - (5103.L*gamma(3.L, (20.L*Pi*y) / 
(3.L)))/(2000.L*pow(Pi, 3)*pow(q, 2)) - (47385.L*gamma(3.L, 
(32.L*Pi*y) / (3.L)))/(16384.L*pow(Pi, 3)*pow(q, 3)) - 
(26973.L*gamma(3.L, (44.L*Pi*y) / (3.L)))/(10648.L*pow(Pi, 3)*pow(q, 
4)) - (31347.L*gamma(3.L, (56.L*Pi*y) / (3.L)))/(10976.L*pow(Pi, 
3)*pow(q, 5)) - (199017.L*gamma(3.L, (68.L*Pi*y) / 
(3.L)))/(78608.L*pow(Pi, 3)*pow(q, 6)) - (372519.L*gamma(3.L, 
(80.L*Pi*y) / (3.L)))/(128000.L*pow(Pi, 3)*pow(q, 7)) - 
(123201.L*gamma(3.L, (92.L*Pi*y) / (3.L)))/(48668.L*pow(Pi, 3)*pow(q, 
8)) - (801171.L*gamma(3.L, (104.L*Pi*y) / (3.L)))/(281216.L*pow(Pi, 
3)*pow(q, 9)) - (987795.L*gamma(3.L, (116.L*Pi*y) / 
(3.L)))/(390224.L*pow(Pi, 3)*pow(q, 10)) - (3033369.L*gamma(3.L, 
(128.L*Pi*y) / (3.L)))/(1048576.L*pow(Pi, 3)*pow(q, 11)) - 
(31347.L*gamma(3.L, (140.L*Pi*y) / (3.L)))/(12250.L*pow(Pi, 3)*pow(q, 
12)) - (1250235.L*gamma(3.L, (152.L*Pi*y) / (3.L)))/(438976.L*pow(Pi, 
3)*pow(q, 13)) - (2791341.L*gamma(3.L, (164.L*Pi*y) / 
(3.L)))/(1102736.L*pow(Pi, 3)*pow(q, 14)) - (1969029.L*gamma(3.L, 
(176.L*Pi*y) / (3.L)))/(681472.L*pow(Pi, 3)*pow(q, 15)) - 
(525609.L*gamma(3.L, (188.L*Pi*y) / (3.L)))/(207646.L*pow(Pi, 
3)*pow(q, 16)) - (11482479.L*gamma(3.L, (200.L*Pi*y) / 
(3.L)))/(4000000.L*pow(Pi, 3)*pow(q, 17)) - (6029559.L*gamma(3.L, 
(212.L*Pi*y) / (3.L)))/(2382032.L*pow(Pi, 3)*pow(q, 18)) - 
(2037555.L*gamma(3.L, (224.L*Pi*y) / (3.L)))/(702464.L*pow(Pi, 
3)*pow(q, 19)) - (4158945.L*gamma(3.L, (236.L*Pi*y) / 
(3.L)))/(1643032.L*pow(Pi, 3)*pow(q, 20)));
 

Y0233=q32*(-729.L/(128.L*pow(Pi, 3)) - (5103.L*q)/(1000.L*pow(Pi, 3)) - 
(47385.L*pow(q, 2))/(8192.L*pow(Pi, 3)) - (26973.L*pow(q, 
3))/(5324.L*pow(Pi, 3)) - (31347.L*pow(q, 4))/(5488.L*pow(Pi, 3)) - 
(199017.L*pow(q, 5))/(39304.L*pow(Pi, 3)) - (372519.L*pow(q, 
6))/(64000.L*pow(Pi, 3)) - (123201.L*pow(q, 7))/(24334.L*pow(Pi, 3)) 
- (801171.L*pow(q, 8))/(140608.L*pow(Pi, 3)) - (987795.L*pow(q, 
9))/(195112.L*pow(Pi, 3)) - (3033369.L*pow(q, 10))/(524288.L*pow(Pi, 
3)) - (31347.L*pow(q, 11))/(6125.L*pow(Pi, 3)) - (1250235.L*pow(q, 
12))/(219488.L*pow(Pi, 3)) - (2791341.L*pow(q, 13))/(551368.L*pow(Pi, 
3)) - (1969029.L*pow(q, 14))/(340736.L*pow(Pi, 3)) - (525609.L*pow(q, 
15))/(103823.L*pow(Pi, 3)) - (11482479.L*pow(q, 
16))/(2000000.L*pow(Pi, 3)) - (6029559.L*pow(q, 
17))/(1191016.L*pow(Pi, 3)) - (2037555.L*pow(q, 
18))/(351232.L*pow(Pi, 3)) - (4158945.L*pow(q, 19))/(821516.L*pow(Pi, 
3)) - (81.L*gamma(3.L, (4.L*Pi*y) / (3.L)))/(32.L*pow(Pi, 3)*q) - 
(5913.L*gamma(3.L, (16.L*Pi*y) / (3.L)))/(2048.L*pow(Pi, 3)*pow(q, 
2)) - (3483.L*gamma(3.L, (28.L*Pi*y) / (3.L)))/(1372.L*pow(Pi, 
3)*pow(q, 3)) - (45927.L*gamma(3.L, (40.L*Pi*y) / 
(3.L)))/(16000.L*pow(Pi, 3)*pow(q, 4)) - (89019.L*gamma(3.L, 
(52.L*Pi*y) / (3.L)))/(35152.L*pow(Pi, 3)*pow(q, 5)) - 
(379161.L*gamma(3.L, (64.L*Pi*y) / (3.L)))/(131072.L*pow(Pi, 
3)*pow(q, 6)) - (138915.L*gamma(3.L, (76.L*Pi*y) / 
(3.L)))/(54872.L*pow(Pi, 3)*pow(q, 7)) - (242757.L*gamma(3.L, 
(88.L*Pi*y) / (3.L)))/(85184.L*pow(Pi, 3)*pow(q, 8)) - 
(1275831.L*gamma(3.L, (100.L*Pi*y) / (3.L)))/(500000.L*pow(Pi, 
3)*pow(q, 9)) - (254259.L*gamma(3.L, (112.L*Pi*y) / 
(3.L)))/(87808.L*pow(Pi, 3)*pow(q, 10)) - (75411.L*gamma(3.L, 
(124.L*Pi*y) / (3.L)))/(29791.L*pow(Pi, 3)*pow(q, 11)) - 
(1791153.L*gamma(3.L, (136.L*Pi*y) / (3.L)))/(628864.L*pow(Pi, 
3)*pow(q, 12)) - (2051487.L*gamma(3.L, (148.L*Pi*y) / 
(3.L)))/(810448.L*pow(Pi, 3)*pow(q, 13)) - (597051.L*gamma(3.L, 
(160.L*Pi*y) / (3.L)))/(204800.L*pow(Pi, 3)*pow(q, 14)) - 
(1610037.L*gamma(3.L, (172.L*Pi*y) / (3.L)))/(636056.L*pow(Pi, 
3)*pow(q, 15)) - (1108809.L*gamma(3.L, (184.L*Pi*y) / 
(3.L)))/(389344.L*pow(Pi, 3)*pow(q, 16)) - (9557433.L*gamma(3.L, 
(196.L*Pi*y) / (3.L)))/(3764768.L*pow(Pi, 3)*pow(q, 17)) - 
(6498387.L*gamma(3.L, (208.L*Pi*y) / (3.L)))/(2249728.L*pow(Pi, 
3)*pow(q, 18)) - (1699299.L*gamma(3.L, (220.L*Pi*y) / 
(3.L)))/(665500.L*pow(Pi, 3)*pow(q, 19)) - (8890155.L*gamma(3.L, 
(232.L*Pi*y) / (3.L)))/(3121792.L*pow(Pi, 3)*pow(q, 20)));









//weight=-1


Y01211=q31*(81.L/(8.L*sqrt(2.L)*pow(Pi, 2)) + 
(1053.L*q)/(128.L*sqrt(2.L)*pow(Pi, 2)) + (2025.L*pow(q, 
2))/(196.L*sqrt(2.L)*pow(Pi, 2)) + (729.L*pow(q, 
3))/(100.L*sqrt(2.L)*pow(Pi, 2)) + (6885.L*pow(q, 
4))/(676.L*sqrt(2.L)*pow(Pi, 2)) + (16605.L*pow(q, 
5))/(2048.L*sqrt(2.L)*pow(Pi, 2)) + (14661.L*pow(q, 
6))/(1444.L*sqrt(2.L)*pow(Pi, 2)) + (3645.L*pow(q, 
7))/(484.L*sqrt(2.L)*pow(Pi, 2)) + (48681.L*pow(q, 
8))/(5000.L*sqrt(2.L)*pow(Pi, 2)) + (26325.L*pow(q, 
9))/(3136.L*sqrt(2.L)*pow(Pi, 2)) + (38961.L*pow(q, 
10))/(3844.L*sqrt(2.L)*pow(Pi, 2)) + (2187.L*pow(q, 
11))/(289.L*sqrt(2.L)*pow(Pi, 2)) + (55485.L*pow(q, 
12))/(5476.L*sqrt(2.L)*pow(Pi, 2)) + (12393.L*pow(q, 
13))/(1600.L*sqrt(2.L)*pow(Pi, 2)) + (74925.L*pow(q, 
14))/(7396.L*sqrt(2.L)*pow(Pi, 2)) + (8019.L*pow(q, 
15))/(1058.L*sqrt(2.L)*pow(Pi, 2)) + (198531.L*pow(q, 
16))/(19208.L*sqrt(2.L)*pow(Pi, 2)) + (6885.L*pow(q, 
17))/(832.L*sqrt(2.L)*pow(Pi, 2)) + (2916.L*sqrt(2.L)*pow(q, 
18))/(605.L*pow(Pi, 2)) + (25515.L*pow(q, 
19))/(3364.L*sqrt(2.L)*pow(Pi, 2)) + (243.L*gamma(2.L, (8.L*Pi*y) / 
(3.L)))/(32.L*sqrt(2.L)*pow(Pi, 2)*q) + (243.L*gamma(2.L, (20.L*Pi*y) 
/ (3.L)))/(25.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 2)) + (4131.L*gamma(2.L, 
(32.L*Pi*y) / (3.L)))/(512.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 3)) + 
(1215.L*gamma(2.L, (44.L*Pi*y) / (3.L)))/(121.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 4)) + (6075.L*gamma(2.L, (56.L*Pi*y) / 
(3.L)))/(784.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 5)) + 
(1458.L*sqrt(2.L)*gamma(2.L, (68.L*Pi*y) / (3.L)))/(289.L*pow(Pi, 
2)*pow(q, 6)) + (3159.L*gamma(2.L, (80.L*Pi*y) / 
(3.L)))/(400.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 7)) + 
(2673.L*sqrt(2.L)*gamma(2.L, (92.L*Pi*y) / (3.L)))/(529.L*pow(Pi, 
2)*pow(q, 8)) + (20655.L*gamma(2.L, (104.L*Pi*y) / 
(3.L)))/(2704.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 9)) + (8505.L*gamma(2.L, 
(116.L*Pi*y) / (3.L)))/(841.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 10)) + 
(66339.L*gamma(2.L, (128.L*Pi*y) / (3.L)))/(8192.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 11)) + (243.L*sqrt(2.L)*gamma(2.L, (140.L*Pi*y) / 
(3.L)))/(49.L*pow(Pi, 2)*pow(q, 12)) + (43983.L*gamma(2.L, 
(152.L*Pi*y) / (3.L)))/(5776.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 13)) + 
(8505.L*sqrt(2.L)*gamma(2.L, (164.L*Pi*y) / (3.L)))/(1681.L*pow(Pi, 
2)*pow(q, 14)) + (15795.L*gamma(2.L, (176.L*Pi*y) / 
(3.L)))/(1936.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 15)) + 
(11178.L*sqrt(2.L)*gamma(2.L, (188.L*Pi*y) / (3.L)))/(2209.L*pow(Pi, 
2)*pow(q, 16)) + (146043.L*gamma(2.L, (200.L*Pi*y) / 
(3.L)))/(20000.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 17)) + 
(28431.L*gamma(2.L, (212.L*Pi*y) / (3.L)))/(2809.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 18)) + (103275.L*gamma(2.L, (224.L*Pi*y) / 
(3.L)))/(12544.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 19)) + 
(35235.L*gamma(2.L, (236.L*Pi*y) / (3.L)))/(3481.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 20)));

Y01212=(-45.L*q)/(8.L*pow(Pi, 2)) - (27.L*pow(q, 2))/(8.L*pow(Pi, 2)) - 
(41.L*pow(q, 3))/(8.L*pow(Pi, 2)) - (585.L*pow(q, 4))/(128.L*pow(Pi, 
2)) - (108.L*pow(q, 5))/(25.L*pow(Pi, 2)) - (15.L*pow(q, 
6))/(4.L*pow(Pi, 2)) - (1125.L*pow(q, 7))/(196.L*pow(Pi, 2)) - 
(459.L*pow(q, 8))/(128.L*pow(Pi, 2)) - (365.L*pow(q, 
9))/(72.L*pow(Pi, 2)) - (81.L*pow(q, 10))/(20.L*pow(Pi, 2)) - 
(540.L*pow(q, 11))/(121.L*pow(Pi, 2)) - (533.L*pow(q, 
12))/(128.L*pow(Pi, 2)) - (3825.L*pow(q, 13))/(676.L*pow(Pi, 2)) - 
(675.L*pow(q, 14))/(196.L*pow(Pi, 2)) - (24.L*pow(q, 
15))/(5.L*pow(Pi, 2)) - (9225.L*pow(q, 16))/(2048.L*pow(Pi, 2)) - 
(1296.L*pow(q, 17))/(289.L*pow(Pi, 2)) - (91.L*pow(q, 
18))/(24.L*pow(Pi, 2)) - (8145.L*pow(q, 19))/(1444.L*pow(Pi, 2)) - 
(351.L*pow(q, 20))/(100.L*pow(Pi, 2)) + pow(y, 2) / (2.L) - 
(81.L*0.7813024129)/(16.L*pow(Pi, 2)) - (9.L*gamma(2.L, 
4.L*Pi*y))/(2.L*pow(Pi, 2)*q) - (135.L*gamma(2.L, 
8.L*Pi*y))/(32.L*pow(Pi, 2)*pow(q, 2)) - (5.L*gamma(2.L, 
12.L*Pi*y))/(pow(Pi, 2)*pow(q, 3)) - (117.L*gamma(2.L, 
16.L*Pi*y))/(32.L*pow(Pi, 2)*pow(q, 4)) - (27.L*gamma(2.L, 
20.L*Pi*y))/(5.L*pow(Pi, 2)*pow(q, 5)) - (123.L*gamma(2.L, 
24.L*Pi*y))/(32.L*pow(Pi, 2)*pow(q, 6)) - (225.L*gamma(2.L, 
28.L*Pi*y))/(49.L*pow(Pi, 2)*pow(q, 7)) - (2295.L*gamma(2.L, 
32.L*Pi*y))/(512.L*pow(Pi, 2)*pow(q, 8)) - (91.L*gamma(2.L, 
36.L*Pi*y))/(18.L*pow(Pi, 2)*pow(q, 9)) - (81.L*gamma(2.L, 
40.L*Pi*y))/(25.L*pow(Pi, 2)*pow(q, 10)) - (675.L*gamma(2.L, 
44.L*Pi*y))/(121.L*pow(Pi, 2)*pow(q, 11)) - (65.L*gamma(2.L, 
48.L*Pi*y))/(16.L*pow(Pi, 2)*pow(q, 12)) - (765.L*gamma(2.L, 
52.L*Pi*y))/(169.L*pow(Pi, 2)*pow(q, 13)) - (3375.L*gamma(2.L, 
56.L*Pi*y))/(784.L*pow(Pi, 2)*pow(q, 14)) - (123.L*gamma(2.L, 
60.L*Pi*y))/(25.L*pow(Pi, 2)*pow(q, 15)) - (1845.L*gamma(2.L, 
64.L*Pi*y))/(512.L*pow(Pi, 2)*pow(q, 16)) - (1620.L*gamma(2.L, 
68.L*Pi*y))/(289.L*pow(Pi, 2)*pow(q, 17)) - (365.L*gamma(2.L, 
72.L*Pi*y))/(96.L*pow(Pi, 2)*pow(q, 18)) - (1629.L*gamma(2.L, 
76.L*Pi*y))/(361.L*pow(Pi, 2)*pow(q, 19)) - (351.L*gamma(2.L, 
80.L*Pi*y))/(80.L*pow(Pi, 2)*pow(q, 20));

Y01231=(9.L*q)/(2.L*pow(Pi, 2)) + (135.L*pow(q, 2))/(32.L*pow(Pi, 2)) + 
(5.L*pow(q, 3))/pow(Pi, 2) + (117.L*pow(q, 4))/(32.L*pow(Pi, 2)) + 
(27.L*pow(q, 5))/(5.L*pow(Pi, 2)) + (123.L*pow(q, 6))/(32.L*pow(Pi, 
2)) + (225.L*pow(q, 7))/(49.L*pow(Pi, 2)) + (2295.L*pow(q, 
8))/(512.L*pow(Pi, 2)) + (91.L*pow(q, 9))/(18.L*pow(Pi, 2)) + 
(81.L*pow(q, 10))/(25.L*pow(Pi, 2)) + (675.L*pow(q, 
11))/(121.L*pow(Pi, 2)) + (65.L*pow(q, 12))/(16.L*pow(Pi, 2)) + 
(765.L*pow(q, 13))/(169.L*pow(Pi, 2)) + (3375.L*pow(q, 
14))/(784.L*pow(Pi, 2)) + (123.L*pow(q, 15))/(25.L*pow(Pi, 2)) + 
(1845.L*pow(q, 16))/(512.L*pow(Pi, 2)) + (1620.L*pow(q, 
17))/(289.L*pow(Pi, 2)) + (365.L*pow(q, 18))/(96.L*pow(Pi, 2)) + 
(1629.L*pow(q, 19))/(361.L*pow(Pi, 2)) + (351.L*pow(q, 
20))/(80.L*pow(Pi, 2)) + pow(y, 2) / (2.L) + (81.L*0.7813024129)/(16.L*pow(Pi, 2)) + (45.L*gamma(2.L, 
4.L*Pi*y))/(8.L*pow(Pi, 2)*q) + (27.L*gamma(2.L, 
8.L*Pi*y))/(8.L*pow(Pi, 2)*pow(q, 2)) + (41.L*gamma(2.L, 
12.L*Pi*y))/(8.L*pow(Pi, 2)*pow(q, 3)) + (585.L*gamma(2.L, 
16.L*Pi*y))/(128.L*pow(Pi, 2)*pow(q, 4)) + (108.L*gamma(2.L, 
20.L*Pi*y))/(25.L*pow(Pi, 2)*pow(q, 5)) + (15.L*gamma(2.L, 
24.L*Pi*y))/(4.L*pow(Pi, 2)*pow(q, 6)) + (1125.L*gamma(2.L, 
28.L*Pi*y))/(196.L*pow(Pi, 2)*pow(q, 7)) + (459.L*gamma(2.L, 
32.L*Pi*y))/(128.L*pow(Pi, 2)*pow(q, 8)) + (365.L*gamma(2.L, 
36.L*Pi*y))/(72.L*pow(Pi, 2)*pow(q, 9)) + (81.L*gamma(2.L, 
40.L*Pi*y))/(20.L*pow(Pi, 2)*pow(q, 10)) + (540.L*gamma(2.L, 
44.L*Pi*y))/(121.L*pow(Pi, 2)*pow(q, 11)) + (533.L*gamma(2.L, 
48.L*Pi*y))/(128.L*pow(Pi, 2)*pow(q, 12)) + (3825.L*gamma(2.L, 
52.L*Pi*y))/(676.L*pow(Pi, 2)*pow(q, 13)) + (675.L*gamma(2.L, 
56.L*Pi*y))/(196.L*pow(Pi, 2)*pow(q, 14)) + (24.L*gamma(2.L, 
60.L*Pi*y))/(5.L*pow(Pi, 2)*pow(q, 15)) + (9225.L*gamma(2.L, 
64.L*Pi*y))/(2048.L*pow(Pi, 2)*pow(q, 16)) + (1296.L*gamma(2.L, 
68.L*Pi*y))/(289.L*pow(Pi, 2)*pow(q, 17)) + (91.L*gamma(2.L, 
72.L*Pi*y))/(24.L*pow(Pi, 2)*pow(q, 18)) + (8145.L*gamma(2.L, 
76.L*Pi*y))/(1444.L*pow(Pi, 2)*pow(q, 19)) + (351.L*gamma(2.L, 
80.L*Pi*y))/(100.L*pow(Pi, 2)*pow(q, 20));

Y01232=q32*(243.L/(32.L*sqrt(2.L)*pow(Pi, 2)) + 
(243.L*q)/(25.L*sqrt(2.L)*pow(Pi, 2)) + (4131.L*pow(q, 
2))/(512.L*sqrt(2.L)*pow(Pi, 2)) + (1215.L*pow(q, 
3))/(121.L*sqrt(2.L)*pow(Pi, 2)) + (6075.L*pow(q, 
4))/(784.L*sqrt(2.L)*pow(Pi, 2)) + (1458.L*sqrt(2.L)*pow(q, 
5))/(289.L*pow(Pi, 2)) + (3159.L*pow(q, 6))/(400.L*sqrt(2.L)*pow(Pi, 
2)) + (2673.L*sqrt(2.L)*pow(q, 7))/(529.L*pow(Pi, 2)) + 
(20655.L*pow(q, 8))/(2704.L*sqrt(2.L)*pow(Pi, 2)) + (8505.L*pow(q, 
9))/(841.L*sqrt(2.L)*pow(Pi, 2)) + (66339.L*pow(q, 
10))/(8192.L*sqrt(2.L)*pow(Pi, 2)) + (243.L*sqrt(2.L)*pow(q, 
11))/(49.L*pow(Pi, 2)) + (43983.L*pow(q, 
12))/(5776.L*sqrt(2.L)*pow(Pi, 2)) + (8505.L*sqrt(2.L)*pow(q, 
13))/(1681.L*pow(Pi, 2)) + (15795.L*pow(q, 
14))/(1936.L*sqrt(2.L)*pow(Pi, 2)) + (11178.L*sqrt(2.L)*pow(q, 
15))/(2209.L*pow(Pi, 2)) + (146043.L*pow(q, 
16))/(20000.L*sqrt(2.L)*pow(Pi, 2)) + (28431.L*pow(q, 
17))/(2809.L*sqrt(2.L)*pow(Pi, 2)) + (103275.L*pow(q, 
18))/(12544.L*sqrt(2.L)*pow(Pi, 2)) + (35235.L*pow(q, 
19))/(3481.L*sqrt(2.L)*pow(Pi, 2)) + (81.L*gamma(2.L, (4.L*Pi*y) / 
(3.L)))/(8.L*sqrt(2.L)*pow(Pi, 2)*q) + (1053.L*gamma(2.L, (16.L*Pi*y) 
/ (3.L)))/(128.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 2)) + (2025.L*gamma(2.L, 
(28.L*Pi*y) / (3.L)))/(196.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 3)) + 
(729.L*gamma(2.L, (40.L*Pi*y) / (3.L)))/(100.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 4)) + (6885.L*gamma(2.L, (52.L*Pi*y) / 
(3.L)))/(676.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 5)) + (16605.L*gamma(2.L, 
(64.L*Pi*y) / (3.L)))/(2048.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 6)) + 
(14661.L*gamma(2.L, (76.L*Pi*y) / (3.L)))/(1444.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 7)) + (3645.L*gamma(2.L, (88.L*Pi*y) / 
(3.L)))/(484.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 8)) + (48681.L*gamma(2.L, 
(100.L*Pi*y) / (3.L)))/(5000.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 9)) + 
(26325.L*gamma(2.L, (112.L*Pi*y) / (3.L)))/(3136.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 10)) + (38961.L*gamma(2.L, (124.L*Pi*y) / 
(3.L)))/(3844.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 11)) + (2187.L*gamma(2.L, 
(136.L*Pi*y) / (3.L)))/(289.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 12)) + 
(55485.L*gamma(2.L, (148.L*Pi*y) / (3.L)))/(5476.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 13)) + (12393.L*gamma(2.L, (160.L*Pi*y) / 
(3.L)))/(1600.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 14)) + 
(74925.L*gamma(2.L, (172.L*Pi*y) / (3.L)))/(7396.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 15)) + (8019.L*gamma(2.L, (184.L*Pi*y) / 
(3.L)))/(1058.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 16)) + 
(198531.L*gamma(2.L, (196.L*Pi*y) / 
(3.L)))/(19208.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 17)) + 
(6885.L*gamma(2.L, (208.L*Pi*y) / (3.L)))/(832.L*sqrt(2.L)*pow(Pi, 
2)*pow(q, 18)) + (2916.L*sqrt(2.L)*gamma(2.L, (220.L*Pi*y) / 
(3.L)))/(605.L*pow(Pi, 2)*pow(q, 19)) + (25515.L*gamma(2.L, 
(232.L*Pi*y) / (3.L)))/(3364.L*sqrt(2.L)*pow(Pi, 2)*pow(q, 20)));











// Weight: w=0


Y0011=1.L;

Y0031=-2.0053522829578814L/(exp(80.L*Pi*y)*pow(q, 20)) - 
1.0051891142646021L/(exp(76.L*Pi*y)*pow(q, 19)) - 
0.15915494309189535L/(exp(72.L*Pi*y)*pow(q, 18)) - 
1.0111019914073351L/(exp(68.L*Pi*y)*pow(q, 17)) - 
1.8501762134432833L/(exp(64.L*Pi*y)*pow(q, 16)) - 
0.38197186342054884L/(exp(60.L*Pi*y)*pow(q, 15)) - 
1.6370222718023522L/(exp(56.L*Pi*y)*pow(q, 14)) - 
1.028385786132247L/(exp(52.L*Pi*y)*pow(q, 13)) - 
0.5570423008216338L/(exp(48.L*Pi*y)*pow(q, 12)) - 
1.041741445692406L/(exp(44.L*Pi*y)*pow(q, 11)) - 
1.7188733853924698L/(exp(40.L*Pi*y)*pow(q, 10)) - 
0.1061032953945969L/(exp(36.L*Pi*y)*pow(q, 9)) - 
1.7904931097838226L/(exp(32.L*Pi*y)*pow(q, 8)) - 
1.091348181201568L/(exp(28.L*Pi*y)*pow(q, 7)) - 
0.477464829275686L/(exp(24.L*Pi*y)*pow(q, 6)) - 
1.1459155902616465L/(exp(20.L*Pi*y)*pow(q, 5)) - 
1.6711269024649011L/(exp(16.L*Pi*y)*pow(q, 4)) - 
0.3183098861837907L/(exp(12.L*Pi*y)*pow(q, 3)) - 
1.432394487827058L/(exp(8.L*Pi*y)*pow(q, 2)) - 
0.954929658551372L/(exp(4.L*Pi*y)*q) - 0.954929658551372L*q - 
1.432394487827058L*pow(q, 2) - 0.3183098861837907L*pow(q, 3) - 
1.6711269024649011L*pow(q, 4) - 1.1459155902616465L*pow(q, 5) - 
0.477464829275686L*pow(q, 6) - 1.091348181201568L*pow(q, 7) - 
1.7904931097838226L*pow(q, 8) - 0.1061032953945969L*pow(q, 9) - 
1.7188733853924698L*pow(q, 10) - 1.041741445692406L*pow(q, 11) - 
0.5570423008216338L*pow(q, 12) - 1.028385786132247L*pow(q, 13) - 
1.6370222718023522L*pow(q, 14) - 0.38197186342054884L*pow(q, 15) - 
1.8501762134432833L*pow(q, 16) - 1.0111019914073351L*pow(q, 17) - 
0.15915494309189535L*pow(q, 18) - 1.0051891142646021L*pow(q, 19) - 
2.0053522829578814L*pow(q, 20) + (-2.471877649503246L + Pi*y)/Pi;

Y0032=q31*((9.L / (2.L) + 270.L/(59.L*exp((236.L*Pi*y) / (3.L))*pow(q, 20)) + 
135.L/(14.L*exp((224.L*Pi*y) / (3.L))*pow(q, 19)) + 
243.L/(53.L*exp((212.L*Pi*y) / (3.L))*pow(q, 18)) + 
837.L/(100.L*exp((200.L*Pi*y) / (3.L))*pow(q, 17)) + 
216.L/(47.L*exp((188.L*Pi*y) / (3.L))*pow(q, 16)) + 
189.L/(22.L*exp((176.L*Pi*y) / (3.L))*pow(q, 15)) + 
189.L/(41.L*exp((164.L*Pi*y) / (3.L))*pow(q, 14)) + 
135.L/(19.L*exp((152.L*Pi*y) / (3.L))*pow(q, 13)) + 
216.L/(35.L*exp((140.L*Pi*y) / (3.L))*pow(q, 12)) + 
567.L/(64.L*exp((128.L*Pi*y) / (3.L))*pow(q, 11)) + 
135.L/(29.L*exp((116.L*Pi*y) / (3.L))*pow(q, 10)) + 
189.L/(26.L*exp((104.L*Pi*y) / (3.L))*pow(q, 9)) + 
108.L/(23.L*exp((92.L*Pi*y) / (3.L))*pow(q, 8)) + 
189.L/(20.L*exp((80.L*Pi*y) / (3.L))*pow(q, 7)) + 
81.L/(17.L*exp((68.L*Pi*y) / (3.L))*pow(q, 6)) + 
54.L/(7.L*exp((56.L*Pi*y) / (3.L))*pow(q, 5)) + 
54.L/(11.L*exp((44.L*Pi*y) / (3.L))*pow(q, 4)) + 
135.L/(16.L*exp((32.L*Pi*y) / (3.L))*pow(q, 3)) + 
27.L/(5.L*exp((20.L*Pi*y) / (3.L))*pow(q, 2)) + 
27.L/(4.L*exp((8.L*Pi*y) / (3.L))*q) + (63.L*q) / (8.L) + 
(36.L*pow(q, 2)) / (7.L) + (81.L*pow(q, 3)) / (10.L) + (63.L*pow(q, 
4)) / (13.L) + (279.L*pow(q, 5)) / (32.L) + (90.L*pow(q, 6)) / (19.L) 
+ (81.L*pow(q, 7)) / (11.L) + (279.L*pow(q, 8)) / (50.L) + 9.L*pow(q, 
9) + (144.L*pow(q, 10)) / (31.L) + (243.L*pow(q, 11)) / (34.L) + 
(171.L*pow(q, 12)) / (37.L) + (81.L*pow(q, 13)) / (8.L) + 
(198.L*pow(q, 14)) / (43.L) + (162.L*pow(q, 15)) / (23.L) + 
(513.L*pow(q, 16)) / (98.L) + (441.L*pow(q, 17)) / (52.L) + 
(324.L*pow(q, 18)) / (55.L) + (405.L*pow(q, 19)) / (58.L))/Pi);

Y0033=q32*((27.L / (4.L) + 405.L/(58.L*exp((232.L*Pi*y) / (3.L))*pow(q, 20)) + 
324.L/(55.L*exp((220.L*Pi*y) / (3.L))*pow(q, 19)) + 
441.L/(52.L*exp((208.L*Pi*y) / (3.L))*pow(q, 18)) + 
513.L/(98.L*exp((196.L*Pi*y) / (3.L))*pow(q, 17)) + 
162.L/(23.L*exp((184.L*Pi*y) / (3.L))*pow(q, 16)) + 
198.L/(43.L*exp((172.L*Pi*y) / (3.L))*pow(q, 15)) + 
81.L/(8.L*exp((160.L*Pi*y) / (3.L))*pow(q, 14)) + 
171.L/(37.L*exp((148.L*Pi*y) / (3.L))*pow(q, 13)) + 
243.L/(34.L*exp((136.L*Pi*y) / (3.L))*pow(q, 12)) + 
144.L/(31.L*exp((124.L*Pi*y) / (3.L))*pow(q, 11)) + 
9.L/(exp((112.L*Pi*y) / (3.L))*pow(q, 10)) + 
279.L/(50.L*exp((100.L*Pi*y) / (3.L))*pow(q, 9)) + 
81.L/(11.L*exp((88.L*Pi*y) / (3.L))*pow(q, 8)) + 
90.L/(19.L*exp((76.L*Pi*y) / (3.L))*pow(q, 7)) + 
279.L/(32.L*exp((64.L*Pi*y) / (3.L))*pow(q, 6)) + 
63.L/(13.L*exp((52.L*Pi*y) / (3.L))*pow(q, 5)) + 
81.L/(10.L*exp((40.L*Pi*y) / (3.L))*pow(q, 4)) + 
36.L/(7.L*exp((28.L*Pi*y) / (3.L))*pow(q, 3)) + 
63.L/(8.L*exp((16.L*Pi*y) / (3.L))*pow(q, 2)) + 
9.L/(2.L*exp((4.L*Pi*y) / (3.L))*q) + (27.L*q) / (5.L) + 
(135.L*pow(q, 2)) / (16.L) + (54.L*pow(q, 3)) / (11.L) + (54.L*pow(q, 
4)) / (7.L) + (81.L*pow(q, 5)) / (17.L) + (189.L*pow(q, 6)) / (20.L) 
+ (108.L*pow(q, 7)) / (23.L) + (189.L*pow(q, 8)) / (26.L) + 
(135.L*pow(q, 9)) / (29.L) + (567.L*pow(q, 10)) / (64.L) + 
(216.L*pow(q, 11)) / (35.L) + (135.L*pow(q, 12)) / (19.L) + 
(189.L*pow(q, 13)) / (41.L) + (189.L*pow(q, 14)) / (22.L) + 
(216.L*pow(q, 15)) / (47.L) + (837.L*pow(q, 16)) / (100.L) + 
(243.L*pow(q, 17)) / (53.L) + (135.L*pow(q, 18)) / (14.L) + 
(270.L*pow(q, 19)) / (59.L))/Pi);

















//weight=1

Y1211=q31*(-3.L*sqrt(2.L)*(1.L + q + 2.L*pow(q, 2) + 2.L*pow(q, 4) + pow(q, 5) 
+ 2.L*pow(q, 6) + pow(q, 8) + 2.L*pow(q, 9) + 2.L*pow(q, 10) + 
2.L*pow(q, 12) + 2.L*pow(q, 14) + 3.L*pow(q, 16) + 2.L*pow(q, 17)));

Y1212=1.L + 6.L*q + 6.L*pow(q, 3) + 6.L*pow(q, 4) + 12.L*pow(q, 7) + 
6.L*pow(q, 9) + 6.L*pow(q, 12) + 12.L*pow(q, 13) + 6.L*pow(q, 16) + 
12.L*pow(q, 19);

Y1231=eulergamma - (6.L*gamma(0.L, 4.L*Pi*y))/q - (6.L*gamma(0.L, 
12.L*Pi*y))/pow(q, 3) - (6.L*gamma(0.L, 16.L*Pi*y))/pow(q, 4) - 
(12.L*gamma(0.L, 28.L*Pi*y))/pow(q, 7) - (6.L*gamma(0.L, 
36.L*Pi*y))/pow(q, 9) - (6.L*gamma(0.L, 48.L*Pi*y))/pow(q, 12) - 
(12.L*gamma(0.L, 52.L*Pi*y))/pow(q, 13) - (6.L*gamma(0.L, 
64.L*Pi*y))/pow(q, 16) - (12.L*gamma(0.L, 76.L*Pi*y))/pow(q, 19) - 
24.L*pow(q, 8)*log(2.L) - 24.L*pow(q, 14)*log(2.L) + log(3.L) - 
12.L*pow(q, 5)*log(5.L) - 12.L*pow(q, 15)*log(5.L) - 12.L*pow(q, 
20)*log(5.L) - 6.L*pow(q, 7)*log(9.L) - 9.L*pow(q, 9)*log(9.L) - 
6.L*pow(q, 13)*log(9.L) - 6.L*pow(q, 19)*log(9.L) - 12.L*pow(q, 
11)*log(11.L) - 12.L*pow(q, 17)*log(17.L) - 3.L*pow(q, 3)*log(81.L) - 
3.L*pow(q, 12)*log(81.L) - q*log(729.L) - pow(q, 4)*log(729.L) - 
pow(q, 16)*log(729.L) - pow(q, 2)*log(4096.L) - pow(q, 6)*log(4096.L) 
- pow(q, 18)*log(4096.L) + log(4.L*Pi) + log(y) - 6.L*log(2.678938535) + 6.L*log(1.354117939);

Y1232=q32*(-((sqrt(2.L)*(3.L*pow(q, 17)*gamma(0.L, (4.L*Pi*y) / (3.L)) + 
3.L*pow(q, 16)*gamma(0.L, (16.L*Pi*y) / (3.L)) + 6.L*pow(q, 
15)*gamma(0.L, (28.L*Pi*y) / (3.L)) + 6.L*pow(q, 13)*gamma(0.L, 
(52.L*Pi*y) / (3.L)) + 3.L*pow(q, 12)*gamma(0.L, (64.L*Pi*y) / (3.L)) 
+ 6.L*pow(q, 11)*gamma(0.L, (76.L*Pi*y) / (3.L)) + 3.L*pow(q, 
9)*gamma(0.L, (100.L*Pi*y) / (3.L)) + 6.L*pow(q, 8)*gamma(0.L, 
(112.L*Pi*y) / (3.L)) + 6.L*pow(q, 7)*gamma(0.L, (124.L*Pi*y) / 
(3.L)) + 6.L*pow(q, 5)*gamma(0.L, (148.L*Pi*y) / (3.L)) + 6.L*pow(q, 
3)*gamma(0.L, (172.L*Pi*y) / (3.L)) + 9.L*q*gamma(0.L, (196.L*Pi*y) / 
(3.L)) + 6.L*gamma(0.L, (208.L*Pi*y) / (3.L)) + 12.L*pow(q, 
20)*log(2.L) + 12.L*pow(q, 22)*log(2.L) + 12.L*pow(q, 26)*log(2.L) + 
18.L*pow(q, 28)*log(2.L) + 12.L*pow(q, 30)*log(2.L) + 24.L*pow(q, 
36)*log(2.L) + 6.L*pow(q, 19)*log(5.L) + 6.L*pow(q, 24)*log(5.L) + 
12.L*pow(q, 29)*log(5.L) + 6.L*pow(q, 21)*log(11.L) + 6.L*pow(q, 
32)*log(11.L) + 6.L*pow(q, 23)*log(17.L) + 6.L*pow(q, 25)*log(23.L) + 
6.L*pow(q, 27)*log(29.L) + 6.L*pow(q, 31)*log(41.L) + 6.L*pow(q, 
33)*log(47.L) + 6.L*pow(q, 35)*log(53.L) + 6.L*pow(q, 37)*log(59.L) + 
pow(q, 18)*log(64.L) + pow(q, 34)*log(64.L)))/pow(q, 18)));












// Weight: w=2


 Y211=-24.L*q - 72.L*pow(q, 2) - 96.L*pow(q, 3) - 168.L*pow(q, 4) - 
144.L*pow(q, 5) - 288.L*pow(q, 6) - 192.L*pow(q, 7) - 360.L*pow(q, 8) 
- 312.L*pow(q, 9) - 432.L*pow(q, 10) - 288.L*pow(q, 11) - 
672.L*pow(q, 12) - 336.L*pow(q, 13) - 576.L*pow(q, 14) - 576.L*pow(q, 
15) - 744.L*pow(q, 16) - 432.L*pow(q, 17) - 936.L*pow(q, 18) - 
480.L*pow(q, 19) - 1008.L*pow(q, 20) + (-3.L + Pi*y)/(Pi*y);


Y231=1.L + 12.L*q + 36.L*pow(q, 2) + 12.L*pow(q, 3) + 84.L*pow(q, 4) + 
72.L*pow(q, 5) + 36.L*pow(q, 6) + 96.L*pow(q, 7) + 180.L*pow(q, 8) + 
12.L*pow(q, 9) + 216.L*pow(q, 10) + 144.L*pow(q, 11) + 84.L*pow(q, 
12) + 168.L*pow(q, 13) + 288.L*pow(q, 14) + 72.L*pow(q, 15) + 
372.L*pow(q, 16) + 216.L*pow(q, 17) + 36.L*pow(q, 18) + 240.L*pow(q, 
19) + 504.L*pow(q, 20);

Y232=q31*(-6.L*(1.L + 7.L*q + 8.L*pow(q, 2) + 18.L*pow(q, 3) + 14.L*pow(q, 4) 
+ 31.L*pow(q, 5) + 20.L*pow(q, 6) + 36.L*pow(q, 7) + 31.L*pow(q, 8) + 
56.L*pow(q, 9) + 32.L*pow(q, 10) + 54.L*pow(q, 11) + 38.L*pow(q, 12) 
+ 90.L*pow(q, 13) + 44.L*pow(q, 14) + 72.L*pow(q, 15) + 57.L*pow(q, 
16) + 98.L*pow(q, 17) + 72.L*pow(q, 18) + 90.L*pow(q, 19)));
 

Y233=q32*(-18.L*(1.L + 2.L*q + 5.L*pow(q, 2) + 4.L*pow(q, 3) + 8.L*pow(q, 4) + 
6.L*pow(q, 5) + 14.L*pow(q, 6) + 8.L*pow(q, 7) + 14.L*pow(q, 8) + 
10.L*pow(q, 9) + 21.L*pow(q, 10) + 16.L*pow(q, 11) + 20.L*pow(q, 12) 
+ 14.L*pow(q, 13) + 28.L*pow(q, 14) + 16.L*pow(q, 15) + 31.L*pow(q, 
16) + 18.L*pow(q, 17) + 40.L*pow(q, 18) + 20.L*pow(q, 19)));














//weight=3


Y3211=q31*(-9.L*sqrt(2.L)*(1.L + 13.L*q + 50.L*pow(q, 2) + 72.L*pow(q, 3) + 
170.L*pow(q, 4) + 205.L*pow(q, 5) + 362.L*pow(q, 6) + 360.L*pow(q, 7) 
+ 601.L*pow(q, 8) + 650.L*pow(q, 9) + 962.L*pow(q, 10) + 864.L*pow(q, 
11) + 1370.L*pow(q, 12) + 1224.L*pow(q, 13) + 1850.L*pow(q, 14) + 
1584.L*pow(q, 15) + 2451.L*pow(q, 16) + 2210.L*pow(q, 17) + 
2880.L*pow(q, 18) + 2520.L*pow(q, 19)));

Y3212=-1.L + 90.L*q + 216.L*pow(q, 2) + 738.L*pow(q, 3) + 1170.L*pow(q, 4) 
+ 1728.L*pow(q, 5) + 2160.L*pow(q, 6) + 4500.L*pow(q, 7) + 
3672.L*pow(q, 8) + 6570.L*pow(q, 9) + 6480.L*pow(q, 10) + 
8640.L*pow(q, 11) + 9594.L*pow(q, 12) + 15300.L*pow(q, 13) + 
10800.L*pow(q, 14) + 17280.L*pow(q, 15) + 18450.L*pow(q, 16) + 
20736.L*pow(q, 17) + 19656.L*pow(q, 18) + 32580.L*pow(q, 19) + 
22464.L*pow(q, 20);

Y3231=sqrt(2.L) + 72.L*sqrt(2.L)*q + 270.L*sqrt(2.L)*pow(q, 2) + 
720.L*sqrt(2.L)*pow(q, 3) + 936.L*sqrt(2.L)*pow(q, 4) + 
2160.L*sqrt(2.L)*pow(q, 5) + 2214.L*sqrt(2.L)*pow(q, 6) + 
3600.L*sqrt(2.L)*pow(q, 7) + 4590.L*sqrt(2.L)*pow(q, 8) + 
6552.L*sqrt(2.L)*pow(q, 9) + 5184.L*sqrt(2.L)*pow(q, 10) + 
10800.L*sqrt(2.L)*pow(q, 11) + 9360.L*sqrt(2.L)*pow(q, 12) + 
12240.L*sqrt(2.L)*pow(q, 13) + 13500.L*sqrt(2.L)*pow(q, 14) + 
17712.L*sqrt(2.L)*pow(q, 15) + 14760.L*sqrt(2.L)*pow(q, 16) + 
25920.L*sqrt(2.L)*pow(q, 17) + 19710.L*sqrt(2.L)*pow(q, 18) + 
26064.L*sqrt(2.L)*pow(q, 19) + 28080.L*sqrt(2.L)*pow(q, 20);

Y3232=q32*(54.L*(1.L + 8.L*q + 17.L*pow(q, 2) + 40.L*pow(q, 3) + 50.L*pow(q, 4) 
+ 96.L*pow(q, 5) + 104.L*pow(q, 6) + 176.L*pow(q, 7) + 170.L*pow(q, 
8) + 280.L*pow(q, 9) + 273.L*pow(q, 10) + 400.L*pow(q, 11) + 
362.L*pow(q, 12) + 560.L*pow(q, 13) + 520.L*pow(q, 14) + 736.L*pow(q, 
15) + 601.L*pow(q, 16) + 936.L*pow(q, 17) + 850.L*pow(q, 18) + 
1160.L*pow(q, 19)));







// Weight: w=4


Y4111=1.L + 240.L*q + 2160.L*pow(q, 2) + 6720.L*pow(q, 3) + 17520.L*pow(q, 
4) + 30240.L*pow(q, 5) + 60480.L*pow(q, 6) + 82560.L*pow(q, 7) + 
140400.L*pow(q, 8) + 181680.L*pow(q, 9) + 272160.L*pow(q, 10) + 
319680.L*pow(q, 11) + 490560.L*pow(q, 12) + 527520.L*pow(q, 13) + 
743040.L*pow(q, 14) + 846720.L*pow(q, 15) + 1123440.L*pow(q, 16) + 
1179360.L*pow(q, 17) + 1635120.L*pow(q, 18) + 1646400.L*pow(q, 19) + 
2207520.L*pow(q, 20);


Y4121=q31*(-12.L*(1.L - 8.L*q + 20.L*pow(q, 2) - 70.L*pow(q, 4) + 64.L*pow(q, 
5) + 56.L*pow(q, 6) - 125.L*pow(q, 8) - 160.L*pow(q, 9) + 
308.L*pow(q, 10) + 110.L*pow(q, 12) - 520.L*pow(q, 14) + 57.L*pow(q, 
16) + 560.L*pow(q, 17)));
 



Y431=1.L - 84.L*q - 756.L*pow(q, 2) - 2028.L*pow(q, 3) - 6132.L*pow(q, 4) 
- 10584.L*pow(q, 5) - 18252.L*pow(q, 6) - 28896.L*pow(q, 7) - 
49140.L*pow(q, 8) - 54516.L*pow(q, 9) - 95256.L*pow(q, 10) - 
111888.L*pow(q, 11) - 148044.L*pow(q, 12) - 184632.L*pow(q, 13) - 
260064.L*pow(q, 14) - 255528.L*pow(q, 15) - 393204.L*pow(q, 16) - 
412776.L*pow(q, 17) - 490644.L*pow(q, 18) - 576240.L*pow(q, 19) - 
772632.L*pow(q, 20);

Y432=q31*(6.L*(1.L + 73.L*q + 344.L*pow(q, 2) + 1134.L*pow(q, 3) + 
2198.L*pow(q, 4) + 4681.L*pow(q, 5) + 6860.L*pow(q, 6) + 
11988.L*pow(q, 7) + 15751.L*pow(q, 8) + 25112.L*pow(q, 9) + 
29792.L*pow(q, 10) + 44226.L*pow(q, 11) + 50654.L*pow(q, 12) + 
73710.L*pow(q, 13) + 79508.L*pow(q, 14) + 109512.L*pow(q, 15) + 
117993.L*pow(q, 16) + 160454.L*pow(q, 17) + 167832.L*pow(q, 18) + 
219510.L*pow(q, 19)));
 

Y433=q32*(54.L*(1.L + 14.L*q + 65.L*pow(q, 2) + 148.L*pow(q, 3) + 344.L*pow(q, 
4) + 546.L*pow(q, 5) + 1022.L*pow(q, 6) + 1352.L*pow(q, 7) + 
2198.L*pow(q, 8) + 2710.L*pow(q, 9) + 4161.L*pow(q, 10) + 
4816.L*pow(q, 11) + 6860.L*pow(q, 12) + 7658.L*pow(q, 13) + 
10804.L*pow(q, 14) + 11536.L*pow(q, 15) + 15751.L*pow(q, 16) + 
16542.L*pow(q, 17) + 22360.L*pow(q, 18) + 22820.L*pow(q, 19)));










//weight=5


Y5211=q31*(-3.L*sqrt(2.L)*(1.L + 241.L*q + 2402.L*pow(q, 2) + 9360.L*pow(q, 3) 
+ 28562.L*pow(q, 4) + 61681.L*pow(q, 5) + 130322.L*pow(q, 6) + 
219600.L*pow(q, 7) + 390001.L*pow(q, 8) + 578882.L*pow(q, 9) + 
923522.L*pow(q, 10) + 1252800.L*pow(q, 11) + 1874162.L*pow(q, 12) + 
2405520.L*pow(q, 13) + 3418802.L*pow(q, 14) + 4197600.L*pow(q, 15) + 
5767203.L*pow(q, 16) + 6883442.L*pow(q, 17) + 9135360.L*pow(q, 18) + 
10609200.L*pow(q, 19)));

Y5212=1.L + 246.L*q + 3600.L*pow(q, 2) + 19686.L*pow(q, 3) + 
59286.L*pow(q, 4) + 149760.L*pow(q, 5) + 295200.L*pow(q, 6) + 
590892.L*pow(q, 7) + 925200.L*pow(q, 8) + 1594326.L*pow(q, 9) + 
2302560.L*pow(q, 10) + 3513600.L*pow(q, 11) + 4744326.L*pow(q, 12) + 
7026252.L*pow(q, 13) + 8647200.L*pow(q, 14) + 12280320.L*pow(q, 15) + 
15173526.L*pow(q, 16) + 20044800.L*pow(q, 17) + 23914800.L*pow(q, 18) 
+ 32059212.L*pow(q, 19) + 36092160.L*pow(q, 20);

Y5221=q32*(-36.L*sqrt(2.L)*(1.L - 7.L*q + 14.L*pow(q, 2) + 4.L*pow(q, 3) - 
28.L*pow(q, 4) - 21.L*pow(q, 5) + 14.L*pow(q, 6) + 188.L*pow(q, 7) - 
112.L*pow(q, 8) - 233.L*pow(q, 9) - 60.L*pow(q, 10) + 196.L*pow(q, 
11) + 560.L*pow(q, 12) - 427.L*pow(q, 13) - 8.L*pow(q, 14) - 
308.L*pow(q, 15) - 257.L*pow(q, 16) + 423.L*pow(q, 17) - 392.L*pow(q, 
18) + 1064.L*pow(q, 19)));

Y5222=q31*(12.L*(1.L - 2.L*q - 28.L*pow(q, 2) + 126.L*pow(q, 3) - 112.L*pow(q, 
4) - 284.L*pow(q, 5) + 560.L*pow(q, 6) - 72.L*pow(q, 7) - 
257.L*pow(q, 8) + 56.L*pow(q, 9) - 364.L*pow(q, 10) + 378.L*pow(q, 
11) - 826.L*pow(q, 12) + 1764.L*pow(q, 13) + 1736.L*pow(q, 14) - 
3384.L*pow(q, 15) - 1617.L*pow(q, 16) + 224.L*pow(q, 17) + 
504.L*pow(q, 18) + 4194.L*pow(q, 19)));

Y5231=-sqrt(2.L) + 240.L*sqrt(2.L)*q + 3690.L*sqrt(2.L)*pow(q, 2) + 
19680.L*sqrt(2.L)*pow(q, 3) + 57840.L*sqrt(2.L)*pow(q, 4) + 
153504.L*sqrt(2.L)*pow(q, 5) + 295290.L*sqrt(2.L)*pow(q, 6) + 
576480.L*sqrt(2.L)*pow(q, 7) + 948330.L*sqrt(2.L)*pow(q, 8) + 
1594320.L*sqrt(2.L)*pow(q, 9) + 2246400.L*sqrt(2.L)*pow(q, 10) + 
3601440.L*sqrt(2.L)*pow(q, 11) + 4742880.L*sqrt(2.L)*pow(q, 12) + 
6854880.L*sqrt(2.L)*pow(q, 13) + 8863380.L*sqrt(2.L)*pow(q, 14) + 
12284064.L*sqrt(2.L)*pow(q, 15) + 14803440.L*sqrt(2.L)*pow(q, 16) + 
20545920.L*sqrt(2.L)*pow(q, 17) + 23914890.L*sqrt(2.L)*pow(q, 18) + 
31277280.L*sqrt(2.L)*pow(q, 19) + 36994464.L*sqrt(2.L)*pow(q, 20);

Y5232=q32*(18.L*(5.L + 208.L*q + 1285.L*pow(q, 2) + 4880.L*pow(q, 3) + 
12010.L*pow(q, 4) + 27840.L*pow(q, 5) + 50128.L*pow(q, 6) + 
93280.L*pow(q, 7) + 142810.L*pow(q, 8) + 235760.L*pow(q, 9) + 
328965.L*pow(q, 10) + 499616.L*pow(q, 11) + 651610.L*pow(q, 12) + 
941920.L*pow(q, 13) + 1176080.L*pow(q, 14) + 1626560.L*pow(q, 15) + 
1950005.L*pow(q, 16) + 2630160.L*pow(q, 17) + 3086570.L*pow(q, 18) + 
4039120.L*pow(q, 19)));





// Weight: w=6

Y611=1.L - 504.L*q - 16632.L*pow(q, 2) - 122976.L*pow(q, 3) - 
532728.L*pow(q, 4) - 1575504.L*pow(q, 5) - 4058208.L*pow(q, 6) - 
8471232.L*pow(q, 7) - 17047800.L*pow(q, 8) - 29883672.L*pow(q, 9) - 
51991632.L*pow(q, 10) - 81170208.L*pow(q, 11) - 129985632.L*pow(q, 
12) - 187132176.L*pow(q, 13) - 279550656.L*pow(q, 14) - 
384422976.L*pow(q, 15) - 545530104.L*pow(q, 16) - 715608432.L*pow(q, 
17) - 986161176.L*pow(q, 18) - 1247954400.L*pow(q, 19) - 
1665307728.L*pow(q, 20);

Y6311=1.L + 252.L*q + 5076.L*pow(q, 2) + 41292.L*pow(q, 3) + 
178884.L*pow(q, 4) + 528552.L*pow(q, 5) + 1333476.L*pow(q, 6) + 
2835936.L*pow(q, 7) + 5727780.L*pow(q, 8) + 9858492.L*pow(q, 9) + 
17422776.L*pow(q, 10) + 27158544.L*pow(q, 11) + 42858324.L*pow(q, 12) 
+ 62773128.L*pow(q, 13) + 93715488.L*pow(q, 14) + 126745992.L*pow(q, 
15) + 182748132.L*pow(q, 16) + 239920056.L*pow(q, 17) + 
325067796.L*pow(q, 18) + 418224240.L*pow(q, 19) + 558154584.L*pow(q, 
20);

Y6312=q31*(-6.L*(1.L + 247.L*q + 3848.L*pow(q, 2) + 23778.L*pow(q, 3) + 
86174.L*pow(q, 4) + 248911.L*pow(q, 5) + 570980.L*pow(q, 6) + 
1229076.L*pow(q, 7) + 2251951.L*pow(q, 8) + 4099736.L*pow(q, 9) + 
6610112.L*pow(q, 10) + 10808694.L*pow(q, 11) + 16000598.L*pow(q, 12) 
+ 24401610.L*pow(q, 13) + 33932444.L*pow(q, 14) + 49019112.L*pow(q, 
15) + 65178777.L*pow(q, 16) + 90569138.L*pow(q, 17) + 
116177832.L*pow(q, 18) + 156178890.L*pow(q, 19)));
 

Y6313=q32*(-18.L*(1.L + 242.L*q + 2645.L*pow(q, 2) + 12244.L*pow(q, 3) + 
42728.L*pow(q, 4) + 109446.L*pow(q, 5) + 254174.L*pow(q, 6) + 
494888.L*pow(q, 7) + 941534.L*pow(q, 8) + 1578970.L*pow(q, 9) + 
2664741.L*pow(q, 10) + 4041616.L*pow(q, 11) + 6286340.L*pow(q, 12) + 
8910254.L*pow(q, 13) + 13094188.L*pow(q, 14) + 17637136.L*pow(q, 15) 
+ 24802351.L*pow(q, 16) + 32177538.L*pow(q, 17) + 43731400.L*pow(q, 
18) + 54989540.L*pow(q, 19)));


Y6321=216.L*q - 1296.L*pow(q, 2) + 1944.L*pow(q, 3) + 864.L*pow(q, 4) + 
1296.L*pow(q, 5) - 11664.L*pow(q, 6) - 8640.L*pow(q, 7) + 
36288.L*pow(q, 8) + 17496.L*pow(q, 9) - 7776.L*pow(q, 10) - 
121824.L*pow(q, 11) + 7776.L*pow(q, 12) + 137808.L*pow(q, 13) + 
51840.L*pow(q, 14) + 11664.L*pow(q, 15) - 245376.L*pow(q, 16) + 
190512.L*pow(q, 17) - 104976.L*pow(q, 18) - 120096.L*pow(q, 19) + 
5184.L*pow(q, 20);

Y6322=q31*(12.L*(-1.L - 4.L*q + 40.L*pow(q, 2) + 36.L*pow(q, 3) - 638.L*pow(q, 
4) + 1136.L*pow(q, 5) + 556.L*pow(q, 6) - 3384.L*pow(q, 7) + 
3089.L*pow(q, 8) + 160.L*pow(q, 9) - 4400.L*pow(q, 10) + 
5292.L*pow(q, 11) + 2410.L*pow(q, 12) - 1008.L*pow(q, 13) - 
9644.L*pow(q, 14) - 5040.L*pow(q, 15) + 15207.L*pow(q, 16) - 
2552.L*pow(q, 17) + 3384.L*pow(q, 18) + 27828.L*pow(q, 19)));
 

Y6323=q32*(72.L*(1.L - q - 28.L*pow(q, 2) + 94.L*pow(q, 3) - 40.L*pow(q, 4) - 
147.L*pow(q, 5) - 4.L*pow(q, 6) + 140.L*pow(q, 7) + 638.L*pow(q, 8) - 
773.L*pow(q, 9) - 240.L*pow(q, 10) + 40.L*pow(q, 11) - 556.L*pow(q, 
12) + 1145.L*pow(q, 13) + 376.L*pow(q, 14) + 3112.L*pow(q, 15) - 
3089.L*pow(q, 16) - 5625.L*pow(q, 17) + 1120.L*pow(q, 18) + 
3014.L*pow(q, 19)));






	//+++++++++++++++++++++++++++++++Lepton+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	complex<long double> me1, me2, me3, me4, me5, me6, me7, me8, me9;
	complex<long double> mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9;
  complex<long double> Mnu1,Mnu2,Mnu3,Mnu4,Mnu5,Mnu6,Mnu7,Mnu8,Mnu9;

cx_lmat vect1 =  {{(sqrt(2.L)*beta*Y3211),(sqrt(2.L)*alpha*Y3231-beta*Y3212),-alpha*Y3232},{-sqrt(2.L)*alpha*Y3232,(-beta*Y3211),(-alpha*Y3231-sqrt(2.L)*beta*Y3212)},{(gammainput*Y0233),gammainput*Y0232,gammainput*Y0231}};

cx_lmat vecnu1 ={{g1 + 2.L*g2*Y0031,-(g2*Y0033),-(g2*Y0032)},{-(g2*Y0033),2.L*g2*Y0032,g1 - g2*Y0031},{-(g2*Y0032),g1 - g2*Y0031,2.L*g2*Y0033}};




  me1 = vect1[0][0];
	me2 = vect1[0][1];
	me3 = vect1[0][2];
	me4 = vect1[1][0];
	me5 = vect1[1][1];
	me6 = vect1[1][2];
	me7 = vect1[2][0];
	me8 = vect1[2][1];
	me9 = vect1[2][2];
       
 





  Mnu1 = vecnu1[0][0];
	Mnu2 = vecnu1[0][1];
	Mnu3 = vecnu1[0][2];
	Mnu4 = vecnu1[1][0];
	Mnu5 = vecnu1[1][1];
	Mnu6 = vecnu1[1][2];
	Mnu7 = vecnu1[2][0];
	Mnu8 = vecnu1[2][1];
	Mnu9 = vecnu1[2][2];
 


	
	complex<long double> Ml[9] = {me1, me2, me3, me4, me5, me6, me7, me8, me9};
	complex<long double> Ml_sq[9];
	Matrixsq(Ml, Ml_sq); 
 //++++++Neutrino mass matrix+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

complex<long double> Mnu[9] = {Mnu1, Mnu2, Mnu3,Mnu4, Mnu5, Mnu6,Mnu7, Mnu8, Mnu9};

	complex<long double> Mnu_sq[9];
	Matrixsq(Mnu, Mnu_sq);
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

	struct Coef lres = extractcoef(Ml_sq, 1.);
	struct Coef nures = extractcoef(Mnu_sq, 0.);

	long double memu, mmutau, Delta2131, s13sql, s12sql, s23sql, m1, m2, m3, msum, mee, overlambda, dcp = 0.;
     
	memu = pow(lres.lambda1 / lres.lambda2, 0.5);
	mmutau = pow(lres.lambda2 / lres.lambda3, 0.5);

	complex<long double> Ul_ct[9] = {conj(lres.v11), conj(lres.v12), conj(lres.v13),
					conj(lres.v21), conj(lres.v22), conj(lres.v23),
					conj(lres.v31), conj(lres.v32), conj(lres.v33)};
	
        // find the phase matrix of "num" 
        // The elements of "U_{nu}^{T}M_{nu}U_{nu}=hat{M_{nu}}*diag(e^{i2phase_1},e^{i2phase_2},e^{i2phase_3})", thus the exact diagonalization-matrix should be U_{nu}*diag(e^{-iphase_1},e^{-iphase_2},e^{-iphase_3})
        //The elements of the matrix "U_{nu}^{T}M_{nu}" is
  
        complex<long double> Unu[9] = {nures.v11, nures.v21, nures.v31, nures.v12, nures.v22, nures.v32, nures.v13, nures.v23, nures.v33};
        complex<long double> UnuT[9] = {nures.v11, nures.v12, nures.v13, nures.v21, nures.v22, nures.v23, nures.v31, nures.v32, nures.v33};
	complex<long double> UTM[9], UTMU[9];
	multimatrix(UnuT, Mnu, UTM);
        multimatrix(UTM, Unu, UTMU);
  
        //The diagonal elements of the matrix "U_{nu}^{T}M_{nu}U_{nu}" is
        //long double phase1 = -(UTMU[0]).Theta()/2, phase2 = -(UTMU[4]).Theta()/2, phase3 = -(UTMU[8]).Theta()/2;
	long double phase1 = -arg(UTMU[0])/2, phase2 = -arg(UTMU[4])/2, phase3 = -arg(UTMU[8])/2;

        complex<long double> Unup[9] = {nures.v11 * exp(i * phase1), nures.v21 * exp(i * phase2), nures.v31 * exp(i * phase3),
					nures.v12 * exp(i * phase1), nures.v22 * exp(i * phase2), nures.v32 * exp(i * phase3),
					nures.v13 * exp(i * phase1), nures.v23 * exp(i * phase2), nures.v33 * exp(i * phase3)};

        complex<long double> Upmns[9];
	multimatrix(Ul_ct, Unup, Upmns);	
        //Construct the PMNS matrix: U_{PMNS}=em^{-1}*num;
        complex<long double> Upmns11 = Upmns[0], Upmns12 = Upmns[1], Upmns13 = Upmns[2], Upmns23 = Upmns[5], Upmns31 = Upmns[6], Upmns33 = Upmns[8];

        //long double lsin13sq,lsin12sq,lsin23sq,ltheta13,ltheta12,ltheta23,ldelta,lalpha1,lbeta,lalpha2,lalpha11,lalpha22 = 0.;
        long double lsin13sq = pow(abs(Upmns13), 2);
        long double lsin12sq = pow(abs(Upmns12), 2) / (1. - pow(abs(Upmns13), 2));
        long double lsin23sq = pow(abs(Upmns23), 2) / (1. - pow(abs(Upmns13), 2));
       	long double ldcp = arg(((Upmns11 * Upmns33 * conj(Upmns13) * conj(Upmns31))/(pow(1.0-lsin12sq, 0.5) * pow(1.0-lsin23sq, 0.5) * (1.0-lsin13sq) * pow(lsin13sq, 0.5)) + pow(1.0-lsin12sq, 0.5) * pow(1.0-lsin23sq, 0.5) * pow(lsin13sq, 0.5))/(pow(lsin12sq, 0.5) * pow(lsin23sq, 0.5)));

        // ltheta13,ltheta12,ltheta23 represent the mixing angles in PMNS matrix
        long double ltheta13 = asin(pow(lsin13sq, 0.5));
        long double ltheta12 = asin(pow(lsin12sq, 0.5));
        long double ltheta23 = asin(pow(lsin23sq, 0.5));

        //Construct kaisq function
        //mbtau = (pow(dres.lambda3, 0.5) / Ml33).Re();
        Delta2131 = (nures.lambda2 - nures.lambda1) / (nures.lambda2 - nures.lambda3);
        //Delta2131 = (nures.lambda2 - nures.lambda1) / (nures.lambda3 - (nures.lambda1 + nures.lambda2) / 2);
 long double overpara = 0.0000742;
	overlambda = overpara / (nures.lambda2 - nures.lambda1);
  m1 = sqrt(overlambda * nures.lambda1);
	m2 = sqrt(overlambda * nures.lambda2);
	m3 = sqrt(overlambda * nures.lambda3);
	msum = m1 + m2 + m3;

	s13sql = lsin13sq;
	s12sql = lsin12sq;
        s23sql = lsin23sq;
       	if(ldcp < 0.)
      {dcp = ldcp/Pi + 2.L;}
        else
         {dcp = ldcp/Pi;}

         long double obs[8] = {memu, mmutau, Delta2131, s13sql, s12sql, s23sql, dcp, msum};

	reschisq[1] = obs[0];
	reschisq[2] = obs[1];
	reschisq[3] = obs[2];
	reschisq[4] = obs[3];
	reschisq[5] = obs[4];
	reschisq[6] = obs[5];
  reschisq[7] = obs[6];
  reschisq[8] = obs[7];

	//NH  ordering****NuFit 6.1 ****************************************************************************************************************************
	long double mean[8] = {0.00473694, 0.05882, 0.03035440998791784, 0.02262, 0.3088, 0.550, 1.52222,0.12};
	long double err_up[8] = {0.00000473694, 0.00005882,0.000628132, 0.00057, 0.0067, 0.013, 0.122222,0.00001};
	long double err_down[8] = {0.00000473694, 0.00005882, 0.000642065, 0.00056, 0.0066, 0.016, 0.138889,100.};

	//tan\beta = 10_fixed_v2********NH  ordering****NuFit 4.1
	//long double mean[15] = {0.00000544127,0.0028213,0.000102026,0.00201939,0.110706,0.22736,0.00349377,0.0401458,69.2133, 0.00473689, 0.0585684, 0.0292326, 0.02237, 0.310, 0.563};
	// long double err_up[15] = {0.00000169179,0.000119545,0.0000114843,0.000119237,0.00290982,0.000727552,0.000125776,0.00342333,3.1146, 0.0002, 0.0045, 0.000882, 0.00066, 0.013, 0.018};
	// long double err_down[15] = {0.00000169179,0.000119545,0.0000114843,0.000119237,0.00290982,0.000727552,0.000125776,0.00342333,3.1146, 0.0002, 0.0045, 0.000882, 0.00065, 0.012, 0.024};


	 //long double mean[14] = {0.00192864, 0.0028213, 0.0505231, 0.018241,0.22736, 0.00349377,0.0401458,69.2133, 0.00473689, 0.0585684, 0.0292326, 0.02237, 0.310, 0.563};
	 //long double err_up[14] = {0.000601677, 0.000119545, 0.0061911, 0.00100515, 0.000727552,0.000125776,0.00342333,3.1146, 0.0000401939, 0.000465432, 0.000882, 0.00066, 0.013, 0.018};
	 //long double err_down[14] = {0.000601677, 0.000119545, 0.0061911, 0.00100515, 0.000727552,0.000125776,0.00342333,3.1146, 0.0000401939, 0.000465432, 0.000882, 0.00065, 0.012, 0.024};

	long double chisq = 0.;
	for (int m = 0; m < 8; m++)
       {
          long double err = obs[m] >= mean[m] ? err_up[m] : err_down[m];
          long double z = (obs[m] - mean[m]) / err;
          chisq += z * z;
       }
	//      }
	reschisq[0] = chisq;
     
      return reschisq;

}      


void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  long double chisq1[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  calculate_chisq(par, chisq1);
  f = chisq1[0];
  return;
}


int Weinberg_operator()
{

   

  TMinuit *ptMinuit = new TMinuit();
  ptMinuit->SetPrintLevel(-1);
  ptMinuit->SetFCN(calc_chi_square);

  Double_t arglist[14];
  Int_t ierflg = 0;
  Int_t flag = 0;

  arglist[0] = 1;
  ptMinuit->mnexcm("SET ERR",arglist,1,ierflg); 

  Double_t vstart[14] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t step[14] = {0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  Int_t icstat = 9; 
  Double_t amin = 9999;
  Double_t edm = 10.;
  
  Double_t errdef;
  Int_t nvpar, nparx, npar;
  
  double fParamVal, fParamErr;
  
  //Double_t min = 50000, min2 = 50000;
  
  Double_t min2 = pow(10.,10.);
  
  Double_t fpar[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Double_t dowpar[5] = {-0.1,  1.,   -1.0,    -1.0,     -2.};
 // Double_t uppar[5] = {0.,   1.1,   0.0,   0.0,  -1.};
 //Double_t dowpar[5] = {-0.5,  0.96,   -10.0,    -5.0,     -5.};
// Double_t uppar[5] = {0.,   2.4,     0.0 ,   0.0,  0.};
Double_t dowpar[6] = {0.3647,  1.10359,   1.0,  0.0,  -30};
 Double_t uppar[6] = {0.36473,   1.10693,    1.2,   0.01,  -29};


  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dist1(dowpar[0], uppar[0]), dist2(dowpar[1], uppar[1]), dist3(dowpar[2], uppar[2]), dist4(dowpar[3], uppar[3]), dist5(dowpar[4], uppar[4]);//,dist6(dowpar[5], uppar[5]) , dist7(dowpar[6], uppar[6]),dist8(dowpar[7], uppar[7]), dist9(dowpar[8], uppar[8]), dist10(dowpar[9], uppar[9]);//, dist11(dowpar[10], uppar[10]), dist12(dowpar[11], uppar[11]);

int parNum = 5;
lvec minParams(parNum);


  for (int m = 1; m <= 100; m++)
  {

	flag = 0;
    //int seed = 100000*m + 55;
    //gRandom->SetSeed( seed );
    
    //Double_t min1 = 50000;
    Double_t min1 = pow(10.,10.);
    //Double_t dmin = 10;

    vstart[0] = dist1(gen);
    vstart[1] = dist2(gen);
    vstart[2] = dist3(gen);
    vstart[3] = dist4(gen);
    vstart[4] = dist5(gen);
    //vstart[5] = dist6(gen);
   // vstart[6] = dist7(gen);
   // vstart[7] = dist8(gen);
  //  vstart[8] = dist9(gen);
   // vstart[9] = dist10(gen);
   // vstart[10] = dist11(gen);
  //  vstart[11] = dist12(gen);
    //vstart[12] = dist13(gen);
    //vstart[13] = dist14(gen);
  int n,flagt=0;     

 //do{
    for(flag = 0; flag < 20; ++flag){

    flagt = flagt + 1;   


    ptMinuit->mnparm(0, "a1", vstart[0], step[0], dowpar[0], uppar[0], ierflg);
    ptMinuit->mnparm(1, "a2", vstart[1], step[1], dowpar[1], uppar[1], ierflg);
    ptMinuit->mnparm(2, "a3", vstart[2], step[2], dowpar[2], uppar[2], ierflg);
    ptMinuit->mnparm(3, "a4", vstart[3], step[3], dowpar[3], uppar[3], ierflg);
    ptMinuit->mnparm(4, "a5", vstart[4], step[4], dowpar[4], uppar[4], ierflg);
    //ptMinuit->mnparm(5, "a6", vstart[5], step[5], dowpar[5], uppar[5], ierflg);
   // ptMinuit->mnparm(6, "a7", vstart[6], step[6], dowpar[6], uppar[6], ierflg);
  //  ptMinuit->mnparm(7, "a8", vstart[7], step[7], dowpar[7], uppar[7], ierflg);
   // ptMinuit->mnparm(8, "a9", vstart[8], step[8], dowpar[8], uppar[8], ierflg);
  //  ptMinuit->mnparm(9, "a10", vstart[9], step[9], dowpar[9], uppar[9], ierflg);
  //  ptMinuit->mnparm(10, "a11", vstart[10], step[10], dowpar[10], uppar[10], ierflg);
   // ptMinuit->mnparm(11, "a12", vstart[11], step[11], dowpar[11], uppar[11], ierflg);
    //ptMinuit->mnparm(12, "a13", vstart[12], step[12], dowpar[12], uppar[12], ierflg);
    //ptMinuit->mnparm(13, "a14", vstart[13], step[13], dowpar[13], uppar[13], ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    ptMinuit->mnexcm("MIGRAD", arglist, 7, ierflg);

    ptMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);


	for (n = 0; n < parNum; ++n)
	{
	  ptMinuit->GetParameter(n,vstart[n],fParamErr);	
	}

    if (amin < min1){
	for(n = 0; n < parNum; ++n){
          minParams[n] = vstart[n];
         }
           min1 = amin;	
	   //flagt = 10;
	}
    
	for(n = 0; n < parNum; ++n){
        //par[n] = minParams[n]*randx(gen);
	vstart[n] = gRandom->Uniform(minParams[n]*0.9, minParams[n]*1.1);

	if(vstart[n]<dowpar[n]){
	  vstart[n]=dowpar[n];
	  }
	else if(vstart[n]>uppar[n]){
          vstart[n]=uppar[n];
          }
       }

    //if(flagt>30 || (flag>2 && min1>50000)){
    if(flagt>20){
         break;
       }

//    dmin = abs(amin - min1);
//    min1 = amin;

  }


    //dmin = abs(amin - min1);
    //min1 = amin;

  //}while(dmin > 0.0001);

	//myfile << m <<": chisq = " << min1 << endl;

   if (min1 < min2)
     {
        min2 = min1;
	
        cout << m <<": chisq = " << min1 << endl;

	fpar[0] = minParams[0];
	fpar[1] = minParams[1];
	fpar[2] = minParams[2];
	fpar[3] = minParams[3];
	fpar[4] = minParams[4];
	//fpar[5] = minParams[5];
  //fpar[6] = minParams[6];
	//fpar[7] = minParams[7];
//	fpar[8] = minParams[8];
	//fpar[9] = minParams[9];
//	fpar[10] = minParams[10];
//	fpar[11] = minParams[11];
	//fpar[12] = minParams[12];
	//fpar[13] = minParams[13];
	
     }
  }
  
  cout << " The scanning is completed." << endl;
 cout << " fchisq = " << min2 << endl;
  cout << setprecision(20);    
  cout << " fpar = {" << fpar[0] << ", " << fpar[1] << ", " << fpar[2] << ", " << fpar[3] <<", " << fpar[4] <<","<< fpar[5] <<"}" << endl;
  //myfile << "         " << fpar[4] << ", " << fpar[5] << "}" << endl;
  
  long double check[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  calculate_chisq(fpar, check);
  
  cout << "Good point as: " << "{" << check[0] << "," << check[1] << "," << check[2] << "," << check[3] << "," << check[4] << "," << check[5] << "," << check[6] << "," << check[7] <<"," << check[8] << "}" << endl;



  
 return 0;
}  


