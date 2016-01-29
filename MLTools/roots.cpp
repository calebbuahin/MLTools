#include "roots.h"
#include <stdlib.h>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <omp.h>

using namespace std;
using namespace Eigen;

//double Roots::c_b41 = 1;
//double Roots::a = 0;
//double Roots::b = 0;
//double Roots::c = 0;
//double Roots::d = 0;
//double Roots::e = 0;
//double Roots::f = 0;
//double Roots::g = 0;
//double Roots::h = 0;
//double Roots::u = 0;
//double Roots::v = 0;
//double Roots::a1 = 0;
//double Roots::a2 = 0;
//double Roots::a3 = 0;
//double Roots::a6 = 0;
//double Roots::a7 = 0;
//double Roots::si = 0;
//double Roots::sr = 0;
//double Roots::lzr = 0;
//double Roots::lzi = 0;
//double Roots::szi = 0;
//double Roots::szr = 0;
//float Roots::are = 0;
//float Roots::eta = 0;
//float Roots::mre = 0;
//int Roots::n = 0;
//int Roots::nn = 0;
//double Roots::k[101];
//double Roots::p[101];
//double Roots::qk[101];
//double Roots::qp[101];
//double Roots::svk[101];

complex<double>* Roots::roots(const std::vector<double> &coeffs , int & size)
{
  int matsz = coeffs.size() - 1;
  complex<double>* vret = new complex<double>[matsz];
  size = matsz;

  MatrixXd companion_mat(matsz,matsz);

  #pragma omp parallel num_threads(8)
  {
    #pragma omp for
    {

      for(int n = 0; n < matsz; n++)
        {
          for(int m = 0; m < matsz; m++)
            {
              if(n == m + 1)
                companion_mat(n,m) = 1.0;
              else
                companion_mat(n,m) = 0.0;

              if(n == 0)
                {
                  companion_mat(n,m) = (-1.0*coeffs[m+1])/(coeffs[0]*1.0);
                }
            }
        }
    }

    //cout << companion_mat << endl;

    MatrixXcd eig = companion_mat.eigenvalues();


    #pragma omp for
    {
      for(int i = 0; i < matsz; i++)
        vret[i] =  eig(i) ;
    }
  }

  return vret;
}

//void Roots::roots(const vector<double>& op,int& degree,double* real,double* imag, bool& fail)
//{
//  /* System generated locals */
//  int i__1;

//  /* Local variables */
//  static float base;
//  static double temp[101];
//  static float cosr, sinr;
//  static int i, j, l;
//  static double t;
//  static float x, infin;
//  static bool zerok;
//  static double aa, bb, cc;
//  static float df, ff;
//  static int jj;
//  static float sc, lo, dx, pt[101], xm;
//  static int nz;
//  static double factor;
//  static float xx, yy, smalno;
//  static int nm1;
//  static float bnd, min_, max_;
//  static int cnt;
//  static float xxx;

//  /* FINDS THE ZEROS OF A REAL POLYNOMIAL                 */
//  /* OP  - DOUBLE PRECISION VECTOR OF COEFFICIENTS IN     */
//  /*       ORDER OF DECREASING POWERS.                    */
//  /* DEGREE   - int DEGREE OF POLYNOMIAL.             */
//  /* ZEROR, ZEROI - OUTPUT DOUBLE PRECISION VECTORS OF    */
//  /*                REAL AND IMAGINARY PARTS OF THE       */
//  /*                ZEROS.                                */
//  /* FAIL  - OUTPUT LOGICAL PARAMETER, TRUE ONLY IF       */
//  /*         LEADING COEFFICIENT IS ZERO OR IF RPOLY      */
//  /*         HAS FOUND FEWER THAN DEGREE ZEROS.           */
//  /*         IN THE LATTER CASE DEGREE IS RESET TO        */
//  /*         THE NUMBER OF ZEROS FOUND.                   */
//  /* TO CHANGE THE SIZE OF POLYNOMIALS WHICH CAN BE       */
//  /* SOLVED, RESET THE DIMENSIONS OF THE ARRAYS IN THE    */
//  /* COMMON AREA AND IN THE FOLLOWING DECLARATIONS.       */
//  /* THE SUBROUTINE USES SINGLE PRECISION CALCULATIONS    */
//  /* FOR SCALING, BOUNDS AND ERROR CALCULATIONS. ALL      */
//  /* CALCULATIONS FOR THE ITERATIONS ARE DONE IN DOUBLE   */
//  /* PRECISION.                                           */
//  /* THE FOLLOWING STATEMENTS SET MACHINE CONSTANTS USED  */
//  /* IN VARIOUS PARTS OF THE PROGRAM. THE MEANING OF THE  */
//  /* FOUR CONSTANTS ARE...                                */
//  /* ETA     THE MAXIMUM RELATIVE REPRESENTATION ERROR    */
//  /*         WHICH CAN BE DESCRIBED AS THE SMALLEST       */
//  /*         POSITIVE FLOATING POINT NUMBER SUCH THAT     */
//  /*         1.D0+ETA IS GREATER THAN 1.                  */
//  /* INFINY  THE LARGEST FLOATING-POINT NUMBER.           */
//  /* SMALNO  THE SMALLEST POSITIVE FLOATING-POINT NUMBER  */
//  /*         IF THE EXPONENT RANGE DIFFERS IN SINGLE AND  */
//  /*         DOUBLE PRECISION THEN SMALNO AND INFIN       */
//  /*         SHOULD INDICATE THE SMALLER RANGE.           */
//  /* BASE    THE BASE OF THE FLOATING-POINT NUMBER        */
//  /*         SYSTEM USED.                                 */
//  /* THE VALUES BELOW CORRESPOND TO THE BURROUGHS B6700   */
//  /* changed for sparc, but these seem better -- awf */

//  base = 2.0f;
//  eta = 2.23e-16f;
//  /* the sun compiler will not compile with the number too large for float */
//#ifdef __SUNPRO_C
//  infin = (float)3.40282346638528860e+38;
//#else
//  infin = (float)1e50; /* on purpose too large to fit in `float' type */
//#endif
//  smalno = (float)1e-45;
//  /* ARE AND MRE REFER TO THE UNIT ERROR IN + AND * */
//  /* RESPECTIVELY. THEY ARE ASSUMED TO BE THE SAME AS */
//  /* ETA. */
//  are = eta;
//  mre = eta;
//  lo = smalno / eta;
//  /* INITIALIZATION OF CONSTANTS FOR SHIFT ROTATION */
//  xx = 0.70710678f;
//  yy = -xx;
//  cosr = -0.069756474f;
//  sinr = 0.99756405f;
//  fail = false;
//  n = degree;
//  nn = n + 1;
//  /* ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO. */
//  if (op[0] == 0.) {
//      fail = true;
//      degree = 0;
//      return;
//    }
//  /* REMOVE THE ZEROS AT THE ORIGIN IF ANY */
//  while (op[nn - 1] == 0.) {
//      j = degree - n;
//      real[j] = 0.;
//      imag[j] = 0.;
//      --nn;
//      --n;
//    }

//  /* MAKE A COPY OF THE COEFFICIENTS */
//  for (i = 0; i < nn; ++i) {
//      p[i] = op[i];
//    }
//  /* START THE ALGORITHM FOR ONE ZERO */
//L40:
//  if (n > 2) {
//      goto L60;
//    }
//  if (n < 1) {
//      return;
//    }
//  /* CALCULATE THE FINAL ZERO OR PAIR OF ZEROS */
//  if (n != 2)
//    {
//      real[degree-1] = -p[1] / p[0];
//      imag[degree-1] = 0.;
//      return;
//    }
//  quad(p, &p[1], &p[2], &real[degree-2], &imag[degree-2], &real[degree-1], &imag[degree-1]);
//  return;
//  /* FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS. */
//L60:
//  max_ = 0.0f;
//  min_ = infin;
//  for (i = 0; i < nn; ++i) {
//      x = (double)abs(p[i]);
//      if (x > max_) {
//          max_ = x;
//        }
//      if (x != 0.0f && x < min_) {
//          min_ = x;
//        }
//    }
//  /* SCALE IF THERE ARE LARGE OR VERY SMALL COEFFICIENTS */
//  /* COMPUTES A SCALE FACTOR TO MULTIPLY THE */
//  /* COEFFICIENTS OF THE POLYNOMIAL. THE SCALING IS DONE */
//  /* TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW */
//  /* INTERFERING WITH THE CONVERGENCE CRITERION. */
//  /* THE FACTOR IS A POWER OF THE BASE */
//  sc = lo / min_;
//  if (sc > 1.0f) {
//      goto L80;
//    }
//  if (max_ < 10.0f) {
//      goto L110;
//    }
//  if (sc == 0.0f) {
//      sc = smalno;
//    }
//  goto L90;
//L80:
//  if (infin / sc < max_) {
//      goto L110;
//    }
//L90:
//  l = (int)(log((double)sc) / log((double)base) + 0.5);
//  factor = (double)base;
//  factor = pow_di(&factor, &l);
//  if (factor == 1.) {
//      goto L110;
//    }
//  for (i = 0; i < nn; ++i) {
//      p[i] *= factor;
//    }
//  /* COMPUTE LOWER BOUND ON MODULI OF ZEROS. */
//L110:
//  for (i = 0; i < nn; ++i) {
//      pt[i] = (float)abs(p[i]);
//    }
//  pt[nn - 1] = -pt[nn - 1];
//  /* COMPUTE UPPER ESTIMATE OF BOUND */
//  x = (float)exp((log(-pt[nn - 1]) - log(pt[0])) / n);
//  if (pt[n - 1] == 0.0f) {
//      goto L130;
//    }
//  /* IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT. */
//  xm = -pt[nn - 1] / pt[n - 1];
//  if (xm < x) {
//      x = xm;
//    }
//  /* CHOP THE INTERVAL (0,X) UNTIL FF .LE. 0 */
//L130:
//  xm = x * 0.1f;
//  ff = pt[0];
//  for (i = 1; i < nn; ++i) {
//      ff = ff * xm + pt[i];
//    }
//  if (ff > 0.0f) {
//      x = xm;
//      goto L130;
//    }
//  dx = x;
//  /* DO NEWTON ITERATION UNTIL X CONVERGES TO TWO */
//  /* DECIMAL PLACES */
//  while (abs(dx/x) > 0.005f) {
//      ff = pt[0];
//      df = ff;
//      for (i = 1; i < n; ++i) {
//          ff = ff * x + pt[i];
//          df = df * x + ff;
//        }
//      ff = ff * x + pt[nn - 1];
//      dx = ff / df;
//      x -= dx;
//    }
//  bnd = x;
//  /* COMPUTE THE DERIVATIVE AS THE INTIAL K POLYNOMIAL */
//  /* AND DO 5 STEPS WITH NO SHIFT */
//  nm1 = n - 1;


//  for (i = 1; i < n; ++i) {
//      k[i] = (nn - i - 1) * p[i] / n;
//    }
//  k[0] = p[0];
//  aa = p[nn - 1];
//  bb = p[n - 1];
//  zerok = k[n - 1] == 0.;
//  for (jj = 1; jj <= 5; ++jj) {
//      cc = k[n - 1];
//      if (zerok) { /* USE UNSCALED FORM OF RECURRENCE */
//          for (i = 0; i < nm1; ++i) {
//              j = nn - i - 2;
//              k[j] = k[j - 1];
//            }
//          k[0] = 0.;
//          zerok = k[n - 1] == 0.;
//        }
//      else { /* USE SCALED FORM OF RECURRENCE */
//          t = -aa / cc;
//          for (i = 0; i < nm1; ++i) {
//              j = nn - i - 2;
//              k[j] = t * k[j - 1] + p[j];
//            }
//          k[0] = p[0];
//          zerok = abs(k[n - 1]) <= abs(bb) * eta * 10.0;
//        }
//    }
//  /* SAVE K FOR RESTARTS WITH NEW SHIFTS */
//  for (i = 0; i < n; ++i) {
//      temp[i] = k[i];
//    }
//  /* LOOP TO SELECT THE QUADRATIC  CORRESPONDING TO EACH */
//  /* NEW SHIFT */
//  for (cnt = 1; cnt <= 20; ++cnt) {
//      /* QUADRATIC CORRESPONDS TO A DOUBLE SHIFT TO A */
//      /* NON-REAL POINT AND ITS COMPLEX CONJUGATE. THE POINT */
//      /* HAS MODULUS BND AND AMPLITUDE ROTATED BY 94 DEGREES */
//      /* FROM THE PREVIOUS SHIFT */
//      xxx = cosr * xx - sinr * yy;
//      yy = sinr * xx + cosr * yy;
//      xx = xxx;
//      sr = bnd * xx;
//      si = bnd * yy;
//      u = sr * -2.;
//      v = bnd;
//      /* SECOND STAGE CALCULATION, FIXED QUADRATIC */
//      i__1 = cnt * 20;
//      fxshfr(&i__1, &nz);
//      if (nz == 0) {
//          goto L260;
//        }
//      /* THE SECOND STAGE JUMPS DIRECTLY TO ONE OF THE THIRD */
//      /* STAGE ITERATIONS AND RETURNS HERE IF SUCCESSFUL. */
//      /* DEFLATE THE POLYNOMIAL, STORE THE ZERO OR ZEROS AND */
//      /* RETURN TO THE MAIN ALGORITHM. */
//      j = degree - n;

//      real[j] = szr;
//      imag[j] = szi;
//      nn -= nz;
//      n = nn - 1;
//      for (i = 0; i < nn; ++i)
//        {
//          p[i] = qp[i];
//        }
//      if (nz == 1)
//        {
//          goto L40;
//        }

//      real[j + 1] = lzr;
//      imag[j + 1] = lzi;

//      goto L40;
//      /* IF THE ITERATION IS UNSUCCESSFUL ANOTHER QUADRATIC */
//      /* IS CHOSEN AFTER RESTORING K */
//L260:
//      for (i = 0; i < n; ++i)
//        {
//          k[i] = temp[i];
//        }
//    }
//  /* RETURN WITH FAILURE IF NO CONVERGENCE WITH 20 */
//  /* SHIFTS */
//  fail = true;
//  degree -= n;
//} /* rpoly_ */


//void Roots::quad(double* a, double* b1, double* c, double* sr, double* si, double* lr, double* li)
//{
//  /* Local variables */
//  static double b, d, e;

//  /* CALCULATE THE ZEROS OF THE QUADRATIC A*Z**2+B1*Z+C. */
//  /* THE QUADRATIC FORMULA, MODIFIED TO AVOID */
//  /* OVERFLOW, IS USED TO FIND THE LARGER ZERO IF THE */
//  /* ZEROS ARE REAL AND BOTH ZEROS ARE COMPLEX. */
//  /* THE SMALLER REAL ZERO IS FOUND DIRECTLY FROM THE */
//  /* PRODUCT OF THE ZEROS C/A. */
//  if (*a != 0.) {
//      goto L20;
//    }
//  *sr = 0.;
//  if (*b1 != 0.) {
//      *sr = -(*c) / *b1;
//    }
//  *lr = 0.;
//L10:
//  *si = 0.;
//  *li = 0.;
//  return;
//L20:
//  if (*c == 0.) {
//      *sr = 0.;
//      *lr = -(*b1) / *a;
//      goto L10;
//    }
//  /* COMPUTE DISCRIMINANT AVOIDING OVERFLOW */
//  b = *b1 / 2.;
//  if (abs(b) >= abs(*c)) {
//      e = 1. - *a / b * (*c / b);
//      d = sqrt(abs(e)) * abs(b);
//    }
//  else {
//      e = *a;
//      if (*c < 0.) {
//          e = -(*a);
//        }
//      e = b * (b / abs(*c)) - e;
//      d = sqrt(abs(e)) * sqrt(abs(*c));
//    }
//  if (e < 0.) {
//      goto L60;
//    }
//  /* REAL ZEROS */
//  if (b >= 0.) {
//      d = -d;
//    }
//  *lr = (-b + d) / *a;
//  *sr = 0.;
//  if (*lr != 0.) {
//      *sr = *c / *lr / *a;
//    }
//  goto L10;
//  /* COMPLEX CONJUGATE ZEROS */
//L60:
//  *sr = -b / *a;
//  *lr = *sr;
//  *si = abs(d / *a);
//  *li = -(*si);
//} /* quad_ */

//double Roots::pow_di(const double *ap, const int *bp)
//{
//  double pow, x;
//  int n;
//  unsigned long u;

//  pow = 1;
//  x = *ap;
//  n = *bp;

//  if(n != 0)
//    {
//      if(n < 0)
//        {
//          n = -n;
//          x = 1/x;
//        }
//      for(u = n; ; )
//        {
//          if(u & 01)
//            pow *= x;
//          if(u >>= 1)
//            x *= x;
//          else
//            break;
//        }
//    }
//  return pow;
//}

//void Roots::fxshfr(int *l2, int *nz)
//{
//  /* Local variables */
//  static int type;
//  static bool stry, vtry;
//  static int i, j, iflag;
//  static double s;
//  static float betas, betav;
//  static bool spass;
//  static bool vpass;
//  static double ui, vi;
//  static float ts, tv, vv;
//  static float ots, otv, tss;
//  static double ss, oss, ovv, svu, svv;
//  static float tvv;

//  /* COMPUTES UP TO  L2  FIXED SHIFT K-POLYNOMIALS, */
//  /* TESTING FOR CONVERGENCE IN THE LINEAR OR QUADRATIC */
//  /* CASE. INITIATES ONE OF THE VARIABLE SHIFT */
//  /* ITERATIONS AND RETURNS WITH THE NUMBER OF ZEROS */
//  /* FOUND. */
//  /* L2 - LIMIT OF FIXED SHIFT STEPS */
//  /* NZ - NUMBER OF ZEROS FOUND */
//  *nz = 0;
//  betav = .25f;
//  betas = .25f;
//  oss = sr;
//  ovv = v;
//  /* EVALUATE POLYNOMIAL BY SYNTHETIC DIVISION */
//  quadsd(&nn, &u, &v, p, qp, &a, &b);
//  calcsc(&type);
//  for (j = 1; j <= *l2; ++j) {
//      /* CALCULATE NEXT K POLYNOMIAL AND ESTIMATE V */
//      nextk(&type);
//      calcsc(&type);
//      newest(&type, &ui, &vi);
//      vv = (float)vi;
//      /* ESTIMATE S */
//      ss = 0.0;
//      if (k[n - 1] != 0.) {
//          ss = -p[nn - 1] / k[n - 1];
//        }
//      tv = 1.0f;
//      ts = 1.0f;
//      if (j == 1 || type == 3) {
//          goto L70;
//        }
//      /* COMPUTE RELATIVE MEASURES OF CONVERGENCE OF S AND V */
//      /* SEQUENCES */
//      if (vv != 0.0f) {
//          tv = (float)abs((vv - ovv) / vv);
//        }
//      if (ss != 0.0) {
//          ts = (float)abs((ss - oss) / ss);
//        }
//      /* IF DECREASING, MULTIPLY TWO MOST RECENT */
//      /* CONVERGENCE MEASURES */
//      tvv = 1.0f;
//      if (tv < otv) {
//          tvv = tv * otv;
//        }
//      tss = 1.0f;
//      if (ts < ots) {
//          tss = ts * ots;
//        }
//      /* COMPARE WITH CONVERGENCE CRITERIA */
//      vpass = tvv < betav;
//      spass = tss < betas;
//      if (! (spass || vpass)) {
//          goto L70;
//        }
//      /* AT LEAST ONE SEQUENCE HAS PASSED THE CONVERGENCE */
//      /* TEST. STORE VARIABLES BEFORE ITERATING */
//      svu = u;
//      svv = v;
//      for (i = 1; i <= n; ++i) {
//          svk[i - 1] = k[i - 1];
//        }
//      s = ss;
//      /* CHOOSE ITERATION ACCORDING TO THE FASTEST */
//      /* CONVERGING SEQUENCE */
//      vtry = false;
//      stry = false;
//      if (spass && (! vpass || tss < tvv)) {
//          goto L40;
//        }
//L20:
//      quadit(&ui, &vi, nz);
//      if (*nz > 0) {
//          return;
//        }
//      /* QUADRATIC ITERATION HAS FAILED. FLAG THAT IT HAS */
//      /* BEEN TRIED AND DECREASE THE CONVERGENCE CRITERION. */
//      vtry = true;
//      betav *= 0.25f;
//      /* TRY LINEAR ITERATION IF IT HAS NOT BEEN TRIED AND */
//      /* THE S SEQUENCE IS CONVERGING */
//      if (stry || ! spass) {
//          goto L50;
//        }
//      for (i = 1; i <= n; ++i) {
//          k[i - 1] = svk[i - 1];
//        }
//L40:
//      realit(&s, nz, &iflag);
//      if (*nz > 0) {
//          return;
//        }
//      /* LINEAR ITERATION HAS FAILED. FLAG THAT IT HAS BEEN */
//      /* TRIED AND DECREASE THE CONVERGENCE CRITERION */
//      stry = true;
//      betas *= 0.25f;
//      /* IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL */
//      /* ZERO ATTEMPT QUADRATIC INTERATION */
//      if (iflag != 0) {
//          ui = -(s + s);
//          vi = s * s;
//          goto L20;
//        }
//      /* RESTORE VARIABLES */
//L50:
//      u = svu;
//      v = svv;
//      for (i = 1; i <= n; ++i) {
//          k[i - 1] = svk[i - 1];
//        }
//      /* TRY QUADRATIC ITERATION IF IT HAS NOT BEEN TRIED */
//      /* AND THE V SEQUENCE IS CONVERGING */
//      if (vpass && ! vtry) {
//          goto L20;
//        }
//      /* RECOMPUTE QP AND SCALAR VALUES TO CONTINUE THE */
//      /* SECOND STAGE */
//      quadsd(&nn, &u, &v, p, qp, &a, &b);
//      calcsc(&type);
//L70:
//      ovv = vv;
//      oss = ss;
//      otv = tv;
//      ots = ts;
//    }
//}

//void Roots::quadsd(int *nn, double *u, double *v, double *p, double *q, double *a, double *b)
//{
//  /* Local variables */
//  static double c;
//  static int i;

//  /* DIVIDES P BY THE QUADRATIC  1,U,V  PLACING THE */
//  /* QUOTIENT IN Q AND THE REMAINDER IN A,B */

//  *b = p[0];
//  q[0] = *b;
//  *a = p[1] - *u * *b;
//  q[1] = *a;
//  for (i = 2; i < *nn; ++i)
//    {
//      c = p[i] - *u * *a - *v * *b;
//      q[i] = c;
//      *b = *a;
//      *a = c;
//    }
//}

//void Roots::calcsc(int *type)
//{
//  /* THIS ROUTINE CALCULATES SCALAR QUANTITIES USED TO */
//  /* COMPUTE THE NEXT K POLYNOMIAL AND NEW ESTIMATES OF */
//  /* THE QUADRATIC COEFFICIENTS. */
//  /* TYPE - INTEGER VARIABLE SET HERE INDICATING HOW THE */
//  /* CALCULATIONS ARE NORMALIZED TO AVOID OVERFLOW */
//  /* SYNTHETIC DIVISION OF K BY THE QUADRATIC 1,U,V */
//  quadsd(&n, &u, &v, k, qk, &c, &d);
//  if (abs(c) > abs(k[n - 1]) * 100.0 * eta)
//    {
//      goto L10;
//    }
//  if (abs(d) > abs(k[n - 2]) * 100.0 * eta)
//    {
//      goto L10;
//    }
//  *type = 3;
//  /* TYPE=3 INDICATES THE QUADRATIC IS ALMOST A FACTOR */
//  /* OF K */
//  return;
//L10:
//  if (abs(d) < abs(c)) {
//      goto L20;
//    }
//  *type = 2;
//  /* TYPE=2 INDICATES THAT ALL FORMULAS ARE DIVIDED BY D */
//  e = a / d;
//  f = c / d;
//  g = u * b;
//  h = v * b;
//  a3 = (a + g) * e + h * (b / d);
//  a1 = b * f - a;
//  a7 = (f + u) * a + h;
//  return;
//L20:
//  *type = 1;
//  /* TYPE=1 INDICATES THAT ALL FORMULAS ARE DIVIDED BY C */
//  e = a / c;
//  f = d / c;
//  g = u * e;
//  h = v * b;
//  a3 = a * e + (h / c + g) * b;
//  a1 = b - a * (d / c);
//  a7 = a + g * d + h * f;
//  return;
//}

//void Roots::newest(int *type, double *uu, double *vv)
//{
//  static double temp, a4, a5, b1, b2, c1, c2, c3, c4;

//  /* COMPUTE NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS */
//  /* USING THE SCALARS COMPUTED IN CALCSC. */
//  /* USE FORMULAS APPROPRIATE TO SETTING OF TYPE. */
//  if (*type == 3) {
//      goto L30;
//    }
//  if (*type == 2) {
//      goto L10;
//    }
//  a4 = a + u * b + h * f;
//  a5 = c + (u + v * f) * d;
//  goto L20;
//L10:
//  a4 = (a + g) * f + h;
//  a5 = (f + u) * c + v * d;
//  /* EVALUATE NEW QUADRATIC COEFFICIENTS. */
//L20:
//  b1 = -k[n - 1] / p[nn - 1];
//  b2 = -(k[n - 2] + b1 * p[n - 1]) / p[nn - 1];
//  c1 = v * b2 * a1;
//  c2 = b1 * a7;
//  c3 = b1 * b1 * a3;
//  c4 = c1 - c2 - c3;
//  temp = a5 + b1 * a4 - c4;
//  if (temp == 0.) {
//      goto L30;
//    }
//  *uu = u - (u * (c3 + c2) + v * (b1 * a1 + b2 * a7)) / temp;
//  *vv = v * (c4 / temp + 1.0);
//  return;
//  /* IF TYPE=3 THE QUADRATIC IS ZEROED */
//L30:
//  *uu = 0.;
//  *vv = 0.;
//  return;
//}

//void Roots::nextk(int *type)
//{
//  /* Local variables */
//  static double temp;
//  static int i;

//  /* COMPUTES THE NEXT K POLYNOMIALS USING SCALARS */
//  /* COMPUTED IN CALCSC */
//  if (*type == 3) {
//      goto L40;
//    }
//  temp = a;
//  if (*type == 1) {
//      temp = b;
//    }
//  if (abs(a1) > abs(temp) * eta * 10.0) {
//      goto L20;
//    }
//  /* IF A1 IS NEARLY ZERO THEN USE A SPECIAL FORM OF THE */
//  /* RECURRENCE */
//  k[0] = 0.;
//  k[1] = -a7 * qp[0];
//  for (i = 3; i <= n; ++i) {
//      k[i - 1] = a3 * qk[i - 3] - a7 * qp[i - 2];
//    }
//  return;
//  /* USE SCALED FORM OF THE RECURRENCE */
//L20:
//  a7 /= a1;
//  a3 /= a1;
//  k[0] = qp[0];
//  k[1] = qp[1] - a7 * qp[0];
//  for (i = 3; i <= n; ++i) {
//      k[i - 1] = a3 * qk[i - 3] - a7 * qp[i - 2] + qp[i - 1];
//    }
//  return;
//  /* USE UNSCALED FORM OF THE RECURRENCE IF TYPE IS 3 */
//L40:
//  k[0] = 0.;
//  k[1] = 0.;
//  for (i = 3; i <= n; ++i) {
//      k[i - 1] = qk[i - 3];
//    }
//}

//void Roots::quadit(double *uu, double *vv, int *nz)
//{
//  /* Local variables */
//  static int type, i, j;
//  static double t;
//  static bool tried;
//  static float ee;
//  static double ui, vi;
//  static float mp, zm;
//  static float relstp, omp;

//  /* VARIABLE-SHIFT K-POLYNOMIAL ITERATION FOR A */
//  /* QUADRATIC FACTOR CONVERGES ONLY IF THE ZEROS ARE */
//  /* EQUIMODULAR OR NEARLY SO. */
//  /* UU,VV - COEFFICIENTS OF STARTING QUADRATIC */
//  /* NZ - NUMBER OF ZERO FOUND */
//  *nz = 0;
//  tried = false;
//  u = *uu;
//  v = *vv;
//  j = 0;
//  /* MAIN LOOP */
//L10:
//  quad(&c_b41, &u, &v, &szr, &szi, &lzr, &lzi);
//  /* RETURN IF ROOTS OF THE QUADRATIC ARE REAL AND NOT */
//  /* CLOSE TO MULTIPLE OR NEARLY EQUAL AND  OF OPPOSITE */
//  /* SIGN */
//  if (abs(abs(szr) - abs(lzr)) > abs(lzr) * .01) {
//      return;
//    }
//  /* EVALUATE POLYNOMIAL BY QUADRATIC SYNTHETIC DIVISION */
//  quadsd(&nn, &u, &v, p, qp, &a, &b);
//  mp = (float)abs(a - szr * b)
//      + (float)abs(szi * b);
//  /* COMPUTE A RIGOROUS  BOUND ON THE ROUNDING ERROR IN */
//  /* EVALUTING P */
//  zm = (float)sqrt(abs(v));
//  ee = (float)abs(qp[0]) * 2.0f;
//  t = -szr * b;
//  for (i = 2; i <= n; ++i) {
//      ee = ee * zm + (float)abs(qp[i - 1]);
//    }
//  ee = ee * zm + (float)abs(a + t);
//  ee = (float)((mre * 5.0 + are * 4.0) * ee
//               - (mre * 5.0 + are * 2.0) * (abs(a + t) + abs(b) * zm)
//               + are * 2.0 * abs(t));
//  /* ITERATION HAS CONVERGED SUFFICIENTLY IF THE */
//  /* POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND */
//  if (mp <= ee * 20.0f) {
//      *nz = 2;
//      return;
//    }
//  /* STOP ITERATION AFTER 20 STEPS */
//  if (++j > 20) {
//      return;
//    }
//  if (j < 2) {
//      goto L50;
//    }
//  if (relstp > 0.01f || mp < omp || tried) {
//      goto L50;
//    }
//  /* A CLUSTER APPEARS TO BE STALLING THE CONVERGENCE. */
//  /* FIVE FIXED SHIFT STEPS ARE TAKEN WITH A U,V CLOSE */
//  /* TO THE CLUSTER */
//  if (relstp < eta) {
//      relstp = eta;
//    }
//  relstp = sqrtf(relstp);
//  u -= u * relstp;
//  v += v * relstp;
//  quadsd(&nn, &u, &v, p, qp, &a, &b);
//  for (i = 1; i <= 5; ++i) {
//      calcsc(&type);
//      nextk(&type);
//    }
//  tried = true;
//  j = 0;
//L50:
//  omp = mp;
//  /* CALCULATE NEXT K POLYNOMIAL AND NEW U AND V */
//  calcsc(&type);
//  nextk(&type);
//  calcsc(&type);
//  newest(&type, &ui, &vi);
//  /* IF VI IS ZERO THE ITERATION IS NOT CONVERGING */
//  if (vi == 0.) {
//      return;
//    }
//  relstp = (float)abs((vi - v) / vi);
//  u = ui;
//  v = vi;
//  goto L10;
//}

//void Roots::realit(double *sss, int *nz, int *iflag)
//{
//  /* Local variables */
//  static int i, j;
//  static double s, t;
//  static float ee, mp, ms;
//  static double kv, pv;
//  static float omp;

//  /* VARIABLE-SHIFT H POLYNOMIAL ITERATION FOR A REAL */
//  /* ZERO. */
//  /* SSS   - STARTING ITERATE */
//  /* NZ    - NUMBER OF ZERO FOUND */
//  /* IFLAG - FLAG TO INDICATE A PAIR OF ZEROS NEAR REAL */
//  /*         AXIS. */
//  *nz = 0;
//  s = *sss;
//  *iflag = 0;
//  j = 0;
//  /* MAIN LOOP */
//L10:
//  pv = p[0];
//  /* EVALUATE P AT S */
//  qp[0] = pv;
//  for (i = 2; i <= nn; ++i) {
//      pv = pv * s + p[i - 1];
//      qp[i - 1] = pv;
//    }
//  mp = (float)abs(pv);
//  /* COMPUTE A RIGOROUS BOUND ON THE ERROR IN EVALUATING */
//  /* P */
//  ms = (float)abs(s);
//  ee = (float)(mre / (are + mre) * abs(qp[0]));
//  for (i = 2; i <= nn; ++i) {
//      ee = ee * ms + (float)abs(qp[i - 1]);
//    }
//  /* ITERATION HAS CONVERGED SUFFICIENTLY IF THE */
//  /* POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND */
//  if (mp <= ((are + mre) * ee - mre * mp) * 20.0f) {
//      *nz = 1;
//      szr = s;
//      szi = 0.;
//      return;
//    }
//  /* STOP ITERATION AFTER 10 STEPS */
//  if (++j > 10) {
//      return;
//    }
//  if (j < 2) {
//      goto L50;
//    }
//  if (abs(t) > abs(s - t) * 0.001 || mp <= omp) {
//      goto L50;
//    }
//  /* A CLUSTER OF ZEROS NEAR THE REAL AXIS HAS BEEN */
//  /* ENCOUNTERED RETURN WITH IFLAG SET TO INITIATE A */
//  /* QUADRATIC ITERATION */
//  *iflag = 1;
//  *sss = s;
//  return;
//  /* RETURN IF THE POLYNOMIAL VALUE HAS INCREASED */
//  /* SIGNIFICANTLY */
//L50:
//  omp = mp;
//  /* COMPUTE T, THE NEXT POLYNOMIAL, AND THE NEW ITERATE */
//  kv = k[0];
//  qk[0] = kv;
//  for (i = 2; i <= n; ++i) {
//      kv = kv * s + k[i - 1];
//      qk[i - 1] = kv;
//    }
//  if (abs(kv) <= abs(k[n - 1]) * 10.0 * eta) {
//      goto L80;
//    }
//  /* USE THE SCALED FORM OF THE RECURRENCE IF THE VALUE */
//  /* OF K AT S IS NONZERO */
//  t = -pv / kv;
//  k[0] = qp[0];
//  for (i = 2; i <= n; ++i) {
//      k[i - 1] = t * qk[i - 2] + qp[i - 1];
//    }
//  goto L100;
//  /* USE UNSCALED FORM */
//L80:
//  k[0] = 0.;
//  for (i = 2; i <= n; ++i) {
//      k[i - 1] = qk[i - 2];
//    }
//L100:
//  kv = k[0];
//  for (i = 2; i <= n; ++i) {
//      kv = kv * s + k[i - 1];
//    }
//  t = 0.;
//  if (abs(kv) > abs(k[n - 1]) * 10.0 * eta) {
//      t = -pv / kv;
//    }
//  s += t;
//  goto L10;
//}
