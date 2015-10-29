#ifndef ROOTS_H
#define ROOTS_H

#include <vector>
#include <complex>

class Roots
{
public:

//  static void roots(const std::vector<double>& op, int& degree, double* real, double* imag, bool& fail);

  static std::complex<double>* roots(const std::vector<double>& coefficients, int & size);

//  static void quad(double* a,double* b1,double* c, double* sr,double* si,double* lr,double* li);

//  static void quadsd(int* nn, double* u, double* v, double* p, double* q, double* a, double* b);

//  static double pow_di(const double* ap, const int* bp);

//  static void fxshfr(int* l2,int* nz);

//  static void calcsc(int* type);

//  static void newest(int *type, double *uu, double *vv);

//  static void nextk(int *type);

//  static void quadit(double *uu, double *vv, int *nz);

//  static void realit(double *sss, int *nz, int *iflag);

//private:
//  static double p[], qp[], k[], qk[], svk[], sr, si, u, v, a, b,
//  c, d, a1, a2, a3, a6, a7, e, f, g, h, szr, szi, lzr, lzi;
//  static float eta, are, mre;
//  static int n, nn;
//  static double c_b41;
};

#endif // ROOTS_H
