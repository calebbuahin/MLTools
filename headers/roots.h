#ifndef ROOTS_H
#define ROOTS_H

#include <vector>
#include <complex>

class Roots
{
public:
  static std::complex<double>* roots(const std::vector<double>& coefficients, int & size);
};

#endif // ROOTS_H
