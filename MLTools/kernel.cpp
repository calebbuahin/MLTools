#include <stdafx.h>
#include "mrvm.h"


Kernel::Kernel(KernelType type , double lengthScale)
  :m_useBias(true), m_polynomailPower(2)
{
  this->m_kernelType = type;
  this->m_lengthScale = lengthScale;
}

Kernel::~Kernel()
{

}

Kernel::KernelType Kernel::type() const
{
  return this->m_kernelType;
}

void Kernel::setKernelType(KernelType type)
{
  this->m_kernelType = type;
}

double Kernel::lengthScale() const
{
  return this->m_lengthScale;
}

void Kernel::setLengthScale(double lengthScale)
{
  this->m_lengthScale = lengthScale;
}

double Kernel::polynomialPower() const
{
  return this->m_polynomailPower;
}

void Kernel::setPolynomialPower(double polynomialPower)
{
  this->m_polynomailPower = polynomialPower;
}


bool Kernel::useBias() const
{
  return this->m_useBias;
}

void Kernel::setUseBias(bool m_useBias)
{
  this->m_useBias = m_useBias;
}

array Kernel::calculateKernel(const array& x1, const array& x2)
{


  switch (m_kernelType)
    {
    case Gaussian:
      return calculateGaussianKernel(x1,x2);
      break;
    case Laplace:
      return calculateLaplaceKernel(x1,x2);
      break;
    case Polynomial:
      return calculatePolynomialKernel(x1,x2);
      break;
    case Spline:
      return calculateSplineKernel(x1,x2);
      break;
    case Cauchy:
      return calculateCauchyKernel(x1,x2);
      break;
    case Cubic:
      return calculateCubicKernel(x1,x2);
      break;
    case Distance:
      return calculateDistanceKernel(x1,x2);
      break;
    case ThinPlateSpline:
      return calculateThinPlateSplineKernel(x1,x2);
      break;
    case Bubble:
      return calculateBubbleKernel(x1,x2);
      break;
    }
}

array Kernel::calculateGaussianKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);
  int colX1 = x1.dims(1);

  int rowX2 = x2.dims(0);
  int colX2 = x2.dims(1);

  double eta = 1/(m_lengthScale*m_lengthScale);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), af::exp( -eta * distanceSquared(x1,x2)/2.0));
    }
  else
    {
      return af::exp( -eta * distanceSquared(x1,x2)/2.0);
    }

}

array Kernel::calculateLaplaceKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);
  
  double eta = 1/(m_lengthScale*m_lengthScale);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), af::exp( -1* af::sqrt(eta * distanceSquared(x1,x2))));
    }
  else
    {
      return af::exp( -1* af::sqrt(eta * distanceSquared(x1,x2)));
    }
}

array Kernel::calculatePolynomialKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1/(m_lengthScale*m_lengthScale);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), af::pow(af::matmul( x1, (eta * x2).T()) + 1, m_polynomailPower));
    }
  else
    {
      return af::pow(af::matmul( x1, (eta * x2).T()) + 1, m_polynomailPower);
    }
}

array Kernel::calculateSplineKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);
  int colX1 = x1.dims(1);

  int rowX2 = x2.dims(0);
  int colX2 = x2.dims(1);

  double eta = 1/(m_lengthScale*m_lengthScale);

  array r2 = eta*distanceSquared(x1,x2);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), 0.5*r2*log(r2+(r2==0)));
    }
  else
    {
      return af::constant(1,rowX1,1), 0.5*r2*log(r2+(r2==0));
    }
}

array Kernel::calculateCauchyKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1/(m_lengthScale*m_lengthScale);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), af::pow(af::matmul( x1, (eta * x2).T()) + 1, m_polynomailPower));
    }
  else
    {
      return af::pow(af::matmul( x1, (eta * x2).T()) + 1, m_polynomailPower);
    }
}

array Kernel::calculateCubicKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1/(m_lengthScale*m_lengthScale);

  array r2 = eta*distanceSquared(x1,x2);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), r2*sqrt(r2));
    }
  else
    {
      return r2*sqrt(r2);
    }
}

array Kernel::calculateDistanceKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1/(m_lengthScale*m_lengthScale);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), sqrt(eta)*sqrt(distanceSquared(x1,x2)));
    }
  else
    {
      return sqrt(eta)*sqrt(distanceSquared(x1,x2));
    }
}

array Kernel::calculateThinPlateSplineKernel(const array& x1, const array& x2)
{

  int rowX1 = x1.dims(0);
  int colX1 = x1.dims(1);

  int rowX2 = x2.dims(0);
  int colX2 = x2.dims(1);

  double eta = 1/(m_lengthScale*m_lengthScale);

  array r2 = eta*distanceSquared(x1,x2);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), 0.5*r2*log(r2+(r2==0)));
    }
  else
    {
      return af::constant(1,rowX1,1), 0.5*r2*log(r2+(r2==0));
    }
}

array Kernel::calculateBubbleKernel(const array& x1, const array& x2)
{
  int rowX1 = x1.dims(0);
  int colX1 = x1.dims(1);

  int rowX2 = x2.dims(0);
  int colX2 = x2.dims(1);

  double eta = 1/(m_lengthScale*m_lengthScale);

  array r2 = eta*distanceSquared(x1,x2);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), r2 < 1);
    }
  else
    {
      return r2 < 1;
    }
}

array Kernel::distanceSquared(const array& x , const array& y)
{
  int nx = x.dims(0);
  int ny = y.dims(0);

  array values = matmul(sum(pow(x,2),1), constant(1,1,ny)) +
                 matmul(constant(1,nx,1), sum(pow(y,2),1).T()) -
                 2*(af::matmul(x,y.T()));

  af_print(values);

  return values;
}
