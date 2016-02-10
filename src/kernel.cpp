#include <stdafx.h>
#include "mrvm.h"

using namespace  af;

Kernel::Kernel(KernelType type , double lengthScale)
  : m_polynomailPower(2),m_useBias(true)
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

af::array Kernel::calculateKernel(const af::array& x1, const af::array& x2)
{
  switch (m_kernelType)
    {
    case Gaussian:
      return calculateGaussianKernel(x1,x2);
      break;
    case Laplace:
      return calculateLaplaceKernel(x1,x2);
      break;
    case HomogeneousPolynomail:
      return calculateHomogeneousPolynomialKernel(x1,x2);
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

af::array Kernel::calculateGaussianKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  float eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array temp = af::exp(-eta * distanceSquared(x1,x2) / 2.0);

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;
}

af::array Kernel::calculateLaplaceKernel(const af::array& x1, const af::array& x2)
{

  int rowX1 = x1.dims(0);
  
  double eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array temp =  af::exp( -1 * af::sqrt(eta * distanceSquared(x1,x2)));

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;
}

af::array Kernel::calculatePolynomialKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array temp =  af::pow(af::matmul(x1,(eta * x2).T()) + 1.0, m_polynomailPower);

  af_print(af::sum( af::sum( af::isNaN(temp) )),1);

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;
}

af::array Kernel::calculateHomogeneousPolynomialKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1/(m_lengthScale*m_lengthScale);

  af::array temp = af::pow(eta * af::matmul(x1, x2.T()) , m_polynomailPower);

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;
}

af::array Kernel::calculateSplineKernel(const af::array& x1, const af::array& x2)
{
  
  af::array x1_ = x1 / m_lengthScale;
  af::array x2_ = x2 / m_lengthScale;

  int rowX1 = x1.dims(0);
  int colX1 = x1.dims(1);

  int rowX2 = x2.dims(0);
//  int colX2 = x2.dims(1);

  af::array k = af::constant(1.0,rowX1,rowX1);

  for(int i = 0 ; i < colX1 ; i++)
    {
      af::array xx = matmul( x1_(span,i) , x2_(span,i).T());
      af::array xx1 = matmul(x1_(span,i) , af::constant(1.0 , 1 , rowX2));
      af::array xx2 = matmul(af::constant(1.0, rowX1,1), x2_(span,i).T());
      af::array minXX = af::min(xx1,xx2);
      af::array temp = 1.0 + xx + xx * minXX - (xx1 +  xx2)/2.0 * af::pow(minXX,2.0) + af::pow(minXX, 3.0)/3.0;
      k = k * temp;
    }

  if(m_useBias)
    {
      return af::join(1,af::constant(1,rowX1,1), k);
    }
  else
    {
      return k;
    }
}

af::array Kernel::calculateCauchyKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array r2 = eta*distanceSquared(x1,x2);

  af::array temp = 1.0 / (1.0 + r2*sqrt(r2));

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;

}

af::array Kernel::calculateCubicKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array r2 = eta*distanceSquared(x1,x2);

  af::array temp = r2*sqrt(r2);

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;

}

af::array Kernel::calculateDistanceKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array temp = sqrt(eta)*sqrt(distanceSquared(x1,x2));

  return m_useBias ? af::join(1, af::constant(1.0,rowX1,1),temp) : temp;

}

af::array Kernel::calculateThinPlateSplineKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1.0/(m_lengthScale*m_lengthScale);

  af::array r2 = eta*distanceSquared(x1,x2);

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), 0.5 * r2 * af::log(r2 + (r2 == 0)));
    }
  else
    {
      return 0.5 * r2 * af::log(r2 + (r2 == 0));
    }
}

af::array Kernel::calculateBubbleKernel(const af::array& x1, const af::array& x2)
{
  int rowX1 = x1.dims(0);

  double eta = 1/(m_lengthScale*m_lengthScale);

  af::array r2 = eta*distanceSquared(x1,x2) ;

  if(m_useBias)
    {
      return af::join(1, af::constant(1,rowX1,1), r2 < 1);
    }
  else
    {
      return r2 < 1;
    }
}

af::array Kernel::distanceSquared(const af::array& x , const af::array& y)
{
  int nx = x.dims(0);
  int ny = y.dims(0);

  af::array values = af::matmul(af::sum(af::pow(x,2.0),1), af::constant(1.0,1,ny)) +
      af::matmul(af::constant(1.0,nx,1), af::sum(af::pow(y,2.0),1).T()) -
      2.0 * (af::matmul(x,y.T()));

  return  values;
}
