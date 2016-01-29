#include "stdafx.h"
#include "mrvm.h"
#include "QDebug"
#include <iostream>
#include "roots.h"
#include <QTime>
#include <omp.h>

using namespace std;
using namespace  af;

bool MRVM::s_gdalRegistered = false;


MRVM::MRVM(const QFileInfo& file)
  :QObject(), m_mode(TrainingAndRegression) , m_converged(false) , m_verbose(true)
{
  if(!s_gdalRegistered)
    {

      GDALAllRegister();
      s_gdalRegistered = true;
    }
  
  m_file = file;
  readProject();
}

MRVM::~MRVM()
{
  qDeleteAll(m_inputItems);
  m_inputItems.clear();
  
  qDeleteAll(m_outputItems);
  m_outputItems.clear();
}

QString MRVM::name() const
{
  return m_name;
}

int MRVM::inputDimension() const
{
  int id = 0;
  
  for(int i = 0 ; i < m_inputItems.length() ; i++)
    {
      id = id + m_inputItems[i]->columnCount();
    }
  
  return id;
}

int MRVM::targetDimension() const
{
  int id = 0;
  
  for(int i = 0 ; i < m_outputItems.length() ; i++)
    {
      id = id + m_outputItems[i]->columnCount();
    }
  
  return id;
}

int MRVM::maxNumberOfIterations() const
{
  return m_maxNumberOfIterations;
}

void MRVM::setMaxNumberOfIterations(int niters)
{
  if(niters> 0)
    {
      this->m_maxNumberOfIterations = niters;
    }
}

bool MRVM::verbose() const
{
  return m_verbose;
}

void MRVM::setVerbose(bool verbose)
{
  this->m_verbose;
}

MRVM::Mode MRVM::mode() const
{
  return m_mode;
}

bool MRVM::converged() const
{
  return m_converged;
}

const Kernel& MRVM::kernel() const
{
  return m_kernel;
}

const QList<MRVMItem*>& MRVM::inputItems() const
{
  return this->m_inputItems;
}

void MRVM::addInputItem(MRVMItem* const inputItem)
{
  if(!m_inputItems.contains(inputItem))
    m_inputItems.append(inputItem);
}

bool MRVM::removeInputItem(MRVMItem* const inputItem)
{
  return m_inputItems.removeOne(inputItem);
}

const QList<MRVMItem*>& MRVM::outputItems() const
{
  return m_outputItems;
}

void MRVM::addOutputItem(MRVMItem* const outputItem)
{
  if(!m_outputItems.contains(outputItem))
    m_outputItems.append(outputItem);
}

bool MRVM::removeOutputItem(MRVMItem* const outputItem)
{
  return m_outputItems.removeOne(outputItem);
}

QString MRVM::matrixOutputFile() const
{
  return this->m_matrixOutputFile;
}

void MRVM::setMatrixOutputFile(const QString& matrixOutputFile)
{
  this->m_matrixOutputFile = matrixOutputFile;
}

const af::array& MRVM::usedRelevantVectors() const
{
  return m_used;
}

const af::array& MRVM::alpha() const
{
  return m_alpha;
}

const af::array& MRVM::invSigma() const
{
  return m_invSigma;
}

const af::array& MRVM::omega() const
{
  return m_omega;
}

const af::array& MRVM::mu() const
{
  return m_Mu;
}

void MRVM::saveProject()
{
  QFile file(m_file.absoluteFilePath());
  
  if(file.open(QIODevice::WriteOnly))
    {
      QXmlStreamWriter writer(&file);
      writer.setAutoFormatting(true);
      
      writer.writeStartDocument();
      writer.writeStartElement("MRVM");
      writer.writeAttribute("name", m_name);

      switch(m_mode)
        {
        case Training:
          writer.writeAttribute("mode", "Training");
          break;
        case Regression:
          writer.writeAttribute("mode", "Regression");
          break;
        case TrainingAndRegression:
          writer.writeAttribute("mode", "TrainingAndRegression");
          break;
        }
      
      writer.writeTextElement("Tolerance", QString::number(m_tolerance));
      
      writer.writeTextElement("MaxNumberOfIterations", QString::number(m_maxNumberOfIterations));
      
      writer.writeTextElement("NumberOfIterations", QString::number(m_numberOfIterations));
      
      writer.writeTextElement("Verbose", m_verbose ? "True" : "False");
      
      writer.writeTextElement("Converged", m_converged ? "True" : "False");

      writer.writeTextElement("MaxAlphaChange", QString::number(m_maxChangeAlpha));

      writer.writeTextElement("MinAlphaChange", QString::number(m_minChangeAlpha));

      af::saveArray("InputMatrix",m_inputMatrix, m_matrixOutputFile.toStdString().c_str());

      af::saveArray("Used",m_used, m_matrixOutputFile.toStdString().c_str(),true);

      af::saveArray("Alpha",m_alpha, m_matrixOutputFile.toStdString().c_str(),true);
      \
      af::saveArray("InverseSigma",m_invSigma, m_matrixOutputFile.toStdString().c_str(),true);

      af::saveArray("Omega",m_omega, m_matrixOutputFile.toStdString().c_str(),true);

      af::saveArray("Mu",m_Mu, m_matrixOutputFile.toStdString().c_str(),true);

      writer.writeTextElement("MatrixOutputFile", m_matrixOutputFile);
      
      writer.writeStartElement("Kernel");
      
      Kernel::KernelType kernelType = m_kernel.type();
      
      switch(kernelType)
        {
        case Kernel::Gaussian:
          writer.writeAttribute("KernelType", "Gaussian");
          break;
        case Kernel::Laplace:
          writer.writeAttribute("KernelType", "Laplace");
          break;
        case Kernel::Polynomial:
          writer.writeAttribute("KernelType", "Polynomial");
          break;
        case Kernel::HomogeneousPolynomail:
          writer.writeAttribute("KernelType", "HomogeneousPolynomail");
          break;
        case Kernel::Spline:
          writer.writeAttribute("KernelType", "Spline");
          break;
        case Kernel::Cauchy:
          writer.writeAttribute("KernelType", "Cauchy");
          break;
        case Kernel::Cubic:
          writer.writeAttribute("KernelType", "Cubic");
          break;
        case Kernel::Distance:
          writer.writeAttribute("KernelType", "Distance");
          break;
        case Kernel::ThinPlateSpline:
          writer.writeAttribute("KernelType", "ThinPlateSpline");
          break;
        case Kernel::Bubble:
          writer.writeAttribute("KernelType", "Bubble");
          break;
        }
      
      writer.writeTextElement("LengthScale", QString::number(m_kernel.lengthScale()));
      
      if(m_kernel.useBias())
        {
          writer.writeTextElement("UseBias", "True");
        }
      else
        {
          writer.writeTextElement("UseBias", "False");
        }
      
      writer.writeTextElement("PolynomialPower",QString::number(m_kernel.polynomialPower()));
      
      writer.writeEndElement();
      
      if(m_inputItems.length() > 0)
        {
          writer.writeStartElement("InputItems");
          
          for(int i = 0 ; i < m_inputItems.length() ; i++)
            {
              MRVMItem* item = m_inputItems[i];
              item->writeXML(writer);
            }
          

          int usedCount = m_used.dims(0);

          if(usedCount)
            {
              //  af_print(m_used);

              af_dtype typess = m_used.type();

              std::cout << typess << endl;

              char* usedData = m_used.host<char>();

              writer.writeStartElement("RelevantVectors");

              for(int i = 0 ; i < usedCount ; i++ )
                {
                  if(usedData[i])
                    writer.writeTextElement("Vector",QString::number(i+1));
                }

              delete[] usedData;

              writer.writeEndElement();
            }


          writer.writeEndElement();
        }


      if(m_outputItems.length() > 0)
        {
          writer.writeStartElement("OutputItems");
          
          for(int i = 0 ; i < m_outputItems.length() ; i++)
            {
              MRVMItem* item = m_outputItems[i];
              item->writeXML(writer);
            }
          
          writer.writeEndElement();
        }
      
      
      writer.writeEndElement();
      writer.writeEndDocument();
      
      file.close();
    }
}

void MRVM::start()
{
  switch (m_mode)
    {
    case Training:
      performTraining();
      break;
    case Regression:
      performRegression();
      break;
    default:
      {
        performTraining();
        performRegression();
      }
      break;
    }
  
  saveProject();
  
}

void MRVM::performTraining()
{
  
  switch (algmode)
    {
    case 0:
      mrvm();
      break;
    default:
      fmrvm();
      break;
    }
}

void MRVM::mrvm()
{

  QTime timer;

  timer.start();

  int row, column;
  
  float* values = getInputMatrix(row,column);
  m_inputMatrix = af::array(row,column, values);
  delete[] values;
  
  //  if(m_verbose)
  //    af_print(m_inputMatrix);
  
  values = getTargetMatrix(row,column);
  m_targetMatrix = af::array(row,column,values);
  delete[] values;
  
  //  if(m_verbose)
  //    af_print(m_targetMatrix);
  
  
  
  af::array phi = m_kernel.calculateKernel(m_inputMatrix,m_inputMatrix);
  
  //  if(m_verbose)
  //    {
  //      af_print(phi);
  //    }
  
  N = m_targetMatrix.dims(0);
  V = m_targetMatrix.dims(1);
  
  assert(phi.dims(0) == N && phi.dims(1) == N+1);
  
  af::array beta, mask, phiSigmaPhi;

#pragma omp parallel num_threads(3)
  {
    m_alpha = af::constant(af::Inf, N+1,1);
    beta = 1.0 / (0.1 * af::var(m_targetMatrix));
    mask = m_alpha < af::NaN;
    phiSigmaPhi = af::constant(0.0,N,N,V);
  }
  
  if(m_verbose)
    af_print(beta);
  
  
  af::array nmask;
  af::array Sigma;
  af::array PhiMask;

  m_converged = false;
  

  m_minChangeAlpha = std::numeric_limits<float>::max();
  m_maxChangeAlpha = std::numeric_limits<float>::min();

  for(m_numberOfIterations = 0 ; m_numberOfIterations < m_maxNumberOfIterations; m_numberOfIterations++)
    {
      af::array sPrime = af::constant(af::NaN, N+1,V);
      af::array qPrime = af::constant(af::NaN, N+1,V);
      
      af::array s = af::constant(af::NaN, N+1,V);
      af::array q = af::constant(af::NaN, N+1,V);
      
      af::array allAlphaInf = af::alltrue(mask == 0);
      
#pragma omp parallel num_threads(8)
      {
#pragma omp for
        {
          for(int i = 0 ; i <  N+1; i++)
            {
              af::array phiSq = af::matmul(phi(span,i).T(),phi(span,i));

              if(af::alltrue<bool>(allAlphaInf > 0))
                {
                  for(int j = 0 ; j < V ; j++)
                    {
                      sPrime(i,j) = beta(j) * phiSq;
                      qPrime(i,j) = beta(j) * matmul(phi(span,i).T(),m_targetMatrix(span,j));
                    }
                }
              else
                {
                  for(int j = 0 ; j < V ; j++)
                    {
                      af::array power = af::pow(beta(j),2.0);
                      float* vvv = power.host<float>();
                      af::array temp  = vvv[0] * matmul(phi(span,i).T(),phiSigmaPhi(span,span,j)) ;
                      delete[] vvv;
                      af::array mop = beta(j) * phiSq;
                      af::array mop1 = matmul(temp,phi(span,i));
                      sPrime(i,j) = mop - mop1;
                      qPrime(i,j) = beta(j) * matmul(phi(span,i).T(),m_targetMatrix(span,j)) - matmul(temp,m_targetMatrix(span,j));
                    }
                }

              if (alltrue<bool>(mask(i) == 0))
                {
                  s(i,span) = sPrime(i,span);
                  q(i,span) = qPrime(i,span);
                }
              else
                {
                  float* tempf = m_alpha(i).host<float>();
                  af::array temp = tempf[0]/(tempf[0]-sPrime(i,span));
                  s(i,span) =  temp * sPrime(i,span);
                  q(i,span) =  temp * qPrime(i,span);
                  delete[] tempf;
                }
            }
        }
      }
      
      af::array deltaL = af::constant(af::NaN, N+1,1);
      af::array task = af::constant(0,N+1,1);
      af::array alphaNew = af::constant(af::NaN, N+1,1);
      
      //      if(m_verbose)
      //        {
      //          af_print(s);
      //          af_print(q);
      //        }
      
      for(int i = 0 ; i < N+1 ; i++)
        {
          af::array polyCoeff = af::constant(0, V+1 , 2*V+1 );
          af::array temp(1,1); temp = V * 1.0f;
          af::array te(1,3);
          
          for(int j = 0;j < V; j++)
            {
              te(0,0) = 1.0f; te(0,1) = 2.0*s(i,j); te(0,2) = af::pow(s(i,j),2.0);
              af::array ff = af::convolve(temp,te,AF_CONV_EXPAND);
              temp = ff;
            }
          
          polyCoeff(0,span) = temp;
          
          for(int j = 0 ; j < V ; j++)
            {
              temp = af::array(1,3);
              temp(0,0) = -1; temp(0,1) = -1 * s(i,j) - af::pow(q(i,j),2); temp(0,2) = 0;
              
              for(int k = 0; k < V; k++)
                {
                  if(k != j)
                    {
                      af::array temp1(1,3); temp1(0,0) = 1; temp1(0,1) = 2*s(i,k); temp1(0,2) = pow(s(i,k),2);
                      af::array ff = convolve(temp,temp1,AF_CONV_EXPAND);
                      temp = ff;
                    }
                }
              polyCoeff(j+1,span) = temp;
            }
          
          af::array sumPolyCoeff = af::sum(polyCoeff);

          if(alltrue<bool>(sumPolyCoeff(0) != 0))
            {
              ASSERT(false, "The order of polynomial should be 2V-1");
            }
          
          int length  = sumPolyCoeff.dims(1);
          
          float* values = sumPolyCoeff.host<float>();
          
          vector<double> sumPolyCoeffs;
          
          bool fnzf = false;
          int k = 0;
          
          for(int ii = 0 ; ii < length ; ii++)
            {
              if(!fnzf)
                {
                  if(fabs(values[ii]) > 10e-16)
                    {
                      fnzf = true;
                      sumPolyCoeffs.push_back(values[ii]);
                      k++;
                    }
                }
              else
                {
                  sumPolyCoeffs.push_back(values[ii]);
                  k++;
                }
            }
          
          delete[] values;
          
          if(k == 0)
            {
              alphaNew(i) = af::Inf;
              QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i), alphaNew(i),sPrime(i,span), qPrime(i,span),s(i,span),q(i,span));
              task(i,span) = t.first;
              deltaL(i) = t.second;
            }
          else
            {
              int size = 0;
              complex<double>* croots = Roots::roots(sumPolyCoeffs, size);
              vector<float> roots;
              af::array posRealRoots;
              
              for(int ff = 0 ; ff < size ; ff++)
                {
                  if(fabs(croots[ff].imag()) < 1e-20 && croots[ff].real() > 0)
                    {
                      roots.push_back(croots[ff].real());
                    }
                }
              
              delete[] croots;
              
              if(roots.size())
                {
                  posRealRoots = af::array(1,roots.size() , & roots[0]);
                }
              

              if(roots.size() == 0)
                {
                  alphaNew(i) = af::Inf;
                  QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i),alphaNew(i),sPrime(i,span),qPrime(i,span),s(i,span),q(i,span));
                  task(i,span) = t.first;
                  deltaL(i) = t.second;
                }
              else if(roots.size() == 1)
                {
                  alphaNew(i) = posRealRoots;
                  QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i),alphaNew(i),sPrime(i,span),qPrime(i,span),s(i,span),q(i,span));
                  task(i,span) = t.first;
                  deltaL(i) = t.second;
                }
              else
                {
                  
                  af::array taskRoots = af::constant(0, roots.size() , 1);
                  af::array deltaRoots = af::constant(af::NaN , roots.size(),1);
                  
                  for(int k = 0; k < roots.size() ; k++)
                    {
                      QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i),posRealRoots(k),sPrime(i,span),qPrime(i,span),s(i,span),q(i,span));
                      taskRoots(k,span) = t.first;
                      deltaRoots(k) = t.second;
                    }
                  
                  af::array vals, indexes;

                  af::max(vals,indexes,deltaRoots);
                  deltaL(i) = vals;
                  alphaNew(i) = posRealRoots(indexes);
                  task(i,span) = taskRoots(indexes,span);
                }
            }
        }
      

      af::array i(1,1);
      
      if(alltrue<bool>(af::isNaN(deltaL)))
        {
          cout << "No further tasks ending iteration" << endl;
          break;
        }
      else if (alltrue<bool>((deltaL * (deltaL != af::NaN))  == -af::Inf))
        {
          cout << "All values of deltaL are -infinity" << endl;
          i = rand() % (N+1) + 1;
        }
      else
        {
          af::array vals(1,1);
          
          vals = -1*af::Inf;

          for(int mm = 0 ; mm < deltaL.dims(0); mm++)
            {
              if(alltrue<bool>(deltaL(mm,0) > vals))
                {
                  vals = deltaL(mm,0);
                  i = mm;
                }
            }

        }

      if(m_numberOfIterations != 0)
        {
#pragma omp parallel num_threads(8)
          {
#pragma omp for
            {
              for(int j = 0 ; j < V; j++)
                {

                  af::array temp = Sigma(span,span,j);
                  af::array diagm;

                  if(temp.dims(0) == 1 && temp.dims(1) == 1)
                    diagm = temp;
                  else
                    diagm = af::diag(Sigma(span,span,j));

                  beta(j) = (N - nmask + matmul(m_alpha(mask).T(), diagm))/af::sum(af::pow(m_targetMatrix(span,j) - matmul(PhiMask, m_Mu(span,j)),2 ));

                }
            }
          }
        }
      
      af::array check = task(i,span);
      
      if(alltrue<bool>(check == 1))
        {
          af::array changeLogAlpha = log(m_alpha(i)/alphaNew(i));

          float* changeLogAlphaf = changeLogAlpha.host<float>();

          m_alpha(i) = alphaNew(i);
          
          if(m_verbose)
            af_print(changeLogAlpha);

          if(fabs(changeLogAlphaf[0]) < m_minChangeAlpha)
            m_minChangeAlpha = fabs(changeLogAlphaf[0]);

          if(fabs(changeLogAlphaf[0]) > m_maxChangeAlpha)
            m_maxChangeAlpha = fabs(changeLogAlphaf[0]) ;


          if(fabs(changeLogAlphaf[0]) < m_tolerance)
            {
              af::array notmask = (mask == 0);

              af::array indexes = af::where(notmask > 0);

              af::array alphaNewEx = alphaNew(indexes);

              if(alltrue<bool>(af::isInf(alphaNewEx)))
                {
                  m_converged = true;
                }
            }

          delete[] changeLogAlphaf;

        }
      else if(alltrue<bool>(check == 2))
        {
          m_alpha(i) = alphaNew(i);
          mask(i) = 1;
          
        }
      else if(alltrue<bool>(check == 3))
        {
          m_alpha(i) = af::Inf;
          mask(i) = 0;
        }
      else
        {
          ASSERT (false,"All tasks are non");
        }
      
      PhiMask = phi(span, mask);
      nmask = af::sum(mask);
      
      uint* vals = nmask.host<uint>();
      int size = vals[0];
      m_invSigma = af::constant(0,size,size,V);
      
      
      Sigma = af::constant(0,size,size,V);
      
      m_Mu = af::constant(0, size,V);
      delete[] vals;
      
      phiSigmaPhi = af::constant(0,N,N,V);
      
      af::array PhiMaskPhiMask = matmul(PhiMask.T(), PhiMask);
      
#pragma omp parallel num_threads(8)
      {
#pragma omp for
        {
          for(int j = 0 ; j < V ; j++)
            {
              af::array almask = m_alpha(mask);

              if(almask.dims(0)  > 1)
                {
                  af::array ddiag = af::diag(almask,0,false);
                  almask = ddiag;
                }

              float* betav = beta(j).host<float>();

              m_invSigma(span,span,j) = betav[0] * PhiMaskPhiMask + almask;

              af::array upper;
              af::cholesky(upper, m_invSigma(span,span,j));

              af::array inVU = af::inverse(upper);

              Sigma(span,span,j) = matmul(inVU,inVU.T());
              m_Mu(span,j) = betav[0]* (matmul(Sigma(span,span,j),matmul(PhiMask.T(), m_targetMatrix(span,j))));
              delete[] betav;

              phiSigmaPhi(span,span,j) = matmul(PhiMask,matmul(Sigma(span,span,j),PhiMask.T()));

            }
        }
      }
      

      cout << "Iteration Number: " <<  (m_numberOfIterations + 1)  << "  Time elapsed: " <<  timer.elapsed() * 1.0/1000.0 << " s"<< endl;

      if(m_converged)
        {
          cout << "Model Converged !" << endl;
          break;
        }
    }

  cout << "Total Number of Iterations : " << m_numberOfIterations + 1 << endl;
  cout << "Training time : " << timer.elapsed() * 1.0/1000.0 << "s" << endl;
  

  m_used = mask;
  
  af::array targetPredict = matmul(phi(span, m_used),m_Mu);
  
  af::array diff = m_targetMatrix - targetPredict;
  
  //  if(m_verbose)
  //    {
  //    //  af_print(targetPredict);
  //    //  af_print(diff);
  //    }
  
  af::array omegaTilde = matmul(diff.T() , diff) / ((N-1)*1.0);
  
  af::array rhat = corrcov(omegaTilde);
  
  af::array dhat = af::constant(0,V,V);
  
  af::array t = 1.0 / af::sqrt(beta);
  
  
  int z = t.dims(0);
  
  gfor(seq f, z)
  {
    dhat(f,f) = t(f);
  }
  
  m_omega = matmul(matmul(dhat,rhat),dhat);
  
  //if(m_verbose)
  //  af_print(m_omega);
  
}

void MRVM::fmrvm()
{
  int row, column;

  float* values = getInputMatrix(row,column);
  m_inputMatrix = af::array(row,column, values);
  delete[] values;

  if(m_verbose)
    af_print(m_inputMatrix);

  values = getTargetMatrix(row,column);
  m_targetMatrix = af::array(row,column,values);
  delete[] values;

  if(m_verbose)
    af_print(m_targetMatrix);

  af::array phi = m_kernel.calculateKernel(m_inputMatrix,m_inputMatrix);

  if(m_verbose)
    af_print(phi);

  N = m_targetMatrix.dims(0);
  V = m_targetMatrix.dims(1);

  assert(phi.dims(0) == N && phi.dims(1) == N+1);

  m_alpha = af::constant(af::Inf, N+1,1);
  af::array beta = 1.0 / (0.1 * af::var(m_targetMatrix));

  if(m_verbose)
    af_print(beta);

  af::array mask = m_alpha < af::NaN;

  af::array nmask;
  af::array Sigma;
  af::array PhiMask;


  m_converged = false;

  af::array phiSigmaPhi = af::constant(0.0,N,N,V);

  m_minChangeAlpha = std::numeric_limits<float>::max();
  m_maxChangeAlpha = std::numeric_limits<float>::min();

  for(m_numberOfIterations = 0 ; m_numberOfIterations < m_maxNumberOfIterations; m_numberOfIterations++)
    {

    }

}

void MRVM::performRegression()
{

  af::array diagf;

  if(m_omega.dims(0) > 1)
    {
      diagf = af::diag(m_omega).T();
    }
  else
    {
      diagf = m_omega;
    }

  af_print(diagf);

  int row, column;
  float* values = getInputMatrix(row,column,false);
  m_inputForecastMatrix = af::array(row,column, values);

  if(m_verbose)
    af_print(m_inputForecastMatrix);

  int rows = diagf.dims(0) ;
  int rowT = rows * row;
  int cols = diagf.dims(1);

  af::array noise(rowT,cols);

  for(int i = 0 ; i < rowT ; i = i + rows)
    {
      af::seq r = seq(i, i + rows - 1);
      noise(r, span) = diagf;
    }

  if(m_verbose)
    af_print(noise);

  af::array varWeight = constant(0.0, row,V);

  af::array phi = m_kernel.calculateKernel(m_inputForecastMatrix , m_inputMatrix);
  if(m_verbose)
    af_print(phi);

  if(m_verbose)
    af_print(m_used);

  af::array indices = af::where( m_used == 1.0);
  if(m_verbose)
    af_print(indices);

  af::array phiUsed = phi(span, indices);
  if(m_verbose)
    af_print(phiUsed);

  m_outputForecastMatrix = matmul( phiUsed , m_Mu);

  if(m_verbose)
    af_print(m_outputForecastMatrix);

  for(int i = 0; i < V ; i++)
    {
      af::array t = af::solve( m_invSigma(span,span,i) ,phiUsed.T());

      af::array f = matmul  ( phi(span,indices) , t);
      af::array ls = af::diag(f);
      varWeight(span,i) = ls;
    }

  m_stdPrediction = af::sqrt(noise + varWeight);

  if(m_verbose)
    af_print(m_stdPrediction);

  writeOutput();

}

bool MRVM::gdalRegistered()
{
  return s_gdalRegistered;
}

void MRVM::readProject()
{
  QFile file(m_file.absoluteFilePath());
  
  if(file.open(QIODevice::ReadOnly))
    {
      QXmlStreamReader xmlReader(& file);
      
      
      while(!xmlReader.isEndDocument() && !xmlReader.hasError())
        {
          xmlReader.readNext();
          
          QXmlStreamReader::TokenType tokenType = xmlReader.tokenType();
          
          if(!xmlReader.name().compare("MRVM", Qt::CaseInsensitive) && !xmlReader.hasError())
            {
              while (!(xmlReader.isEndElement() && !xmlReader.name().compare("MRVM", Qt::CaseInsensitive)) && !xmlReader.hasError())
                {
                  switch(tokenType)
                    {
                    case QXmlStreamReader::StartElement:
                      {
                        if(!xmlReader.name().compare("MRVM", Qt::CaseInsensitive))
                          {
                            QXmlStreamAttributes attributes = xmlReader.attributes();
                            
                            for(QXmlStreamAttributes::iterator it = attributes.begin() ; it != attributes.end() ; it++)
                              {
                                QXmlStreamAttribute attribute = *it;
                                
                                if(!attribute.name().compare("name", Qt::CaseInsensitive))
                                  {
                                    m_name = attribute.value().toString();
                                  }
                                else if(!attribute.name().compare("mode", Qt::CaseInsensitive))
                                  {
                                    QString mode = attribute.value().toString();
                                    
                                    if(!mode.compare(mode,"Training" , Qt::CaseInsensitive))
                                      {
                                        m_mode = Training;
                                      }
                                    else if(!mode.compare(mode,"Regression" , Qt::CaseInsensitive))
                                      {
                                        m_mode = Regression;
                                      }
                                    else if(!mode.compare(mode,"TrainingAndRegression" , Qt::CaseInsensitive))
                                      {
                                        m_mode = TrainingAndRegression;
                                      }
                                  }
                              }
                          }
                        else if(!xmlReader.name().compare("Tolerance", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();
                            bool ok = false;
                            double v = tol.toDouble(&ok);
                            if(ok)
                              {
                                m_tolerance = v;
                              }
                          }
                        else if(!xmlReader.name().compare("MaxNumberOfIterations", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();
                            bool ok = false;
                            double v = tol.toDouble(&ok);
                            if(ok)
                              {
                                m_maxNumberOfIterations = v;
                              }
                          }
                        else if(!xmlReader.name().compare("NumberOfIterations", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();
                            bool ok = false;
                            double v = tol.toDouble(&ok);
                            if(ok)
                              {
                                m_numberOfIterations = v;
                              }
                          }
                        else if(!xmlReader.name().compare("Converged", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();

                            if(!tol.compare("true" , Qt::CaseInsensitive))
                              m_converged = true;
                            else
                              m_converged = false;
                          }
                        else if(!xmlReader.name().compare("MaxAlphaChange", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();
                            bool ok = false;
                            double v = tol.toDouble(&ok);
                            if(ok)
                              {
                                m_maxChangeAlpha = v;
                              }
                          }
                        else if(!xmlReader.name().compare("MinAlphaChange", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();
                            bool ok = false;
                            double v = tol.toDouble(&ok);
                            if(ok)
                              {
                                m_minChangeAlpha = v;
                              }
                          }
                        else if(!xmlReader.name().compare("Verbose", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();

                            if(!tol.compare("true" , Qt::CaseInsensitive))
                              m_verbose = true;
                            else
                              m_verbose = false;
                          }
                        else if(!xmlReader.name().compare("MatrixOutputFile", Qt::CaseInsensitive))
                          {
                            m_matrixOutputFile = xmlReader.readElementText();
                            QFileInfo fileInfo(m_matrixOutputFile);

                            if(fileInfo.exists())
                              {

                                if(af::readArrayCheck(m_matrixOutputFile.toStdString().c_str() , "InputMatrix") != -1)
                                  m_inputMatrix = readArray(m_matrixOutputFile.toStdString().c_str() , "InputMatrix");

                                if(af::readArrayCheck(m_matrixOutputFile.toStdString().c_str() , "Used") != -1)
                                  m_used = readArray(m_matrixOutputFile.toStdString().c_str() , "Used");

                                if(af::readArrayCheck(m_matrixOutputFile.toStdString().c_str() , "Alpha") != -1)
                                  m_alpha = readArray(m_matrixOutputFile.toStdString().c_str() , "Alpha");

                                if(af::readArrayCheck(m_matrixOutputFile.toStdString().c_str() , "InverseSigma") != -1)
                                  m_invSigma = readArray(m_matrixOutputFile.toStdString().c_str() , "InverseSigma");

                                if(af::readArrayCheck(m_matrixOutputFile.toStdString().c_str() , "Omega") != -1)
                                  m_omega = readArray(m_matrixOutputFile.toStdString().c_str() , "Omega");

                                if(af::readArrayCheck(m_matrixOutputFile.toStdString().c_str() , "Mu") != -1)
                                  m_Mu = readArray(m_matrixOutputFile.toStdString().c_str() , "Mu");
                              }
                          }
                        else if(!xmlReader.name().compare("Kernel", Qt::CaseInsensitive))
                          {
                            while(!(xmlReader.isEndElement() && !xmlReader.name().compare("Kernel", Qt::CaseInsensitive)) && !xmlReader.hasError())
                              {
                                QXmlStreamAttributes attributes = xmlReader.attributes();
                                
                                for(QXmlStreamAttributes::iterator it = attributes.begin() ; it != attributes.end() ; it++)
                                  {
                                    QXmlStreamAttribute attribute = *it;
                                    
                                    if(!attribute.name().compare("KernelType", Qt::CaseInsensitive))
                                      {
                                        QString kernelType = attribute.value().toString();
                                        
                                        if(!kernelType.compare("Gaussian", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Gaussian);
                                          }
                                        else if(!kernelType.compare("Laplace", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Laplace);
                                          }
                                        else if(!kernelType.compare("Polynomial", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Polynomial);
                                          }
                                        else if(!kernelType.compare("HomogeneousPolynomail", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::HomogeneousPolynomail);
                                          }
                                        else if(!kernelType.compare("Spline", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Spline);
                                          }
                                        else if(!kernelType.compare("Cubic", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Cubic);
                                          }
                                        else if(!kernelType.compare("Distance", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Distance);
                                          }
                                        else if(!kernelType.compare("ThinPlateSpline", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::ThinPlateSpline);
                                          }
                                        else if(!kernelType.compare("Bubble", Qt::CaseInsensitive))
                                          {
                                            m_kernel.setKernelType(Kernel::Bubble);
                                          }
                                      }
                                  }
                                
                                if(!xmlReader.name().compare("LengthScale",Qt::CaseInsensitive))
                                  {
                                    QString lengthScale = xmlReader.readElementText();
                                    m_kernel.setLengthScale(lengthScale.toDouble());
                                  }
                                else if(!xmlReader.name().compare("UseBias",Qt::CaseInsensitive))
                                  {
                                    QString useBias = xmlReader.readElementText();
                                    
                                    if(!useBias.compare("true" , Qt::CaseInsensitive))
                                      m_kernel.setUseBias(true);
                                    else
                                      m_kernel.setUseBias(false);
                                  }
                                else if(!xmlReader.name().compare("PolynomialPower",Qt::CaseInsensitive))
                                  {
                                    QString power = xmlReader.readElementText();
                                    
                                    m_kernel.setPolynomialPower(power.toDouble());
                                  }
                                
                                xmlReader.readNext();
                              }
                          }
                        else if(!xmlReader.name().compare("InputItems", Qt::CaseInsensitive))
                          {
                            cout << "Input Items\n" << endl;
                            
                            while(!(xmlReader.isEndElement() && !xmlReader.name().compare("InputItems", Qt::CaseInsensitive)) && !xmlReader.hasError())
                              {
                                if(!xmlReader.name().compare("MRVMItem",Qt::CaseInsensitive))
                                  {
                                    MRVMItem* item = readMRVMItem(MRVMItem::Input, xmlReader);
                                    
                                    cout << item->toString().toStdString() << endl;
                                    
                                    if(item)
                                      m_inputItems.append(item);
                                  }
                                
                                xmlReader.readNext();
                              }
                          }
                        else if(!xmlReader.name().compare("OutputItems", Qt::CaseInsensitive))
                          {
                            cout << "Output Items\n" << endl;
                            

                            while (!(xmlReader.isEndElement()
                                     && !xmlReader.name().compare("OutputItems", Qt::CaseInsensitive))
                                   && !xmlReader.hasError())
                              {
                                if(!xmlReader.name().compare("MRVMItem",Qt::CaseInsensitive))
                                  {
                                    MRVMItem* item = readMRVMItem(MRVMItem::Output, xmlReader);
                                    
                                    cout << item->toString().toStdString() << endl;
                                    
                                    if(item)
                                      m_outputItems.append(item);
                                  }
                                
                                xmlReader.readNext();
                              }
                          }
                        
                      }
                      break;
                    }
                  xmlReader.readNext();
                }
            }
        }
      
      ASSERT(!xmlReader.hasError(), xmlReader.errorString().toStdString().c_str());
      
      qaqc();
      
      file.flush();
      file.close();
    }
}

MRVMItem* MRVM::readMRVMItem(MRVMItem::IOType iotype, QXmlStreamReader& xmlReader)
{
  MRVMItem* item = NULL;
  
  if(!xmlReader.name().compare("MRVMItem" , Qt::CaseInsensitive))
    {
      QXmlStreamAttributes attributes =  xmlReader.attributes();
      QString type, name;
      
      for(QXmlStreamAttributes::iterator it = attributes.begin() ;
          it != attributes.end(); it ++)
        {
          QXmlStreamAttribute attribute  = *it;
          
          if(!attribute.name().compare("type",Qt::CaseInsensitive))
            {
              type = attribute.value().toString();
            }
          else if(!attribute.name().compare("name",Qt::CaseInsensitive))
            {
              name = attribute.value().toString();
            }
        }
      
      if(!type.compare("RealMRVMItem", Qt::CaseInsensitive))
        {
          item = new RealMRVMItem(iotype, name);
          item->readXML(xmlReader);
        }
      else if(!type.compare("RealArrayMRVMItem", Qt::CaseInsensitive))
        {
          item = new RealArrayMRVMItem(iotype, name);
          item->readXML(xmlReader);
        }
      else if(!type.compare("CategoricalMRVMItem", Qt::CaseInsensitive))
        {
          item = new CategoricalMRVMItem(iotype, name);
          item->readXML(xmlReader);
        }
      else if(!type.compare("RealRaster", Qt::CaseInsensitive))
        {
          item = new RealRaster(iotype, name);
          item->readXML(xmlReader);
        }
      else if(!type.compare("CategoricalRaster", Qt::CaseInsensitive))
        {
          //item = new CategoricalRaster(iotype, name);
           //item->readXML(xmlReader);
        }
    }
  
  return item;
}

void MRVM::qaqc()
{
  
  if(m_inputItems.count())
    {
      
      int numberOfTraining = m_inputItems[0]->numTrainingValues();
      int numberOfForecastValues = m_inputItems[0]->numForecastValues();
      
      for(QList<MRVMItem*>::iterator it = m_inputItems.begin() ; it != m_inputItems.end() ; it++)
        {
          MRVMItem* item = (*it);
          
          ASSERT(item->numTrainingValues() == numberOfTraining, "Number of training values are not equal for all items" );
          ASSERT(item->numForecastValues() == numberOfForecastValues, "Number of forecast values are not equal for all items" );
        }
      
      for(QList<MRVMItem*>::iterator it = m_outputItems.begin() ; it != m_outputItems.end() ; it++)
        {
          MRVMItem* item = (*it);
          ASSERT(item->numTrainingValues() == numberOfTraining, "Number of training values are not equal for all items" );
          
          const QList<QString> forecastValuesAsString  = item->forecastValuesAsString();
          
          if(forecastValuesAsString.count() != numberOfForecastValues)
            {
              QList<QString> forecast;
              
              for(int m = 0 ; m < numberOfForecastValues ; m++)
                {
                  forecast.append("");
                }
              
              item->setForecastValuesAsString(forecast);
            }
          
          const QList<QString> forecastUncertaintyAsString  = item->forecastUncertaintyValuesAsString();
          
          if(forecastUncertaintyAsString.count() != numberOfForecastValues)
            {
              
              QList<QString> forecast;
              
              for(int m = 0 ; m < numberOfForecastValues ; m++)
                {
                  forecast.append("");
                }
              
              item->setForecastUncertaintyValueAsString(forecast);
            }
        }
    }
}

float* MRVM::getInputMatrix(int& row, int& column, bool training)
{
  column = inputDimension();
  
  float* values = NULL;
  
  if(training)
    {
      row =  m_inputItems[0]->numTrainingValues();
      
      values  = new float[row*column];
      
      for(int r = 0 ; r < row ; r++)
        {
          
          int currentC = 0;
          
          for(int iitem = 0 ; iitem < m_inputItems.length() ; iitem++)
            {
              MRVMItem* item = m_inputItems[iitem];
              
              int ilength = item->columnCount();
              
              float* ivalues = item->trainingValues(r);
              
              for(int c = 0; c < ilength; c++)
                {
                  values[r + currentC*row] = ivalues[c];
                  
                  currentC ++;
                }
              
              delete[] ivalues;
            }
        }
    }
  else
    {
      row =  m_inputItems[0]->numForecastValues();
      
      values  = new float[row*column];
      
      for(int r = 0 ; r < row ; r++)
        {
          
          int currentC = 0;
          
          for(int iitem = 0 ; iitem < m_inputItems.length() ; iitem++)
            {
              MRVMItem* item = m_inputItems[iitem];
              
              int ilength = item->columnCount();
              
              float* ivalues = item->forecastValues(r);
              
              for(int c = 0; c < ilength; c++)
                {
                  values[r + currentC*row] = ivalues[c];
                  
                  currentC ++;
                }
              
              delete[] ivalues;
            }
        }
    }
  
  
  return values;
}

float* MRVM::getTargetMatrix(int& row, int& column, bool training)
{
  column = targetDimension();
  
  float* values = NULL;
  
  if(training)
    {
      row = m_outputItems[0]->numTrainingValues();
      
      values = new float[row*column];
      
      for(int r = 0 ; r < row ; r++)
        {
          
          int currentC = 0;
          
          for(int iitem = 0 ; iitem < m_outputItems.length() ; iitem++)
            {
              MRVMItem* item = m_outputItems[iitem];
              
              int ilength = item->columnCount();
              float* ivalues = item->trainingValues(r);
              
              for(int c = 0; c < ilength; c++)
                {
                  values[r + currentC*row] = ivalues[c];
                  
                  currentC ++;
                }
              
              delete[] ivalues;
            }
        }
    }
  else
    {
      row = m_outputItems[0]->numForecastValues();
      
      values = new float[row*column];
      
      for(int r = 0 ; r < row ; r++)
        {
          
          int currentC = 0;
          
          for(int iitem = 0 ; iitem < m_outputItems.length() ; iitem++)
            {
              MRVMItem* item = m_outputItems[iitem];
              
              int ilength = item->columnCount();
              float* ivalues = item->forecastValues(r);
              
              for(int c = 0; c < ilength; c++)
                {
                  values[r + currentC*row] = ivalues[c];
                  
                  currentC ++;
                }
              
              delete[] ivalues;
            }
        }
    }
  
  return values;
}

void MRVM::writeOutput()
{
  int r = m_inputForecastMatrix.dims(0);
  
  ASSERT(m_inputForecastMatrix.dims(0) == m_outputForecastMatrix.dims(0), "Output dimensions are not equal to input dimensions");
  
  for(int i = 0; i < r; i++)
    {
      int c = 0;
      
      for(int j = 0 ; j < m_outputItems.count() ; j++)
        {
          MRVMItem* item = m_outputItems[j];
          int columnCount = item->columnCount();

          seq jc = seq(c , c + columnCount - 1);

          af::array rvalues = m_outputForecastMatrix(i, jc);

          af::array uvalues = m_stdPrediction(i, jc);

          float* rvaluesf = rvalues.host<float>();
          item->setForecastValues(i,rvaluesf);
          delete[] rvaluesf;

          float* uvaluesf = uvalues.host<float>();
          item->setForecastUncertaintyValues(i, uvaluesf);
          delete[] uvaluesf;

          c = c + columnCount;
        }
    }


}

QPair<int,double> MRVM::calculateDeltaL(const af::array &mask, const af::array &alpha, const af::array &alphaNew, const af::array &sPrime, const af::array &qPrime, const af::array &s, const af::array &q)
{
  int task = 0;
  af::array deltaL;

  if(!alltrue<bool>(isInf(alphaNew)))
    {
      if(alltrue<bool>(mask > 0))
        {
          task = 1;
          af::array tempCalc = (1.0 / alphaNew) - (1.0 / alpha);
          //af_print(tempCalc);

          float* tCalc = tempCalc.host<float>();

          af::array qPrimeSq = af::pow(qPrime,2);
          // af_print(qPrimeSq);

          af::array temp1 = qPrimeSq / (sPrime + 1.0/tCalc[0]);
          // af_print(temp1);


          af::array temp2 = af::log(1+ sPrime*tCalc[0]);
          // af_print(temp2);


          delete[] tCalc;

          deltaL = af::sum(temp1 - temp2);
        }
      else
        {
          task = 2;

          float* values = alphaNew.host<float>();
          af::array tempCalc = s + values[0];

          // af_print(tempCalc);


          af::array alphaTemp(tempCalc.dims(0) , tempCalc.dims(1));
          alphaTemp = values[0];

          delete[] values;


          deltaL =af::sum(af::pow(q,2)/tempCalc + af::log(alphaTemp/tempCalc));
        }
    }
  else if(alltrue<bool>(mask > 0))
    {
      task =  3;

      float* a_alpha = alpha.host<float>();

      af::array qPrimeSq = af::pow(qPrime,2);
      deltaL = af::sum(qPrimeSq / (sPrime - a_alpha[0]) - af::log(1 - sPrime/a_alpha[0]));

      delete[] a_alpha;
    }
  else
    {
      return  QPair<int,double>(0,af::NaN);
    }

  if(af::alltrue<bool>(af::imag(deltaL) > 0))
    {
      return QPair<int,double>(task,-af::Inf);
    }
  else
    {
      float* f = deltaL.host<float>();
      double v = f[0];
      delete[] f;
      return QPair<int,double>(task,v);
    }
}

af::array MRVM::corrcov(const af::array &cov)
{
  int m = cov.dims(0);
  int n = cov.dims(1);
  
  assert(m == n);
  
  af::array output(m,m);
  
  //  gfor(seq i , m)
  for(int i = 0 ; i < m ; i++)
    {
      for(int j = 0; j < n ; j++)
        {
          output(i,j) = cov(i,j) / af::sqrt(cov(i,i) * cov(j,j));
        }
    }
  
  return output;
}
