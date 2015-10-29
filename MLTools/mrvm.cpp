#include "stdafx.h"
#include "mrvm.h"
#include "QDebug"
#include <iostream>
#include "roots.h"

using namespace std;

bool MRVM::gdalRegistered = false;


MRVM::MRVM(const QFileInfo& file)
  :QObject(), m_mode(TrainingAndRegression)
{
  if(!gdalRegistered)
    {
      GDALAllRegister();
      gdalRegistered = true;
    }
  
  m_file = file;
  readFile();
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

int MRVM::trainingItemsCount() const
{
  return m_trainingItemsCount;
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

MRVM::Mode MRVM::mode() const
{
  return m_mode;
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

void MRVM::save()
{
  QFile file(m_file.absoluteFilePath());

  if(file.open(QIODevice::WriteOnly | QIODevice::Truncate))
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

      writer.writeTextElement("AverageNumberOfIterations", QString::number(m_averageNumberOfIterations));

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
              writeMRVMItem(item, writer);
            }

          writer.writeEndElement();
        }

      if(m_outputItems.length() > 0)
        {
          writer.writeStartElement("OutputItems");

          for(int i = 0 ; i < m_outputItems.length() ; i++)
            {
              MRVMItem* item = m_outputItems[i];
              writeMRVMItem(item, writer);
            }

          writer.writeEndElement();
        }


      if(m_inputItems.count() > 0 && m_trainingItemsCount <= m_inputItems[0]->values().length())
        {
          writer.writeStartElement("InputValues");

          for (int i= 0 ; i < m_trainingItemsCount ; i++)
            {
              QString value = m_inputItems[0]->values()[i];

              for(int j =1 ; j < m_inputItems.length() ; j++)
                {
                  value = value + "," + m_inputItems[j]->values()[i];
                }

              writer.writeTextElement("Value", value);
            }

          writer.writeEndElement();

          writer.writeStartElement("TargetValues");

          for (int i= 0 ; i < m_trainingItemsCount ; i++)
            {
              QString value = m_outputItems[0]->values()[i];

              for(int j =1 ; j < m_outputItems.length() ; j++)
                {
                  value = value + "," + m_outputItems[j]->values()[i];
                }

              writer.writeTextElement("Value", value);
            }

          writer.writeEndElement();


          if(m_inputItems[0]->values().length() - m_trainingItemsCount > 0)
            {
              writer.writeStartElement("ForecastInputValues");

              int max = m_inputItems[0]->values().length();

              for (int i= m_trainingItemsCount  ; i < max ; i++)
                {
                  QString value = m_inputItems[0]->values()[i];

                  for(int j =1 ; j < m_inputItems.length() ; j++)
                    {
                      value = value + "," + m_inputItems[j]->values()[i];
                    }

                  writer.writeTextElement("Value", value);
                }

              writer.writeEndElement();

              writer.writeStartElement("ForecastOutputValues");

              max = m_outputItems[0]->values().length();

              for (int i= m_trainingItemsCount ; i < max  ; i++)
                {
                  QString value = m_outputItems[0]->values()[i];

                  for(int j =1 ; j < m_outputItems.length() ; j++)
                    {
                      value = value + "," + m_outputItems[j]->values()[i];
                    }

                  writer.writeTextElement("Value", value);
                }

              writer.writeEndElement();
            }
        }

      writer.writeEndElement();
      writer.writeEndDocument();

      file.flush();
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
}

void MRVM::performTraining()
{

  switch (algmode)
    {
    case 0:
      mrvm1();
      break;
    case 1:
      mrvm2();
      break;
    default:
      fmrvm1();
      break;
    }
}

void MRVM::mrvm1()
{

}

void MRVM::mrvm2()
{
  int row, column;

  float* values = getInputMatrix(row,column);
  m_inputMatrix = array(row,column, values);
  delete[] values;
  //af_print(m_inputMatrix);

  values = getTargetMatrix(row,column);
  m_targetMatrix = array(row,column,values );
  delete[] values;
  //af_print(m_targetMatrix);

  array phi = m_kernel.calculateKernel(m_inputMatrix,m_inputMatrix);
  //af_print(phi);

  int N = m_targetMatrix.dims(0);
  int V = m_targetMatrix.dims(1);

  assert(phi.dims(0) == N && phi.dims(1) == N+1);

  m_alpha = af::constant(af::Inf, N+1,1);
  //af_print(m_alpha);

  array beta = 1.0 / (0.1 * af::var(m_targetMatrix));
  // af_print(beta);

  array mask = m_alpha < af::NaN;

  array nmask;
  array Sigma;
  array PhiMask;
  array Mu;
  //array invSigma;

  //af_print(mask);

  bool converged = false;
  array phiSigmaPhi = af::constant(0.0,N,N,V);
  //  af_print(phiSigmaPhi);
  int iterCount = 0;
  
  for(iterCount = 0 ; iterCount < m_maxNumberOfIterations; iterCount++)
    {
      array sPrime = af::constant(af::NaN, N+1,1);
      array qPrime = af::constant(af::NaN, N+1,V);

      array s = af::constant(af::NaN, N+1,V);
      array q = af::constant(af::NaN, N+1,V);

      //      af_print(mask);
      array allAlphaInf = af::allTrue(mask == 0);
      //      af_print(allAlphaInf);

      for(int i = 0 ; i <  N+1; i++)
        {
          array phiSq = matmul(phi(span,i).T(),phi(span,i));
          //          af_print(phiSq);

          if(af::allTrue<bool>(allAlphaInf))
            {
              for(int j = 0 ; j < V ; j++)
                {
                  sPrime(i,j) = beta(j) * phiSq;
                  qPrime(i,j) = beta(j) * matmul(phi(span,i).T(),m_targetMatrix(span,j));
                }
            }
          else
            {
              //  af_print(phiSigmaPhi(span,span));

              for(int j = 0 ; j < V ; j++)
                {
                  array power = af::pow(beta(j),2);
                  //                  af_print(power);

                  float* vvv = power.host<float>();
                  //                  af_print(power);

                  //  array ftemp =  matmul(phi(span,i).T(),phiSigmaPhi(span,span,j)) ;
                  //                  af_print(ftemp);

                  array temp  = vvv[0] * matmul(phi(span,i).T(),phiSigmaPhi(span,span,j)) ;
                  delete[] vvv;
                  sPrime(i,j) = beta(j) * phiSq - matmul(temp,phi(span,i));
                  qPrime(i,j) = beta(j) * matmul(phi(span,i).T(),m_targetMatrix(span,j)) - matmul(temp,m_targetMatrix(span,j));
                }
            }

          if (!allTrue<bool>(mask))
            {
              s(i,span) = sPrime(i,span);
              q(i,span) = qPrime(i,span);
            }
          else
            {
              array temp = m_alpha(i)/(m_alpha(i)-sPrime(i,span));
              s(i,span) = temp * sPrime(i,span);
              q(i,span) = temp * qPrime(i,span);
            }
        }

      //af_print(sPrime);
      // af_print(qPrime);

      // af_print(s);
      //af_print(q);

      array deltaL = af::constant(af::NaN, N+1,1);
      // af_print(deltaL);

      array task = af::constant(0,N+1,1);
      //  af_print(task);

      array alphaNew = af::constant(af::NaN, N+1,1);

      for(int i = 0 ; i < N+1 ; i++)
        {
          array polyCoeff = af::constant(0, V+1 , 2*V+1 );
          array temp(1,1); temp = V;

          for(int j = 0;j < V; j++)
            {
              // af_print(temp);
              array te(1,3); te(0,0) = 1; te(0,1) = 2*s(i,j); te(0,2) = af::pow(s(i,j),2);
              //  af_print(te);
              array ff = af::convolve(temp,te,AF_CONV_EXPAND);
              temp = ff;
              //af_print(temp);
            }

          //af_print(polyCoeff);
          polyCoeff(0,span) = temp;
          // af_print(polyCoeff);

          for(int j = 0 ; j < V ; j++)
            {
              temp = af::array(1,3);
              temp(0,0) = -1; temp(0,1) = -1 * s(i,j) - af::pow(q(i,j),2); temp(0,2) = 0;

              for(int k = 0; k < V; k++)
                {
                  if(k != j)
                    {
                      array temp1(1,3); temp1(0,0) = 1; temp1(0,1) = 2*s(i,k); temp1(0,2) = pow(s(i,k),2);
                      // af_print(temp1);
                      //af_print(temp);

                      array ff = convolve(temp,temp1,AF_CONV_EXPAND);
                      temp = ff;
                      //af_print(temp);
                    }
                }
              polyCoeff(j+1,span) = temp;
            }

          //af_print(polyCoeff);
          array sumPolyCoeff = sum(polyCoeff);
          //af_print(sumPolyCoeff);

          if(allTrue<bool>(sumPolyCoeff(0) != 0))
            {
              cout << "The order of polynomial should be 2V-1" << endl;
            }

          int length  = sumPolyCoeff.dims(1);

          //cout << "sum dims " << length <<  endl;

          float* values = sumPolyCoeff.host<float>();

          vector<double> sumPolyCoeffs;

          // cout << "coeffs fixed ";

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
                      //cout << values[ii]  << ",";
                    }
                }
              else
                {
                  sumPolyCoeffs.push_back(values[ii]);
                  // cout << values[ii]  << ",";
                  k++;
                }
            }

          //cout << endl;

          delete[] values;

          if(k == 0)
            {
              alphaNew(i) = af::Inf;
              QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i), alphaNew(i),sPrime(i,span), qPrime(i,span),s(i,span),q(i,span));
              task(i,span) = t.first;
              deltaL(i) = t.second;
              // af_print(task);
              //af_print(deltaL);

            }
          else
            {
              //af_print(sumPolyCoeff);
              int size = 0;
              complex<double>* croots = Roots::roots(sumPolyCoeffs, size);
              vector<float> roots;
              array posRealRoots;

              for(int ff = 0 ; ff < size ; ff++)
                {
                  // cout << croots[ff] << endl;

                  if(fabs(croots[ff].imag()) < 10e-16 && croots[ff].real() > 0)
                    {
                      roots.push_back(croots[ff].real());
                    }
                }

              delete[] croots;

              if(roots.size())
                {
                  posRealRoots = array(1,roots.size() , & roots[0]);
                  //af_print(posRealRoots);
                }


              if(roots.size() == 0)
                {
                  alphaNew(i) = af::Inf;
                  QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i),alphaNew(i),sPrime(i,span),qPrime(i,span),s(i,span),q(i,span));

                  task(i,span) = t.first;
                  deltaL(i) = t.second;
                  //af_print(task);
                  //af_print(deltaL);
                }
              else if(roots.size() == 1)
                {
                  alphaNew(i) = posRealRoots;
                  //af_print(alphaNew(i));

                  QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i),alphaNew(i),sPrime(i,span),qPrime(i,span),s(i,span),q(i,span));
                  task(i,span) = t.first;
                  deltaL(i) = t.second;
                  //af_print(task);
                  // af_print(deltaL);

                }
              else
                {

                  array taskRoots = af::constant(0, roots.size() , 1);
                  //af_print(taskRoots);

                  array deltaRoots = af::constant(af::NaN , roots.size(),1);
                  //af_print(deltaRoots);


                  for(int k = 0; k < roots.size() ; k++)
                    {
                      QPair<int,double> t = calculateDeltaL(mask(i),m_alpha(i),posRealRoots(i),sPrime(i,span),qPrime(i,span),s(i,span),q(i,span));
                      taskRoots(k,span) = t.first;
                      deltaRoots(k) = t.second;
                    }

                  // af_print(taskRoots); af_print(deltaRoots);

                  array vals, indexes;
                  af::max(vals,indexes,deltaRoots);

                  // af_print(vals);
                  // af_print(indexes);

                  //af_print(posRealRoots);
                  //af_print(posRealRoots(indexes));

                  //af_print(task);
                  //af_print(task(indexes,span));

                  alphaNew(i) = posRealRoots(indexes);
                  task(i,span) = taskRoots(indexes,span);
                }
            }
        }

      array i(1,1);

      //af_print(deltaL);

      if(allTrue<bool>(af::isNaN(deltaL)))
        {
          cout << "No further tasks ending iteration" << endl;
          break;
        }
      else if (allTrue<bool>((deltaL * (deltaL != af::NaN))  == -af::Inf))
        {
          i = rand() % (N+1) + 1;
        }
      else
        {
          array vals(1,1);

          vals = -1*af::Inf;

          for(int mm = 0 ; mm < deltaL.dims(0); mm++)
            {
              if(allTrue<bool>(deltaL(mm,0) > vals))
                {
                  vals = deltaL(mm,0);
                  i = mm;
                }
            }

          //af_print(deltaL);

          // af::max(vals,i,deltaL);
          //af_print(vals);
          //af_print(i);
        }

      if(iterCount != 0)
        {
          for(int j = 0 ; j < V; j++)
            {

              array temp = Sigma(span,span,j);
              array diagm;

              if(temp.dims(0) == 1 && temp.dims(1) == 1)
                diagm = temp;
              else
                diagm = af::diag(Sigma(span,span,j));

              //af_print(diagm);

              beta(j) = (N - nmask + matmul(m_alpha(mask).T(), diagm))/af::sum(af::pow(m_targetMatrix(span,j) - matmul(PhiMask, Mu(span,j)),2 ));

              //af_print(beta);
              // af_print(beta(j));
            }
        }

      // af_print(task);
      // af_print(i);

      array check = task(i,span);

      // af_print(task);
      // af_print(check);

      if(allTrue<bool>(check == 1))
        {
          array changeLogAlpha = log(m_alpha(i)/alphaNew(i));
          m_alpha(i) = alphaNew(i);

          // af_print(m_alpha);
          // af_print(alphaNew);
          // af_print(abs(m_alpha - alphaNew));


          if(allTrue<bool>(abs(changeLogAlpha) < m_tolerance))
            {
              if(allTrue<bool>(af::isInf(alphaNew(mask == 0))))
                {
                  converged = true;
                }
            }
        }
      else if(allTrue<bool>(check == 2))
        {
          m_alpha(i) = alphaNew(i);
          mask(i) = 1;

        }
      else if(allTrue<bool>(check == 3))
        {
          m_alpha(i) = af::Inf;
          mask(i) = 0;
        }
      else
        {
          //throw exception
          return ;
        }

      PhiMask = phi(span, mask);
      nmask = af::sum(mask);
      // af_print(phi)
      // af_print(PhiMask);
      // af_print(mask) ; af_print(nmask);

      // dtype ttype = nmask.type();
      //cout << ttype << endl;

      uint* vals = nmask.host<uint>();
      int size = vals[0];
      invSigma = af::constant(0,size,size,V); //invSigma = 0;
      //af_print(invSigma);

      Sigma = af::constant(0,size,size,V);

      Mu = af::constant(0, size,V); //Mu = 0;//  af::constant(0,(int)vals(0),V);
      delete[] vals;

      phiSigmaPhi = af::constant(0,N,N,V);

      array PhiMaskPhiMask = matmul(PhiMask.T(), PhiMask);

      for(int j = 0 ; j < V ; j++)
        {
          // af_print(mask);
          array almask = m_alpha(mask);

          //af_print(m_alpha);
          // af_print(almask);

          if(almask.dims(0)  > 1)
            {
              array ddiag = af::diag(almask,0,false);
              //    af_print(ddiag);
              almask = ddiag;
            }

          //af_print(PhiMaskPhiMask);
          //cout << "size size " << size << endl;

          // af_print(beta(j));

          float* betav = beta(j).host<float>();

          invSigma(span,span,j) = betav[0] * PhiMaskPhiMask + almask;


          array upper;
          int success = af::cholesky(upper, invSigma(span,span,j));

          //af_print(upper);

          array inVU = af::inverse(upper);

          Sigma(span,span,j) = matmul(inVU,inVU.T());
          Mu(span,j) = betav[0]* (matmul(Sigma(span,span,j),matmul(PhiMask.T(), m_targetMatrix(span,j))));
          delete[] betav;

          //af_print(Mu);

          phiSigmaPhi(span,span,j) = matmul(PhiMask,matmul(Sigma(span,span,j),PhiMask.T()));
        }

      iterCount++;
      cout << "Iteration Number : " << iterCount << endl;

      if(converged)
        break;
    }

  cout << "Total Number of Iterations : " << iterCount << endl;

  m_used = mask;

  //af_print(m_used);
  // af_print(phi);

  //array(phi(span, m_used));

  array targetPredict = matmul(phi(span, m_used),Mu);

  // af_print(m_targetMatrix);
  array diff = m_targetMatrix - targetPredict;
  af_print(targetPredict);
  af_print(diff);

  array omegaC = matmul(diff.T() , diff);
  af_print(omegaC);

  array omegaTilde = matmul(diff.T() , diff) / ((N-1)*1.0);
  af_print(omegaTilde);

  array rhat = corrcov(omegaTilde);
  af_print(rhat);

  array dhat = af::constant(0,V,V);
  af_print(dhat);

  array t = 1.0 / af::sqrt(beta);
  af_print(t);

  int z = t.dims(0);

  gfor(seq f, z)
  {
    dhat(f,f) = t(f);
  }

  af_print(dhat);

  m_omega = (matmul(dhat,rhat),dhat);

  af_print(m_omega);
}

QPair<int,double> MRVM::calculateDeltaL(const array &mask, const array &alpha, const array &alphaNew, const array &sPrime, const array &qPrime, const array &s, const array &q)
{
  af_print(alphaNew);
  af_print(isInf(alphaNew));
  af_print(mask);
  int task = 0;
  array deltaL;

  if(!allTrue<bool>(isInf(alphaNew)))
    {
      if(allTrue<bool>(mask))
        {
          task = 1;
          array tempCalc = (1.0 / alphaNew) - (1.0 / alpha);
          array qPrimeSq = af::pow(qPrime,2);
          deltaL = af::sum(qPrimeSq / (sPrime + 1.0/tempCalc )
                           - af::log(1+ sPrime*tempCalc));
        }
      else
        {
          task = 2;
          af_print(alphaNew) ; af_print(s);

          float* values = alphaNew.host<float>();
          array tempCalc = s + values[0];
          delete[] values;

          array alphaTemp(tempCalc.dims(0) , tempCalc.dims(1));

          gfor(array mm , alphaTemp.dims(1))
          {
            for(int ll = 0 ; ll < alphaTemp.dims(0) ; ll++)
              {
                alphaTemp(ll,mm) = alphaNew;
              }
          }

          af_print(alphaTemp);
          deltaL =af::sum(af::pow(q,2)/tempCalc + af::log(alphaTemp/tempCalc));
        }
    }
  else if(allTrue<bool>(mask))
    {
      task =  3;
      array qPrimeSq = af::pow(qPrime,2);
      deltaL = af::sum(qPrimeSq / (sPrime - alpha) - af::log(1 - sPrime/alpha));
    }
  else
    {
      return  QPair<int,double>(0,af::NaN);
    }

  if(af::allTrue<bool>(af::imag(deltaL) > 0))
    {
      return QPair<int,double>(task,-af::Inf);
    }
  else
    {
      af_print(deltaL);
      float* f = deltaL.host<float>();
      double v = f[0];
      delete[] f;
      return QPair<int,double>(task,v);
    }
}

void MRVM::fmrvm1()
{
  int row, column;

  float* values = getInputMatrix(row,column);
  m_inputMatrix = array(row,column, values);
  delete[] values;
  af_print(m_inputMatrix);

  values = getTargetMatrix(row,column);
  m_targetMatrix = array(row,column,values);
  delete[] values;
  af_print(m_targetMatrix);

  array phi = m_kernel.calculateKernel(m_inputMatrix,m_inputMatrix);

  int N = m_targetMatrix.dims(0);
  int V = m_targetMatrix.dims(1);

  assert(phi.dims(0) == N && phi.dims(1) == N+1);
  
  //  array omega = 0.1 * af::cov(m_targetMatrix);
  //  m_alpha = af::constant(af::Inf , N+1, 1);
  //  array mask = m_alpha < af::Inf;
  
  
}

void MRVM::performRegression()
{

}

void MRVM::readFile()
{
  QFile file(m_file.absoluteFilePath());
  
  if(file.open(QIODevice::ReadOnly))
    {
      QXmlStreamReader xmlReader(& file);
      m_trainingItemsCount = 0;

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
                        else if(!xmlReader.name().compare("AverageNumberOfIterations", Qt::CaseInsensitive))
                          {
                            QString tol = xmlReader.readElementText();
                            bool ok = false;
                            double v = tol.toDouble(&ok);
                            if(ok)
                              {
                                m_averageNumberOfIterations = v;
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
                            while(!(xmlReader.isEndElement() && !xmlReader.name().compare("InputItems", Qt::CaseInsensitive)) && !xmlReader.hasError())
                              {
                                if(!xmlReader.name().compare("MRVMItem",Qt::CaseInsensitive))
                                  {
                                    MRVMItem* item = readMRVMItem(xmlReader);

                                    if(item)
                                      m_inputItems.append(item);
                                  }

                                xmlReader.readNext();
                              }
                          }
                        else if(!xmlReader.name().compare("OutputItems", Qt::CaseInsensitive))
                          {
                            while (!(xmlReader.isEndElement()
                                     && !xmlReader.name().compare("OutputItems", Qt::CaseInsensitive))
                                   && !xmlReader.hasError())
                              {
                                if(!xmlReader.name().compare("MRVMItem",Qt::CaseInsensitive))
                                  {
                                    MRVMItem* item = readMRVMItem(xmlReader);

                                    if(item)
                                      m_outputItems.append(item);
                                  }

                                xmlReader.readNext();

                              }
                          }
                        else if(!xmlReader.name().compare("InputValues", Qt::CaseInsensitive))
                          {
                            while (!(xmlReader.isEndElement()
                                     && !xmlReader.name().compare("InputValues", Qt::CaseInsensitive))
                                   && !xmlReader.hasError())
                              {

                                if(!xmlReader.name().compare("Value",Qt::CaseInsensitive))
                                  {
                                    QString valuesAsString = xmlReader.readElementText();

                                    QStringList valuesSplit = valuesAsString.split(",");

                                    assert(valuesSplit.length() == m_inputItems.count());

                                    for(int i = 0 ; i < m_inputItems.length() ;i++)
                                      {
                                        m_inputItems[i]->addValue(valuesSplit[i]);
                                      }

                                    m_trainingItemsCount ++;
                                  }
                                xmlReader.readNext();
                              }
                          }
                        else if(!xmlReader.name().compare("TargetValues", Qt::CaseInsensitive))
                          {
                            while (!(xmlReader.isEndElement()
                                     && !xmlReader.name().compare("TargetValues", Qt::CaseInsensitive))
                                   && !xmlReader.hasError())
                              {

                                if(!xmlReader.name().compare("Value",Qt::CaseInsensitive))
                                  {
                                    QString valuesAsString = xmlReader.readElementText();

                                    QStringList valuesSplit = valuesAsString.split(",");

                                    assert(valuesSplit.length() == m_outputItems.count());

                                    for(int i = 0 ; i < m_outputItems.length() ;i++)
                                      {
                                        m_outputItems[i]->addValue(valuesSplit[i]);
                                      }
                                  }

                                xmlReader.readNext();
                              }
                          }
                        else if(!xmlReader.name().compare("ForecastInputValues", Qt::CaseInsensitive))
                          {
                            while (!(xmlReader.isEndElement()
                                     && !xmlReader.name().compare("ForecastInputValues", Qt::CaseInsensitive))
                                   && !xmlReader.hasError())
                              {

                                if(!xmlReader.name().compare("Value",Qt::CaseInsensitive))
                                  {
                                    QString valuesAsString = xmlReader.readElementText();

                                    QStringList valuesSplit = valuesAsString.split(",");

                                    assert(valuesSplit.length() == m_inputItems.count());

                                    for(int i = 0 ; i < m_inputItems.length() ;i++)
                                      {
                                        m_inputItems[i]->addValue(valuesSplit[i]);
                                      }
                                  }

                                xmlReader.readNext();
                              }
                          }
                        else if(!xmlReader.name().compare("ForecastOutputValues", Qt::CaseInsensitive))
                          {
                            while (!(xmlReader.isEndElement()
                                     && !xmlReader.name().compare("ForecastOutputValues", Qt::CaseInsensitive))
                                   && !xmlReader.hasError())
                              {

                                if(!xmlReader.name().compare("Value",Qt::CaseInsensitive))
                                  {
                                    QString valuesAsString = xmlReader.readElementText();

                                    QStringList valuesSplit = valuesAsString.split(",");

                                    assert(valuesSplit.length() == m_outputItems.count());

                                    for(int i = 0 ; i < m_outputItems.length() ;i++)
                                      {
                                        m_outputItems[i]->addValue(valuesSplit[i]);
                                      }
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
      if(xmlReader.hasError())
        {
          cout << xmlReader.errorString().toStdString().c_str() << endl;
        }
      
      file.flush();
      file.close();
    }
}

MRVMItem* MRVM::readMRVMItem( QXmlStreamReader& xmlReader)
{
  MRVMItem* item = NULL;

  while (!(xmlReader.isEndElement() && !xmlReader.name().compare("MRVMItem" , Qt::
                                                                 CaseInsensitive)) && !xmlReader.hasError())
    {

      QXmlStreamAttributes attributes =  xmlReader.attributes();

      QString type, name;
      
      if(!xmlReader.name().compare("MRVMItem" , Qt::CaseInsensitive))
        {
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
              item = new RealMRVMItem(name);
            }
          else if(!type.compare("RealArrayMRVMItem", Qt::CaseInsensitive))
            {
              item = new RealArrayMRVMItem(name);
            }
          else if(!type.compare("CategoricalMRVMItem", Qt::CaseInsensitive))
            {
              item = new CategoricalMRVMItem(name);
            }
          else if(!type.compare("RealRaster", Qt::CaseInsensitive))
            {
              item = new RealRaster(name);
            }
          else if(!type.compare("CategoricalRaster", Qt::CaseInsensitive))
            {
              item = new CategoricalRaster(name);
            }
        }
      else if(!xmlReader.name().compare("Properties" , Qt::CaseInsensitive) && item)
        {
          QMap<QString,QString> properties;
          
          while (!(xmlReader.isEndElement() && !xmlReader.name().compare("Properties" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {

              if(!xmlReader.name().compare("KeyValue", Qt::CaseInsensitive))
                {
                  QXmlStreamAttributes propAttributes =  xmlReader.attributes();
                  QString key("");
                  QString value = xmlReader.readElementText();


                  for(QXmlStreamAttributes::iterator pit = propAttributes.begin() ;
                      pit != propAttributes.end(); pit ++)
                    {
                      QXmlStreamAttribute attribute = *pit;

                      if(!attribute.name().compare("Key",Qt::CaseInsensitive))
                        {
                          key = attribute.value().toString();
                          break;
                        }
                    }

                  if(!key.isNull()  && !key.isEmpty() &&  ! value.isNull() && !value.isEmpty())
                    {
                      properties[key] = value;
                    }

                }

              xmlReader.readNext() ;
            }

          item->setProperties(properties);
        }
      else if(!xmlReader.name().compare("Categories" , Qt::CaseInsensitive) && item)
        {
          QMap<int,QString> categories;

          while (!(xmlReader.isEndElement() && !xmlReader.name().compare("Categories" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {

              if(!xmlReader.name().compare("Category", Qt::CaseInsensitive))
                {
                  QString  name("");

                  QXmlStreamAttributes propAttributes =  xmlReader.attributes();

                  for(QXmlStreamAttributes::iterator pit = propAttributes.begin() ;
                      pit != propAttributes.end(); pit ++)
                    {
                      QXmlStreamAttribute attribute = *pit;

                      if(!attribute.name().compare("name",Qt::CaseInsensitive))
                        {
                          name = attribute.value().toString();
                          break;
                        }
                    }

                  QString index = xmlReader.readElementText();

                  if(!name.isNull()  && !name.isEmpty()
                     &&  ! index.isNull() && !index.isEmpty()
                     )
                    {
                      bool ok = false;

                      int ind = index.toInt(&ok);

                      if(ok)
                        {
                          categories[ind] = name;
                        }
                    }

                }

              xmlReader.readNext();
            }

          CategoricalMRVMItem* citem = static_cast<CategoricalMRVMItem*>(item);

          if(citem)
            {
              citem->setCategories(categories);
            }
        }

      xmlReader.readNext();
    }
  
  return item;
}

void MRVM::writeMRVMItem(MRVMItem* const item,  QXmlStreamWriter& xmlWriter)
{

  xmlWriter.writeStartElement("MRVMItem");
  xmlWriter.writeAttribute("name", item->name());
  xmlWriter.writeAttribute("type", item->type());

  if(item->properties().count() > 0)
    {
      xmlWriter.writeStartElement("Properties");

      QMap<QString,QString> properties = item->properties();

      for(QMap<QString,QString>::iterator it = properties.begin()
          ; it != properties.end()
          ; it ++)
        {
          xmlWriter.writeStartElement("KeyValue");
          xmlWriter.writeAttribute("Key",it.key());
          xmlWriter.writeCharacters(it.value());
          xmlWriter.writeEndElement();
        }

      xmlWriter.writeEndElement();
    }

  CategoricalMRVMItem* catTtem = dynamic_cast<CategoricalMRVMItem*>(item);

  if (catTtem)
    {
      QMap<int,QString> categories = catTtem->categories();

      if(categories.count() > 0)
        {
          xmlWriter.writeStartElement("Categories");

          for(QMap<int,QString>::iterator it = categories.begin() ;
              it != categories.end()
              ; it++)
            {
              xmlWriter.writeStartElement("Category");
              xmlWriter.writeAttribute("name",it.value());
              xmlWriter.writeCharacters(QString::number( it.key()));
              xmlWriter.writeEndElement();
            }

          xmlWriter.writeEndElement();
        }
    }

  xmlWriter.writeEndElement();
}

float* MRVM::getInputMatrix(int& row, int& column)
{
  column = inputDimension();
  row = m_trainingItemsCount;


  float* values = new float[row*column];

  for(int r = 0 ; r < row ; r++)
    {

      int currentC = 0;

      for(int iitem = 0 ; iitem < m_inputItems.length() ; iitem++)
        {
          MRVMItem* item = m_inputItems[iitem];

          int ilength = item->columnCount();
          double* ivalues = item->values(r);

          for(int c = 0; c < ilength; c++)
            {
              values[r + currentC*row] = ivalues[c];

              currentC ++;
            }

          delete[] ivalues;
        }
    }

  return values;
}

float* MRVM::getTargetMatrix(int& row, int& column)
{
  column = targetDimension();
  row = m_trainingItemsCount;

  float* values = new float[row*column];

  for(int r = 0 ; r < row ; r++)
    {

      int currentC = 0;

      for(int iitem = 0 ; iitem < m_outputItems.length() ; iitem++)
        {
          MRVMItem* item = m_outputItems[iitem];

          int ilength = item->columnCount();
          double* ivalues = item->values(r);

          for(int c = 0; c < ilength; c++)
            {
              values[r + currentC*row] = ivalues[c];

              currentC ++;
            }

          delete[] ivalues;
        }
    }

  return values;
}

void MRVM::setTargetMatrixRegressionValues(float *&values)
{
  
}

array MRVM::corrcov(const array &cov)
{
  int m = cov.dims(0);
  int n = cov.dims(1);

  assert(m == n);

  array output(m,m);

//  gfor(seq i , m)
  for(int i = 0 ; i < m ; i++)
  {
    for(int j = 0; j < n ; j++)
      {
        output(i,j) = cov(i,j) / af::sqrt(cov(i,i) + cov(j,j));
      }
  }

  return output;
}
