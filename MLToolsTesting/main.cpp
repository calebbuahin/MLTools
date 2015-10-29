#include <QtCore/QCoreApplication>
#include "mrvm.h"
#include <iostream>
#include <cstdio>
#include <QFileInfo>
#include <arrayfire.h>
#include <roots.h>

using namespace std;
using namespace af;

int main(int argc, char *argv[])
{
  // QCoreApplication a(argc, argv);
  //  Kernel k;

  //  float x[] = {0.098492,0.738511 ,0.438356 ,0.791304 ,0.826749,
  //               0.487803, 0.781121,0.206121 ,0.348345, 0.629994,
  //               0.581608,0.102269,0.295694,0.070442, 0.972727 ,
  //               0.810519,  0.877167, 0.593647, 0.644570, 0.133485,
  //               0.598017, 0.770663, 0.850891, 0.395262,0.821617 };

  //  float y[] = {0.098492,0.738511 ,0.438356 ,0.791304 ,0.826749,
  //               0.487803, 0.781121,0.206121 ,0.348345, 0.629994,
  //               0.581608,0.102269,0.295694,0.070442, 0.972727 ,
  //               0.810519,  0.877167, 0.593647, 0.644570, 0.133485,
  //               0.598017, 0.770663, 0.850891, 0.395262,0.821617 };

  //  array f(5,5,x);
  //  array t(5,5,y);

  //  float* fdas = t.host<float>();
  //  delete[] fdas;

  //  af_print(f);
  //  array ff = f < 0.5;
  //  af_print(ff);

  //  k.setLengthScale(1);
  //  array c = k.calculateGaussianKernel(f,t);
  //  af_print(c);

  //  float v[] = {-3,-1,6};
  //  array temp (1,3,v);
  //  af_print(temp);

  //  array t(1,3);
  //  t = 0;
  //  af_print(t);

  //  array tempR = af::root(temp,0);
  //  af_print(tempR);


  //  double  p[] = {0,-2,-4, 6};
  //  double  r[3]; //= {0,0};
  //  double  i[3];// = {0,0};
  //  bool failed = true;
  //  int count = 3;

  //  MRVM::roots(p,count,r,i,failed);

  //  if(!failed)
  //    {
  //      for(int m = 0 ; m < count; m++)
  //        cout<< r[m] << "+" << i[m] << endl;

  //    }

  //  array u (1,3);
  //  array v (1,2);

  //  u(0,0) = 1;
  //  u(0,1) = 0;
  //  u(0,2) = 1;

  //  v(0,0) = 2;
  //  v(0,1) = 7;

  //  array c = af::convolve(u,v, AF_CONV_EXPAND);
  //  af_print(c)
  //  array A(1,3);
  //  A(0,0) = 1; A(0,1) = 1; A(0,2) = 3;
  //  A(1,0) = 1; A(1,1) = 1; A(1,2) = 3;
  //  A(2,0) = 1; A(2,1) = 1; A(2,2) = 3;

  //  bool B = allTrue<bool>(A);
  //  af_print(A);
  //  array B = af::sum(A);
  //  cout << B << endl;

  //  array f = af::randn(5,5);
  //  af_print(f);
  //  af_print(f < 0);
  //  af_print(f (f < 0));

//  vector<double> c(4);
//  c[0] = -2.625e+05; c[1] = -9.329e+07;
//  c[2]= -8.279e+09; c[3] = 2e+09 ;

// // c[0] = 1 ; c[1]= 0; c[2]= 0; c[3]= 0; c[4]= -1.0;

//  vector<complex<double> > values = Roots::roots(c);

//  for(int i =0 ; i < values.size() ; i++)
//    {
//      complex<double> dd = values[i];
//      cout << dd << endl;
//    }


//  array temp = af::randu(10,1,f32);
//  af_print(temp);


//  array ind = af::corrcoef<double>(temp,temp);
//  af_print(ind);

//  af_print(temp(ind));

//  array test = af::constant(5,10);
//  temp(ind) = test;

//  af_print(temp);

  if (argc > 1)
    {
      const char* filePath = argv[1];
      QFileInfo inputFile = QFileInfo(QString(filePath));

      if (inputFile.exists())
        {
          MRVM mRVMModel(inputFile);
          mRVMModel.start();
          mRVMModel.save();

        }
      else
        {
          cout << "Enter valid input file path" << endl;
        }
    }

  return 0;
}
