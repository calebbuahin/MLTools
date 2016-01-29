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

//  float values[] = {1,4,6,3,2,1,2,5,4};
//  af::array f = af::array(3,3,values);
//  af_print(f);
//  af_print(sum(f));

  if (argc > 1)
    {
      const char* filePath = argv[1];
      QFileInfo inputFile = QFileInfo(QString(filePath));

      if (inputFile.exists())
        {
          MRVM mRVMModel(inputFile);
          mRVMModel.start();


        }
      else
        {
          cout << "Enter valid input file path" << endl;
        }
    }

  return 0;
}
