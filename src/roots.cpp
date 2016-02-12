#include <include/stdafx.h>
#include <include/roots.h>
#include <stdlib.h>
#include <math.h>
#include <src/Eigen/Core>
#include <src/Eigen/Eigenvalues>
#include <iostream>
#include <omp.h>

using namespace std;
using namespace Eigen;


complex<double>* Roots::roots(const std::vector<double> &coeffs , int & size)
{
    int matsz = coeffs.size() - 1;
    complex<double>* vret = new complex<double>[matsz];
    size = matsz;

    MatrixXd companion_mat(matsz,matsz);

#pragma omp parallel for
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

    MatrixXcd eig = companion_mat.eigenvalues();


#pragma omp parallel for
    for(int i = 0; i < matsz; i++)
        vret[i] =  eig(i) ;

    return vret;
}

