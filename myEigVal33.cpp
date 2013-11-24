#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;


int myMonic3(double*, complex<double>*, bool* reality, bool& allreal);

int myEigVal33(int &dim, double** Mpp, complex<double>* EigVals){


  int ireturn=0; //we live in hope!


  double a=Mpp[0][0]; double b=Mpp[0][1]; double c=Mpp[0][2];
  double d=Mpp[1][0]; double e=Mpp[1][1]; double f=Mpp[1][2]; 
  double p=Mpp[2][0]; double q=Mpp[2][1]; double r=Mpp[2][2]; 

        double A,A1,A2,A3;

        A=a*(e*r-f*q)-b*(d*r-f*p)+c*(d*q-e*p);

        A1=a*r-c*p;
        A2=e*r-f*q;
        A3=a*e-d*b;


        double coeffs[3];

        coeffs[0]=-A;
        coeffs[1]=(A1+A2+A3);
        coeffs[2]=-(a+e+r);


        int monic_return;
        bool allreal;
        bool reality[3];

        monic_return=myMonic3(coeffs,EigVals,reality,allreal);
        if(monic_return != 0){
           cout << "monic_return non zero in EigVal33\n";
        }
        
     return ireturn;
}



