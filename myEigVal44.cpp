#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;


int myMonic4(double*, complex<double>*, bool* reality, bool& allreal);

int myEigVal44(int &dim, double** Mpp, complex<double>* EigVals){


  int ireturn=0; //we live in hope!

  double a=Mpp[0][0]; double b=Mpp[0][1]; double c=Mpp[0][2];  double d=Mpp[0][3];
  double e=Mpp[1][0]; double f=Mpp[1][1]; double g=Mpp[1][2];  double h=Mpp[1][3];
  double p=Mpp[2][0]; double q=Mpp[2][1]; double r=Mpp[2][2];  double s=Mpp[2][3];
  double t=Mpp[3][0]; double u=Mpp[3][1]; double v=Mpp[3][2];  double w=Mpp[3][3];

        double D1,D2,D3,D4,D5,D6,D7,D8,D9;

        D1=a*f-e*b;  D2=a*r-c*p;  D3=a*w-d*t;
        D4=f*r-g*q;  D5=f*w-h*u;  D6=r*w-v*s;
        D7=e*w-t*h;  D8=p*w-s*t;  D9=h*p-e*s;

        double A1,A2,A3,A4,alpha,beta,gamma;

        A1=f*D6-g*(q*w-u*s)+h*(q*v-r*u);
        A2=a*D6-c*(p*w-s*t)+d*(p*v-r*t);
        A3=t*(b*h-f*d)-u*(a*h-e*d)+w*D1;
        A4=p*(b*g-f*c)-q*(a*g-c*e)+r*D1;

        alpha=r*D7-g*D8+v*D9;
        beta =q*D7-f*D8+u*D9;
        gamma=e*(q*v-r*u)-f*(p*v-r*t)+g*(p*u-q*t);

        double A;

        A=a*A1-b*alpha+c*beta-d*gamma;

        double coeffs[4];

        coeffs[0]=A;
        coeffs[1]=-(A1+A2+A3+A4);
        coeffs[2]=D1+D2+D3+D4+D5+D6;
        coeffs[3]=-(a+f+r+w);

        int monic_return;
        bool allreal;
        bool reality[4];

        monic_return=myMonic4(coeffs,EigVals,reality,allreal);
        
      
        


     return ireturn;
}



