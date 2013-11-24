#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;


int myMonic3(double*, complex<double>*, bool* reality, bool& allreal);


int myEigen33(int &dim, double** Mpp, complex<double>* EigVals, complex<double>** EigVecs, bool* reality, bool &allreal){


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

        reality[0]=false; reality[1]=false; reality[2]=false; 
        allreal=false;

        int monic_return;

        monic_return=myMonic3(coeffs,EigVals,reality,allreal);
        if(monic_return !=0){
          cout << "monic_return non zero in Eigen33\n";
        }

        double zero=0.0e0;
        double one=1.0e0;

//      eigenvectors (x1,x2,x3)
//      OK, this is bad, but for now we "fix" x3 to be 1.0+0.0 i
//      regardless and "hope for the best" Need to come back and 
//      do a "proper job" here.

//      The way things are fixed up, we can only have 
//      real real real
//      real comp comp



        double normalise;

        if(!reality[0]){
          cout << "Can't have complex first root in myEigen33!\n"; exit(0);
        }

        double diag1R,diag2R;
        double ex,why;
        double DET;

        complex<double>diag1Z,diag2Z;
        complex<double>exZ,whyZ;
        complex<double>DETZ;

        diag1R=a-real(EigVals[0]);
        diag2R=e-real(EigVals[0]);

        // first root;

        ex=b*f-diag2R*c;
        why=d*c-diag1R*f;

        DET=diag1R*diag2R-b*d;
        ex=ex/DET;
        why=why/DET;

        normalise=sqrt(one+ex*ex+why*why);


  //      normalise=one;

        EigVecs[0][0]=complex<double>(ex/normalise,zero);
        EigVecs[0][1]=complex<double>(why/normalise,zero);
        EigVecs[0][2]=complex<double>(one/normalise,zero); 

        if(reality[1]){  
             //then all are real
             if(!allreal){
                cout << "contradiction in myEigen33\n"; exit(0);
             }

        diag1R=a-real(EigVals[1]);
        diag2R=e-real(EigVals[1]);

        ex=-c*diag2R+b*f;
        why=d*c-f*diag1R;

        DET=diag1R*diag2R-b*d;
        ex=ex/DET;
        why=why/DET;

        normalise=sqrt(one+ex*ex+why*why);

        EigVecs[1][0]=complex<double>(ex/normalise,zero);
        EigVecs[1][1]=complex<double>(why/normalise,zero);
        EigVecs[1][2]=complex<double>(one/normalise,zero); 

        diag1R=a-real(EigVals[2]);
        diag2R=e-real(EigVals[2]);

        ex=-c*diag2R+b*f;
        why=d*c-f*diag1R;

        DET=diag1R*diag2R-b*d;
        ex=ex/DET;
        why=why/DET;

        normalise=sqrt(one+ex*ex+why*why);

        EigVecs[2][0]=complex<double>(ex/normalise,zero);
        EigVecs[2][1]=complex<double>(why/normalise,zero);
        EigVecs[2][2]=complex<double>(one/normalise,zero); 

        }

        if(!reality[1]){
          //complex conjugate pair

        diag1Z=a-EigVals[1];
        diag2Z=e-EigVals[1];

        exZ=-c*diag2Z+b*f;
        whyZ=d*c-f*diag1Z;

        DETZ=diag1Z*diag2Z-b*d;
        exZ=exZ/DETZ;
        whyZ=whyZ/DETZ;

        normalise=sqrt( abs(one+exZ*conj(exZ)+whyZ*conj(whyZ)) );

        EigVecs[1][0]=exZ/normalise;
        EigVecs[1][1]=whyZ/normalise;
        EigVecs[1][2]=complex<double>(one/normalise,zero); 

        EigVecs[2][0]=conj(EigVecs[1][0]);
        EigVecs[2][1]=conj(EigVecs[1][1]);
        EigVecs[2][2]=conj(EigVecs[1][2]);

        }



     return ireturn;
}



