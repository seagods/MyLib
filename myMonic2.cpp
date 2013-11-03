#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;

int myMonic2(double* realcoeffs, complex<double>* roots, bool* reality){

     int ireturn=0;
     double realpart[2],imagpart[2];

     double discrim;
     double b,c,zero;

     zero=0.0;

     b=realcoeffs[1];
     c=realcoeffs[0];
     // x^2+bx+c=0;

     discrim=b*b-4.0*c;

     realpart[0]=-b/2.0;
     realpart[1]=-b/2.0;

     if(discrim>=zero){
       realpart[0]=realpart[0]+sqrt(discrim)/2.0;
       realpart[1]=realpart[1]-sqrt(discrim)/2.0;
       imagpart[0]=zero;
       imagpart[1]=zero;
       reality[0]=true;
       reality[1]=true;
     }
     else{
       imagpart[0]=sqrt(-discrim)/2.0;
       imagpart[1]=-sqrt(-discrim)/2.0;
       reality[0]=false;
       reality[1]=false;
     }

     roots[0]=complex<double>(realpart[0],imagpart[0]);
     roots[1]=complex<double>(realpart[1],imagpart[1]);

     bool polish=true;
     if(polish){

//   just one Newton Raphson Step

     if(reality[0]){
        double xval=real(roots[0]);
        double xvalsq=xval*xval;
        double fval=xvalsq+realcoeffs[1]*xval+realcoeffs[0];
        double deriv=2.0*xval+realcoeffs[1];
        double newxval;
        if(abs(deriv)>1.0e-16){
          newxval=xval-fval/deriv;
          roots[0]=complex<double>(newxval,zero);
        }
        else{
           cout << "Trouble in Newton-Raphson in myMonic2\n"; exit(0);
        }
        xval=real(roots[1]);
        xvalsq=xval*xval;
        fval=xvalsq+realcoeffs[1]*xval+realcoeffs[0];
        deriv=2.0*xval+realcoeffs[1];
        newxval;
        if(abs(deriv)>1.0e-16){
          newxval=xval-fval/deriv;
          roots[1]=complex<double>(newxval,zero);
        }
        else{
           cout << "Trouble in Newton-Raphson in myMonic2\n"; exit(0);
        }
      }
      else{

        complex<double> xval=roots[0];
        complex<double> xvalsq=xval*xval;
        complex<double> fval=xvalsq+realcoeffs[1]*xval+realcoeffs[0];
        complex<double> deriv=2.0*xval+realcoeffs[1];
        complex<double> newxval;
        if(abs(deriv)>1.0e-16){
          newxval=xval-fval/deriv;
          roots[0]=newxval;
          roots[1]=conj(newxval);
        }

      }

      } //endif for polish

     return ireturn;
}



