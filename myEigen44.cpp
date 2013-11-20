#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;


int myMonic4(double*, complex<double>*, bool* reality, bool& allreal);

int myGaussElim(int &dimless, double** Areal, double* breal, double* xreal);
int myGaussElimZ(int &dimless, complex<double>** Acomp,  complex<double>* bcomp, complex<double>* xcomp);

int myEigen44(int &dim, double** Mpp, complex<double>* EigVals, complex<double>** EigVecs, bool* reality, bool &allreal){


  int ireturn=0; //we live in hope!
  int ireturnReal=0;
  int ireturnComp=0;



  double a=Mpp[0][0]; double b=Mpp[0][1]; double c=Mpp[0][2];  double d=Mpp[0][3];
  double e=Mpp[1][0]; double f=Mpp[1][1]; double g=Mpp[1][2];  double h=Mpp[1][3];
  double p=Mpp[2][0]; double q=Mpp[2][1]; double r=Mpp[2][2];  double s=Mpp[2][3];
  double t=Mpp[3][0]; double u=Mpp[3][1]; double v=Mpp[3][2];  double w=Mpp[3][3];


 // cout << "ARSE!" << endl; exit(0);


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

        reality[0]=false; reality[1]=false; reality[2]=false; reality[3]=false;
        allreal=false;

        int monic_return;
    //    bool allreal;
    //    bool reality[4];

        monic_return=myMonic4(coeffs,EigVals,reality,allreal);

        double zero=0.0e0;
        double one=1.0e0;
        int dimless=dim-1;  // dimension of system of linear equations


        bool somereal=false;
        if(reality[0])somereal=true;
        if(reality[2])somereal=true;

        double** Areal=new double*[dimless];
        double* breal=new double[dimless];
        double* xreal=new double[dimless];

        for(int i=0; i<dimless; i++){
           Areal[i]=new double[dimless]; 
        }

        complex<double>** Acomp=new complex<double>*[dimless];
        complex<double>*  bcomp=new complex<double>[dimless];
        complex<double>*  xcomp=new complex<double>[dimless];
     
        for(int i=0; i<dimless; i++){
           Acomp[i]=new complex<double>[dimless]; 
        }


//      OK, this is bad, but for now we "fix" x4 to be 1.0+0.0 i
//      regardless and "hope for the best" Need to come back and 
//      do a "proper job" here.

//      The way things are fixed up, we can only have 
//      real real real real
//      real real comp comp
//      comp comp real real
//      comp comp comp comp.

//      We CANNNOT possibly encounter 
//      real comp comp real
//      real comp real comp 
//      etc.
//      so we can safely assume one of four possibilities
        double normalise;

/**************   Vectors one  and two *******************************/

        if(reality[0]){
        for(int i=0; i<dimless; i++){
          for(int j=0; j<dimless; j++){ 
             Areal[i][j]=Mpp[i][j];       
        }
             Areal[i][i]=Areal[i][i]-real(EigVals[0]);
             breal[i]=-Mpp[i][dimless];
        }
             ireturnReal=myGaussElim(dimless, Areal, breal, xreal);
             if(ireturnReal != 0){
               cout << "myEigen44 failed in myGaussElim\n"; return ireturnReal;
             }

             normalise=sqrt(one+xreal[0]*xreal[0]+xreal[1]*xreal[1]+xreal[2]*xreal[2]);

             EigVecs[0][0]=complex<double>(xreal[0]/normalise,zero);
             EigVecs[0][1]=complex<double>(xreal[1]/normalise,zero);
             EigVecs[0][2]=complex<double>(xreal[2]/normalise,zero);
             EigVecs[0][3]=complex<double>(one/normalise,zero); 

        // The second eigenvalue must also be real!
        for(int i=0; i<dimless; i++){
          for(int j=0; j<dimless; j++){ 
             Areal[i][j]=Mpp[i][j];       
        }
             Areal[i][i]=Areal[i][i]-real(EigVals[1]);
             breal[i]=-Mpp[i][dimless];
        }
             ireturnReal=myGaussElim(dimless, Areal, breal, xreal);
             if(ireturnReal != 0){
               cout << "myEigen44 failed in myGaussElim\n"; return ireturnReal;
             }

             normalise=sqrt(one+xreal[0]*xreal[0]+xreal[1]*xreal[1]+xreal[2]*xreal[2]);

             EigVecs[1][0]=complex<double>(xreal[0]/normalise,zero);
             EigVecs[1][1]=complex<double>(xreal[1]/normalise,zero);
             EigVecs[1][2]=complex<double>(xreal[2]/normalise,zero);
             EigVecs[1][3]=complex<double>(one/normalise,zero); 

        }

        if(!reality[0]){
//      complex conjugate pair
        for(int i=0; i<dimless; i++){
          for(int j=0; j<dimless; j++){ 
             Acomp[i][j]=complex<double>(Mpp[i][j],zero);       
        }
             Acomp[i][i]=Acomp[i][i]-EigVals[0];
             bcomp[i]=-complex<double>(Mpp[i][dimless],zero);
        }
             ireturnComp=myGaussElimZ(dimless, Acomp, bcomp, xcomp);
             if(ireturnComp != 0){
               cout << "myEigen44 failed in myGaussElimZ\n"; return ireturnComp;
             }
             double normalise;
             normalise=sqrt( one+real(xcomp[0]*conj(xcomp[0]))
                                +real(xcomp[1]*conj(xcomp[1]))
                                +real(xcomp[2]*conj(xcomp[2]))  
                           );

             EigVecs[0][0]=xcomp[0]/normalise; //eigenvectors stored as ROWS in eigvecs
             EigVecs[0][1]=xcomp[1]/normalise;
             EigVecs[0][2]=xcomp[2]/normalise;
             EigVecs[0][3]=complex<double>(one/normalise,zero);
             EigVecs[1][0]=conj(EigVecs[0][0]);
             EigVecs[1][1]=conj(EigVecs[0][1]);
             EigVecs[1][2]=conj(EigVecs[0][2]);
             EigVecs[1][3]=conj(EigVecs[0][3]);
     
        }

/**************   Vectors three and four *******************************/

       if(reality[2]){
        for(int i=0; i<dimless; i++){
          for(int j=0; j<dimless; j++){ 
             Areal[i][j]=Mpp[i][j];       
        }
             Areal[i][i]=Areal[i][i]-real(EigVals[2]);
             breal[i]=-Mpp[i][dimless];
        }
             ireturnReal=myGaussElim(dimless, Areal, breal, xreal);
             if(ireturnReal != 0){
               cout << "myEigen44 failed in myGaussElim\n"; return ireturnReal;
             }

             normalise=sqrt(one+xreal[0]*xreal[0]+xreal[1]*xreal[1]+xreal[2]*xreal[2]);

             EigVecs[2][0]=complex<double>(xreal[0]/normalise,zero);
             EigVecs[2][1]=complex<double>(xreal[1]/normalise,zero);
             EigVecs[2][2]=complex<double>(xreal[2]/normalise,zero);
             EigVecs[2][3]=complex<double>(one/normalise,zero); 

        // The second eigenvalue must also be real!
        for(int i=0; i<dimless; i++){
          for(int j=0; j<dimless; j++){ 
             Areal[i][j]=Mpp[i][j];       
        }
             Areal[i][i]=Areal[i][i]-real(EigVals[3]);
             breal[i]=-Mpp[i][dimless];
        }
             ireturnReal=myGaussElim(dimless, Areal, breal, xreal);
             if(ireturnReal != 0){
               cout << "myEigen44 failed in myGaussElim\n"; return ireturnReal;
             }

             normalise=sqrt(one+xreal[0]*xreal[0]+xreal[1]*xreal[1]+xreal[2]*xreal[2]);

             EigVecs[3][0]=complex<double>(xreal[0]/normalise,zero);
             EigVecs[3][1]=complex<double>(xreal[1]/normalise,zero);
             EigVecs[3][2]=complex<double>(xreal[2]/normalise,zero);
             EigVecs[3][3]=complex<double>(one/normalise,zero); 

        }



        if(!reality[2]){
//      complex conjugate pair
        for(int i=0; i<dimless; i++){
          for(int j=0; j<dimless; j++){ 
             Acomp[i][j]=complex<double>(Mpp[i][j],zero);       
        }
             Acomp[i][i]=Acomp[i][i]-EigVals[2];
             bcomp[i]=-complex<double>(Mpp[i][dimless],zero);
        }
             ireturnComp=myGaussElimZ(dimless, Acomp, bcomp, xcomp);
             if(ireturnComp != 0){
               cout << "myEigen44 failed in myGaussElimZ\n"; return ireturnComp;
             }
             double normalise;
             normalise=sqrt( one+real(xcomp[0]*conj(xcomp[0]))
                                +real(xcomp[1]*conj(xcomp[1]))
                                +real(xcomp[2]*conj(xcomp[2]))  
                           );

             EigVecs[2][0]=xcomp[0]/normalise; //eigenvectors stored as ROWS in eigvecs
             EigVecs[2][1]=xcomp[1]/normalise;
             EigVecs[2][2]=xcomp[2]/normalise;
             EigVecs[2][3]=complex<double>(one/normalise,zero);
             EigVecs[3][0]=conj(EigVecs[2][0]);
             EigVecs[3][1]=conj(EigVecs[2][1]);
             EigVecs[3][2]=conj(EigVecs[2][2]);
             EigVecs[3][3]=conj(EigVecs[2][3]);
     
        }




        if(somereal){
        for(int i=0; i<dimless;i++){
           delete[] Areal[i];
        }
        delete[] Areal;
        }
        if(!allreal){
        for(int i=0; i<dimless;i++){
           delete[] Acomp[i];
        }
        delete[] Acomp;
        }

        delete[] bcomp;
        delete[] xcomp;
        delete[] breal;
        delete[] xreal;




     return ireturn;
}



