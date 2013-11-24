#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>

using namespace::std;
// Given N roots of a monic (X^N+a X^N-1+b X^N-2+...+C) polynomial we can recover all the coefficients
// X[N]  = N roots
// coeffs[N+1] = coefficients ALTERNATING IN SIGN.  i.e coeffs[0]=1, coeffs[1]=-a, coeffs[2]=+b etc.
// set verbose=false unless you want to check out the results of picking n objects from N in Pickem
int myPickem(int** picks, int& m, int& N, int & ifacs);

int myRecoverMonic(double* X, double* coeffs, int &N){

 //   bool verbose=false;

    int ireturn;
    ireturn=0;

    double mult=1.0;
   
    coeffs[0]=1.0;   //coeff of x^N=1.0
    coeffs[N]=1.0;  // coeff of X^0 is to be product of all roots
    for(int i=0; i<N; i++){
       coeffs[N]=coeffs[N]*X[i];
    }

    int ifacs;

    int Nend;   // we only have to do half the work
    Nend=N/2;   // the second part of the list is the first half - except the numbers are
                // to be regarded as missing instead of present 
    int iremain=N%2;


    for(int i=0; i<Nend; i++){  //begin main loop over i

       ifacs=1;
       int nobj=i+1;
       // pick i+1 out of N objects
        
       for(int j=0; j<nobj; j++){
      //    simple to prove j+1 always goes exactly...
            ifacs=ifacs*(N-j)/(j+1);
       }
    //   if(verbose)cout << "i=" << i <<  "  ifacs=" << ifacs << endl;

   //    if(verbose)cout <<" declare picks " <<  ifacs << "  times " << nobj << endl;

       int** picks=new int*[ifacs];


       for(int k=0; k< ifacs; k++){
           picks[k]=new int[nobj];
       }

        // we only want <= m/2
       //if >m/2 have what we picked as missing from rather than present!
       ireturn=myPickem(picks, nobj, N, ifacs);

       bool mask[ifacs][N];
       for(int k=0; k< ifacs; k++){
         for(int l=0; l<N; l++){
            mask[k][l]=true;
         }
       } 
         coeffs[nobj]=0;
    //   if(verbose)cout << " for n=" << nobj << " I can pick\n";
       for(int k=0; k< ifacs; k++){
         mult=1.0;
         for(int l=0; l<nobj; l++){
           //  if(verbose)cout << picks[k][l] << "  ";
             mask[k][picks[k][l]]=false;
             mult=mult*X[picks[k][l]];
         }
       //  if(verbose)cout <<   endl;
         coeffs[nobj]=coeffs[nobj]+mult;
       } 

       if(i<Nend-1+iremain){
         coeffs[N-nobj]=0;
    //   if(verbose)cout << " for n=" << N-nobj << " I can pick\n";
       for(int k=0; k< ifacs; k++){
         mult=1.0;
         for(int l=0; l<N; l++){
            if(mask[k][l]){
                   mult=mult*X[l];
            }
         }
      //   if(verbose)cout << endl;
          coeffs[N-nobj]=coeffs[N-nobj]+mult;
       }
       } // end i loop



       for(int k=0; k< ifacs; k++){
          delete[] picks[k];
       } 
       delete[] picks; 
       
   
     }  


    return ireturn;

}


