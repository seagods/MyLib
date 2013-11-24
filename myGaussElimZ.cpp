#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>
using namespace std;

int myGaussElimZ(int &dim, complex<double>** Mpp, complex<double>* rhsV, complex<double>* lhsV){

	int SwapVec[dim]; //initialise with [0,1,2,3,4,5,...,dim]
	int Whereis[dim]; //initialise with copy of SwapVec

	for(int i=0;i<dim;i++){
                SwapVec[i]=i; // keep tabs on swaps
                Whereis[i]=i; // recover swaps
        }

	int colPOS, LOW, HIGH, whereLOW, whereHIGH;  

        double zero=0.0e0;
        complex<double> ZERO;
        ZERO=complex<double>(zero,zero);

	complex<double> divide, multiply, factor;

        double largest;
        double largestinrow;
        int ireturn;

        ireturn=0; //we live in hope!


	int klarge;
        klarge=0;

        double absval;
        double Tol=1e-14;

        bool pivot=true;


        double scalevec[dim];
        double avalue;
        for(int i=0;i<dim;i++){
        largestinrow=0.0;
        for(int j=0;j<dim;j++){
             avalue=abs(Mpp[i][j]);
             if(avalue > largestinrow)largestinrow=avalue;
        }
        scalevec[i]=largestinrow;



           if(scalevec[i] < Tol){
             cout << "row " << i << "'s largest element is << " << largestinrow << endl;
             cout <<" either matrix is singular, or problem needs re-scaling\n";
             cout << "return with value of 6\n"; 
             ireturn=6;
             return ireturn;
                       }
        }
       
	for(int i=0; i<dim; i++){
                if(pivot){
		largest=0.0; 

                if(i<dim-1){

		     for(int k=i; k<dim; k++){
			absval=abs(Mpp[Whereis[k]][i])/scalevec[Whereis[k]];
			if(absval>largest){
			       largest=absval;
			       klarge=Whereis[k];
			}}
                      if(fabs(largest)<Tol)return(1);

		if( Whereis[klarge] != i){  
                        colPOS=Whereis[i];  //pos in swap of row i
                                            // klarge is pos in swap of largest val 
		        LOW=SwapVec[colPOS]; 
                        HIGH=SwapVec[klarge];

			SwapVec[colPOS]=HIGH; 
                        SwapVec[klarge]=LOW;

			whereHIGH=Whereis[HIGH];
			whereLOW=Whereis[LOW];

			Whereis[HIGH]=whereLOW;
			Whereis[LOW]=whereHIGH; //test
		} // Whereis[klarge] !=i

                }//i < dim-1

                }// endif pivot
                       
               int index1, index2;

                index1=Whereis[i];
                divide=Mpp[index1][i];

	        for( int j=i; j<dim; j++){
                  index2=Whereis[j];
		  multiply=Mpp[index2][i]; 

                  if(index1 != index2){ 

                      factor=multiply/divide;
           
                      for(int p=i; p<dim; p++){
                        if(p==i){
                            Mpp[index2][p]=0.0;
                        }
                        else{
                           Mpp[index2][p]-=Mpp[index1][p]*factor;
                        }
               
                      } //end for p

                      rhsV[index2]-=rhsV[index1]*factor;     
                        

	         } // end for j    


                 }  //index1 not equal to index2
           
	} //end loop for i
    




       //Back Subst.!
	complex<double> sum;
	int irow;

	for(int i=0; i<dim; i++){

	      irow=dim-1-i;
              sum=ZERO; 
	      for(int j=dim-i; j<dim; j++){
		sum=sum+Mpp[Whereis[irow]][j]*lhsV[j]; 
	      }
              lhsV[irow]=(rhsV[Whereis[irow]]-sum)/Mpp[Whereis[irow]][irow];
	}


     return ireturn;
}



