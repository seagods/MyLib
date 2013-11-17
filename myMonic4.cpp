#include <stdio.h>  //standard input output
#include <iostream> // input output (cin if(verbose)cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;

typedef complex<double> dcmplx;


//Declare function for solving resolvent cubic
int myMonic3(double*, complex<double>*, bool* , bool& );

int myMonic4(double* realcoeffs, complex<double>* roots, bool* reality, bool& allreal){

     int ireturn=0;

     double a0,a1,a2,a3;

     a0=realcoeffs[0];
     a1=realcoeffs[1];
     a2=realcoeffs[2];
     a3=realcoeffs[3];

     double tolval=1.0e-10;  //tol for root polish
     double tolvalreal=1.0e-15; //tol for deciding reality
     int iqmax=10; //maximum number for Newton Raphson

     reality[0]=false; reality[1]=false;
     reality[2]=false; reality[3]=false;

     double p,q,r; //coefficients of depressed quartic x^4+px^2+qx+r=0
                   //shift=a_3/4

     double a3sq, a3sqdiv4,a3sqdiv8;

     a3sq=a3*a3;
     a3sqdiv4=a3sq/4.0;
     a3sqdiv8=a3sq/8.0;   

//   coefficients of depressed quartic  

     p=a2-3.0*a3sqdiv8;
     q=a1-a2*a3/2.0+a3*a3sqdiv8;
     r=a0-a1*a3/4.0+a3sqdiv8*(a2/2.0-0.75*a3sqdiv8);

      //resolvant cubic for depressed quartic

     double cubic_coeffsDep[3];
     dcmplx cubic_rootsDep[3];

     bool cubic_allrealDep=false;
     bool cubicRealDep[3];
     cubicRealDep[0]=false; cubicRealDep[1]=false; cubicRealDep[2]=false;

     cubic_coeffsDep[0]=4.0*r*p-q*q;
     cubic_coeffsDep[1]=-4.0*r;
     cubic_coeffsDep[2]=-p;

     int icube=0;

     icube=myMonic3(cubic_coeffsDep, cubic_rootsDep, cubicRealDep, cubic_allrealDep);


     if(icube !=0)ireturn=icube;

     dcmplx Uval,UmP,rootUmP;

     Uval=dcmplx(0.0,0.0);
     UmP=dcmplx(0.0,0.0);
     rootUmP=dcmplx(0.0,0.0);

     int iump=0;
     if(cubic_allrealDep){
       double max_ump=0.0;

       if(abs(cubic_rootsDep[0]-p) > max_ump)iump=0;
       if(abs(cubic_rootsDep[1]-p) > max_ump)iump=1;
       if(abs(cubic_rootsDep[2]-p) > max_ump)iump=2;
     }
     
     
     Uval=cubic_rootsDep[iump];
     UmP=Uval-p;   
     rootUmP=sqrt(UmP);

     dcmplx bee1,bee2,cee1,cee2;
 //    C++ constuctor class initialises to (0.0,0.0) anyway!
 //    bee1=dcmplx(0.0,0.0);
 //    bee2=dcmplx(0.0,0.0);
 //    cee1=dcmplx(0.0,0.0);
 //    cee2=dcmplx(0.0,0.0);


     bee1=rootUmP;
     cee1=Uval/2.0-q/2.0/rootUmP;

     bee2=-rootUmP;
     cee2=Uval/2.0+q/2.0/rootUmP;

     roots[0]=(-bee1+sqrt( bee1*bee1-4.0*cee1))/2.0-a3/4.0;
     roots[1]=(-bee1-sqrt( bee1*bee1-4.0*cee1))/2.0-a3/4.0;
     roots[2]=(-bee2+sqrt( bee2*bee2-4.0*cee2))/2.0-a3/4.0;
     roots[3]=(-bee2-sqrt( bee2*bee2-4.0*cee2))/2.0-a3/4.0;

//   Use Newton-Raphson to "polish" the roots;

     if(abs(imag(roots[0]))<tolvalreal){
                reality[0]=true; reality[1]=true;
     }
     if(abs(imag(roots[2]))<tolvalreal){
                reality[2]=true; reality[3]=true;
     }


     bool polish;

     dcmplx testfun,testfundiff;

     double coeffsdiff[4];

     coeffsdiff[3]=4.0;
     coeffsdiff[2]=3.0*realcoeffs[3];
     coeffsdiff[1]=2.0*realcoeffs[2];
     coeffsdiff[0]=realcoeffs[1];


     dcmplx root,nextroot;
     dcmplx cube,square,powfour;

     int iq=0;
     polish=true;
     root=roots[0];

     while(polish){ 

       square=root*root;
       cube=root*square;
       powfour=root*cube;

       testfun=powfour+realcoeffs[3]*cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;



       if(polish){

         //  cout <<"qth root[0]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;


           testfundiff=coeffsdiff[3]*cube+coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic4" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic4\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[0]=root;



     if(!reality[0]){
        roots[1]=conj(roots[0]);
     }
     else{ //might need to polish root[1]


     iq=0;
     polish=true;
     root=roots[1];

     while(polish){ 

       square=root*root;
       cube=root*square;
       powfour=root*cube;

       testfun=powfour+realcoeffs[3]*cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;



       if(polish){

     //      cout <<"qth root[1]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;

           testfundiff=coeffsdiff[3]*cube+coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic4" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic4\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[1]=root;

     }  // end if else reality[0]



     iq=0;
     polish=true;
     root=roots[2];

     while(polish){ 

       square=root*root;
       cube=root*square;
       powfour=root*cube;

       testfun=powfour+realcoeffs[3]*cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;

       if(polish){

     //      cout <<"qth root[2]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;

           testfundiff=coeffsdiff[3]*cube+coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic4" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic4\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[2]=root;


     if(!reality[2]){
        roots[3]=conj(roots[2]);
     }
     else{ //might need to polish root[1]

     iq=0;
     polish=true;
     root=roots[3];

     while(polish){ 

       square=root*root;
       cube=root*square;
       powfour=root*cube;

       testfun=powfour+realcoeffs[3]*cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;



       if(polish){

       //    cout <<"qth root[3]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;

           testfundiff=coeffsdiff[3]*cube+coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic4" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic4\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[3]=root;

     }  // end if else reality[2]





     return ireturn;
}



