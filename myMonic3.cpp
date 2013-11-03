#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision

#include <complex>

using namespace std;

typedef complex<double> dcmplx;

int myMonic3(double* realcoeffs, complex<double>* roots, bool* reality, bool& allreal){

     int ireturn=0;

     double a0,a1,a2;

     double q,r,discrim;

     double root3=1.7320508075688772;
     double root3over2=0.86602540378443860;

     a0=realcoeffs[0];
     a1=realcoeffs[1];
     a2=realcoeffs[2];

     double tolval=1.0e-10;  //tol for root polish
     double tolvalreal=1.0e-15; //tol for deciding reality
     int iqmax=10; //maximum number for Newton Raphson

     double a2sq=a2*a2;

     q=a1/3.0-a2sq/9.0;
     r=(a1*a2-3.0*a0)/6.0-a2*a2sq/27.0;

     discrim=q*q*q+r*r;

     if(discrim <=0.0)allreal=true;



     double Magnit;
     double MagnitCroot; //cube root of magnitude;
     double phase1,phase2;



     if(allreal){

     complex<double>qr1,s1,qr2,s2;

 //    cout << "Discriminant=" << discrim << "  q="  << q << "  r=" << r << endl;

 //        cout << "All Roots Real \n";

         reality[0]=true;
         reality[1]=true;
         reality[2]=true;

         double imp;
          
         imp=sqrt(-discrim);  //discrim known to be negative

         qr1=complex<double>(r,imp);
         qr2=complex<double>(r,-imp);

   //     cout << " qr1=" << qr1 << "  qr2=" << qr2 << endl;

         Magnit=abs(qr1);

         //phase angle is between -pi and pi
         //(complex arg function uses atan2)
         phase1=arg(qr1);
         phase2=-phase1;

     //    cout <<"Magnitude="<< Magnit << "  phases are " << phase1 << "  and  " << phase2 << endl;

         MagnitCroot=exp(log(Magnit)/3.0);

         phase1=phase1/3.0;
         phase2=phase2/3.0;

    //     cout << "Cube root of magitude=" << MagnitCroot <<
      //           "   " <<   pow(MagnitCroot,3.0) << endl;

         s1=polar(MagnitCroot,phase1);
         s2=polar(MagnitCroot,phase2);

   //     cout << " s1=" << s1 << "  s2=" << s2 << endl;

         //since phase1=-phase2, s1+s2 and Real(s1)=Real(s2)
         //  s1+s2=2.*real(s1)

         double realsum=real(s1+s2);

         double realval1=realsum-a2/3.0;
         //s1+s2 is just 2*Re(s1)
         double realval2=-realsum/2.0-a2/3.0;
         double imagval=root3over2*imag(s1-s2);

    //     cout << "rv1=" << realval1 << "  rv2=" << realval2 << "  imv=" << imagval << endl;

         roots[0]=complex<double>( realval1 , 0.0);
         roots[1]=complex<double>( realval2+imagval, 0.0);
         roots[2]=complex<double>( realval2-imagval, 0.0);
      }

      if(!allreal){


         reality[0]=true;
         reality[1]=false;
         reality[2]=false;

        double qr1,s1,qr2,s2;
/*
        cout << "Discriminant=" << discrim << "  q="  << q << "  r=" << r << 
                 "   q^3 +r^2=" << q*q*q+r*r <<  endl;

        cout << "Root Of Discriminant=" << sqrt(discrim) <<  endl;

        cout << "Only One  Root Real \n";
*/

        double rootdiscrim;

        rootdiscrim=sqrt(discrim);  

//        cout <<"root discrim=" << rootdiscrim << "  " << 100.0/3.0/sqrt(3.0) << endl;

                              // discrim known to be positive
                             // imp not imaginary part 
                             // if allreal false

         qr1=r+rootdiscrim;  //quadratic roots one and two 
         qr2=r-rootdiscrim;

 //       cout << "roots of quadratic are " << qr1  << "  and  " << qr2 <<endl;

         if(qr1>0){
           s1=exp(log(qr1)/3.0);
         }
         else{ 
           s1=-exp(log(-qr1)/3.0);
         }

         if(qr2>0){
           s2=exp(log(qr2)/3.0);
         }
         else{ 
           s2=-exp(log(-qr2)/3.0);
         }
/*
         cout << "s1 and s2 are " << s1 << "  and " << s2 << endl; 
         cout << "s1^3 and s2^3 are " << s1*s1*s1 << "  and " << s2*s2*s2 << endl; 
         cout << "s1 + s2 =" << s1+s2 << endl; 
         cout << "s1 - s2 =" << s1-s2 << endl; 
*/

         double sum,diff;
         sum=s1+s2;
         diff=s1-s2;


//       all start off as -a2/3
//         cout << " -a2/3=" << -a2/3.0 << endl;
         roots[0]=complex<double>(-a2/3.0,0.0);
         roots[1]=roots[0];
         roots[2]=roots[0];


         roots[0]=roots[0]+(complex<double>(sum,0.0));
         roots[1]=roots[1]-(complex<double>(sum/2.0,0.0))
                   +(complex<double>(diff,0.0))*(complex<double>(0.0,root3over2));
         roots[2]=roots[2]-(complex<double>(sum/2.0,0.0))
                   -(complex<double>(diff,0.0))*(complex<double>(0.0,root3over2));


      }



  //   bool MyArse=true;

 //    if(MyArse){
 //    cout << "start polish\n";


     bool polish;

     dcmplx testfun,testfundiff;

     double coeffsdiff[3];

     coeffsdiff[2]=3.0;
     coeffsdiff[1]=2.0*realcoeffs[2];
     coeffsdiff[0]=realcoeffs[1];

     dcmplx root,nextroot;
     dcmplx cube,square;

     int iq=0;
     polish=true;
     root=roots[0];

     while(polish){ 

       square=root*root;
       cube=root*square;

       testfun=cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;



       if(polish){

       //    cout <<"qth root[0]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;


           testfundiff=coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic3" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic3\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[0]=root;


     iq=0;
     polish=true;
     root=roots[1];

     while(polish){ 

       square=root*root;
       cube=root*square;

       testfun=cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;

       if(polish){

  //         cout <<"qth root[2]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;

           testfundiff=coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic3" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic3\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[1]=root;


     if(!reality[1]){
        roots[2]=conj(roots[1]);
     }
     else{ //might need to polish root[2]

     iq=0;
     polish=true;
     root=roots[2];

     while(polish){ 

       square=root*root;
       cube=root*square;

       testfun=cube+realcoeffs[2]*square+realcoeffs[1]*root+realcoeffs[0];

       if(abs(testfun)<tolval)polish=false;



       if(polish){

      //     cout <<"qth root[3]=" << root << "  q=" << iq << "  f(z)=" << testfun <<  endl;

           testfundiff=coeffsdiff[2]*square+coeffsdiff[1]*root+coeffsdiff[0];

           if(abs(testfundiff)>tolval){
               nextroot=root-testfun/testfundiff; iq++;
               if(iq>iqmax){
                 cout << " Stuck in root polish in myMonic3" << endl; exit(0);
               }
               root=nextroot;
           }
           else{
              cout << "disaster in polishing root in myMonic3\n";
              exit(0);
           }
           
        
       }  // end if polish
     }  //end while polish
     roots[2]=root;

     }  // end if else reality[1]

 //    } //endif MyArse


     return ireturn;
}



