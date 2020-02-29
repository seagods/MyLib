#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>

#include <gmp.h>
#include <gmpxx.h>
//We are using C++ interface for GMP  ---  no need for mpz_clear(x) mpz_ini(x) etc

int base=10;
int defprec=256;

mpz_class ZeroZ("0");
mpz_class OneZ("1");

mpq_class ZeroQ("0/1");
mpq_class OneQ("1/1");

mpf_class ZeroF("0.0");
mpf_class OneF("1.0");

int maxintbits=defprec; //if denom of properfrac more than 256 bits jump to floating point

using namespace::std;

//  non member functions (NEVER defined in class body) Here if any
// classes with friend functions must be forward declared, and so must the friend function! Here if any


class GMPfrac;
//begin class body
class  GMPfrac{
   public:

      //frac=a+b/c=B/c
       
      GMPfrac(); //constructor defalult = 0+0/1=0/1.
      GMPfrac(mpz_class&, mpz_class&, bool);  //construct from B/c and give it a sign
      GMPfrac(mpz_class&, mpz_class&, mpz_class&, bool);  //construct from a+b/c and give it a sign

      GMPfrac(const GMPfrac&);  //copy constructor by reference
      GMPfrac(const GMPfrac*);  //copy constructor by pointer

      void* operator new(size_t); 
      void* operator new(size_t, mpz_class&, mpz_class&, bool); 
      void* operator new(size_t, mpz_class&, mpz_class&, mpz_class&, bool); 
      void* operator new(size_t, mpz_class&, mpq_class&, bool); 
      void* operator new(size_t, mpq_class&, bool); 

      void  operator delete(void*);

      GMPfrac operator=(const GMPfrac&);
      GMPfrac operator=(const GMPfrac*);

      GMPfrac operator+(const GMPfrac&);
      GMPfrac operator+(const GMPfrac*);

      GMPfrac operator-(const GMPfrac&);
      GMPfrac operator-(const GMPfrac*);

      GMPfrac operator*(const GMPfrac&);
      GMPfrac operator*(const GMPfrac*);

      GMPfrac operator/(const GMPfrac&);
      GMPfrac operator/(const GMPfrac*);


//    more member functions
      void Show() const;

      bool GetSign(){ return sign; }
      void SetSign(bool sgn){sign=sgn; }

      bool GetUseFrac(){return usefrac;}
      
      mpf_class ExtractFloat(); // return float(a+b/c) or floatval

      void JumpToFloat();  //a+b/c -->  float value

      void Negate(){ if(sign==true){sign=false;} else {sign=true; }  }

      void SetFrac(mpz_class&, mpz_class&, mpz_class&, bool); //(a+b/c)
      void SetFrac(mpz_class&, mpz_class&, bool); //(B/c)

      void SetFrac(mpz_class&, mpq_class&, bool);  //a+gmp rational
      void SetFrac(mpq_class&, bool);  //gmp rational


      mpz_class* GetWhole(){return whole;}
      mpq_class* GetProper(){return properfrac;}

      bool NotZero();
   
     ~GMPfrac();


  private:
// BEST assign  pointers here, they will contain convincing looking garbage stack address otherwise
// them to
   mpz_class *whole=NULL; mpq_class *properfrac=NULL; mpf_class *floatval=NULL; 

   bool usefrac; bool sign;  //if  usefrac false use floating point

};  //end class body


#include "GMPfracConstr.h"      //GMPfrac Constuctors
#include "GMPfracOpsRef.h"      //GMPfrac operators rhs = reference
#include "GMPfracOpsPtr.h"      //GMPfrac operators rhs = pointer
#include "GMPfracNew.h"      //GMPfrac ooverloaded new and delete

void GMPfrac::Show()const { 

/*

     cout << "GMP frac Show() function\n";

     cout <<" usefrac and sign are " << this->usefrac << " and  "  << this->sign << endl;
     cout << "whole ptr=" << this->whole  << "    ";
     if(whole != NULL){cout << "whole part=" <<*whole << endl;} else { cout << endl; }
     cout << "properfrac ptr=" << this->whole  << "    ";
     if(properfrac != NULL){ cout << "properfrac part=" <<*properfrac << endl;} else { cout << endl;}

     cout << "floatval ptr=" << this->floatval  << "    ";
     if(floatval != NULL){cout << "float value=" <<*floatval << endl;} else{ cout << endl;}
*/

     if(this->usefrac){
             cout << "**********************************************\n";

             cout << "Show() fration GMPfrac=" << sign << "  (" << *whole  <<  "  +  " << *properfrac <<  ")" << endl;
             cout << "**********************************************\n";
     }
     else{
             cout << "Show() Floating Point GMPfrac="  << *floatval  <<  endl;

     }
     
     
}

mpf_class GMPfrac::ExtractFloat(){

   mpf_class returnval=0;

   if(usefrac){
      if(whole==NULL || properfrac==NULL){cout << "invalid fraction in ExtractFloat -exit now\n"; exit(1);}
      returnval=mpf_class(*whole)+mpf_class(*properfrac);
      if(!sign)returnval=-returnval;
   }
   else{
      if(floatval==NULL){cout << "invalid float in ExtractFloat -exit now\n"; exit(1);}
      returnval=*floatval;
   }
   return returnval;  
}

void GMPfrac::JumpToFloat(){
  if(!usefrac){cout <<"Invalid Jump to Float -exit now!\n";  exit(1); }  
     if(whole==NULL || properfrac==NULL){cout << "invalid fraction in JumpToFloat -exit now\n"; exit(1); }
     if(floatval!=NULL){cout << "Pre existing float in JumpToFloat ---exit now!\n";  exit(1); }

    floatval=new mpf_class;

    *floatval=mpf_class(*whole)+mpf_class(*properfrac);
    if(!sign){ *floatval=-*floatval; }

    delete whole; delete properfrac;
    whole=NULL; properfrac=NULL;
    usefrac=false;

}

void GMPfrac::SetFrac(mpz_class &aval, mpz_class &bval, mpz_class &cval, bool sgn){

  if(cval==0){cout <<"division by zero in constructor from a+b/c - exit now!\n"; exit(1);}

  if(aval<0){cout <<"all signs must be +ve in SetFrac from a+b/c - exit now!\n"; exit(1);}
  if(bval<0){cout <<"all signs must be +ve in SetFrac from a+b/c - exit now!\n"; exit(1);}
  if(cval<0){cout <<"all signs must be +ve in SetFrac from a+b/c - exit now!\n"; exit(1);}

  if(bval>=cval){
      mpz_class addwhole=bval/cval;
      aval=aval+addwhole;
      bval=bval%cval;
     }

  if(usefrac){  

     if(whole==NULL){cout << "(usefrac=true) whole doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(properfrac==NULL){cout <<"usefrac=true) properfrac doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(floatval!=NULL){cout <<"fusefrac=true) floatval exists in SetFrac ---exit now!\n"; exit(1);}

     mpq_class bslashc;
     mpz_class& numer=bslashc.get_num();  mpz_class& denom=bslashc.get_den();
     numer=bval; denom=cval;
     *whole=aval;  *properfrac=bslashc; sign=sgn;
     properfrac->canonicalize();
  }
  else{ //setfrac will override the floating point version

     if(whole!=NULL){cout <<"(usefrac=false) whole exists in SetFrac ---exit now! \n"; exit(1);}
     if(properfrac!=NULL){cout <<"(usefrac=false) properfrac  exists in SetFrac ---exit now!\n"; exit(1);}
     if(floatval==NULL){cout <<"(usefrac=false) floatval doesn't exist in SetFrac ---exit now!\n"; exit(1);}

     whole=new mpz_class; properfrac=new mpq_class;

     mpq_class bslashc;
     mpz_class& numer=bslashc.get_num();  mpz_class& denom=bslashc.get_den();
     numer=bval; denom=cval;

     *whole=aval;  *properfrac=bslashc; sign=sgn;
     properfrac->canonicalize();

     usefrac=true;
     delete floatval;
     floatval=NULL;
  }

}

void GMPfrac::SetFrac(mpz_class &Bval, mpz_class &cval, bool sgn){

  if(cval==0){cout <<"division by zero in constructor from B/c - exit now!\n"; exit(1);}
  if(Bval<0){cout <<"all signs must be +ve in SetFrac from B/c - exit now!\n"; exit(1);}
  if(cval<0){cout <<"all signs must be +ve in SetFrac from B/c - exit now!\n"; exit(1);}

  *whole=ZeroZ;

  mpz_class aval=ZeroZ;
  mpz_class bval=ZeroZ;
  aval=Bval/cval;
  bval=Bval%cval;

  if(usefrac){  

     if(whole==NULL){cout << "(usefrac=true) whole doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(properfrac==NULL){cout <<"usefrac=true) properfrac doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(floatval!=NULL){cout <<"fusefrac=true) floatval exists in SetFrac ---exit now!\n"; exit(1);}

     mpq_class bslashc;
     mpz_class& numer=bslashc.get_num();  mpz_class& denom=bslashc.get_den();
     numer=bval; denom=cval;
     *whole=aval;  *properfrac=bslashc; sign=sgn;
     properfrac->canonicalize();
  }
  else{ //setfrac will override the floating point version

     if(whole!=NULL){cout <<"(usefrac=false) whole exists in SetFrac ---exit now! \n"; exit(1);}
     if(properfrac!=NULL){cout <<"(usefrac=false) properfrac  exists in SetFrac ---exit now!\n"; exit(1);}
     if(floatval==NULL){cout <<"(usefrac=false) floatval doesn't exist in SetFrac ---exit now!\n"; exit(1);}

     whole=new mpz_class; properfrac=new mpq_class;

     mpq_class bslashc;
     mpz_class& numer=bslashc.get_num();  mpz_class& denom=bslashc.get_den();
     numer=bval; denom=cval;

     *whole=aval;  *properfrac=bslashc; sign=sgn;
     properfrac->canonicalize();

     usefrac=true;
     delete floatval;
     floatval=NULL;
  }

}



void GMPfrac::SetFrac(mpz_class &setwhole, mpq_class &setproperfrac, bool sgn){

     mpz_class& numer=setproperfrac.get_num();  mpz_class& denom=setproperfrac.get_den();

     if(numer>=denom){
       numer=numer%denom;
       setwhole=setwhole+numer/denom;
       setproperfrac.canonicalize();
     }

  if(usefrac){  

     if(whole==NULL){cout << "(usefrac=true) whole doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(properfrac==NULL){cout <<"usefrac=true) properfrac doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(floatval!=NULL){cout <<"fusefrac=true) floatval exists in SetFrac ---exit now!\n"; exit(1);}


     *whole=setwhole;
     *properfrac=setproperfrac;
     sign=sgn;
  }
  else{  //setfrac will override the floating point version

     if(whole!=NULL){cout <<"(usefrac=false) whole exists in SetFrac ---exit now! \n"; exit(1);}
     if(properfrac!=NULL){cout <<"(usefrac=false) properfrac  exists in SetFrac ---exit now!\n"; exit(1);}
     if(floatval==NULL){cout <<"(usefrac=false) floatval doesn't exist in SetFrac ---exit now!\n"; exit(1);}

     whole=new mpz_class; properfrac=new mpq_class;

     *whole=setwhole;
     *properfrac=setproperfrac;
     sign=sgn;

     usefrac=true;
     delete floatval;
     floatval=NULL;
  }

}



void GMPfrac::SetFrac(mpq_class &setproperfrac, bool sgn){

     mpz_class& numer=setproperfrac.get_num();  mpz_class& denom=setproperfrac.get_den();

     mpz_class setwhole=ZeroZ;

     if(numer>=denom){
       numer=numer%denom;
       setwhole=setwhole+numer/denom;
       setproperfrac.canonicalize();
     }

  if(usefrac){  

     if(whole==NULL){cout << "(usefrac=true) whole doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(properfrac==NULL){cout <<"usefrac=true) properfrac doesn't exist in SetFrac ---exit now!\n"; exit(1);}
     if(floatval!=NULL){cout <<"fusefrac=true) floatval exists in SetFrac ---exit now!\n"; exit(1);}


     *whole=setwhole;
     *properfrac=setproperfrac;
     sign=sgn;
  }
  else{  //setfrac will override the floating point version

     if(whole!=NULL){cout <<"(usefrac=false) whole exists in SetFrac ---exit now! \n"; exit(1);}
     if(properfrac!=NULL){cout <<"(usefrac=false) properfrac  exists in SetFrac ---exit now!\n"; exit(1);}
     if(floatval==NULL){cout <<"(usefrac=false) floatval doesn't exist in SetFrac ---exit now!\n"; exit(1);}

     whole=new mpz_class; properfrac=new mpq_class;

     *whole=setwhole;
     *properfrac=setproperfrac;
     sign=sgn;

     usefrac=true;
     delete floatval;
     floatval=NULL;
  }

}

bool GMPfrac::NotZero(){


  bool NotZ=true;

  if(usefrac){

    if(whole==NULL || properfrac==NULL){ cout << "Bad input to NotZero  -exit now!\n";  exit(1); }
    mpz_class a, b, c;

    mpq_class F;

    a=*whole;

    mpz_class& numer=properfrac->get_num();
    mpz_class& denom=properfrac->get_den();

    b=numer;
    c=denom;

    cout << "a=" << a << endl;
    cout << "b=" << b << endl;
    cout << "c=" << c << "   ZeroZ=" << ZeroZ << endl;
    cout << "b/c=" << *properfrac << endl;

    if(c==ZeroZ){ "c=0 in NotZero  -exit now!\n"; exit(1); }


    if(a==ZeroZ &&  b==ZeroZ)NotZ=false;

    return NotZ;

  }

   //must be float

   if(floatval==NULL){ cout << "Bad input to NotZero  -exit now!\n";  exit(1); }

   if(*floatval==ZeroF)NotZ=false;

   return NotZ;

  

}







