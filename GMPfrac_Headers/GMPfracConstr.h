//default constuctor
GMPfrac::GMPfrac(){

  cout << "Constructing default GMPfrac\n";
  usefrac=true;  //that is -- use fraction a+b/c not floating point

  whole=new mpz_class;  properfrac=new mpq_class;  floatval=NULL;
  //whole and properfrac have have basic memory assigned but NOT initialised
  //floatval is not initialised and has no memory assigned to it
  

 if(whole==NULL || properfrac==NULL){
         cout << "error in GMPfrac constructor -exit now!\n"; exit(1); }

  *whole=0; *properfrac="0/1"; 
  sign=true; //that is sign means positive  

}                                                  
//end default constructor

//default destructor
GMPfrac::~GMPfrac(){   

    cout << "destructor call for GMPfrac\n";

    if(usefrac){ // a fraction
       if(whole==NULL || properfrac==NULL){ 
               cout << "error 1 in GMPfrac destructor -exit now!\n"; exit(1); }
       delete whole; delete properfrac;
       whole=NULL; properfrac=NULL; 
    }
    else{ // a floating point value has had memory allocated to it 
      if(floatval==NULL){ 
                 cout << "error 2 in GMPfrac destructor -exit now!\n"; exit(1); }
      delete floatval;
      floatval=NULL;
      
    }
    //No safety belts 
    //all routines MUST ensure either/or for fraction and floatval
    // and MUST ensure BOTH whole and properfrac exist (or don't).   
}
//end default destructor




GMPfrac::GMPfrac(mpz_class &Bval ,mpz_class &cval, bool signval){           //constructor from B/c and sign

     if(cval==0){cout <<"division by zero in constructor from B/c - exit now!\n"; exit(1);}
     if(Bval<0){cout <<"all signs must be +ve in onstructor from B/c - exit now!\n"; exit(1);}
     if(cval<0){cout <<"all signs must be +ve in constructor from B/c - exit now!\n"; exit(1);}

     usefrac=true;
     sign=signval;

     whole=new mpz_class; properfrac=new mpq_class;

     cout << "Constructing from B/c \n";
     cout <<"Bval=" << Bval << endl;
     cout <<"cval=" << cval << endl;

     *whole=Bval/cval;
    
     mpz_class bval=Bval%cval;

     mpz_class& numer=properfrac->get_num();
     mpz_class& denom=properfrac->get_den();

     numer=bval; denom=cval;
     cout <<"fraction before canonicalize=" << *whole << " + " << *properfrac << endl;

     properfrac->canonicalize();

     cout <<"fraction=" << *whole << " + " << *properfrac << endl;
    
}                                                  
// end construct from B/c
GMPfrac::GMPfrac(mpz_class &aval, mpz_class &bval ,mpz_class &cval, bool signval){    //constructor from a+b/c and sign


     if(cval==0){cout <<"division by zero in constructor from a+b/c - exit now!\n"; exit(1);}

     if(aval<0){cout <<"all signs must be +ve in onstructor from a+b/c - exit now!\n"; exit(1);}
     if(bval<0){cout <<"all signs must be +ve in onstructor from a+b/c - exit now!\n"; exit(1);}
     if(cval<0){cout <<"all signs must be +ve in constructor from a+b/c - exit now!\n"; exit(1);}

     usefrac=true;
     sign=signval;

     whole=new mpz_class; properfrac=new mpq_class;

     cout << "Constructing from a+b/c\n";

     cout <<"aval=" << aval << endl;
     cout <<"bval=" << bval << endl;
     cout <<"cval=" << cval << endl;

     *whole=aval;
     mpz_class newbval=bval;

     mpz_class addwhole=bval/cval;
     if(addwhole>0){
        *whole+=addwhole;
        newbval=bval%cval;
     }

     mpz_class& numer=properfrac->get_num();
     mpz_class& denom=properfrac->get_den();

     numer=newbval; denom=cval;

     cout <<"fraction before canonicalize=" << *whole << " + " << *properfrac << endl;

     properfrac->canonicalize();

     cout << "fraction=" << *whole << " + " << *properfrac << endl;
    
}                                                   // end construct from a+b/c

GMPfrac::GMPfrac(const GMPfrac& rhs){  //copy constructors



    usefrac=rhs.usefrac;
    sign=rhs.sign;

    if(usefrac){
       if(rhs.whole==NULL || rhs.properfrac==NULL){"error 1 in copy from ref to GMPfrac- exit now!\n"; exit(1);}
        whole=new mpz_class; properfrac=new mpq_class;
       *whole=*(rhs.whole);  *properfrac=*(rhs.properfrac);
    }
    else{
       cout <<"copy 1  " << rhs.floatval << endl;
       if(rhs.floatval==NULL ){"error 2 in copy from ref to GMPfrac- exit now!\n"; exit(1);}
       floatval=new mpf_class;
       *floatval=*(rhs.floatval);
        }
    }                                                  // end copy constructor given a GMPfrac reference

GMPfrac::GMPfrac(const GMPfrac* rhs){  //copy constructors

   sign=rhs->sign;
   usefrac=rhs->usefrac;

   if(rhs->usefrac){
        if(rhs->whole==NULL || rhs->properfrac==NULL){"error 1 in copy from ptr to GMPfrac- exit now!\n"; exit(1);}
        whole=new mpz_class; properfrac=new mpq_class;
        *whole=*(rhs->whole);  *properfrac=*(rhs->properfrac);
   }
   else{
      if(rhs->floatval==NULL){"error 2 in copy from ref to GMPfrac- exit now!\n"; exit(1);}
       floatval=new mpf_class;
      *floatval=*(rhs->floatval); 
   }                         
}                         // end copy constructor given a GMPfrac ptr








