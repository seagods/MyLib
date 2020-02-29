//Begin Overloaded arithmetic operators

/**********************************************************************************/
/********************* Begin = operators  ****************************************/
/*********************************************************************************/

//assignment operators, lhs exists and must be initialised beforhand call to assignments
GMPfrac GMPfrac::operator=(const GMPfrac& rhs){

  //can't just use copy constructor for Temp GMPfrac because of const rhs
  // remember  const member functions cannot change any member variables as well


  bool rhsfrac=rhs.usefrac;
  bool rhssgn=rhs.sign;

  cout <<"this whole and whole are " << this->whole << "  and  "  <<  whole << endl;
  cout <<"this proper and proper are " << this->properfrac << "  and  "  << properfrac << endl;

  if(this->usefrac && rhsfrac){
      if(this->whole==NULL ||  this->properfrac==NULL || this->floatval !=NULL){
         cout <<" error in GMPfrac assignment 1a -exit now!\n"; exit(1);  }
      if(rhs.whole==NULL ||  rhs.properfrac==NULL || rhs.floatval != NULL){
         cout <<" error in GMPfrac assignment 1b-exit now!\n"; exit(1);  }

       *(this->whole)=*(rhs.whole); 
       *(this->properfrac)=*(rhs.properfrac);
       this->sign=rhs.sign;
       this->usefrac=rhs.usefrac;
  }



  if(this->usefrac && !rhsfrac){
      if(this->whole==NULL ||  this->properfrac==NULL || this->floatval !=NULL){
         cout <<" error in GMPfrac assignment 2a -exit now!\n"; exit(1);  }
      if(rhs.whole!=NULL ||  rhs.properfrac!=NULL || rhs.floatval==NULL){
         cout <<" error in GMPfrac assignment 2b-exit now!\n"; exit(1);  }

      delete this->whole; delete this->properfrac;
      this->whole=NULL; this->properfrac=NULL;

      this->usefrac=false;
      this->floatval=new mpf_class;

      *(this->floatval)=*(rhs.floatval);

      this->sign=rhs.sign;
      this->usefrac=rhs.usefrac;
  }

  if(!(this->usefrac) && rhsfrac){
      if(this->whole!=NULL ||  this->properfrac!=NULL || this->floatval==NULL){
         cout <<" error in GMPfrac assignment 3a -exit now!\n"; exit(1);  }
      if(rhs.whole==NULL ||  rhs.properfrac==NULL || rhs.floatval!=NULL){
         cout <<" error in GMPfrac assignment 3b-exit now!\n"; exit(1);  }

      delete floatval;
      floatval=NULL; 

      this->whole=new mpz_class; this->properfrac=new mpq_class;
      *(this->whole)=*(rhs.whole);  *(this->properfrac)=*(rhs.properfrac);

      this->sign=rhs.sign;
      this->usefrac=rhs.usefrac;

  }
  if(!(this->usefrac) && !rhsfrac){
      if(this->whole!=NULL ||  this->properfrac!=NULL || this->floatval==NULL){
         cout <<" error in GMPfrac assignment 4a -exit now!\n"; exit(1);  }
      if(rhs.whole!=NULL ||  rhs.properfrac!=NULL || rhs.floatval==NULL){
         cout <<" error in GMPfrac assignment 4b-exit now!\n"; exit(1);  }


      *(this->floatval)=*(rhs.floatval);

      this->sign=rhs.sign;
      this->usefrac=rhs.usefrac;

  }


return *this;
}


/**********************************************************************/
/*****               begin + operator ref rhs     *********************/
/**********************************************************************/

GMPfrac GMPfrac::operator+(const GMPfrac& rhs){

  bool usef1, usef2, sgn1, sgn2;

  usef1=this->usefrac; usef2=rhs.usefrac;
  sgn1=this->sign;  sgn2=rhs.sign;

  if(usef1 && usef2){  //we are adding two fractions
     if(this->whole==NULL || rhs.whole==NULL){
         cout <<" error in GMPfrac operator 1-exit now!\n"; exit(1);  }
     if(this->properfrac==NULL || rhs.properfrac==NULL){
         cout <<" error in GMPfrac operator 2-exit now!\n"; exit(1);  }

/***************************************************************************************/
/*************************   START ADDITION    *****************************************/
/***************************************************************************************/



     if(sgn1==sgn2){

       GMPfrac Result(this);  //copy to temporary fraction
       //Result has new whole and properfrac and they are initialised

       mpz_class result_whole;  mpq_class result_properfrac;

       result_whole=*(this->whole)+*(rhs.whole);
       result_properfrac=*(this->properfrac)+*(rhs.properfrac);
         

       mpz_class& numer=result_properfrac.get_num();
       mpz_class& denom=result_properfrac.get_den();

/*
       denom=denom*denom;  //test jump to float
       denom=denom*denom;
       denom=denom*denom;
*/

       if(numer>=denom){ numer=numer-denom; result_whole+=1; result_properfrac.canonicalize(); }


       mpz_t dee1;  mpz_init(dee1);  //need an mpz_t denom value  for getting the number of "limbs" in denom
       mpz_set(dee1, denom.get_mpz_t() );  //sets d1 from denom using get_mpz()

       cout << "dee1 is " << mpz_size(dee1) << " 64 bit words (limbs) long\n"<< endl;

       *(Result.whole)=result_whole;
       *(Result.properfrac)=result_properfrac;
       Result.sign=sgn1;

       if(64*mpz_size(dee1) <= maxintbits){
          Result.floatval=NULL;
          Result.usefrac=true;  
          return Result;
       }
       else{
          cout <<"jump to float mpf_class value\n";  
          Result.JumpToFloat();
          return Result;
       }


       mpz_clear(dee1);
      

    } //end sign1==sgn2


    if( (sgn1==true) && (sgn2==false) ){
       GMPfrac Result(this); 
       GMPfrac RHS(rhs);
       RHS.Negate(); //need to subtract positive copy of rhs
       Result=Result-RHS;
       return Result;
    }
    if( (sgn1==false) && (sgn2==true) ){
       GMPfrac Result(*this); 
       Result.Negate(); 
       Result=Result-rhs;
       Result.Negate();
       return Result;
    }

   //    mpz_clear(dee1);


  } //end of both usefracs true

  cout <<" cannot get here in operator+  \n";   exit(1);
}


/***********************************************************************/
/*****              begin - operator  ref rhs     *********************/
/**********************************************************************/
GMPfrac GMPfrac::operator-(const GMPfrac& rhs){

  bool usef1, usef2, sgn1, sgn2;

  bool negwhole=false; bool  negproper=false;  bool negfrac=false; bool zerowhole=false;




  usef1=this->usefrac; usef2=rhs.usefrac;
  sgn1=this->sign;  sgn2=rhs.sign;

  if(usef1 && usef2){  //we are adding two fractions
     if(this->whole==NULL || rhs.whole==NULL){
         cout <<" error in GMPfrac operator 1-exit now!\n"; exit(1);  }
     if(this->properfrac==NULL || rhs.properfrac==NULL){
         cout <<" error in GMPfrac operator 2-exit now!\n"; exit(1);  }
     //can begin to add two fractions 

     if(sgn1==sgn2){

/***************************************************************************************/
/*************************   START SUBTRACTION *****************************************/
/***************************************************************************************/

       GMPfrac Result(this);  //copy to temporary fraction

       //Result has new whole and properfrac and they are initialised

       mpz_class result_whole;  mpq_class result_properfrac;

       result_whole=*(this->whole)-*(rhs.whole);
       result_properfrac=*(this->properfrac)-*(rhs.properfrac);


       cout <<"SUBTRACT------------------------------------------------\n";
        cout << *(rhs.properfrac) << endl;
        cout << "FROM\n";
        cout << *(this->properfrac) << endl;
        cout << "GET\n";
        cout << result_properfrac << endl;

        if(result_properfrac <ZeroQ )cout << "result properfrac LESS THAN ZERO\n";

       cout <<"SUBTRACT------------------------------------------------\n";

         
       if(result_whole < ZeroZ)negwhole=true;
       if(result_whole == ZeroZ)zerowhole=true;
       if(result_properfrac < ZeroQ)negproper=true;

       //if both !negwhole and !negpeoper, negfrac already false, do nothing.

       if(negwhole && negproper){  //eg -5-1/4 ->  -(+5+1/4)
           negfrac=true;
           result_whole=ZeroZ-result_whole;
           result_properfrac=ZeroQ-result_properfrac;
       }

       if(zerowhole && negproper){  //  eg -1/4 -> -(0+1/4)
           negfrac=true;
           result_properfrac=ZeroQ-result_properfrac;
       }
       

       if(negwhole && !negproper){ // eg -5+1/4=-(+4+3/4)
           negfrac=true;
           result_whole=OneZ+result_whole;
           result_whole=-result_whole;
           result_properfrac=OneQ-result_properfrac; //automatic canonicalize
       }
       //last case
       if(!negwhole && negproper){
           if(result_whole >= OneZ){  //eg 5-1/4 -> +(4+3/4)
             result_whole=result_whole-OneZ;  
             result_properfrac=OneQ+result_properfrac;
           }
           else{
             negfrac=true;  //whole can only be zero
             result_properfrac==ZeroQ-result_properfrac;
           }
             
       }

       *(Result.whole)=result_whole;
       *(Result.properfrac)=result_properfrac;


        if(!sgn1){Result.sign=negfrac; } else {Result.sign=!negfrac; }

       mpz_class& numer=result_properfrac.get_num();
       mpz_class& denom=result_properfrac.get_den();

/*
       denom=denom*denom;  //test jump to float
       denom=denom*denom;
       denom=denom*denom;
*/

       if(numer>=denom){ numer=numer-denom; result_whole+=1; }


       mpz_t dee1;  mpz_init(dee1);  //need an mpz_t denom value  for getting the number of "limbs" in denom
       mpz_set(dee1, denom.get_mpz_t() );  //sets dee1 from denom using get_mpz()

       cout << "dee1 is " << mpz_size(dee1) << " 64 bit words (limbs) long\n"<< endl;

       *(Result.whole)=result_whole;
       *(Result.properfrac)=result_properfrac;


       if(64*mpz_size(dee1) <= maxintbits){

          Result.floatval=NULL;
     //     Result.sign=sgn1;  //not for subtraction it aint!
          Result.usefrac=true;        
          return Result;
       }
       else{
             cout <<"jump to float mpf_class value\n";  
             Result.JumpToFloat();
             return Result;
       }

       mpz_clear(dee1);
      
    } //end sign1==sgn2

// cases sign1!=sign2

   if( (sgn1==true) && (sgn2==false) ){
      GMPfrac Result(this);
      GMPfrac RHS(rhs);
      RHS.Negate();
      Result=Result+RHS;
      return Result;
   }


   if( (sgn1==false) && (sgn2==true) ){
      GMPfrac Result(this);
      Result.Negate();
      Result=Result+rhs;
      Result.Negate();
      return Result;
   }

 } //end of both usefracs true


  cout <<" cannot get here in operator-  \n";   exit(1);
}



/***********************************************************************/
/*****              begin * operator  ref rhs     *********************/
/**********************************************************************/
GMPfrac GMPfrac::operator*(const GMPfrac& rhs){

/*
   (a1+b1/c1)*(a2+b2/c2);
   avoid large numerator potential at cost of speed
   cost is outweighed if large numerator * large numerator
   is very expensive!
*/

  bool usef1, usef2, sgn1, sgn2;

  usef1=this->usefrac; usef2=rhs.usefrac;
  sgn1=this->sign;  sgn2=rhs.sign;



  if(usef1 && usef2){  //we are multiplying two fractions
     if(this->whole==NULL || rhs.whole==NULL){
         cout <<" error in GMPfrac operator 1-exit now!\n"; exit(1);  }
     if(this->properfrac==NULL || rhs.properfrac==NULL){
         cout <<" error in GMPfrac operator 2-exit now!\n"; exit(1);  }
     //can begin to multply two fractions 


/***************************************************************************************/
/*************************   START MULTIPLICATION **************************************/
/***************************************************************************************/

       GMPfrac Result;  //default0=0+0/1

       mpz_class result_whole=*( Result.GetWhole() );  
       mpq_class result_properfrac=*( Result.GetProper() );


       mpz_class& numer=result_properfrac.get_num();
       mpz_class& denom=result_properfrac.get_den();

       mpz_class a1,a2;
       a1=*(this->whole);
       a2=*(rhs.whole);

       mpq_class prop1=*(this->properfrac);
       mpq_class prop2=*(rhs.properfrac);

       result_whole=a1*a2;  // added to later
       result_properfrac=prop1*prop2;  // added to later

      cout << "Whole starts off as " << result_whole << endl;

       mpq_class add1, add2, add3;

       add1=mpq_class(a1)*prop2;  //automatic conanical form
       add2=mpq_class(a2)*prop1;

       mpz_class adder1=ZeroZ; mpz_class adder2=ZeroZ;  //whole part=a1*a2+adder1+adder2

       //refs to numerator and denominator of add1 and add2
       mpz_class& nadd1=add1.get_num(); mpz_class& dadd1=add1.get_den();
       mpz_class& nadd2=add2.get_num(); mpz_class& dadd2=add2.get_den();

       mpz_class num1,num2,den1,den2;
       num1=nadd1; num2=nadd2; den1=dadd1, den2=dadd2;

//     is declartion of num1, num2, etc really needed?

       cout <<"adders start as " << adder1 <<"  " << adder2 << endl;

//       if(num1>den1){
       if(nadd1>dadd1){
           adder1=nadd1/dadd1;
           cout << "FIRST PART " << num1 <<  "  " << adder1 << "  " << den1 << endl;
           nadd1=nadd1-adder1*dadd1;
    //       nadd1=num1;
      //     dadd1=den1;
           add1.canonicalize();
       }

//       if(num2>den2){
        if(nadd2>dadd2){
           adder2=nadd2/dadd2;

           cout << "SECOND PART " << num2 <<  "  " << adder2 << "  " << den2 << endl;
     //      num2=num2-adder2*den2;
           nadd2=nadd2-adder2*dadd2;
       //    nadd2=num2;
         //  dadd2=den2;
           add2.canonicalize();
       }

       cout <<"adders now " << adder1 <<"  " << adder2 << endl;


       result_whole=result_whole+adder1+adder2;

       cout << "result whole now " << result_whole << endl;
       cout << "add1 and add2 are " << add1 << "   and   " << add2 << endl;

       result_properfrac=result_properfrac+add1+add2;

//     addition of 2 proper fractions to original (prop1*prop2) is anything up to 3, 


 //      mpz_class& nres=result_properfrac.get_num(); 
   //    mpz_class& dres=result_properfrac.get_den();

        //nres and dres are just numer and denom!

     //  mpz_class numerRes=nres; mpz_class denomRes=dres;
      //just use nemer and denom

       if(numer>denom){
            mpz_class addmore=numer/denom;
            result_whole=result_whole+addmore;
            numer=numer-addmore*denom;
        //    nres=numerRes;
          //  dres=denomRes;
            result_properfrac.canonicalize();
       }


       mpz_t dee1;  mpz_init(dee1);  //need an mpz_t denom value  for getting the number of "limbs" in denom
       mpz_set(dee1, denom.get_mpz_t() );  //sets dee1 from denom using get_mpz()

       cout << "dee1 is " << mpz_size(dee1) << " 64 bit words (limbs) long\n"<< endl;



       bool res_sign=true;
       if(sgn1!=sgn2)res_sign=false;

          *(Result.whole)=result_whole;
          *(Result.properfrac)=result_properfrac;
          Result.sign=res_sign;


       

     //  Result.SetFrac(result_whole, result_properfrac, res_sign);

       return Result;

       mpz_clear(dee1);

 } //end of both usefracs true


  cout <<" cannot get here in operator* \n";   exit(1);


}

/***********************************************************************/
/*****              begin / operator  ref rhs     *********************/
/********** ************************************************************/
GMPfrac GMPfrac::operator/(const GMPfrac& rhs){

/*
   (a1+b1/c1)/(a2+b2/c2)=(a1*c1_b1)/(a2+c2*b2)   *     c2/c1
   avoid large numerator potential at cost of speed
   cost is outweighed if large numerator * large numerator
   is very expensive!
*/

  bool usef1, usef2, sgn1, sgn2;

  usef1=this->usefrac; usef2=rhs.usefrac;
  sgn1=this->sign;  sgn2=rhs.sign;



  if(usef1 && usef2){  //we are multiplying two fractions
     if(this->whole==NULL || rhs.whole==NULL){
         cout <<" error in GMPfrac operator 1-exit now!\n"; exit(1);  }
     if(this->properfrac==NULL || rhs.properfrac==NULL){
         cout <<" error in GMPfrac operator 2-exit now!\n"; exit(1);  }
     //can begin to multply two fractions 

     //leave jump to frac to call on * operator for now
//     mpz_t dee1;  mpz_init(dee1);  //need an mpz_t denom value  for getting the number of "limbs" in denom

/***************************************************************************************/
/*************************   START DIVISION **************************************/
/***************************************************************************************/

       GMPfrac Result;  //default0=0+0/1

       mpz_class result_whole=*( Result.GetWhole() );  
       mpq_class result_properfrac=*( Result.GetProper() );


       mpz_class a1,a2;
       a1=*(this->whole);
       a2=*(rhs.whole);

       mpq_class prop1=*(this->properfrac);
       mpq_class prop2=*(rhs.properfrac);

  /*     mpz_class& nfrac1=prop1.get_num(); 
       mpz_class& nfrac2=prop2.get_num();
       mpz_class& dfrac1=prop1.get_den(); 
       mpz_class& dfrac2=prop2.get_den();
*/
       mpz_class& b1=prop1.get_num(); 
       mpz_class& b2=prop2.get_num();
       mpz_class& c1=prop1.get_den(); 
       mpz_class& c2=prop2.get_den();
    //   b1=nfrac1; c1=dfrac1; b2=nfrac2, c2=dfrac2;

       mpq_class part1, part2;
      
       mpz_class& np1=part1.get_num();
       mpz_class& np2=part2.get_num();
       mpz_class& dp1=part1.get_den();
       mpz_class& dp2=part2.get_den();


//     just use np1, dp1 etc
    //   mpz_class nump1=a1*c1+b1;  mpz_class denp1=a2*c2+b2;
     //  mpz_class nump2=c2;        mpz_class denp2=c1;
       np1=a1*c1+b1;  dp1=a2*c2+b2;
       np2=c2;        dp2=c1;

       part1.canonicalize();
       part1.canonicalize();

   //    np1=nump1; dp1=denp1; part1.canonicalize();
  //     np2=nump2; dp2=denp2; part2.canonicalize();

   //    nump1=np1; denp1=dp1; nump2=np2; denp2=dp2;  //(canonicalise may have changed them)

       GMPfrac Mult1(np1,dp1,sgn1); GMPfrac Mult2(np2,dp2,sgn2);

       Result=Mult1*Mult2;

       return Result;
   //    mpz_clear(dee1);
 } //end of both usefracs true


  cout <<" cannot get here in operator* \n";   exit(1);


}






