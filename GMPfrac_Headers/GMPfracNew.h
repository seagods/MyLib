//         THIS WAY  ::operator new and ::operator delete

void* GMPfrac::operator new(size_t size){
      void* FracStore= ::operator new(size);
      cout << "size of FracStore=" << size << endl;
      cout <<"Fraction Stored at " << FracStore << endl;
      return FracStore;
}


void  GMPfrac::operator delete(void*  FracStore){
      cout << "freeing up FracStore memory at " << FracStore << endl;
      ::operator delete(FracStore);
}


void* GMPfrac::operator new(size_t size, mpz_class& Bval, mpz_class& cval, bool signval){
      void* FracStore= ::operator new(size);
      cout << "size of FracStore=" << size << endl;
      cout <<"Fraction Stored at " << FracStore << endl;

      GMPfrac* FracPtr = NULL;
      FracPtr=static_cast< GMPfrac*>(FracStore);
     //need FracPtr cast to GMPfrac* in order to put data into FracStore memory!
     //compile now knows what usefrac and sign are....

     FracPtr->usefrac=true;

     FracPtr->SetFrac(Bval,cval,signval);
     
    return FracStore;

}
void* GMPfrac::operator new(size_t size, mpz_class& aval, mpz_class& bval, mpz_class& cval, bool signval){
      void* FracStore= ::operator new(size);
      cout << "size of FracStore=" << size << endl;
      cout <<"Fraction Stored at " << FracStore << endl;

      GMPfrac* FracPtr = NULL;
      FracPtr=static_cast< GMPfrac*>(FracStore);
     //need FracPtr cast to GMPfrac* in order to put data into FracStore memory!
     //compile now knows what usefrac and sign are....

     FracPtr->usefrac=true;

     FracPtr->SetFrac(aval, bval, cval,signval);
     
    return FracStore;

}
void* GMPfrac::operator new(size_t size, mpz_class& aval, mpq_class& setproperfrac, bool signval){
      void* FracStore= ::operator new(size);
      cout << "size of FracStore=" << size << endl;
      cout <<"Fraction Stored at " << FracStore << endl;

      GMPfrac* FracPtr = NULL;
      FracPtr=static_cast< GMPfrac*>(FracStore);
     //need FracPtr cast to GMPfrac* in order to put data into FracStore memory!
     //compile now knows what usefrac and sign are....

     FracPtr->usefrac=true;

     FracPtr->SetFrac(aval, setproperfrac ,signval);
     
    return FracStore;

}

void* GMPfrac::operator new(size_t size, mpq_class& setproperfrac, bool signval){
      void* FracStore= ::operator new(size);
      cout << "size of FracStore=" << size << endl;
      cout <<"Fraction Stored at " << FracStore << endl;

      GMPfrac* FracPtr = NULL;
      FracPtr=static_cast< GMPfrac*>(FracStore);
     //need FracPtr cast to GMPfrac* in order to put data into FracStore memory!
     //compile now knows what usefrac and sign are....

     FracPtr->usefrac=true;

     FracPtr->SetFrac( setproperfrac ,signval);
     
    return FracStore;

}




