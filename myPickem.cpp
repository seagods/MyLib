int myPickem(int** picks, int& m, int & N, int & ifacs){

    //m is the nuber of objects to picked out of N obects

    int ireturn;
    ireturn=0;

    if(m==1){
    for(int i1=0; i1<ifacs; i1++){
         picks[i1][0]=i1;
    }
    return ireturn;
    }

    //we are picking more than one object out of N

    int nextN=N-1;
    int nextm=m-1;
    int kount=0;
    int start=0;
    int lastkount;

    while(kount<ifacs){

       for(int j=0; j<nextm; j++){

       lastkount=kount;

       int nextNend=nextN/2;
       int nextiremain=nextNend%2;
       nextNend=nextNend+nextiremain;

       int nextifacs=1;

       for(int k=0; k<nextm; k++){
            nextifacs=nextifacs*(nextN-k)/(k+1);
       }

       int** nextpicks=new int*[nextifacs];

       for(int k=0; k< nextifacs; k++){
           nextpicks[k]=new int[nextm];
       }


       ireturn=myPickem(nextpicks, nextm, nextN, nextifacs);
    
       for(int k=0; k < nextifacs; k++){
 
          picks[k+lastkount][0]=start;
          
          for(int l=0; l<nextm; l++){
             picks[k+lastkount][l+1]=nextpicks[k][l]+start+1;
          } 
          kount++;     
       }
       start++;

       nextN--;

       lastkount=kount;

       for(int k=0; k< nextifacs; k++){
          delete[] nextpicks[k];
       } 
       delete[] nextpicks; 
       
       }  //end j loop

   } // end while kount    

   return ireturn;
}
