Regarding splpmu.f, splpmn.f,  splpmu.f, and  splpmn.f
      
These are part of the slatec library (all available from NETLIB).
and are all subsidiary to splp --- a slatex linear programming
routine. The originals do not compile under gfortran due to ancient
FORTRAN ASSIGN statements. I have hacked these four routines
 and replaced the ASSIGNS with a set of GOTO statements.
I keep them copies of them here. The are not actually compiled into 
 libmylib.


Programs actually compiled into libmylib:
 
1:  programs written by me
2:  my modifications of programs written by other people

My Eigen33 and 44  and mymonic 3 and 4 use analytical solutions for cubic and quartic equations.
 (corrects error in Abramowitz and Stegun on solution of quartics).


**************************************************************************************
1:     programs written by me
**************************************************************************************
       myGaussElim, myMonic2, myMonic3, myMonic4.
1.1:   myGaussElim - Gaussian elimination with scaled pivoting, matrix
       is arbitrary N by N. Gaussian elimination program: fairly well
       tested, but use linpack or lapack instead!
       It's here to be developed into something else for an RT application:
       In this application the matrix elements will depend on the RHS!
       So myGaussElim assumes we have only one right hand side, and we don't
       bother storing multipliers or returning permutation vector.  We have 
       actually overwritten matrix(i,j) with zero values where the upper
       triangular is zero --- not necessary, but makes things more obvious at
       this stage. 

1.1.a: myGaussElimZ - Same as myGaussElim, except the coefficient matrix, the
       r.h.s, and the solution vector are all C++ complex numbers (#include <complex>)

1.2:   myMonic2.
       Solves a monic quadratic polynomial with one step Newton-Raphson.
       Trivial I know

1.3    myMonic3: Solves a monic cubic polynomial using Abramowitz and Stegun's
       analytical formula (Page.17 of Handbook of Mathematical Functions) with
       Newton-Raphson. (Last step can be necessary under certain circumstances)
       Needs more work, hope it will "do" for now. 
  
1.4    myMonic4: Analytical solution a monic quartic polynomial: This time we
       DO NOT USE Abramowitz and Stegun's formula (or Wolfram's website formula
       either). Again, we  use a Newton-Raphson polish at end.
       Needs more work, hope it will "do" for now.  

1.5    myEigVal44: Eigenvectors of a general real 4 by 4 matrix. Uses determinants
       to calculate coefficients for monic polynomial and calculates analytical solution
       using myMonic4.

1.6    myEigen44: Eigenvalues and eigenvectors of a general real 4 by 4 matrix. Uses determinants
       to calculate coefficients for monic polynomial and calculates analytical solution
       using myMonic4. Uses myGaussElim and myGaussElimZ to calculate real or complex eigenvectors.
       NOTE: Unless you KNOW IN ADVANCE the eigenvalues are distinct all non-zero, and in each case
       the eigenvector (x1,x2,x3,x4) always has x4 sufficiently distinct from zero.
       Regard it as a toy routine, there are no safety belts!
       Use LAPACK for serious stuff!

1.7    myEigVal33: Similar to myEigVal44, but for a general 3 by 3 matrix, and uses myMonic3

1.8    myEigen33: Similar to myEigen44, but for a general 3 by 3 matrix, and uses myMonic3

1.9    myPickem: stores lists of m choices out N objects in an array

1.9    myRecoverMonic: Calculates the coefficients of a monic polynomial from a vector
       containing N real roots (easily change over to complex - not bothered yet)
       

**************************************************************************************
2:  my modifications of programs written by other people
**************************************************************************************
    mydbesy.f //minor hack of slatec routine to get to compile with gfortran


**************************************************************************************
3:  programs  written by other people - no modification
**************************************************************************************
_______________________________________________________________________________
Problems and Fixes
1:  Monic4 had bug, fixed by using abs in test for roots reality
