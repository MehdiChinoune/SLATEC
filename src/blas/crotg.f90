!DECK CROTG
SUBROUTINE CROTG(Ca,Cb,C,S)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CROTG
  !***PURPOSE  Construct a Givens transformation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B10
  !***TYPE      COMPLEX (SROTG-S, DROTG-D, CROTG-C)
  !***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
  !             LINEAR ALGEBRA, VECTOR
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !    Complex Givens transformation
  !
  !    Construct the Givens transformation
  !
  !             (C    S)
  !       G  =  (      ),  C**2 + ABS(S)**2 =1,
  !             (-S   C)
  !
  !    which zeros the second entry of the complex 2-vector (CA,CB)**T
  !
  !    The quantity CA/ABS(CA)*NORM(CA,CB) overwrites CA in storage.
  !
  !    Input:
  !        CA (Complex)
  !        CB (Complex)
  !
  !    Output:
  !        CA (Complex)      CA/ABS(CA)*NORM(CA,CB)
  !        C  (Real)
  !        S  (Complex)
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CROTG
  COMPLEX Ca, Cb, S
  REAL C
  REAL norm, scale
  COMPLEX alpha
  !***FIRST EXECUTABLE STATEMENT  CROTG
  IF ( ABS(Ca)==0.0 ) THEN
    C = 0.0
    S = (1.0,0.0)
    Ca = Cb
  ELSE
    scale = ABS(Ca) + ABS(Cb)
    norm = scale*SQRT((ABS(Ca/scale))**2+(ABS(Cb/scale))**2)
    alpha = Ca/ABS(Ca)
    C = ABS(Ca)/norm
    S = alpha*CONJG(Cb)/norm
    Ca = alpha*norm
  ENDIF
END SUBROUTINE CROTG
