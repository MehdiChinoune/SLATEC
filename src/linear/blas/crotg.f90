!** CROTG
SUBROUTINE CROTG(Ca,Cb,C,S)
  !>
  !  Construct a Givens transformation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B10
  !***
  ! **Type:**      COMPLEX (SROTG-S, DROTG-D, CROTG-C)
  !***
  ! **Keywords:**  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
  !             LINEAR ALGEBRA, VECTOR
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX(SP) Ca, Cb, S
  REAL(SP) C
  REAL(SP) norm, scalee
  COMPLEX(SP) alpha
  !* FIRST EXECUTABLE STATEMENT  CROTG
  IF ( ABS(Ca)==0.0 ) THEN
    C = 0.0
    S = (1.0,0.0)
    Ca = Cb
  ELSE
    scalee = ABS(Ca) + ABS(Cb)
    norm = scalee*SQRT((ABS(Ca/scalee))**2+(ABS(Cb/scalee))**2)
    alpha = Ca/ABS(Ca)
    C = ABS(Ca)/norm
    S = alpha*CONJG(Cb)/norm
    Ca = alpha*norm
  END IF
END SUBROUTINE CROTG
