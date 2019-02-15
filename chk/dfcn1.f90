!DECK DFCN1
SUBROUTINE DFCN1(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DFCN1
  !***PURPOSE  Subsidiary to DNLS1Q.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FCN1-S, DFCN1-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine which evaluates the function for test program
  !   used in quick check of DNLS1E.
  !
  !   Numerical approximation of Jacobian is used.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added.  (WRB)
  !***END PROLOGUE  DFCN1
  !     .. Scalar Arguments ..
  REAL(8) :: Fjac
  INTEGER Iflag, Ldfjac, M, N
  !     .. Array Arguments ..
  REAL(8) :: Fvec(*), X(*)
  !     .. Local Scalars ..
  REAL(8) :: temp, two
  INTEGER i
  !     .. Intrinsic Functions ..
  INTRINSIC EXP
  !     .. Data statements ..
  DATA two/2.0D0/
  !***FIRST EXECUTABLE STATEMENT  DFCN1
  IF ( Iflag/=1 ) RETURN
  DO i = 1, M
    temp = i
    Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
  ENDDO
END SUBROUTINE DFCN1
