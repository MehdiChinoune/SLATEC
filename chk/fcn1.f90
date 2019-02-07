!*==FCN1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FCN1
SUBROUTINE FCN1(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
  IMPLICIT NONE
  !*--FCN15
  !***BEGIN PROLOGUE  FCN1
  !***PURPOSE  Subsidiary to SNLS1Q.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FCN1-S, DFCN1-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine which evaluates the function for test program
  !   used in quick check of SNLS1E.
  !
  !   Numerical approximation of Jacobian is used.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added.  (WRB)
  !***END PROLOGUE  FCN1
  !     .. Scalar Arguments ..
  REAL Fjac
  INTEGER Iflag, Ldfjac, M, N
  !     .. Array Arguments ..
  REAL Fvec(*), X(*)
  !     .. Local Scalars ..
  REAL temp, two
  INTEGER i
  !     .. Intrinsic Functions ..
  INTRINSIC EXP
  !     .. Data statements ..
  DATA two/2.0E0/
  !***FIRST EXECUTABLE STATEMENT  FCN1
  IF ( Iflag/=1 ) RETURN
  DO i = 1, M
    temp = i
    Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
  ENDDO
END SUBROUTINE FCN1
