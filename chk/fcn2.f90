!*==FCN2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FCN2
SUBROUTINE FCN2(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
  IMPLICIT NONE
  !*--FCN25
  !***BEGIN PROLOGUE  FCN2
  !***PURPOSE  Subsidiary to SNLS1Q.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FCN2-S, DFCN2-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine to evaluate function and full Jacobian for test
  !   problem in quick check of SNLS1E.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added and code polished.
  !           (WRB)
  !***END PROLOGUE  FCN2
  !     .. Scalar Arguments ..
  INTEGER Iflag, Ldfjac, M, N
  !     .. Array Arguments ..
  REAL Fjac(Ldfjac,*), Fvec(*), X(*)
  !     .. Local Scalars ..
  REAL temp, two
  INTEGER i
  !     .. Intrinsic Functions ..
  INTRINSIC EXP
  !     .. Data statements ..
  DATA two/2.0E0/
  !***FIRST EXECUTABLE STATEMENT  FCN2
  IF ( Iflag==0 ) RETURN
  !
  !     Should we evaluate functions or Jacobian?
  !
  IF ( Iflag==1 ) THEN
    !
    !       Evaluate functions.
    !
    DO i = 1, M
      temp = i
      Fvec(i) = two + two*temp - EXP(temp*X(1)) - EXP(temp*X(2))
    ENDDO
  ELSE
    !
    !       Evaluate Jacobian.
    !
    IF ( Iflag/=2 ) RETURN
    DO i = 1, M
      temp = i
      Fjac(i,1) = -temp*EXP(temp*X(1))
      Fjac(i,2) = -temp*EXP(temp*X(2))
    ENDDO
  ENDIF
END SUBROUTINE FCN2
