!DECK DFCN2
SUBROUTINE DFCN2(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DFCN2
  !***PURPOSE  Subsidiary to DNLS1Q.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FCN2-S, DFCN2-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine to evaluate function and full Jacobian for test
  !   problem in quick check of DNLS1E.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added and code polished.
  !           (WRB)
  !***END PROLOGUE  DFCN2
  !     .. Scalar Arguments ..
  INTEGER Iflag, Ldfjac, M, N
  !     .. Array Arguments ..
  REAL(8) :: Fjac(Ldfjac,*), Fvec(*), X(*)
  !     .. Local Scalars ..
  REAL(8) :: temp, two
  INTEGER i
  !     .. Intrinsic Functions ..
  INTRINSIC EXP
  !     .. Data statements ..
  DATA two/2.0D0/
  !***FIRST EXECUTABLE STATEMENT  DFCN2
  IF ( Iflag==0 ) RETURN
  !
  !     Should we evaluate function or Jacobian?
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
END SUBROUTINE DFCN2
