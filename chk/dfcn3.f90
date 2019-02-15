!DECK DFCN3
SUBROUTINE DFCN3(Iflag,M,N,X,Fvec,Fjrow,Nrow)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DFCN3
  !***PURPOSE  Subsidiary to DNLS1Q.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FCN3-S, DFCN3-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine to evaluate the Jacobian, one row at a time, for
  !   test problem used in quick check of DNLS1E.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added and code polished.
  !           (WRB)
  !***END PROLOGUE  DFCN3
  !     .. Scalar Arguments ..
  INTEGER Iflag, M, N, Nrow
  !     .. Array Arguments ..
  REAL(8) :: Fjrow(*), Fvec(*), X(*)
  !     .. Local Scalars ..
  REAL(8) :: temp, two
  INTEGER i
  !     .. Intrinsic Functions ..
  INTRINSIC EXP
  !     .. Data statements ..
  DATA two/2.0D0/
  !***FIRST EXECUTABLE STATEMENT  DFCN3
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
    !       Evaluate one row of Jacobian.
    !
    IF ( Iflag/=3 ) RETURN
    temp = Nrow
    Fjrow(1) = -temp*EXP(temp*X(1))
    Fjrow(2) = -temp*EXP(temp*X(2))
  ENDIF
END SUBROUTINE DFCN3
