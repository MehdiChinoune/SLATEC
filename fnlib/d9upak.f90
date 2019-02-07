!*==D9UPAK.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK D9UPAK
SUBROUTINE D9UPAK(X,Y,N)
  IMPLICIT NONE
  !*--D9UPAK5
  !*** Start of declarations inserted by SPAG
  INTEGER N
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  D9UPAK
  !***PURPOSE  Unpack a floating point number X so that X = Y*2**N.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  A6B
  !***TYPE      DOUBLE PRECISION (R9UPAK-S, D9UPAK-D)
  !***KEYWORDS  FNLIB, UNPACK
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !   Unpack a floating point number X so that X = Y*2.0**N, where
  !   0.5 .LE. ABS(Y) .LT. 1.0.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   780701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900820  Corrected code to find Y between 0.5 and 1.0 rather than
  !           between 0.05 and 1.0.  (WRB)
  !***END PROLOGUE  D9UPAK
  DOUBLE PRECISION X , Y , absx
  !***FIRST EXECUTABLE STATEMENT  D9UPAK
  absx = ABS(X)
  N = 0
  IF ( X==0.0D0 ) THEN
    !
    Y = SIGN(absx,X)
  ELSE
    !
    DO WHILE ( absx<0.5D0 )
      N = N - 1
      absx = absx*2.0D0
    ENDDO
    !
    DO WHILE ( absx>=1.0D0 )
      N = N + 1
      absx = absx*0.5D0
    ENDDO
    Y = SIGN(absx,X)
  ENDIF
  !
END SUBROUTINE D9UPAK
