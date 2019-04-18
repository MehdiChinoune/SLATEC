!** D9UPAK
SUBROUTINE D9UPAK(X,Y,N)
  !>
  !  Unpack a floating point number X so that X = Y*2**N.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  A6B
  !***
  ! **Type:**      DOUBLE PRECISION (R9UPAK-S, D9UPAK-D)
  !***
  ! **Keywords:**  FNLIB, UNPACK
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   Unpack a floating point number X so that X = Y*2.0**N, where
  !   0.5 .LE. ABS(Y) .LT. 1.0.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   780701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900820  Corrected code to find Y between 0.5 and 1.0 rather than
  !           between 0.05 and 1.0.  (WRB)

  INTEGER N
  REAL(8) :: X, Y, absx
  !* FIRST EXECUTABLE STATEMENT  D9UPAK
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
    END DO
    !
    DO WHILE ( absx>=1.0D0 )
      N = N + 1
      absx = absx*0.5D0
    END DO
    Y = SIGN(absx,X)
  END IF
  !
END SUBROUTINE D9UPAK
