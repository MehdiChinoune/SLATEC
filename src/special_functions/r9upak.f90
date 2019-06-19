!** R9UPAK
SUBROUTINE R9UPAK(X,Y,N)
  !> Unpack a floating point number X so that X = Y*2**N.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  A6B
  !***
  ! **Type:**      SINGLE PRECISION (R9UPAK-S, D9UPAK-D)
  !***
  ! **Keywords:**  FNLIB, UNPACK
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   Unpack a floating point number X so that X = Y*2.0**N, where
  !   0.5 <= ABS(Y) < 1.0.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   780701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  !* FIRST EXECUTABLE STATEMENT  R9UPAK
  REAL(SP) :: absx, X, Y
  INTEGER :: N
  absx = ABS(X)
  N = 0
  IF( X==0.0E0 ) THEN
    !
    Y = SIGN(absx,X)
  ELSE
    !
    DO WHILE( absx<0.5E0 )
      N = N - 1
      absx = absx*2.0E0
    END DO
    !
    DO WHILE( absx>=1.0E0 )
      N = N + 1
      absx = absx*0.5E0
    END DO
    Y = SIGN(absx,X)
  END IF
  !
END SUBROUTINE R9UPAK
