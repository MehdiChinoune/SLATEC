!** ZEXP
SUBROUTINE ZEXP(Ar,Ai,Br,Bi)
  !> Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZEXP-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX(DP) :: B
  REAL(DP) :: Ar, Ai, Br, Bi
  !* FIRST EXECUTABLE STATEMENT  ZEXP
  B = EXP( CMPLX( Ar, Ai, DP ) )
  Br = REAL(B,DP)
  Bi = AIMAG(B)

END SUBROUTINE ZEXP