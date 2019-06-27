!** ZLOG
SUBROUTINE ZLOG(Ar,Ai,Br,Bi)
  !> Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZLOG-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION COMPLEX LOGARITHM B=LOG(A)
  !     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  ZABS

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX(DP) :: B
  REAL(DP) :: Ar, Ai, Br, Bi
  !* FIRST EXECUTABLE STATEMENT  ZLOG

  B = LOG( CMPLX( Ar, Ai, DP ) )
  Br = REAL( B, DP )
  Bi = AIMAG( B )

  RETURN
END SUBROUTINE ZLOG
