!** BSRH
REAL FUNCTION BSRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BCRH-S, BSRH-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  BLKTRI
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    CBLKT

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL A, Bh, C, CNV, dx, EPS, F, Sgn, x, xl, Xll, xr, Xrr
  INTEGER IK, Iz, K, NCMplx, NM, NPP
  DIMENSION A(*), C(*), Bh(*)
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !* FIRST EXECUTABLE STATEMENT  BSRH
  xl = Xll
  xr = Xrr
  dx = .5*ABS(xr-xl)
  100  x = .5*(xl+xr)
  IF ( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF ( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    BSRH = .5*(xl+xr)
    RETURN
  ELSE
    xr = x
  ENDIF
  dx = .5*dx
  IF ( dx>CNV ) GOTO 100
  BSRH = .5*(xl+xr)
  RETURN
END FUNCTION BSRH
