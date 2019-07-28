!** BSRH
REAL(SP) PURE FUNCTION BSRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  !> Subsidiary to BLKTRI
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
  USE CBLKT, ONLY : cnv_com
  !
  INTERFACE
    REAL(SP) PURE FUNCTION F(X,Iz,C,A,Bh)
      IMPORT SP
      INTEGER, INTENT(IN) :: Iz
      REAL(SP), INTENT(IN) :: X, A(Iz), Bh(Iz), C(Iz)
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Iz
  REAL(SP), INTENT(IN) :: A(Iz), Bh(Iz), C(Iz)
  REAL(SP), INTENT(IN) :: Sgn, Xll, Xrr
  !
  REAL(SP) :: dx, x, xl, xr
  !* FIRST EXECUTABLE STATEMENT  BSRH
  xl = Xll
  xr = Xrr
  dx = 0.5_SP*ABS(xr-xl)
  100  x = 0.5_SP*(xl+xr)
  IF( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    BSRH = 0.5_SP*(xl+xr)
    RETURN
  ELSE
    xr = x
  END IF
  dx = 0.5_SP*dx
  IF( dx>cnv_com ) GOTO 100
  BSRH = 0.5_SP*(xl+xr)
  !
  RETURN
END FUNCTION BSRH