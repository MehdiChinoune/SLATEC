!** DEFER
SUBROUTINE DEFER(COFX,COFY,Idmn,Usol,Grhs)
  !>
  !  Subsidiary to SEPELI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (DEFER-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine first approximates the truncation error given by
  !     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY where
  !     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 on the interior and
  !     at the boundaries if periodic (here UXXX,UXXXX are the third
  !     and fourth partial derivatives of U with respect to X).
  !     TX is of the form AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
  !     at X=A or X=B if the boundary condition there is mixed.
  !     TX=0.0 along specified boundaries.  TY has symmetric form
  !     in Y with X,AFUN(X),BFUN(X) replaced by Y,DFUN(Y),EFUN(Y).
  !     The second order solution in USOL is used to approximate
  !     (via second order finite differencing) the truncation error
  !     and the result is added to the right hand side in GRHS
  !     and then transferred to USOL to be used as a new right
  !     hand side when calling BLKTRI for a fourth order solution.
  !
  !***
  ! **See also:**  SEPELI
  !***
  ! **Routines called:**  DX, DY
  !***
  ! COMMON BLOCKS    SPLPCM

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE SPLPCM, ONLY : L, AIT, CIT, DLX, DLY, IS, JS, K, KSWx, KSWy, MS, NS
  INTERFACE
    SUBROUTINE COFX(X,A,B,C)
      REAL :: X, A, B, C
    END SUBROUTINE COFX
    SUBROUTINE COFY(Y,D,E,F)
      REAL :: Y, D, E, F
    END SUBROUTINE COFY
  END INTERFACE
  INTEGER :: Idmn
  REAL :: Grhs(Idmn,NS), Usol(Idmn,NS)
  INTEGER :: i, j
  REAL :: ai, bi, ci, dj, ej, fj, tx, ty, uxxx, uxxxx, uyyy, uyyyy, xi, yj
  !* FIRST EXECUTABLE STATEMENT  DEFER
  DO j = JS, NS
    yj = CIT + (j-1)*DLY
    CALL COFY(yj,dj,ej,fj)
    DO i = IS, MS
      xi = AIT + (i-1)*DLX
      CALL COFX(xi,ai,bi,ci)
      !
      !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
      !
      CALL DX(Usol,Idmn,i,j,uxxx,uxxxx)
      CALL DY(Usol,Idmn,i,j,uyyy,uyyyy)
      tx = ai*uxxxx/12.0 + bi*uxxx/6.0
      ty = dj*uyyyy/12.0 + ej*uyyy/6.0
      !
      !     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
      !
      IF ( .NOT.(KSWx==1.OR.(i>1.AND.i<K)) )&
        tx = ai/3.0*(uxxxx/4.0+uxxx/DLX)
      IF ( .NOT.(KSWy==1.OR.(j>1.AND.j<L)) )&
        ty = dj/3.0*(uyyyy/4.0+uyyy/DLY)
      Grhs(i,j) = Grhs(i,j) + DLX**2*tx + DLY**2*ty
    END DO
  END DO
  !
  !     RESET THE RIGHT HAND SIDE IN USOL
  !
  DO i = IS, MS
    DO j = JS, NS
      Usol(i,j) = Grhs(i,j)
    END DO
  END DO
END SUBROUTINE DEFER
