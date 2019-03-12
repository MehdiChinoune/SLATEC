!DECK DEFE4
SUBROUTINE DEFE4(COFX,Idmn,Usol,Grhs)
  IMPLICIT NONE
  REAL ai, AIT, bi, BIT, ci, CIT, DIT, DLX, DLX4, DLY, DLY4, &
    Grhs, TDLx3, TDLy3, tx, ty, Usol, uxxx, uxxxx, uyyy
  REAL uyyyy, xi
  INTEGER i, Idmn, IS, j, JS, K, KSWx, KSWy, L, MIT, MS, NIT, NS
  !***BEGIN PROLOGUE  DEFE4
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SEPX4
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (DEFE4-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
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
  !***SEE ALSO  SEPX4
  !***ROUTINES CALLED  DX4, DY4
  !***COMMON BLOCKS    SPL4
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  DEFE4
  !
  COMMON /SPL4  / KSWx, KSWy, K, L, AIT, BIT, CIT, DIT, MIT, NIT, &
    IS, MS, JS, NS, DLX, DLY, TDLx3, TDLy3, DLX4, &
    DLY4
  DIMENSION Grhs(Idmn,*), Usol(Idmn,*)
  EXTERNAL COFX
  !***FIRST EXECUTABLE STATEMENT  DEFE4
  DO i = IS, MS
    xi = AIT + (i-1)*DLX
    CALL COFX(xi,ai,bi,ci)
    DO j = JS, NS
      !
      !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
      !
      CALL DX4(Usol,Idmn,i,j,uxxx,uxxxx)
      CALL DY4(Usol,Idmn,i,j,uyyy,uyyyy)
      tx = ai*uxxxx/12.0 + bi*uxxx/6.0
      ty = uyyyy/12.0
      !
      !     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
      !
      IF ( .NOT.(KSWx==1.OR.(i>1.AND.i<K)) )&
        tx = ai/3.0*(uxxxx/4.0+uxxx/DLX)
      IF ( .NOT.(KSWy==1.OR.(j>1.AND.j<L)) ) ty = (uyyyy/4.0+uyyy/DLY)/3.0
      Grhs(i,j) = Grhs(i,j) + DLY**2*(DLX**2*tx+DLY**2*ty)
    ENDDO
  ENDDO
  !
  !     RESET THE RIGHT HAND SIDE IN USOL
  !
  DO i = IS, MS
    DO j = JS, NS
      Usol(i,j) = Grhs(i,j)
    ENDDO
  ENDDO
END SUBROUTINE DEFE4
