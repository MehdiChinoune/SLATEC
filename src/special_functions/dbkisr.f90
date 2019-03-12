!DECK DBKISR
SUBROUTINE DBKISR(X,N,Sum,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DBKISR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBSKIN
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BKISR-S, DBKISR-D)
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     DBKISR computes repeated integrals of the K0 Bessel function
  !     by the series for N=0,1, and 2.
  !
  !***SEE ALSO  DBSKIN
  !***ROUTINES CALLED  D1MACH, DPSIXN
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DBKISR
  INTEGER i, Ierr, k, kk, kkn, k1, N, np
  REAL(8) :: ak, atol, bk, c, fk, fn, hx, hxs, pol, pr, &
    Sum, tkp, tol, trm, X, xln
  REAL(8) :: DPSIXN, D1MACH
  DIMENSION c(2)
  SAVE c
  !
  DATA c(1), c(2)/1.57079632679489662D+00, 1.0D0/
  !***FIRST EXECUTABLE STATEMENT  DBKISR
  Ierr = 0
  tol = MAX(D1MACH(4),1.0D-18)
  IF ( X>=tol ) THEN
    pr = 1.0D0
    pol = 0.0D0
    IF ( N/=0 ) THEN
      DO i = 1, N
        pol = -pol*X + c(i)
        pr = pr*X/i
      ENDDO
    ENDIF
    hx = X*0.5D0
    hxs = hx*hx
    xln = LOG(hx)
    np = N + 1
    tkp = 3.0D0
    fk = 2.0D0
    fn = N
    bk = 4.0D0
    ak = 2.0D0/((fn+1.0D0)*(fn+2.0D0))
    Sum = ak*(DPSIXN(N+3)-DPSIXN(3)+DPSIXN(2)-xln)
    atol = Sum*tol*0.75D0
    DO k = 2, 20
      ak = ak*(hxs/bk)*((tkp+1.0D0)/(tkp+fn+1.0D0))*(tkp/(tkp+fn))
      k1 = k + 1
      kk = k1 + k
      kkn = kk + N
      trm = (DPSIXN(k1)+DPSIXN(kkn)-DPSIXN(kk)-xln)*ak
      Sum = Sum + trm
      IF ( ABS(trm)<=atol ) GOTO 100
      tkp = tkp + 2.0D0
      bk = bk + tkp
      fk = fk + 1.0D0
    ENDDO
    Ierr = 2
    RETURN
    !-----------------------------------------------------------------------
    !     SMALL X CASE, X.LT.WORD TOLERANCE
    !-----------------------------------------------------------------------
  ELSEIF ( N>0 ) THEN
    Sum = c(N)
    RETURN
  ELSE
    hx = X*0.5D0
    Sum = DPSIXN(1) - LOG(hx)
    RETURN
  ENDIF
  100  Sum = (Sum*hxs+DPSIXN(np)-xln)*pr
  IF ( N==1 ) Sum = -Sum
  Sum = pol + Sum
  RETURN
END SUBROUTINE DBKISR
