!** BKISR
SUBROUTINE BKISR(X,N,Summ,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BKISR-S, DBKISR-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     BKISR computes repeated integrals of the K0 Bessel function
  !     by the series for N=0,1, and 2.
  !
  !***
  ! **See also:**  BSKIN
  !***
  ! **Routines called:**  PSIXN, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)

  INTEGER i, Ierr, k, kk, kkn, k1, N, np
  REAL ak, atol, bk, fk, fn, hx, hxs, pol, pr, Summ, tkp, tol, trm, X, xln
  REAL PSIXN, R1MACH
  !
  REAL, PARAMETER :: c(2) = [ 1.57079632679489662E+00, 1.0E0 ]
  !* FIRST EXECUTABLE STATEMENT  BKISR
  Ierr = 0
  tol = MAX(R1MACH(4),1.0E-18)
  IF ( X>=tol ) THEN
    pr = 1.0E0
    pol = 0.0E0
    IF ( N/=0 ) THEN
      DO i = 1, N
        pol = -pol*X + c(i)
        pr = pr*X/i
      END DO
    END IF
    hx = X*0.5E0
    hxs = hx*hx
    xln = LOG(hx)
    np = N + 1
    tkp = 3.0E0
    fk = 2.0E0
    fn = N
    bk = 4.0E0
    ak = 2.0E0/((fn+1.0E0)*(fn+2.0E0))
    Summ = ak*(PSIXN(N+3)-PSIXN(3)+PSIXN(2)-xln)
    atol = Summ*tol*0.75E0
    DO k = 2, 20
      ak = ak*(hxs/bk)*((tkp+1.0E0)/(tkp+fn+1.0E0))*(tkp/(tkp+fn))
      k1 = k + 1
      kk = k1 + k
      kkn = kk + N
      trm = (PSIXN(k1)+PSIXN(kkn)-PSIXN(kk)-xln)*ak
      Summ = Summ + trm
      IF ( ABS(trm)<=atol ) GOTO 100
      tkp = tkp + 2.0E0
      bk = bk + tkp
      fk = fk + 1.0E0
    END DO
    Ierr = 2
    RETURN
    !-----------------------------------------------------------------------
    !     SMALL X CASE, X.LT.WORD TOLERANCE
    !-----------------------------------------------------------------------
  ELSEIF ( N>0 ) THEN
    Summ = c(N)
    RETURN
  ELSE
    hx = X*0.5E0
    Summ = PSIXN(1) - LOG(hx)
    RETURN
  END IF
  100  Summ = (Summ*hxs+PSIXN(np)-xln)*pr
  IF ( N==1 ) Summ = -Summ
  Summ = pol + Summ
  RETURN
END SUBROUTINE BKISR
