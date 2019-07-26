!** BKISR
ELEMENTAL SUBROUTINE BKISR(X,N,Summ,Ierr)
  !> Subsidiary to BSKIN
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
  USE service, ONLY : eps_sp
  !
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(OUT) :: Ierr
  REAL(SP), INTENT(IN) :: X
  REAL(SP), INTENT(OUT) :: Summ
  INTEGER :: i, k, kk, kkn, k1, np
  REAL(SP) :: ak, atol, bk, fk, fn, hx, hxs, pol, pr, tkp, tol, trm, xln
  !
  REAL(SP), PARAMETER :: c(2) = [ 1.57079632679489662_SP, 1._SP ]
  !* FIRST EXECUTABLE STATEMENT  BKISR
  Ierr = 0
  tol = MAX(eps_sp,1.0E-18_SP)
  IF( X>=tol ) THEN
    pr = 1._SP
    pol = 0._SP
    IF( N/=0 ) THEN
      DO i = 1, N
        pol = -pol*X + c(i)
        pr = pr*X/i
      END DO
    END IF
    hx = X*0.5_SP
    hxs = hx*hx
    xln = LOG(hx)
    np = N + 1
    tkp = 3._SP
    fk = 2._SP
    fn = N
    bk = 4._SP
    ak = 2._SP/((fn+1._SP)*(fn+2._SP))
    Summ = ak*(PSIXN(N+3)-PSIXN(3)+PSIXN(2)-xln)
    atol = Summ*tol*0.75_SP
    DO k = 2, 20
      ak = ak*(hxs/bk)*((tkp+1._SP)/(tkp+fn+1._SP))*(tkp/(tkp+fn))
      k1 = k + 1
      kk = k1 + k
      kkn = kk + N
      trm = (PSIXN(k1)+PSIXN(kkn)-PSIXN(kk)-xln)*ak
      Summ = Summ + trm
      IF( ABS(trm)<=atol ) GOTO 100
      tkp = tkp + 2._SP
      bk = bk + tkp
      fk = fk + 1._SP
    END DO
    Ierr = 2
    RETURN
    !-----------------------------------------------------------------------
    !     SMALL X CASE, X<WORD TOLERANCE
    !-----------------------------------------------------------------------
  ELSEIF( N>0 ) THEN
    Summ = c(N)
    RETURN
  ELSE
    hx = X*0.5_SP
    Summ = PSIXN(1) - LOG(hx)
    RETURN
  END IF
  100  Summ = (Summ*hxs+PSIXN(np)-xln)*pr
  IF( N==1 ) Summ = -Summ
  Summ = pol + Summ

  RETURN
END SUBROUTINE BKISR