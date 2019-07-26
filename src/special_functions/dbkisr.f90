!** DBKISR
ELEMENTAL SUBROUTINE DBKISR(X,N,Summ,Ierr)
  !> Subsidiary to DBSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (BKISR-S, DBKISR-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     DBKISR computes repeated integrals of the K0 Bessel function
  !     by the series for N=0,1, and 2.
  !
  !***
  ! **See also:**  DBSKIN
  !***
  ! **Routines called:**  D1MACH, DPSIXN

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE service, ONLY : eps_dp
  !
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(OUT) :: Ierr
  REAL(DP), INTENT(IN) :: X
  REAL(DP), INTENT(OUT) :: Summ
  !
  INTEGER :: i, k, kk, kkn, k1, np
  REAL(DP) :: ak, atol, bk, fk, fn, hx, hxs, pol, pr, tkp, tol, trm, xln
  !
  REAL(DP), PARAMETER :: c(2) = [ 1.57079632679489662_DP, 1._DP ]
  !* FIRST EXECUTABLE STATEMENT  DBKISR
  Ierr = 0
  tol = MAX(eps_dp,1.E-18_DP)
  IF( X>=tol ) THEN
    pr = 1._DP
    pol = 0._DP
    IF( N/=0 ) THEN
      DO i = 1, N
        pol = -pol*X + c(i)
        pr = pr*X/i
      END DO
    END IF
    hx = X*0.5_DP
    hxs = hx*hx
    xln = LOG(hx)
    np = N + 1
    tkp = 3._DP
    fk = 2._DP
    fn = N
    bk = 4._DP
    ak = 2._DP/((fn+1._DP)*(fn+2._DP))
    Summ = ak*(DPSIXN(N+3)-DPSIXN(3)+DPSIXN(2)-xln)
    atol = Summ*tol*0.75_DP
    DO k = 2, 20
      ak = ak*(hxs/bk)*((tkp+1._DP)/(tkp+fn+1._DP))*(tkp/(tkp+fn))
      k1 = k + 1
      kk = k1 + k
      kkn = kk + N
      trm = (DPSIXN(k1)+DPSIXN(kkn)-DPSIXN(kk)-xln)*ak
      Summ = Summ + trm
      IF( ABS(trm)<=atol ) GOTO 100
      tkp = tkp + 2._DP
      bk = bk + tkp
      fk = fk + 1._DP
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
    hx = X*0.5_DP
    Summ = DPSIXN(1) - LOG(hx)
    RETURN
  END IF
  100  Summ = (Summ*hxs+DPSIXN(np)-xln)*pr
  IF( N==1 ) Summ = -Summ
  Summ = pol + Summ

  RETURN
END SUBROUTINE DBKISR