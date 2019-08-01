!** ZMLRI
PURE SUBROUTINE ZMLRI(Z,Fnu,Kode,N,Y,Nz,Tol)
  !> Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CMLRI-A, ZMLRI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z)>=0.0 BY THE
  !     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
  USE service, ONLY : tiny_dp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, iaz, ifnu, inu, itime, k, kk, km, m
  COMPLEX(DP) :: ck, cnorm, pt, p1, p2, rz, summ
  REAL(DP) :: ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, rho, &
    rho2, scle, tfnf, tst, x
  REAL(DP), PARAMETER :: sqrt_tiny = SQRT( tiny_dp )
  !* FIRST EXECUTABLE STATEMENT  CMLRI
  scle = tiny_dp/Tol
  Nz = 0
  az = ABS(Z)
  x = REAL(Z,DP)
  iaz = INT( az )
  ifnu = INT( Fnu )
  inu = ifnu + N - 1
  at = iaz + 1._DP
  ck = CMPLX(at,0._DP,DP)/Z
  rz = (2._DP,0._DP)/Z
  p1 = (0._DP,0._DP)
  p2 = (1._DP,0._DP)
  ack = (at+1._DP)/az
  rho = ack + SQRT(ack*ack-1._DP)
  rho2 = rho*rho
  tst = (rho2+rho2)/((rho2-1._DP)*(rho-1._DP))
  tst = tst/Tol
  !-----------------------------------------------------------------------
  !     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
  !-----------------------------------------------------------------------
  ak = at
  DO i = 1, 80
    pt = p2
    p2 = p1 - ck*p2
    p1 = pt
    ck = ck + rz
    ap = ABS(p2)
    IF( ap>tst*ak*ak ) GOTO 100
    ak = ak + 1._DP
  END DO
  Nz = -2
  RETURN
  100  i = i + 1
  k = 0
  IF( inu>=iaz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
    !-----------------------------------------------------------------------
    p1 = (0._DP,0._DP)
    p2 = (1._DP,0._DP)
    at = inu + 1._DP
    ck = CMPLX(at,0._DP,DP)/Z
    ack = at/az
    tst = SQRT(ack/Tol)
    itime = 1
    DO k = 1, 80
      pt = p2
      p2 = p1 - ck*p2
      p1 = pt
      ck = ck + rz
      ap = ABS(p2)
      IF( ap>=tst ) THEN
        IF( itime==2 ) GOTO 200
        ack = ABS(ck)
        flam = ack + SQRT(ack*ack-1._DP)
        fkap = ap/ABS(p1)
        rho = MIN(flam,fkap)
        tst = tst*SQRT(rho/(rho*rho-1._DP))
        itime = 2
      END IF
    END DO
    Nz = -2
    RETURN
  END IF
  !-----------------------------------------------------------------------
  !     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
  !-----------------------------------------------------------------------
  200  k = k + 1
  kk = MAX(i+iaz,k+inu)
  fkk = kk
  p1 = (0._DP,0._DP)
  !-----------------------------------------------------------------------
  !     SCALE P2 AND SUM BY SCLE
  !-----------------------------------------------------------------------
  p2 = CMPLX(scle,0._DP,DP)
  fnf = Fnu - ifnu
  tfnf = fnf + fnf
  bk = LOG_GAMMA(fkk+tfnf+1._DP) - LOG_GAMMA(fkk+1._DP)- LOG_GAMMA(tfnf+1._DP)
  bk = EXP(bk)
  summ = (0._DP,0._DP)
  km = kk - inu
  DO i = 1, km
    pt = p2
    p2 = p1 + CMPLX(fkk+fnf,0._DP,DP)*rz*p2
    p1 = pt
    ak = 1._DP - tfnf/(fkk+tfnf)
    ack = bk*ak
    summ = summ + CMPLX(ack+bk,0._DP,DP)*p1
    bk = ack
    fkk = fkk - 1._DP
  END DO
  Y(N) = p2
  IF( N/=1 ) THEN
    DO i = 2, N
      pt = p2
      p2 = p1 + CMPLX(fkk+fnf,0._DP,DP)*rz*p2
      p1 = pt
      ak = 1._DP - tfnf/(fkk+tfnf)
      ack = bk*ak
      summ = summ + CMPLX(ack+bk,0._DP,DP)*p1
      bk = ack
      fkk = fkk - 1._DP
      m = N - i + 1
      Y(m) = p2
    END DO
  END IF
  IF( ifnu>0 ) THEN
    DO i = 1, ifnu
      pt = p2
      p2 = p1 + CMPLX(fkk+fnf,0._DP,DP)*rz*p2
      p1 = pt
      ak = 1._DP - tfnf/(fkk+tfnf)
      ack = bk*ak
      summ = summ + CMPLX(ack+bk,0._DP,DP)*p1
      bk = ack
      fkk = fkk - 1._DP
    END DO
  END IF
  pt = Z
  IF( Kode==2 ) pt = pt - CMPLX(x,0._DP,DP)
  p1 = -CMPLX(fnf,0._DP,DP)*LOG(rz) + pt
  ap = LOG_GAMMA(1._DP+fnf)
  pt = p1 - CMPLX(ap,0._DP,DP)
  !-----------------------------------------------------------------------
  !     THE DIVISION EXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
  !     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
  !-----------------------------------------------------------------------
  p2 = p2 + summ
  ap = ABS(p2)
  IF( ap<tiny_dp )  ap = ABS(p2/sqrt_tiny) * sqrt_tiny
  p1 = CMPLX(1._DP/ap,0._DP,DP)
  ck = EXP(pt)*p1
  pt = CONJG(p2)*p1
  cnorm = ck*pt
  Y = Y*cnorm
  !
  RETURN
END SUBROUTINE ZMLRI