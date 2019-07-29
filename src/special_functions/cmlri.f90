!** CMLRI
PURE SUBROUTINE CMLRI(Z,Fnu,Kode,N,Y,Nz,Tol)
  !> Subsidiary to CBESI and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CMLRI-A, ZMLRI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z)>=0.0 BY THE
  !     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
  !
  !***
  ! **See also:**  CBESI, CBESK
  !***
  ! **Routines called:**  GAMLN, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_sp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Fnu, Tol
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, iaz, ifnu, inu, itime, k, kk, km, m
  COMPLEX(SP) :: ck, cnorm, pt, p1, p2, rz, summ
  REAL(SP) :: ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, rho, &
    rho2, scle, tfnf, tst, x
  !* FIRST EXECUTABLE STATEMENT  CMLRI
  scle = 1.E+3_SP*tiny_sp/Tol
  Nz = 0
  az = ABS(Z)
  x = REAL(Z,SP)
  iaz = INT( az )
  ifnu = INT( Fnu )
  inu = ifnu + N - 1
  at = iaz + 1._SP
  ck = CMPLX(at,0._SP,SP)/Z
  rz = (2._SP,0._SP)/Z
  p1 = (0._SP,0._SP)
  p2 = (1._SP,0._SP)
  ack = (at+1._SP)/az
  rho = ack + SQRT(ack*ack-1._SP)
  rho2 = rho*rho
  tst = (rho2+rho2)/((rho2-1._SP)*(rho-1._SP))
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
    ak = ak + 1._SP
  END DO
  Nz = -2
  RETURN
  100  i = i + 1
  k = 0
  IF( inu>=iaz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
    !-----------------------------------------------------------------------
    p1 = (0._SP,0._SP)
    p2 = (1._SP,0._SP)
    at = inu + 1._SP
    ck = CMPLX(at,0._SP,SP)/Z
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
        flam = ack + SQRT(ack*ack-1._SP)
        fkap = ap/ABS(p1)
        rho = MIN(flam,fkap)
        tst = tst*SQRT(rho/(rho*rho-1._SP))
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
  p1 = (0._SP,0._SP)
  !-----------------------------------------------------------------------
  !     SCALE P2 AND SUM BY SCLE
  !-----------------------------------------------------------------------
  p2 = CMPLX(scle,0._SP,SP)
  fnf = Fnu - ifnu
  tfnf = fnf + fnf
  bk = LOG_GAMMA(fkk+tfnf+1._SP) - LOG_GAMMA(fkk+1._SP)- LOG_GAMMA(tfnf+1._SP)
  bk = EXP(bk)
  summ = (0._SP,0._SP)
  km = kk - inu
  DO i = 1, km
    pt = p2
    p2 = p1 + CMPLX(fkk+fnf,0._SP,SP)*rz*p2
    p1 = pt
    ak = 1._SP - tfnf/(fkk+tfnf)
    ack = bk*ak
    summ = summ + CMPLX(ack+bk,0._SP,SP)*p1
    bk = ack
    fkk = fkk - 1._SP
  END DO
  Y(N) = p2
  IF( N/=1 ) THEN
    DO i = 2, N
      pt = p2
      p2 = p1 + CMPLX(fkk+fnf,0._SP,SP)*rz*p2
      p1 = pt
      ak = 1._SP - tfnf/(fkk+tfnf)
      ack = bk*ak
      summ = summ + CMPLX(ack+bk,0._SP,SP)*p1
      bk = ack
      fkk = fkk - 1._SP
      m = N - i + 1
      Y(m) = p2
    END DO
  END IF
  IF( ifnu>0 ) THEN
    DO i = 1, ifnu
      pt = p2
      p2 = p1 + CMPLX(fkk+fnf,0._SP,SP)*rz*p2
      p1 = pt
      ak = 1._SP - tfnf/(fkk+tfnf)
      ack = bk*ak
      summ = summ + CMPLX(ack+bk,0._SP,SP)*p1
      bk = ack
      fkk = fkk - 1._SP
    END DO
  END IF
  pt = Z
  IF( Kode==2 ) pt = pt - CMPLX(x,0._SP,SP)
  p1 = -CMPLX(fnf,0._SP,SP)*LOG(rz) + pt
  ap = LOG_GAMMA(1._SP+fnf)
  pt = p1 - CMPLX(ap,0._SP,SP)
  !-----------------------------------------------------------------------
  !     THE DIVISION EXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
  !     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
  !-----------------------------------------------------------------------
  p2 = p2 + summ
  ap = ABS(p2)
  p1 = CMPLX(1._SP/ap,0._SP,SP)
  ck = EXP(pt)*p1
  pt = CONJG(p2)*p1
  cnorm = ck*pt
  Y = Y*cnorm
  !
  RETURN
END SUBROUTINE CMLRI