!** CMLRI
SUBROUTINE CMLRI(Z,Fnu,Kode,N,Y,Nz,Tol)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBESI and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CMLRI-A, ZMLRI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
  !     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
  !
  !***
  ! **See also:**  CBESI, CBESK
  !***
  ! **Routines called:**  GAMLN, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER i, iaz, idum, ifnu, inu, itime, k, kk, km, Kode, m, N, Nz
  COMPLEX ck, cnorm, pt, p1, p2, rz, summ, Y(N), Z
  REAL ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, Fnu, rho, &
    rho2, scle, tfnf, Tol, tst, x, GAMLN, R1MACH
  COMPLEX, PARAMETER :: czero = (0.0E0,0.0E0), cone = (1.0E0,0.0E0), ctwo = (2.0E0,0.0E0)
  scle = 1.0E+3*R1MACH(1)/Tol
  !* FIRST EXECUTABLE STATEMENT  CMLRI
  Nz = 0
  az = ABS(Z)
  x = REAL(Z)
  iaz = INT( az )
  ifnu = INT( Fnu )
  inu = ifnu + N - 1
  at = iaz + 1.0E0
  ck = CMPLX(at,0.0E0)/Z
  rz = ctwo/Z
  p1 = czero
  p2 = cone
  ack = (at+1.0E0)/az
  rho = ack + SQRT(ack*ack-1.0E0)
  rho2 = rho*rho
  tst = (rho2+rho2)/((rho2-1.0E0)*(rho-1.0E0))
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
    IF ( ap>tst*ak*ak ) GOTO 100
    ak = ak + 1.0E0
  END DO
  Nz = -2
  RETURN
  100  i = i + 1
  k = 0
  IF ( inu>=iaz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
    !-----------------------------------------------------------------------
    p1 = czero
    p2 = cone
    at = inu + 1.0E0
    ck = CMPLX(at,0.0E0)/Z
    ack = at/az
    tst = SQRT(ack/Tol)
    itime = 1
    DO k = 1, 80
      pt = p2
      p2 = p1 - ck*p2
      p1 = pt
      ck = ck + rz
      ap = ABS(p2)
      IF ( ap>=tst ) THEN
        IF ( itime==2 ) GOTO 200
        ack = ABS(ck)
        flam = ack + SQRT(ack*ack-1.0E0)
        fkap = ap/ABS(p1)
        rho = MIN(flam,fkap)
        tst = tst*SQRT(rho/(rho*rho-1.0E0))
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
  p1 = czero
  !-----------------------------------------------------------------------
  !     SCALE P2 AND SUM BY SCLE
  !-----------------------------------------------------------------------
  p2 = CMPLX(scle,0.0E0)
  fnf = Fnu - ifnu
  tfnf = fnf + fnf
  bk = GAMLN(fkk+tfnf+1.0E0,idum) - GAMLN(fkk+1.0E0,idum)- GAMLN(tfnf+1.0E0,idum)
  bk = EXP(bk)
  summ = czero
  km = kk - inu
  DO i = 1, km
    pt = p2
    p2 = p1 + CMPLX(fkk+fnf,0.0E0)*rz*p2
    p1 = pt
    ak = 1.0E0 - tfnf/(fkk+tfnf)
    ack = bk*ak
    summ = summ + CMPLX(ack+bk,0.0E0)*p1
    bk = ack
    fkk = fkk - 1.0E0
  END DO
  Y(N) = p2
  IF ( N/=1 ) THEN
    DO i = 2, N
      pt = p2
      p2 = p1 + CMPLX(fkk+fnf,0.0E0)*rz*p2
      p1 = pt
      ak = 1.0E0 - tfnf/(fkk+tfnf)
      ack = bk*ak
      summ = summ + CMPLX(ack+bk,0.0E0)*p1
      bk = ack
      fkk = fkk - 1.0E0
      m = N - i + 1
      Y(m) = p2
    END DO
  END IF
  IF ( ifnu>0 ) THEN
    DO i = 1, ifnu
      pt = p2
      p2 = p1 + CMPLX(fkk+fnf,0.0E0)*rz*p2
      p1 = pt
      ak = 1.0E0 - tfnf/(fkk+tfnf)
      ack = bk*ak
      summ = summ + CMPLX(ack+bk,0.0E0)*p1
      bk = ack
      fkk = fkk - 1.0E0
    END DO
  END IF
  pt = Z
  IF ( Kode==2 ) pt = pt - CMPLX(x,0.0E0)
  p1 = -CMPLX(fnf,0.0E0)*CLOG(rz) + pt
  ap = GAMLN(1.0E0+fnf,idum)
  pt = p1 - CMPLX(ap,0.0E0)
  !-----------------------------------------------------------------------
  !     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
  !     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
  !-----------------------------------------------------------------------
  p2 = p2 + summ
  ap = ABS(p2)
  p1 = CMPLX(1.0E0/ap,0.0E0)
  ck = CEXP(pt)*p1
  pt = CONJG(p2)*p1
  cnorm = ck*pt
  DO i = 1, N
    Y(i) = Y(i)*cnorm
  END DO
  RETURN
END SUBROUTINE CMLRI
