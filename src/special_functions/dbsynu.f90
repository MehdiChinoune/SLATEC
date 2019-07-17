!** DBSYNU
PURE SUBROUTINE DBSYNU(X,Fnu,N,Y)
  !> Subsidiary to DBESY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (BESYNU-S, DBSYNU-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract  **** A DOUBLE PRECISION routine ****
  !         DBSYNU computes N member sequences of Y Bessel functions
  !         Y_{FNU+I-1}(X), I=1,N for non-negative orders FNU and
  !         positive X. Equations of the references are implemented on
  !         small orders DNU for Y_{DNU}(X) and Y_{DNU+1}(X).
  !         Forward recursion with the three term recursion relation
  !         generates higher orders FNU+I-1, I=1,...,N.
  !
  !         To start the recursion FNU is normalized to the interval
  !         -0.5<=DNU<0.5. A special form of the power series is
  !         implemented on 0<X<=X1 while the Miller algorithm for the
  !         K Bessel function in terms of the confluent hypergeometric
  !         function U(FNU+0.5,2*FNU+1,I*X) is implemented on X1<X<=X
  !         Here I is the complex number SQRT(-1.).
  !         For X>X2, the asymptotic expansion for large X is used.
  !         When FNU is a half odd integer, a special formula for
  !         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
  !
  !         The maximum number of significant digits obtainable
  !         is the smaller of 14 and the number of digits carried in
  !         DOUBLE PRECISION arithmetic.
  !
  !         DBSYNU assumes that a significant digit SINH function is
  !         available.
  !
  !     Description of Arguments
  !
  !         INPUT
  !           X      - X>0.0D0
  !           FNU    - Order of initial Y function, FNU>=0.0D0
  !           N      - Number of members of the sequence, N>=1
  !
  !         OUTPUT
  !           Y      - A vector whose first N components contain values
  !                    for the sequence Y(I)=Y_{FNU+I-1}, I=1,N.
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow - a fatal error
  !
  !***
  ! **See also:**  DBESY
  !***
  ! **References:**  N. M. Temme, On the numerical evaluation of the ordinary
  !                 Bessel function of the second kind, Journal of
  !                 Computational Physics 21, (1976), pp. 343-350.
  !               N. M. Temme, On the numerical evaluation of the modified
  !                 Bessel function of the third kind, Journal of
  !                 Computational Physics 19, (1975), pp. 324-337.
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   900328  Added TYPE section.  (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : D1MACH
  !
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(IN) :: Fnu, X
  REAL(DP), INTENT(OUT) :: Y(N)
  INTEGER :: i, inu, j, k, kk, nn
  REAL(DP) :: a(120), ak, arg, a1, a2, bk, cb(120), cbk, cck, ck, coef, cpt, cp1, cp2, &
    cs, cs1, cs2, cx, dnu, dnu2, etest, etx, f, fc, fhs, fk, fks, flrx, fmu, fn, &
    fx, g, g1, g2, p, pt, q, rb(120), rbk, rck, relb, rpt, rp1, rp2, rs, &
    rs1, rs2, rx, s, sa, sb, smu, ss, st, s1, s2, tb, tm, tol, t1, t2
  REAL(DP), PARAMETER ::  x1 = 3._DP, x2 = 20._DP
  REAL(DP), PARAMETER ::  pi = 3.14159265358979E+00_DP, rthpi = 7.97884560802865E-01_DP
  REAL(DP), PARAMETER ::  hpi = 1.57079632679490_DP
  REAL(DP), PARAMETER :: cc(8) = [ 5.77215664901533E-01_DP, -4.20026350340952E-02_DP, &
    -4.21977345555443E-02_DP, 7.21894324666300E-03_DP, -2.15241674114900E-04_DP, &
    -2.01348547807000E-05_DP, 1.13302723200000E-06_DP, 6.11609500000000E-09_DP ]
  !* FIRST EXECUTABLE STATEMENT  DBSYNU
  ak = D1MACH(3)
  tol = MAX(ak,1.E-15_DP)
  s2 = 0._DP
  IF( X<=0._DP ) THEN
    ERROR STOP 'DBSYNU : X < 0'
  ELSEIF( Fnu<0._DP ) THEN
    ERROR STOP 'DBSYNU : FNU < 0'
  ELSEIF( N<1 ) THEN
    ERROR STOP 'DBSYNU : N < 1'
  ELSE
    rx = 2._DP/X
    inu = INT(Fnu+0.5_DP)
    dnu = Fnu - inu
    IF( ABS(dnu)==0.5_DP ) THEN
      !
      !     FNU=HALF ODD INTEGER CASE
      !
      coef = rthpi/SQRT(X)
      s1 = coef*SIN(X)
      s2 = -coef*COS(X)
    ELSE
      dnu2 = 0._DP
      IF( ABS(dnu)>=tol ) dnu2 = dnu*dnu
      IF( X>x1 ) THEN
        coef = rthpi/SQRT(X)
        IF( X>x2 ) THEN
          !
          !     ASYMPTOTIC EXPANSION FOR LARGE X, X>X2
          !
          nn = 2
          IF( inu==0 .AND. N==1 ) nn = 1
          dnu2 = dnu + dnu
          fmu = 0._DP
          IF( ABS(dnu2)>=tol ) fmu = dnu2*dnu2
          arg = X - hpi*(dnu+0.5_DP)
          sa = SIN(arg)
          sb = COS(arg)
          etx = 8._DP*X
          DO k = 1, nn
            s1 = s2
            t2 = (fmu-1._DP)/etx
            ss = t2
            relb = tol*ABS(t2)
            t1 = etx
            s = 1._DP
            fn = 1._DP
            ak = 0._DP
            DO j = 1, 13
              t1 = t1 + etx
              ak = ak + 8._DP
              fn = fn + ak
              t2 = -t2*(fmu-fn)/t1
              s = s + t2
              t1 = t1 + etx
              ak = ak + 8._DP
              fn = fn + ak
              t2 = t2*(fmu-fn)/t1
              ss = ss + t2
              IF( ABS(t2)<=relb ) EXIT
            END DO
            s2 = coef*(s*sa+ss*sb)
            fmu = fmu + 8._DP*dnu + 4._DP
            tb = sa
            sa = -sb
            sb = tb
          END DO
          IF( nn<=1 ) THEN
            s1 = s2
            GOTO 100
          END IF
        ELSE
          !
          !     MILLER ALGORITHM FOR X1<X<=X2
          !
          etest = COS(pi*dnu)/(pi*X*tol)
          fks = 1._DP
          fhs = 0.25_DP
          fk = 0._DP
          rck = 2._DP
          cck = X + X
          rp1 = 0._DP
          cp1 = 0._DP
          rp2 = 1._DP
          cp2 = 0._DP
          k = 0
          DO
            k = k + 1
            fk = fk + 1._DP
            ak = (fhs-dnu2)/(fks+fk)
            pt = fk + 1._DP
            rbk = rck/pt
            cbk = cck/pt
            rpt = rp2
            cpt = cp2
            rp2 = rbk*rpt - cbk*cpt - ak*rp1
            cp2 = cbk*rpt + rbk*cpt - ak*cp1
            rp1 = rpt
            cp1 = cpt
            rb(k) = rbk
            cb(k) = cbk
            a(k) = ak
            rck = rck + 2._DP
            fks = fks + fk + fk + 1._DP
            fhs = fhs + fk + fk
            pt = MAX(ABS(rp1),ABS(cp1))
            fc = (rp1/pt)**2 + (cp1/pt)**2
            pt = pt*SQRT(fc)*fk
            IF( etest<=pt ) THEN
              kk = k
              rs = 1._DP
              cs = 0._DP
              rp1 = 0._DP
              cp1 = 0._DP
              rp2 = 1._DP
              cp2 = 0._DP
              DO i = 1, k
                rpt = rp2
                cpt = cp2
                rp2 = (rb(kk)*rpt-cb(kk)*cpt-rp1)/a(kk)
                cp2 = (cb(kk)*rpt+rb(kk)*cpt-cp1)/a(kk)
                rp1 = rpt
                cp1 = cpt
                rs = rs + rp2
                cs = cs + cp2
                kk = kk - 1
              END DO
              pt = MAX(ABS(rs),ABS(cs))
              fc = (rs/pt)**2 + (cs/pt)**2
              pt = pt*SQRT(fc)
              rs1 = (rp2*(rs/pt)+cp2*(cs/pt))/pt
              cs1 = (cp2*(rs/pt)-rp2*(cs/pt))/pt
              fc = hpi*(dnu-0.5_DP) - X
              p = COS(fc)
              q = SIN(fc)
              s1 = (cs1*q-rs1*p)*coef
              IF( inu>0 .OR. N>1 ) THEN
                pt = MAX(ABS(rp2),ABS(cp2))
                fc = (rp2/pt)**2 + (cp2/pt)**2
                pt = pt*SQRT(fc)
                rpt = dnu + 0.5_DP - (rp1*(rp2/pt)+cp1*(cp2/pt))/pt
                cpt = X - (cp1*(rp2/pt)-rp1*(cp2/pt))/pt
                cs2 = cs1*cpt - rs1*rpt
                rs2 = rpt*cs1 + rs1*cpt
                s2 = (rs2*q+cs2*p)*coef/X
                EXIT
              ELSE
                Y(1) = s1
                RETURN
              END IF
            END IF
          END DO
        END IF
      ELSE
        !
        !     SERIES FOR X<=X1
        !
        a1 = 1._DP - dnu
        a2 = 1._DP + dnu
        t1 = 1._DP/GAMMA(a1)
        t2 = 1._DP/GAMMA(a2)
        IF( ABS(dnu)>0.1_DP ) THEN
          g1 = (t1-t2)/dnu
        ELSE
          !     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
          s = cc(1)
          ak = 1._DP
          DO k = 2, 8
            ak = ak*dnu2
            tm = cc(k)*ak
            s = s + tm
            IF( ABS(tm)<tol ) EXIT
          END DO
          g1 = -(s+s)
        END IF
        g2 = t1 + t2
        smu = 1._DP
        fc = 1._DP/pi
        flrx = LOG(rx)
        fmu = dnu*flrx
        tm = 0._DP
        IF( dnu/=0._DP ) THEN
          tm = SIN(dnu*hpi)/dnu
          tm = (dnu+dnu)*tm*tm
          fc = dnu/SIN(dnu*pi)
          IF( fmu/=0._DP ) smu = SINH(fmu)/fmu
        END IF
        f = fc*(g1*COSH(fmu)+g2*flrx*smu)
        fx = EXP(fmu)
        p = fc*t1*fx
        q = fc*t2/fx
        g = f + tm*q
        ak = 1._DP
        ck = 1._DP
        bk = 1._DP
        s1 = g
        s2 = p
        IF( inu>0 .OR. N>1 ) THEN
          IF( X>=tol ) THEN
            cx = X*X*0.25_DP
            DO
              f = (ak*f+p+q)/(bk-dnu2)
              p = p/(ak-dnu)
              q = q/(ak+dnu)
              g = f + tm*q
              ck = -ck*cx/ak
              t1 = ck*g
              s1 = s1 + t1
              t2 = ck*(p-ak*g)
              s2 = s2 + t2
              bk = bk + ak + ak + 1._DP
              ak = ak + 1._DP
              s = ABS(t1)/(1._DP+ABS(s1)) + ABS(t2)/(1._DP+ABS(s2))
              IF( s<=tol ) EXIT
            END DO
          END IF
          s2 = -s2*rx
          s1 = -s1
        ELSE
          IF( X>=tol ) THEN
            cx = X*X*0.25_DP
            DO
              f = (ak*f+p+q)/(bk-dnu2)
              p = p/(ak-dnu)
              q = q/(ak+dnu)
              g = f + tm*q
              ck = -ck*cx/ak
              t1 = ck*g
              s1 = s1 + t1
              bk = bk + ak + ak + 1._DP
              ak = ak + 1._DP
              s = ABS(t1)/(1._DP+ABS(s1))
              IF( s<=tol ) EXIT
            END DO
          END IF
          Y(1) = -s1
          RETURN
        END IF
      END IF
    END IF
    !
    !     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
    !
    ck = (dnu+dnu+2._DP)/X
    IF( N==1 ) inu = inu - 1
    IF( inu>0 ) THEN
      DO i = 1, inu
        st = s2
        s2 = ck*s2 - s1
        s1 = st
        ck = ck + rx
      END DO
      IF( N==1 ) s1 = s2
    ELSEIF( N<=1 ) THEN
      s1 = s2
    END IF
  END IF
  100  Y(1) = s1
  IF( N==1 ) RETURN
  Y(2) = s2
  IF( N==2 ) RETURN
  DO i = 3, N
    Y(i) = ck*Y(i-1) - Y(i-2)
    ck = ck + rx
  END DO

  RETURN
END SUBROUTINE DBSYNU