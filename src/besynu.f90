!*==BESYNU.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BESYNU
      SUBROUTINE BESYNU(X,Fnu,N,Y)
      IMPLICIT NONE
!*--BESYNU5
!***BEGIN PROLOGUE  BESYNU
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BESY
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BESYNU-S, DBSYNU-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         BESYNU computes N member sequences of Y Bessel functions
!         Y/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
!         positive X. Equations of the references are implemented on
!         small orders DNU for Y/SUB(DNU)/(X) and Y/SUB(DNU+1)/(X).
!         Forward recursion with the three term recursion relation
!         generates higher orders FNU+I-1, I=1,...,N.
!
!         To start the recursion FNU is normalized to the interval
!         -0.5.LE.DNU.LT.0.5. A special form of the power series is
!         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the
!         K Bessel function in terms of the confluent hypergeometric
!         function U(FNU+0.5,2*FNU+1,I*X) is implemented on X1.LT.X.LE.X
!         Here I is the complex number SQRT(-1.).
!         For X.GT.X2, the asymptotic expansion for large X is used.
!         When FNU is a half odd integer, a special formula for
!         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
!
!         BESYNU assumes that a significant digit SINH(X) function is
!         available.
!
!     Description of Arguments
!
!         Input
!           X      - X.GT.0.0E0
!           FNU    - Order of initial Y function, FNU.GE.0.0E0
!           N      - Number of members of the sequence, N.GE.1
!
!         Output
!           Y      - A vector whose first N components contain values
!                    for the sequence Y(I)=Y/SUB(FNU+I-1), I=1,N.
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!
!***SEE ALSO  BESY
!***REFERENCES  N. M. Temme, On the numerical evaluation of the ordinary
!                 Bessel function of the second kind, Journal of
!                 Computational Physics 21, (1976), pp. 343-350.
!               N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  GAMMA, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BESYNU
!
      INTEGER i , inu , j , k , kk , N , nn
      REAL a , ak , arg , a1 , a2 , bk , cb , cbk , cc , cck , ck , coef , cpt , 
     &     cp1 , cp2 , cs , cs1 , cs2 , cx , dnu , dnu2 , etest , etx , f , fc , 
     &     fhs , fk , fks , flrx , fmu , fn , Fnu , fx , g , g1 , g2 , hpi , p , 
     &     pi , pt , q , rb , rbk , rck , relb , rpt , rp1 , rp2 , rs , rs1 , 
     &     rs2 , rthpi , rx , s , sa , sb , smu , ss , st , s1 , s2 , tb , tm , 
     &     tol , t1 , t2 , X , x1 , x2 , Y
      DIMENSION a(120) , rb(120) , cb(120) , Y(*) , cc(8)
      REAL GAMMA , R1MACH
      EXTERNAL GAMMA
      SAVE x1 , x2 , pi , rthpi , hpi , cc
      DATA x1 , x2/3.0E0 , 20.0E0/
      DATA pi , rthpi/3.14159265358979E+00 , 7.97884560802865E-01/
      DATA hpi/1.57079632679490E+00/
      DATA cc(1) , cc(2) , cc(3) , cc(4) , cc(5) , cc(6) , cc(7) , 
     &     cc(8)/5.77215664901533E-01 , -4.20026350340952E-02 , 
     &     -4.21977345555443E-02 , 7.21894324666300E-03 , 
     &     -2.15241674114900E-04 , -2.01348547807000E-05 , 
     &     1.13302723200000E-06 , 6.11609500000000E-09/
!***FIRST EXECUTABLE STATEMENT  BESYNU
      ak = R1MACH(3)
      tol = MAX(ak,1.0E-15)
      IF ( X<=0.0E0 ) THEN
!
!
        CALL XERMSG('SLATEC','BESYNU','X NOT GREATER THAN ZERO',2,1)
        RETURN
      ELSEIF ( Fnu<0.0E0 ) THEN
        CALL XERMSG('SLATEC','BESYNU','FNU NOT ZERO OR POSITIVE',2,1)
        RETURN
      ELSEIF ( N<1 ) THEN
        CALL XERMSG('SLATEC','BESYNU','N NOT GREATER THAN 0',2,1)
        GOTO 99999
      ELSE
        rx = 2.0E0/X
        inu = INT(Fnu+0.5E0)
        dnu = Fnu - inu
        IF ( ABS(dnu)==0.5E0 ) THEN
!
!     FNU=HALF ODD INTEGER CASE
!
          coef = rthpi/SQRT(X)
          s1 = coef*SIN(X)
          s2 = -coef*COS(X)
        ELSE
          dnu2 = 0.0E0
          IF ( ABS(dnu)>=tol ) dnu2 = dnu*dnu
          IF ( X>x1 ) THEN
            coef = rthpi/SQRT(X)
            IF ( X>x2 ) THEN
!
!     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2
!
              nn = 2
              IF ( inu==0.AND.N==1 ) nn = 1
              dnu2 = dnu + dnu
              fmu = 0.0E0
              IF ( ABS(dnu2)>=tol ) fmu = dnu2*dnu2
              arg = X - hpi*(dnu+0.5E0)
              sa = SIN(arg)
              sb = COS(arg)
              etx = 8.0E0*X
              DO k = 1 , nn
                s1 = s2
                t2 = (fmu-1.0E0)/etx
                ss = t2
                relb = tol*ABS(t2)
                t1 = etx
                s = 1.0E0
                fn = 1.0E0
                ak = 0.0E0
                DO j = 1 , 13
                  t1 = t1 + etx
                  ak = ak + 8.0E0
                  fn = fn + ak
                  t2 = -t2*(fmu-fn)/t1
                  s = s + t2
                  t1 = t1 + etx
                  ak = ak + 8.0E0
                  fn = fn + ak
                  t2 = t2*(fmu-fn)/t1
                  ss = ss + t2
                  IF ( ABS(t2)<=relb ) EXIT
                ENDDO
                s2 = coef*(s*sa+ss*sb)
                fmu = fmu + 8.0E0*dnu + 4.0E0
                tb = sa
                sa = -sb
                sb = tb
              ENDDO
              IF ( nn<=1 ) THEN
                s1 = s2
                GOTO 100
              ENDIF
            ELSE
!
!     MILLER ALGORITHM FOR X1.LT.X.LE.X2
!
              etest = COS(pi*dnu)/(pi*X*tol)
              fks = 1.0E0
              fhs = 0.25E0
              fk = 0.0E0
              rck = 2.0E0
              cck = X + X
              rp1 = 0.0E0
              cp1 = 0.0E0
              rp2 = 1.0E0
              cp2 = 0.0E0
              k = 0
              DO
                k = k + 1
                fk = fk + 1.0E0
                ak = (fhs-dnu2)/(fks+fk)
                pt = fk + 1.0E0
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
                rck = rck + 2.0E0
                fks = fks + fk + fk + 1.0E0
                fhs = fhs + fk + fk
                pt = MAX(ABS(rp1),ABS(cp1))
                fc = (rp1/pt)**2 + (cp1/pt)**2
                pt = pt*SQRT(fc)*fk
                IF ( etest<=pt ) THEN
                  kk = k
                  rs = 1.0E0
                  cs = 0.0E0
                  rp1 = 0.0E0
                  cp1 = 0.0E0
                  rp2 = 1.0E0
                  cp2 = 0.0E0
                  DO i = 1 , k
                    rpt = rp2
                    cpt = cp2
                    rp2 = (rb(kk)*rpt-cb(kk)*cpt-rp1)/a(kk)
                    cp2 = (cb(kk)*rpt+rb(kk)*cpt-cp1)/a(kk)
                    rp1 = rpt
                    cp1 = cpt
                    rs = rs + rp2
                    cs = cs + cp2
                    kk = kk - 1
                  ENDDO
                  pt = MAX(ABS(rs),ABS(cs))
                  fc = (rs/pt)**2 + (cs/pt)**2
                  pt = pt*SQRT(fc)
                  rs1 = (rp2*(rs/pt)+cp2*(cs/pt))/pt
                  cs1 = (cp2*(rs/pt)-rp2*(cs/pt))/pt
                  fc = hpi*(dnu-0.5E0) - X
                  p = COS(fc)
                  q = SIN(fc)
                  s1 = (cs1*q-rs1*p)*coef
                  IF ( inu>0.OR.N>1 ) THEN
                    pt = MAX(ABS(rp2),ABS(cp2))
                    fc = (rp2/pt)**2 + (cp2/pt)**2
                    pt = pt*SQRT(fc)
                    rpt = dnu + 0.5E0 - (rp1*(rp2/pt)+cp1*(cp2/pt))/pt
                    cpt = X - (cp1*(rp2/pt)-rp1*(cp2/pt))/pt
                    cs2 = cs1*cpt - rs1*rpt
                    rs2 = rpt*cs1 + rs1*cpt
                    s2 = (rs2*q+cs2*p)*coef/X
                    EXIT
                  ELSE
                    Y(1) = s1
                    RETURN
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ELSE
!
!     SERIES FOR X.LE.X1
!
            a1 = 1.0E0 - dnu
            a2 = 1.0E0 + dnu
            t1 = 1.0E0/GAMMA(a1)
            t2 = 1.0E0/GAMMA(a2)
            IF ( ABS(dnu)>0.1E0 ) THEN
              g1 = (t1-t2)/dnu
            ELSE
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
              s = cc(1)
              ak = 1.0E0
              DO k = 2 , 8
                ak = ak*dnu2
                tm = cc(k)*ak
                s = s + tm
                IF ( ABS(tm)<tol ) EXIT
              ENDDO
              g1 = -(s+s)
            ENDIF
            g2 = t1 + t2
            smu = 1.0E0
            fc = 1.0E0/pi
            flrx = LOG(rx)
            fmu = dnu*flrx
            tm = 0.0E0
            IF ( dnu/=0.0E0 ) THEN
              tm = SIN(dnu*hpi)/dnu
              tm = (dnu+dnu)*tm*tm
              fc = dnu/SIN(dnu*pi)
              IF ( fmu/=0.0E0 ) smu = SINH(fmu)/fmu
            ENDIF
            f = fc*(g1*COSH(fmu)+g2*flrx*smu)
            fx = EXP(fmu)
            p = fc*t1*fx
            q = fc*t2/fx
            g = f + tm*q
            ak = 1.0E0
            ck = 1.0E0
            bk = 1.0E0
            s1 = g
            s2 = p
            IF ( inu>0.OR.N>1 ) THEN
              IF ( X>=tol ) THEN
                cx = X*X*0.25E0
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
                  bk = bk + ak + ak + 1.0E0
                  ak = ak + 1.0E0
                  s = ABS(t1)/(1.0E0+ABS(s1)) + ABS(t2)/(1.0E0+ABS(s2))
                  IF ( s<=tol ) EXIT
                ENDDO
              ENDIF
              s2 = -s2*rx
              s1 = -s1
            ELSE
              IF ( X>=tol ) THEN
                cx = X*X*0.25E0
                DO
                  f = (ak*f+p+q)/(bk-dnu2)
                  p = p/(ak-dnu)
                  q = q/(ak+dnu)
                  g = f + tm*q
                  ck = -ck*cx/ak
                  t1 = ck*g
                  s1 = s1 + t1
                  bk = bk + ak + ak + 1.0E0
                  ak = ak + 1.0E0
                  s = ABS(t1)/(1.0E0+ABS(s1))
                  IF ( s<=tol ) EXIT
                ENDDO
              ENDIF
              Y(1) = -s1
              RETURN
            ENDIF
          ENDIF
        ENDIF
!
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
!
        ck = (dnu+dnu+2.0E0)/X
        IF ( N==1 ) inu = inu - 1
        IF ( inu>0 ) THEN
          DO i = 1 , inu
            st = s2
            s2 = ck*s2 - s1
            s1 = st
            ck = ck + rx
          ENDDO
          IF ( N==1 ) s1 = s2
        ELSEIF ( N<=1 ) THEN
          s1 = s2
        ENDIF
      ENDIF
 100  Y(1) = s1
      IF ( N==1 ) RETURN
      Y(2) = s2
      IF ( N==2 ) RETURN
      DO i = 3 , N
        Y(i) = ck*Y(i-1) - Y(i-2)
        ck = ck + rx
      ENDDO
      RETURN
99999 END SUBROUTINE BESYNU
