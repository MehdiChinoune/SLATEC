!** DBSKNU
PURE SUBROUTINE DBSKNU(X,Fnu,Kode,N,Y,Nz)
  !> Subsidiary to DBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (BESKNU-S, DBSKNU-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract  **** A DOUBLE PRECISION routine ****
  !         DBSKNU computes N member sequences of K Bessel functions
  !         K_{FNU+I-1}(X), I=1,N for non-negative orders FNU and
  !         positive X. Equations of the references are implemented on
  !         small orders DNU for K_{DNU}(X) and K_{DNU+1}(X).
  !         Forward recursion with the three term recursion relation
  !         generates higher orders FNU+I-1, I=1,...,N. The parameter
  !         KODE permits K_{FNU+I-1}(X) values or scaled values
  !         EXP(X)*K_{FNU+I-1}(X), I=1,N to be returned.
  !
  !         To start the recursion FNU is normalized to the interval
  !         -0.5<=DNU<0.5. A special form of the power series is
  !         implemented on 0<X<=X1 while the Miller algorithm for the
  !         K Bessel function in terms of the confluent hypergeometric
  !         function U(FNU+0.5,2*FNU+1,X) is implemented on X1<X<=X2.
  !         For X>X2, the asymptotic expansion for large X is used.
  !         When FNU is a half odd integer, a special formula for
  !         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
  !
  !         The maximum number of significant digits obtainable
  !         is the smaller of 14 and the number of digits carried in
  !         DOUBLE PRECISION arithmetic.
  !
  !         DBSKNU assumes that a significant digit SINH function is
  !         available.
  !
  !     Description of Arguments
  !
  !         INPUT      X,FNU are DOUBLE PRECISION
  !           X      - X>0.0D0
  !           FNU    - Order of initial K function, FNU>=0.0D0
  !           N      - Number of members of the sequence, N>=1
  !           KODE   - A parameter to indicate the scaling option
  !                    KODE= 1  returns
  !                             Y(I)=       K_{FNU+I-1}(X)
  !                                  I=1,...,N
  !                        = 2  returns
  !                             Y(I)=EXP(X)*K_{FNU+I-1}(X)
  !                                  I=1,...,N
  !
  !         OUTPUT     Y is DOUBLE PRECISION
  !           Y      - A vector whose first N components contain values
  !                    for the sequence
  !                    Y(I)=       K_{FNU+I-1}(X), I=1,...,N or
  !                    Y(I)=EXP(X)*K_{FNU+I-1}(X), I=1,...,N
  !                    depending on KODE
  !           NZ     - Number of components set to zero due to
  !                    underflow,
  !                    NZ= 0  , normal return
  !                    NZ/=0, first NZ components of Y set to zero
  !                              due to underflow, Y(I)=0.0D0,I=1,...,NZ
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow - a fatal error
  !         Underflow with KODE=1 - a non-fatal error (NZ/=0)
  !
  !***
  ! **See also:**  DBESK
  !***
  ! **References:**  N. M. Temme, On the numerical evaluation of the modified
  !                 Bessel function of the third kind, Journal of
  !                 Computational Physics 19, (1975), pp. 324-337.
  !***
  ! **Routines called:**  D1MACH, I1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : D1MACH, I1MACH
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Fnu, X
  REAL(DP), INTENT(OUT) :: Y(N)
  INTEGER :: i, iflag, inu, j, k, kk, koded, nn
  REAL(DP) :: a(160), ak, a1, a2, b(160), bk, ck, coef, cx, dk, dnu, dnu2, elim, &
    etest, ex, f, fc, fhs, fk, fks, flrx, fmu, g1, g2, p, pt, p1, p2, q, &
    rx, s, smu, sqk, st, s1, s2, tm, tol, t1, t2
  REAL(DP), PARAMETER ::  x1 = 2._DP, x2 = 17._DP
  REAL(DP), PARAMETER ::  pi = 3.14159265358979E+00_DP, rthpi = 1.25331413731550_DP
  REAL(DP), PARAMETER :: cc(8) = [ 5.77215664901533E-01_DP, -4.20026350340952E-02_DP, &
    -4.21977345555443E-02_DP, 7.21894324666300E-03_DP, -2.15241674114900E-04_DP, &
    -2.01348547807000E-05_DP, 1.13302723200000E-06_DP, 6.11609500000000E-09_DP ]
  !* FIRST EXECUTABLE STATEMENT  DBSKNU
  kk = -I1MACH(15)
  elim = 2.303_DP*(kk*D1MACH(5)-3._DP)
  ak = D1MACH(3)
  tol = MAX(ak,1.E-15_DP)
  IF( X<=0._DP ) THEN
    ERROR STOP 'DBSKNU : X <= 0'
  ELSEIF( Fnu<0._DP ) THEN
    ERROR STOP 'DBSKNU : FNU < 0'
  ELSEIF( Kode<1 .OR. Kode>2 ) THEN
    ERROR STOP 'DBSKNU : KODE NOT 1 OR 2'
  ELSEIF( N<1 ) THEN
    ERROR STOP 'DBSKNU : N NOT GREATER THAN 0'
  ELSE
    Nz = 0
    iflag = 0
    koded = Kode
    rx = 2._DP/X
    inu = INT(Fnu+0.5_DP)
    dnu = Fnu - inu
    IF( ABS(dnu)/=0.5_DP ) THEN
      dnu2 = 0._DP
      IF( ABS(dnu)>=tol ) dnu2 = dnu*dnu
      IF( X<=x1 ) THEN
        !
        !     SERIES FOR X<=X1
        !
        a1 = 1._DP - dnu
        a2 = 1._DP + dnu
        t1 = 1._DP/GAMMA(a1)
        t2 = 1._DP/GAMMA(a2)
        IF( ABS(dnu)>0.1_DP ) THEN
          g1 = (t1-t2)/(dnu+dnu)
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
          g1 = -s
        END IF
        g2 = (t1+t2)*0.5_DP
        smu = 1._DP
        fc = 1._DP
        flrx = LOG(rx)
        fmu = dnu*flrx
        IF( dnu/=0._DP ) THEN
          fc = dnu*pi
          fc = fc/SIN(fc)
          IF( fmu/=0._DP ) smu = SINH(fmu)/fmu
        END IF
        f = fc*(g1*COSH(fmu)+g2*flrx*smu)
        fc = EXP(fmu)
        p = 0.5_DP*fc/t2
        q = 0.5_DP/(fc*t1)
        ak = 1._DP
        ck = 1._DP
        bk = 1._DP
        s1 = f
        s2 = p
        IF( inu>0 .OR. N>1 ) THEN
          IF( X>=tol ) THEN
            cx = X*X*0.25_DP
            DO
              f = (ak*f+p+q)/(bk-dnu2)
              p = p/(ak-dnu)
              q = q/(ak+dnu)
              ck = ck*cx/ak
              t1 = ck*f
              s1 = s1 + t1
              t2 = ck*(p-ak*f)
              s2 = s2 + t2
              bk = bk + ak + ak + 1._DP
              ak = ak + 1._DP
              s = ABS(t1)/(1._DP+ABS(s1)) + ABS(t2)/(1._DP+ABS(s2))
              IF( s<=tol ) EXIT
            END DO
          END IF
          s2 = s2*rx
          IF( koded/=1 ) THEN
            f = EXP(X)
            s1 = s1*f
            s2 = s2*f
          END IF
          GOTO 20
        ELSE
          IF( X>=tol ) THEN
            cx = X*X*0.25_DP
            DO
              f = (ak*f+p+q)/(bk-dnu2)
              p = p/(ak-dnu)
              q = q/(ak+dnu)
              ck = ck*cx/ak
              t1 = ck*f
              s1 = s1 + t1
              bk = bk + ak + ak + 1._DP
              ak = ak + 1._DP
              s = ABS(t1)/(1._DP+ABS(s1))
              IF( s<=tol ) EXIT
            END DO
          END IF
          Y(1) = s1
          IF( koded==1 ) RETURN
          Y(1) = s1*EXP(X)
          RETURN
        END IF
      END IF
    END IF
    DO
      coef = rthpi/SQRT(X)
      IF( koded==2 ) EXIT
      IF( X>elim ) THEN
        koded = 2
        iflag = 1
      ELSE
        coef = coef*EXP(-X)
        EXIT
      END IF
    END DO
    IF( ABS(dnu)==0.5_DP ) THEN
      !
      !     FNU=HALF ODD INTEGER CASE
      !
      s1 = coef
      s2 = coef
    ELSEIF( X>x2 ) THEN
      !
      !     ASYMPTOTIC EXPANSION FOR LARGE X, X>X2
      !
      !     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
      !     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
      !     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
      !     RECURSION
      nn = 2
      IF( inu==0 .AND. N==1 ) nn = 1
      dnu2 = dnu + dnu
      fmu = 0._DP
      IF( ABS(dnu2)>=tol ) fmu = dnu2*dnu2
      ex = X*8._DP
      s2 = 0._DP
      DO k = 1, nn
        s1 = s2
        s = 1._DP
        ak = 0._DP
        ck = 1._DP
        sqk = 1._DP
        dk = ex
        DO j = 1, 30
          ck = ck*(fmu-sqk)/dk
          s = s + ck
          dk = dk + ex
          ak = ak + 8._DP
          sqk = sqk + ak
          IF( ABS(ck)<tol ) EXIT
        END DO
        s2 = s*coef
        fmu = fmu + 8._DP*dnu + 4._DP
      END DO
      IF( nn<=1 ) THEN
        s1 = s2
        GOTO 50
      END IF
    ELSE
      !
      !     MILLER ALGORITHM FOR X1<X<=X2
      !
      etest = COS(pi*dnu)/(pi*X*tol)
      fks = 1._DP
      fhs = 0.25_DP
      fk = 0._DP
      ck = X + X + 2._DP
      p1 = 0._DP
      p2 = 1._DP
      k = 0
      DO
        k = k + 1
        fk = fk + 1._DP
        ak = (fhs-dnu2)/(fks+fk)
        bk = ck/(fk+1._DP)
        pt = p2
        p2 = bk*p2 - ak*p1
        p1 = pt
        a(k) = ak
        b(k) = bk
        ck = ck + 2._DP
        fks = fks + fk + fk + 1._DP
        fhs = fhs + fk + fk
        IF( etest<=fk*p1 ) THEN
          kk = k
          s = 1._DP
          p1 = 0._DP
          p2 = 1._DP
          DO i = 1, k
            pt = p2
            p2 = (b(kk)*p2-p1)/a(kk)
            p1 = pt
            s = s + p2
            kk = kk - 1
          END DO
          s1 = coef*(p2/s)
          IF( inu<=0 .AND. N<=1 ) GOTO 50
          s2 = s1*(X+dnu+0.5_DP-p1/p2)/X
          EXIT
        END IF
      END DO
    END IF
    !
    !     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
    !
    20  ck = (dnu+dnu+2._DP)/X
    IF( N==1 ) inu = inu - 1
    IF( inu>0 ) THEN
      DO i = 1, inu
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        ck = ck + rx
      END DO
      IF( N==1 ) s1 = s2
    ELSEIF( N<=1 ) THEN
      s1 = s2
    END IF
    50 CONTINUE
    IF( iflag==1 ) THEN
      !     IFLAG=1 CASES
      s = -X + LOG(s1)
      Y(1) = 0._DP
      Nz = 1
      IF( s>=-elim ) THEN
        Y(1) = EXP(s)
        Nz = 0
      END IF
      IF( N==1 ) RETURN
      s = -X + LOG(s2)
      Y(2) = 0._DP
      Nz = Nz + 1
      IF( s>=-elim ) THEN
        Nz = Nz - 1
        Y(2) = EXP(s)
      END IF
      IF( N==2 ) RETURN
      kk = 2
      IF( Nz>=2 ) THEN
        DO i = 3, N
          kk = i
          st = s2
          s2 = ck*s2 + s1
          s1 = st
          ck = ck + rx
          s = -X + LOG(s2)
          Nz = Nz + 1
          Y(i) = 0._DP
          IF( s>=-elim ) THEN
            Y(i) = EXP(s)
            Nz = Nz - 1
            GOTO 100
          END IF
        END DO
        RETURN
      END IF
    ELSE
      Y(1) = s1
      IF( N==1 ) RETURN
      Y(2) = s2
      IF( N==2 ) RETURN
      DO i = 3, N
        Y(i) = ck*Y(i-1) + Y(i-2)
        ck = ck + rx
      END DO
      RETURN
    END IF
  END IF
  100 CONTINUE
  IF( kk==N ) RETURN
  s2 = s2*ck + s1
  ck = ck + rx
  kk = kk + 1
  Y(kk) = EXP(-X+LOG(s2))
  IF( kk==N ) RETURN
  kk = kk + 1
  DO i = kk, N
    Y(i) = ck*Y(i-1) + Y(i-2)
    ck = ck + rx
  END DO

  RETURN
END SUBROUTINE DBSKNU