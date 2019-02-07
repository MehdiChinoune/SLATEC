!*==DBSKNU.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBSKNU
SUBROUTINE DBSKNU(X,Fnu,Kode,N,Y,Nz)
  IMPLICIT NONE
  !*--DBSKNU5
  !***BEGIN PROLOGUE  DBSKNU
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBESK
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BESKNU-S, DBSKNU-D)
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract  **** A DOUBLE PRECISION routine ****
  !         DBSKNU computes N member sequences of K Bessel functions
  !         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
  !         positive X. Equations of the references are implemented on
  !         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X).
  !         Forward recursion with the three term recursion relation
  !         generates higher orders FNU+I-1, I=1,...,N. The parameter
  !         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values
  !         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned.
  !
  !         To start the recursion FNU is normalized to the interval
  !         -0.5.LE.DNU.LT.0.5. A special form of the power series is
  !         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the
  !         K Bessel function in terms of the confluent hypergeometric
  !         function U(FNU+0.5,2*FNU+1,X) is implemented on X1.LT.X.LE.X2.
  !         For X.GT.X2, the asymptotic expansion for large X is used.
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
  !           X      - X.GT.0.0D0
  !           FNU    - Order of initial K function, FNU.GE.0.0D0
  !           N      - Number of members of the sequence, N.GE.1
  !           KODE   - A parameter to indicate the scaling option
  !                    KODE= 1  returns
  !                             Y(I)=       K/SUB(FNU+I-1)/(X)
  !                                  I=1,...,N
  !                        = 2  returns
  !                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X)
  !                                  I=1,...,N
  !
  !         OUTPUT     Y is DOUBLE PRECISION
  !           Y      - A vector whose first N components contain values
  !                    for the sequence
  !                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or
  !                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N
  !                    depending on KODE
  !           NZ     - Number of components set to zero due to
  !                    underflow,
  !                    NZ= 0   , normal return
  !                    NZ.NE.0 , first NZ components of Y set to zero
  !                              due to underflow, Y(I)=0.0D0,I=1,...,NZ
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow - a fatal error
  !         Underflow with KODE=1 - a non-fatal error (NZ.NE.0)
  !
  !***SEE ALSO  DBESK
  !***REFERENCES  N. M. Temme, On the numerical evaluation of the modified
  !                 Bessel function of the third kind, Journal of
  !                 Computational Physics 19, (1975), pp. 324-337.
  !***ROUTINES CALLED  D1MACH, DGAMMA, I1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
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
  !***END PROLOGUE  DBSKNU
  !
  INTEGER i , iflag , inu , j , k , kk , Kode , koded , N , nn , Nz
  INTEGER I1MACH
  REAL(8) :: a , ak , a1 , a2 , b , bk , cc , ck , coef , cx , dk , &
    dnu , dnu2 , elim , etest , ex , f , fc , fhs , fk , &
    fks , flrx , fmu , Fnu , g1 , g2 , p , pi , pt , p1 , &
    p2 , q , rthpi , rx , s , smu , sqk , st , s1 , s2 , tm , &
    tol , t1 , t2 , X , x1 , x2 , Y
  DIMENSION a(160) , b(160) , Y(*) , cc(8)
  REAL(8) :: DGAMMA , D1MACH
  EXTERNAL DGAMMA
  SAVE x1 , x2 , pi , rthpi , cc
  DATA x1 , x2/2.0D0 , 17.0D0/
  DATA pi , rthpi/3.14159265358979D+00 , 1.25331413731550D+00/
  DATA cc(1) , cc(2) , cc(3) , cc(4) , cc(5) , cc(6) , cc(7) , &
    cc(8)/5.77215664901533D-01 , -4.20026350340952D-02 , &
    -4.21977345555443D-02 , 7.21894324666300D-03 , &
    -2.15241674114900D-04 , -2.01348547807000D-05 , &
    1.13302723200000D-06 , 6.11609500000000D-09/
  !***FIRST EXECUTABLE STATEMENT  DBSKNU
  kk = -I1MACH(15)
  elim = 2.303D0*(kk*D1MACH(5)-3.0D0)
  ak = D1MACH(3)
  tol = MAX(ak,1.0D-15)
  IF ( X<=0.0D0 ) THEN
    !
    !
    CALL XERMSG('SLATEC','DBSKNU','X NOT GREATER THAN ZERO',2,1)
    RETURN
  ELSEIF ( Fnu<0.0D0 ) THEN
    CALL XERMSG('SLATEC','DBSKNU','FNU NOT ZERO OR POSITIVE',2,1)
    RETURN
  ELSE
    IF ( Kode<1.OR.Kode>2 ) THEN
      CALL XERMSG('SLATEC','DBSKNU','KODE NOT 1 OR 2',2,1)
      RETURN
    ELSE
      IF ( N<1 ) THEN
        CALL XERMSG('SLATEC','DBSKNU','N NOT GREATER THAN 0',2,1)
        GOTO 99999
      ELSE
        Nz = 0
        iflag = 0
        koded = Kode
        rx = 2.0D0/X
        inu = INT(Fnu+0.5D0)
        dnu = Fnu - inu
        IF ( ABS(dnu)/=0.5D0 ) THEN
          dnu2 = 0.0D0
          IF ( ABS(dnu)>=tol ) dnu2 = dnu*dnu
          IF ( X<=x1 ) THEN
            !
            !     SERIES FOR X.LE.X1
            !
            a1 = 1.0D0 - dnu
            a2 = 1.0D0 + dnu
            t1 = 1.0D0/DGAMMA(a1)
            t2 = 1.0D0/DGAMMA(a2)
            IF ( ABS(dnu)>0.1D0 ) THEN
              g1 = (t1-t2)/(dnu+dnu)
            ELSE
              !     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
              s = cc(1)
              ak = 1.0D0
              DO k = 2 , 8
                ak = ak*dnu2
                tm = cc(k)*ak
                s = s + tm
                IF ( ABS(tm)<tol ) EXIT
              ENDDO
              g1 = -s
            ENDIF
            g2 = (t1+t2)*0.5D0
            smu = 1.0D0
            fc = 1.0D0
            flrx = LOG(rx)
            fmu = dnu*flrx
            IF ( dnu/=0.0D0 ) THEN
              fc = dnu*pi
              fc = fc/SIN(fc)
              IF ( fmu/=0.0D0 ) smu = SINH(fmu)/fmu
            ENDIF
            f = fc*(g1*COSH(fmu)+g2*flrx*smu)
            fc = EXP(fmu)
            p = 0.5D0*fc/t2
            q = 0.5D0/(fc*t1)
            ak = 1.0D0
            ck = 1.0D0
            bk = 1.0D0
            s1 = f
            s2 = p
            IF ( inu>0.OR.N>1 ) THEN
              IF ( X>=tol ) THEN
                cx = X*X*0.25D0
                DO
                  f = (ak*f+p+q)/(bk-dnu2)
                  p = p/(ak-dnu)
                  q = q/(ak+dnu)
                  ck = ck*cx/ak
                  t1 = ck*f
                  s1 = s1 + t1
                  t2 = ck*(p-ak*f)
                  s2 = s2 + t2
                  bk = bk + ak + ak + 1.0D0
                  ak = ak + 1.0D0
                  s = ABS(t1)/(1.0D0+ABS(s1)) + ABS(t2)/(1.0D0+ABS(s2))
                  IF ( s<=tol ) EXIT
                ENDDO
              ENDIF
              s2 = s2*rx
              IF ( koded/=1 ) THEN
                f = EXP(X)
                s1 = s1*f
                s2 = s2*f
              ENDIF
              GOTO 20
            ELSE
              IF ( X>=tol ) THEN
                cx = X*X*0.25D0
                DO
                  f = (ak*f+p+q)/(bk-dnu2)
                  p = p/(ak-dnu)
                  q = q/(ak+dnu)
                  ck = ck*cx/ak
                  t1 = ck*f
                  s1 = s1 + t1
                  bk = bk + ak + ak + 1.0D0
                  ak = ak + 1.0D0
                  s = ABS(t1)/(1.0D0+ABS(s1))
                  IF ( s<=tol ) EXIT
                ENDDO
              ENDIF
              Y(1) = s1
              IF ( koded==1 ) RETURN
              Y(1) = s1*EXP(X)
              RETURN
            ENDIF
          ENDIF
        ENDIF
        DO
          coef = rthpi/SQRT(X)
          IF ( koded==2 ) EXIT
          IF ( X>elim ) THEN
            koded = 2
            iflag = 1
          ELSE
            coef = coef*EXP(-X)
            EXIT
          ENDIF
        ENDDO
        IF ( ABS(dnu)==0.5D0 ) THEN
          !
          !     FNU=HALF ODD INTEGER CASE
          !
          s1 = coef
          s2 = coef
        ELSEIF ( X>x2 ) THEN
          !
          !     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2
          !
          !     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
          !     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
          !     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
          !     RECURSION
          nn = 2
          IF ( inu==0.AND.N==1 ) nn = 1
          dnu2 = dnu + dnu
          fmu = 0.0D0
          IF ( ABS(dnu2)>=tol ) fmu = dnu2*dnu2
          ex = X*8.0D0
          s2 = 0.0D0
          DO k = 1 , nn
            s1 = s2
            s = 1.0D0
            ak = 0.0D0
            ck = 1.0D0
            sqk = 1.0D0
            dk = ex
            DO j = 1 , 30
              ck = ck*(fmu-sqk)/dk
              s = s + ck
              dk = dk + ex
              ak = ak + 8.0D0
              sqk = sqk + ak
              IF ( ABS(ck)<tol ) EXIT
            ENDDO
            s2 = s*coef
            fmu = fmu + 8.0D0*dnu + 4.0D0
          ENDDO
          IF ( nn<=1 ) THEN
            s1 = s2
            GOTO 50
          ENDIF
        ELSE
          !
          !     MILLER ALGORITHM FOR X1.LT.X.LE.X2
          !
          etest = COS(pi*dnu)/(pi*X*tol)
          fks = 1.0D0
          fhs = 0.25D0
          fk = 0.0D0
          ck = X + X + 2.0D0
          p1 = 0.0D0
          p2 = 1.0D0
          k = 0
          DO
            k = k + 1
            fk = fk + 1.0D0
            ak = (fhs-dnu2)/(fks+fk)
            bk = ck/(fk+1.0D0)
            pt = p2
            p2 = bk*p2 - ak*p1
            p1 = pt
            a(k) = ak
            b(k) = bk
            ck = ck + 2.0D0
            fks = fks + fk + fk + 1.0D0
            fhs = fhs + fk + fk
            IF ( etest<=fk*p1 ) THEN
              kk = k
              s = 1.0D0
              p1 = 0.0D0
              p2 = 1.0D0
              DO i = 1 , k
                pt = p2
                p2 = (b(kk)*p2-p1)/a(kk)
                p1 = pt
                s = s + p2
                kk = kk - 1
              ENDDO
              s1 = coef*(p2/s)
              IF ( inu<=0.AND.N<=1 ) GOTO 50
              s2 = s1*(X+dnu+0.5D0-p1/p2)/X
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      !
      !     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
      !
      20       ck = (dnu+dnu+2.0D0)/X
      IF ( N==1 ) inu = inu - 1
      IF ( inu>0 ) THEN
        DO i = 1 , inu
          st = s2
          s2 = ck*s2 + s1
          s1 = st
          ck = ck + rx
        ENDDO
        IF ( N==1 ) s1 = s2
      ELSEIF ( N<=1 ) THEN
        s1 = s2
      ENDIF
    ENDIF
    50     IF ( iflag==1 ) THEN
    !     IFLAG=1 CASES
    s = -X + LOG(s1)
    Y(1) = 0.0D0
    Nz = 1
    IF ( s>=-elim ) THEN
      Y(1) = EXP(s)
      Nz = 0
    ENDIF
    IF ( N==1 ) RETURN
    s = -X + LOG(s2)
    Y(2) = 0.0D0
    Nz = Nz + 1
    IF ( s>=-elim ) THEN
      Nz = Nz - 1
      Y(2) = EXP(s)
    ENDIF
    IF ( N==2 ) RETURN
    kk = 2
    IF ( Nz>=2 ) THEN
      DO i = 3 , N
        kk = i
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        ck = ck + rx
        s = -X + LOG(s2)
        Nz = Nz + 1
        Y(i) = 0.0D0
        IF ( s>=-elim ) THEN
          Y(i) = EXP(s)
          Nz = Nz - 1
          GOTO 100
        ENDIF
      ENDDO
      RETURN
    ENDIF
  ELSE
    Y(1) = s1
    IF ( N==1 ) RETURN
    Y(2) = s2
    IF ( N==2 ) RETURN
    DO i = 3 , N
      Y(i) = ck*Y(i-1) + Y(i-2)
      ck = ck + rx
    ENDDO
    RETURN
  ENDIF
ENDIF
100  IF ( kk==N ) RETURN
s2 = s2*ck + s1
ck = ck + rx
kk = kk + 1
Y(kk) = EXP(-X+LOG(s2))
IF ( kk==N ) RETURN
kk = kk + 1
DO i = kk , N
  Y(i) = ck*Y(i-1) + Y(i-2)
  ck = ck + rx
ENDDO
RETURN
99999 END SUBROUTINE DBSKNU
