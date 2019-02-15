!DECK CBIRY
SUBROUTINE CBIRY(Z,Id,Kode,Bi,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CBIRY
  !***PURPOSE  Compute the Airy function Bi(z) or its derivative dBi/dz
  !            for complex argument z.  A scaling option is available
  !            to help avoid overflow.
  !***LIBRARY   SLATEC
  !***CATEGORY  C10D
  !***TYPE      COMPLEX (CBIRY-C, ZBIRY-C)
  !***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
  !             BESSEL FUNCTION OF ORDER TWO THIRDS
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !         On KODE=1, CBIRY computes the complex Airy function Bi(z)
  !         or its derivative dBi/dz on ID=0 or ID=1 respectively.
  !         On KODE=2, a scaling option exp(abs(Re(zeta)))*Bi(z) or
  !         exp(abs(Re(zeta)))*dBi/dz is provided to remove the
  !         exponential behavior in both the left and right half planes
  !         where zeta=(2/3)*z**(3/2).
  !
  !         The Airy functions Bi(z) and dBi/dz are analytic in the
  !         whole z-plane, and the scaling option does not destroy this
  !         property.
  !
  !         Input
  !           Z      - Argument of type COMPLEX
  !           ID     - Order of derivative, ID=0 or ID=1
  !           KODE   - A parameter to indicate the scaling option
  !                    KODE=1  returns
  !                            BI=Bi(z)  on ID=0
  !                            BI=dBi/dz on ID=1
  !                            at z=Z
  !                        =2  returns
  !                            BI=exp(abs(Re(zeta)))*Bi(z)  on ID=0
  !                            BI=exp(abs(Re(zeta)))*dBi/dz on ID=1
  !                            at z=Z where zeta=(2/3)*z**(3/2)
  !
  !         Output
  !           BI     - Result of type COMPLEX
  !           IERR   - Error flag
  !                    IERR=0  Normal return     - COMPUTATION COMPLETED
  !                    IERR=1  Input error       - NO COMPUTATION
  !                    IERR=2  Overflow          - NO COMPUTATION
  !                            (Re(Z) too large with KODE=1)
  !                    IERR=3  Precision warning - COMPUTATION COMPLETED
  !                            (Result has less than half precision)
  !                    IERR=4  Precision error   - NO COMPUTATION
  !                            (Result has no precision)
  !                    IERR=5  Algorithmic error - NO COMPUTATION
  !                            (Termination condition not met)
  !
  ! *Long Description:
  !
  !         Bi(z) and dBi/dz are computed from I Bessel functions by
  !
  !                Bi(z) =  c*sqrt(z)*( I(-1/3,zeta) + I(1/3,zeta) )
  !               dBi/dz =  c*   z   *( I(-2/3,zeta) + I(2/3,zeta) )
  !                    c =  1/sqrt(3)
  !                 zeta =  (2/3)*z**(3/2)
  !
  !         when abs(z)>1 and from power series when abs(z)<=1.
  !
  !         In most complex variable computation, one must evaluate ele-
  !         mentary functions.  When the magnitude of Z is large, losses
  !         of significance by argument reduction occur.  Consequently, if
  !         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR),
  !         then losses exceeding half precision are likely and an error
  !         flag IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.
  !         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then
  !         all significance is lost and IERR=4.  In order to use the INT
  !         function, ZETA must be further restricted not to exceed
  !         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA
  !         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2,
  !         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single
  !         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision.
  !         This makes U2 limiting is single precision and U3 limiting
  !         in double precision.  This means that the magnitude of Z
  !         cannot exceed approximately 3.4E+4 in single precision and
  !         2.1E+6 in double precision.  This also means that one can
  !         expect to retain, in the worst cases on 32-bit machines,
  !         no digits in single precision and only 6 digits in double
  !         precision.
  !
  !         The approximate relative error in the magnitude of a complex
  !         Bessel function can be expressed as P*10**S where P=MAX(UNIT
  !         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
  !         sents the increase in error due to argument reduction in the
  !         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
  !         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
  !         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
  !         have only absolute accuracy.  This is most likely to occur
  !         when one component (in magnitude) is larger than the other by
  !         several orders of magnitude.  If one component is 10**K larger
  !         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
  !         0) significant digits; or, stated another way, when K exceeds
  !         the exponent of P, no significant digits remain in the smaller
  !         component.  However, the phase angle retains absolute accuracy
  !         because, in complex arithmetic with precision P, the smaller
  !         component will not (as a rule) decrease below P times the
  !         magnitude of the larger component. In these extreme cases,
  !         the principal phase angle is on the order of +P, -P, PI/2-P,
  !         or -PI/2+P.
  !
  !***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
  !                 matical Functions, National Bureau of Standards
  !                 Applied Mathematics Series 55, U. S. Department
  !                 of Commerce, Tenth Printing (1972) or later.
  !               2. D. E. Amos, Computation of Bessel Functions of
  !                 Complex Argument and Large Order, Report SAND83-0643,
  !                 Sandia National Laboratories, Albuquerque, NM, May
  !                 1983.
  !               3. D. E. Amos, A Subroutine Package for Bessel Functions
  !                 of a Complex Argument and Nonnegative Order, Report
  !                 SAND85-1018, Sandia National Laboratory, Albuquerque,
  !                 NM, May 1985.
  !               4. D. E. Amos, A portable package for Bessel functions
  !                 of a complex argument and nonnegative order, ACM
  !                 Transactions on Mathematical Software, 12 (September
  !                 1986), pp. 265-273.
  !
  !***ROUTINES CALLED  CBINU, I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890801  REVISION DATE from Version 3.2
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   920128  Category corrected.  (WRB)
  !   920811  Prologue revised.  (DWL)
  !***END PROLOGUE  CBIRY
  COMPLEX Bi, cone, csq, cy, s1, s2, trm1, trm2, Z, zta, z3
  REAL aa, ad, ak, alim, atrm, az, az3, bb, bk, ck, coef, c1, &
    c2, dig, dk, d1, d2, elim, fid, fmr, fnu, fnul, pi, rl, &
    r1m5, sfac, tol, tth, zi, zr, z3i, z3r, R1MACH
  INTEGER Id, Ierr, k, Kode, k1, k2, nz, I1MACH
  DIMENSION cy(2)
  DATA tth, c1, c2, coef, pi/6.66666666666666667E-01, &
    6.14926627446000736E-01, 4.48288357353826359E-01, &
    5.77350269189625765E-01, 3.14159265358979324E+00/
  DATA cone/(1.0E0,0.0E0)/
  !***FIRST EXECUTABLE STATEMENT  CBIRY
  Ierr = 0
  nz = 0
  IF ( Id<0.OR.Id>1 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  az = ABS(Z)
  tol = MAX(R1MACH(4),1.0E-18)
  fid = Id
  IF ( az>1.0E0 ) THEN
    !-----------------------------------------------------------------------
    !     CASE FOR ABS(Z).GT.1.0
    !-----------------------------------------------------------------------
    fnu = (1.0E0+fid)/3.0E0
    !-----------------------------------------------------------------------
    !     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    !     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    !     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    !     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    !     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    !     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    !     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    !     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    !     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
    !-----------------------------------------------------------------------
    k1 = I1MACH(12)
    k2 = I1MACH(13)
    r1m5 = R1MACH(5)
    k = MIN(ABS(k1),ABS(k2))
    elim = 2.303E0*(k*r1m5-3.0E0)
    k1 = I1MACH(11) - 1
    aa = r1m5*k1
    dig = MIN(aa,18.0E0)
    aa = aa*2.303E0
    alim = elim + MAX(-aa,-41.45E0)
    rl = 1.2E0*dig + 3.0E0
    fnul = 10.0E0 + 6.0E0*(dig-3.0E0)
    !-----------------------------------------------------------------------
    !     TEST FOR RANGE
    !-----------------------------------------------------------------------
    aa = 0.5E0/tol
    bb = I1MACH(9)*0.5E0
    aa = MIN(aa,bb)
    aa = aa**tth
    IF ( az>aa ) THEN
      Ierr = 4
      nz = 0
      GOTO 99999
    ELSE
      aa = SQRT(aa)
      IF ( az>aa ) Ierr = 3
      csq = CSQRT(Z)
      zta = Z*csq*CMPLX(tth,0.0E0)
      !-----------------------------------------------------------------------
      !     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
      !-----------------------------------------------------------------------
      sfac = 1.0E0
      zi = AIMAG(Z)
      zr = REAL(Z)
      ak = AIMAG(zta)
      IF ( zr<0.0E0 ) THEN
        bk = REAL(zta)
        ck = -ABS(bk)
        zta = CMPLX(ck,ak)
      ENDIF
      IF ( zi==0.0E0.AND.zr<=0.0E0 ) zta = CMPLX(0.0E0,ak)
      aa = REAL(zta)
      IF ( Kode/=2 ) THEN
        !-----------------------------------------------------------------------
        !     OVERFLOW TEST
        !-----------------------------------------------------------------------
        bb = ABS(aa)
        IF ( bb>=alim ) THEN
          bb = bb + 0.25E0*ALOG(az)
          sfac = tol
          IF ( bb>elim ) GOTO 50
        ENDIF
      ENDIF
      fmr = 0.0E0
      IF ( aa<0.0E0.OR.zr<=0.0E0 ) THEN
        fmr = pi
        IF ( zi<0.0E0 ) fmr = -pi
        zta = -zta
      ENDIF
      !-----------------------------------------------------------------------
      !     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
      !     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBINU
      !-----------------------------------------------------------------------
      CALL CBINU(zta,fnu,Kode,1,cy,nz,rl,fnul,tol,elim,alim)
      IF ( nz>=0 ) THEN
        aa = fmr*fnu
        z3 = CMPLX(sfac,0.0E0)
        s1 = cy(1)*CMPLX(COS(aa),SIN(aa))*z3
        fnu = (2.0E0-fid)/3.0E0
        CALL CBINU(zta,fnu,Kode,2,cy,nz,rl,fnul,tol,elim,alim)
        cy(1) = cy(1)*z3
        cy(2) = cy(2)*z3
        !-----------------------------------------------------------------------
        !     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
        !-----------------------------------------------------------------------
        s2 = cy(1)*CMPLX(fnu+fnu,0.0E0)/zta + cy(2)
        aa = fmr*(fnu-1.0E0)
        s1 = (s1+s2*CMPLX(COS(aa),SIN(aa)))*CMPLX(coef,0.0E0)
        IF ( Id==1 ) THEN
          s1 = Z*s1
          Bi = s1*CMPLX(1.0E0/sfac,0.0E0)
          RETURN
        ELSE
          s1 = csq*s1
          Bi = s1*CMPLX(1.0E0/sfac,0.0E0)
          RETURN
        ENDIF
      ELSEIF ( nz/=(-1) ) THEN
        GOTO 100
      ENDIF
    ENDIF
    50     nz = 0
    Ierr = 2
    RETURN
  ELSE
    !-----------------------------------------------------------------------
    !     POWER SERIES FOR ABS(Z).LE.1.
    !-----------------------------------------------------------------------
    s1 = cone
    s2 = cone
    IF ( az<tol ) THEN
      aa = c1*(1.0E0-fid) + fid*c2
      Bi = CMPLX(aa,0.0E0)
      RETURN
    ELSE
      aa = az*az
      IF ( aa>=tol/az ) THEN
        trm1 = cone
        trm2 = cone
        atrm = 1.0E0
        z3 = Z*Z*Z
        az3 = az*aa
        ak = 2.0E0 + fid
        bk = 3.0E0 - fid - fid
        ck = 4.0E0 - fid
        dk = 3.0E0 + fid + fid
        d1 = ak*dk
        d2 = bk*ck
        ad = MIN(d1,d2)
        ak = 24.0E0 + 9.0E0*fid
        bk = 30.0E0 - 9.0E0*fid
        z3r = REAL(z3)
        z3i = AIMAG(z3)
        DO k = 1, 25
          trm1 = trm1*CMPLX(z3r/d1,z3i/d1)
          s1 = s1 + trm1
          trm2 = trm2*CMPLX(z3r/d2,z3i/d2)
          s2 = s2 + trm2
          atrm = atrm*az3/ad
          d1 = d1 + ak
          d2 = d2 + bk
          ad = MIN(d1,d2)
          IF ( atrm<tol*ad ) EXIT
          ak = ak + 18.0E0
          bk = bk + 18.0E0
        ENDDO
      ENDIF
      IF ( Id==1 ) THEN
        Bi = s2*CMPLX(c2,0.0E0)
        IF ( az>tol ) Bi = Bi + Z*Z*s1*CMPLX(c1/(1.0E0+fid),0.0E0)
        IF ( Kode==1 ) RETURN
        zta = Z*CSQRT(Z)*CMPLX(tth,0.0E0)
        aa = REAL(zta)
        aa = -ABS(aa)
        Bi = Bi*CMPLX(EXP(aa),0.0E0)
        RETURN
      ELSE
        Bi = s1*CMPLX(c1,0.0E0) + Z*s2*CMPLX(c2,0.0E0)
        IF ( Kode==1 ) RETURN
        zta = Z*CSQRT(Z)*CMPLX(tth,0.0E0)
        aa = REAL(zta)
        aa = -ABS(aa)
        Bi = Bi*CMPLX(EXP(aa),0.0E0)
        RETURN
      ENDIF
    ENDIF
  ENDIF
  100  nz = 0
  Ierr = 5
  RETURN
  99999 CONTINUE
  END SUBROUTINE CBIRY
