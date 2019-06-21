!** ZBIRY
SUBROUTINE ZBIRY(Zr,Zi,Id,Kode,Bir,Bii,Ierr)
  !> Compute the Airy function Bi(z) or its derivative dBi/dz
  !            for complex argument z.  A scaling option is available
  !            to help avoid overflow.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      COMPLEX (CBIRY-C, ZBIRY-C)
  !***
  ! **Keywords:**  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
  !             BESSEL FUNCTION OF ORDER TWO THIRDS
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !                      ***A DOUBLE PRECISION ROUTINE***
  !         On KODE=1, ZBIRY computes the complex Airy function Bi(z)
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
  !           ZR     - DOUBLE PRECISION real part of argument Z
  !           ZI     - DOUBLE PRECISION imag part of argument Z
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
  !           BIR    - DOUBLE PRECISION real part of result
  !           BII    - DOUBLE PRECISION imag part of result
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
  !- Long Description:
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
  !         flag IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is
  !         double precision unit roundoff limited to 18 digits precision.
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
  !***
  ! **References:**  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
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
  !***
  ! **Routines called:**  D1MACH, I1MACH, ZABS, ZBINU, ZDIV, ZSQRT

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890801  REVISION DATE from Version 3.2
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   920128  Category corrected.  (WRB)
  !   920811  Prologue revised.  (DWL)
  !   930122  Added ZSQRT to EXTERNAL statement.  (RWC)
  USE service, ONLY : D1MACH, I1MACH
  !     COMPLEX BI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
  REAL(DP) :: aa, ad, ak, alim, atrm, az, az3, bb, Bii, Bir, bk, cc, ck, &
    csqi, csqr, cyi(2), cyr(2), dig, dk, d1, d2, eaa, elim, fid, fmr, &
    fnu, fnul, rl, r1m5, sfac, sti, str, s1i, s1r, s2i, s2r, tol, trm1i, &
    trm1r, trm2i, trm2r, Zi, Zr, ztai, ztar, z3i, z3r
  INTEGER :: Id, Ierr, k, Kode, k1, k2, nz
  REAL(DP), PARAMETER :: tth = 6.66666666666666667E-01_DP, c1 = 6.14926627446000736E-01_DP, &
    c2 = 4.48288357353826359E-01_DP, coef = 5.77350269189625765E-01_DP, &
    pi = 3.14159265358979324_DP
  REAL(DP), PARAMETER :: coner = 1._DP, conei = 0._DP
  !* FIRST EXECUTABLE STATEMENT  ZBIRY
  Ierr = 0
  nz = 0
  IF( Id<0 .OR. Id>1 ) Ierr = 1
  IF( Kode<1 .OR. Kode>2 ) Ierr = 1
  IF( Ierr/=0 ) RETURN
  az = ZABS(Zr,Zi)
  tol = MAX(D1MACH(4),1.E-18_DP)
  fid = Id
  IF( az>1._SP ) THEN
    !-----------------------------------------------------------------------
    !     CASE FOR ABS(Z)>1.0
    !-----------------------------------------------------------------------
    fnu = (1._DP+fid)/3._DP
    !-----------------------------------------------------------------------
    !     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    !     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    !     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    !     EXP(-ELIM)<EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    !     EXP(ELIM)>EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    !     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    !     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    !     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    !     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
    !-----------------------------------------------------------------------
    k1 = I1MACH(15)
    k2 = I1MACH(16)
    r1m5 = D1MACH(5)
    k = MIN(ABS(k1),ABS(k2))
    elim = 2.303_DP*(k*r1m5-3._DP)
    k1 = I1MACH(14) - 1
    aa = r1m5*k1
    dig = MIN(aa,18._DP)
    aa = aa*2.303_DP
    alim = elim + MAX(-aa,-41.45_DP)
    rl = 1.2_DP*dig + 3._DP
    fnul = 10._DP + 6._DP*(dig-3._DP)
    !-----------------------------------------------------------------------
    !     TEST FOR RANGE
    !-----------------------------------------------------------------------
    aa = 0.5_DP/tol
    bb = I1MACH(9)*0.5_DP
    aa = MIN(aa,bb)
    aa = aa**tth
    IF( az>aa ) THEN
      Ierr = 4
      nz = 0
      RETURN
    ELSE
      aa = SQRT(aa)
      IF( az>aa ) Ierr = 3
      CALL ZSQRT(Zr,Zi,csqr,csqi)
      ztar = tth*(Zr*csqr-Zi*csqi)
      ztai = tth*(Zr*csqi+Zi*csqr)
      !-----------------------------------------------------------------------
      !     RE(ZTA)<=0 WHEN RE(Z)<0, ESPECIALLY WHEN IM(Z) IS SMALL
      !-----------------------------------------------------------------------
      sfac = 1._DP
      ak = ztai
      IF( Zr<0._DP ) THEN
        bk = ztar
        ck = -ABS(bk)
        ztar = ck
        ztai = ak
      END IF
      IF( Zi==0._DP .AND. Zr<=0._DP ) THEN
        ztar = 0._DP
        ztai = ak
      END IF
      aa = ztar
      IF( Kode/=2 ) THEN
        !-----------------------------------------------------------------------
        !     OVERFLOW TEST
        !-----------------------------------------------------------------------
        bb = ABS(aa)
        IF( bb>=alim ) THEN
          bb = bb + 0.25_DP*LOG(az)
          sfac = tol
          IF( bb>elim ) GOTO 50
        END IF
      END IF
      fmr = 0._DP
      IF( aa<0._DP .OR. Zr<=0._DP ) THEN
        fmr = pi
        IF( Zi<0._DP ) fmr = -pi
        ztar = -ztar
        ztai = -ztai
      END IF
      !-----------------------------------------------------------------------
      !     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
      !     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBESI
      !-----------------------------------------------------------------------
      CALL ZBINU(ztar,ztai,fnu,Kode,1,cyr,cyi,nz,rl,fnul,tol,elim,alim)
      IF( nz>=0 ) THEN
        aa = fmr*fnu
        z3r = sfac
        str = COS(aa)
        sti = SIN(aa)
        s1r = (str*cyr(1)-sti*cyi(1))*z3r
        s1i = (str*cyi(1)+sti*cyr(1))*z3r
        fnu = (2._DP-fid)/3._DP
        CALL ZBINU(ztar,ztai,fnu,Kode,2,cyr,cyi,nz,rl,fnul,tol,elim,alim)
        cyr(1) = cyr(1)*z3r
        cyi(1) = cyi(1)*z3r
        cyr(2) = cyr(2)*z3r
        cyi(2) = cyi(2)*z3r
        !-----------------------------------------------------------------------
        !     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
        !-----------------------------------------------------------------------
        CALL ZDIV(cyr(1),cyi(1),ztar,ztai,str,sti)
        s2r = (fnu+fnu)*str + cyr(2)
        s2i = (fnu+fnu)*sti + cyi(2)
        aa = fmr*(fnu-1._DP)
        str = COS(aa)
        sti = SIN(aa)
        s1r = coef*(s1r+s2r*str-s2i*sti)
        s1i = coef*(s1i+s2r*sti+s2i*str)
        IF( Id==1 ) THEN
          str = Zr*s1r - Zi*s1i
          s1i = Zr*s1i + Zi*s1r
          s1r = str
          Bir = s1r/sfac
          Bii = s1i/sfac
          RETURN
        ELSE
          str = csqr*s1r - csqi*s1i
          s1i = csqr*s1i + csqi*s1r
          s1r = str
          Bir = s1r/sfac
          Bii = s1i/sfac
          RETURN
        END IF
      ELSEIF( nz/=(-1) ) THEN
        GOTO 100
      END IF
    END IF
    50  Ierr = 2
    nz = 0
    RETURN
  ELSE
    !-----------------------------------------------------------------------
    !     POWER SERIES FOR ABS(Z)<=1.
    !-----------------------------------------------------------------------
    s1r = coner
    s1i = conei
    s2r = coner
    s2i = conei
    IF( az<tol ) THEN
      aa = c1*(1._DP-fid) + fid*c2
      Bir = aa
      Bii = 0._DP
      RETURN
    ELSE
      aa = az*az
      IF( aa>=tol/az ) THEN
        trm1r = coner
        trm1i = conei
        trm2r = coner
        trm2i = conei
        atrm = 1._DP
        str = Zr*Zr - Zi*Zi
        sti = Zr*Zi + Zi*Zr
        z3r = str*Zr - sti*Zi
        z3i = str*Zi + sti*Zr
        az3 = az*aa
        ak = 2._DP + fid
        bk = 3._DP - fid - fid
        ck = 4._DP - fid
        dk = 3._DP + fid + fid
        d1 = ak*dk
        d2 = bk*ck
        ad = MIN(d1,d2)
        ak = 24._DP + 9._DP*fid
        bk = 30._DP - 9._DP*fid
        DO k = 1, 25
          str = (trm1r*z3r-trm1i*z3i)/d1
          trm1i = (trm1r*z3i+trm1i*z3r)/d1
          trm1r = str
          s1r = s1r + trm1r
          s1i = s1i + trm1i
          str = (trm2r*z3r-trm2i*z3i)/d2
          trm2i = (trm2r*z3i+trm2i*z3r)/d2
          trm2r = str
          s2r = s2r + trm2r
          s2i = s2i + trm2i
          atrm = atrm*az3/ad
          d1 = d1 + ak
          d2 = d2 + bk
          ad = MIN(d1,d2)
          IF( atrm<tol*ad ) EXIT
          ak = ak + 18._DP
          bk = bk + 18._DP
        END DO
      END IF
      IF( Id==1 ) THEN
        Bir = s2r*c2
        Bii = s2i*c2
        IF( az>tol ) THEN
          cc = c1/(1._DP+fid)
          str = s1r*Zr - s1i*Zi
          sti = s1r*Zi + s1i*Zr
          Bir = Bir + cc*(str*Zr-sti*Zi)
          Bii = Bii + cc*(str*Zi+sti*Zr)
        END IF
        IF( Kode==1 ) RETURN
        CALL ZSQRT(Zr,Zi,str,sti)
        ztar = tth*(Zr*str-Zi*sti)
        ztai = tth*(Zr*sti+Zi*str)
        aa = ztar
        aa = -ABS(aa)
        eaa = EXP(aa)
        Bir = Bir*eaa
        Bii = Bii*eaa
        RETURN
      ELSE
        Bir = c1*s1r + c2*(Zr*s2r-Zi*s2i)
        Bii = c1*s1i + c2*(Zr*s2i+Zi*s2r)
        IF( Kode==1 ) RETURN
        CALL ZSQRT(Zr,Zi,str,sti)
        ztar = tth*(Zr*str-Zi*sti)
        ztai = tth*(Zr*sti+Zi*str)
        aa = ztar
        aa = -ABS(aa)
        eaa = EXP(aa)
        Bir = Bir*eaa
        Bii = Bii*eaa
        RETURN
      END IF
    END IF
  END IF
  100  nz = 0
  Ierr = 5
  RETURN
END SUBROUTINE ZBIRY
