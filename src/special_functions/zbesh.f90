!** ZBESH
SUBROUTINE ZBESH(Zr,Zi,Fnu,Kode,M,N,Cyr,Cyi,Nz,Ierr)
  !>
  !***
  !  Compute a sequence of the Hankel functions H(m,a,z)
  !            for superscript m=1 or 2, real nonnegative orders a=b,
  !            b+1,... where b>0, and nonzero complex argument z.  A
  !            scaling option is available to help avoid overflow.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10A4
  !***
  ! **Type:**      COMPLEX (CBESH-C, ZBESH-C)
  !***
  ! **Keywords:**  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
  !             BESSEL FUNCTIONS OF THE THIRD KIND, H BESSEL FUNCTIONS,
  !             HANKEL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !                      ***A DOUBLE PRECISION ROUTINE***
  !         On KODE=1, ZBESH computes an N member sequence of complex
  !         Hankel (Bessel) functions CY(L)=H(M,FNU+L-1,Z) for super-
  !         script M=1 or 2, real nonnegative orders FNU+L-1, L=1,...,
  !         N, and complex nonzero Z in the cut plane -pi<arg(Z)<=pi
  !         where Z=ZR+i*ZI.  On KODE=2, CBESH returns the scaled
  !         functions
  !
  !            CY(L) = H(M,FNU+L-1,Z)*exp(-(3-2*M)*Z*i),  i**2=-1
  !
  !         which removes the exponential behavior in both the upper
  !         and lower half planes.  Definitions and notation are found
  !         in the NBS Handbook of Mathematical Functions (Ref. 1).
  !
  !         Input
  !           ZR     - DOUBLE PRECISION real part of nonzero argument Z
  !           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
  !           FNU    - DOUBLE PRECISION initial order, FNU>=0
  !           KODE   - A parameter to indicate the scaling option
  !                    KODE=1  returns
  !                            CY(L)=H(M,FNU+L-1,Z), L=1,...,N
  !                        =2  returns
  !                            CY(L)=H(M,FNU+L-1,Z)*exp(-(3-2M)*Z*i),
  !                            L=1,...,N
  !           M      - Superscript of Hankel function, M=1 or 2
  !           N      - Number of terms in the sequence, N>=1
  !
  !         Output
  !           CYR    - DOUBLE PRECISION real part of result vector
  !           CYI    - DOUBLE PRECISION imag part of result vector
  !           NZ     - Number of underflows set to zero
  !                    NZ=0    Normal return
  !                    NZ>0    CY(L)=0 for NZ values of L (if M=1 and
  !                            Im(Z)>0 or if M=2 and Im(Z)<0, then
  !                            CY(L)=0 for L=1,...,NZ; in the com-
  !                            plementary half planes, the underflows
  !                            may not be in an uninterrupted sequence)
  !           IERR   - Error flag
  !                    IERR=0  Normal return     - COMPUTATION COMPLETED
  !                    IERR=1  Input error       - NO COMPUTATION
  !                    IERR=2  Overflow          - NO COMPUTATION
  !                            (abs(Z) too small and/or FNU+N-1
  !                            too large)
  !                    IERR=3  Precision warning - COMPUTATION COMPLETED
  !                            (Result has half precision or less
  !                            because abs(Z) or FNU+N-1 is large)
  !                    IERR=4  Precision error   - NO COMPUTATION
  !                            (Result has no precision because
  !                            abs(Z) or FNU+N-1 is too large)
  !                    IERR=5  Algorithmic error - NO COMPUTATION
  !                            (Termination condition not met)
  !
  !- Long Description:
  !
  !         The computation is carried out by the formula
  !
  !            H(m,a,z) = (1/t)*exp(-a*t)*K(a,z*exp(-t))
  !                   t = (3-2*m)*i*pi/2
  !
  !         where the K Bessel function is computed as described in the
  !         prologue to CBESK.
  !
  !         Exponential decay of H(m,a,z) occurs in the upper half z
  !         plane for m=1 and the lower half z plane for m=2.  Exponential
  !         growth occurs in the complementary half planes.  Scaling
  !         by exp(-(3-2*m)*z*i) removes the exponential behavior in the
  !         whole z plane as z goes to infinity.
  !
  !         For negative orders, the formula
  !
  !            H(m,-a,z) = H(m,a,z)*exp((3-2*m)*a*pi*i)
  !
  !         can be used.
  !
  !         In most complex variable computation, one must evaluate ele-
  !         mentary functions.  When the magnitude of Z or FNU+N-1 is
  !         large, losses of significance by argument reduction occur.
  !         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
  !         losses exceeding half precision are likely and an error flag
  !         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
  !         precision unit roundoff limited to 18 digits precision.  Also,
  !         if either is larger than U2=0.5/UR, then all significance is
  !         lost and IERR=4.  In order to use the INT function, arguments
  !         must be further restricted not to exceed the largest machine
  !         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
  !         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
  !         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
  !         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
  !         makes U2 limiting in single precision and U3 limiting in
  !         double precision.  This means that one can expect to retain,
  !         in the worst cases on IEEE machines, no digits in single pre-
  !         cision and only 6 digits in double precision.  Similar con-
  !         siderations hold for other machines.
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
  !         magnitude of the larger component.  In these extreme cases,
  !         the principal phase angle is on the order of +P, -P, PI/2-P,
  !         or -PI/2+P.
  !
  !***
  ! **References:**  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
  !                 matical Functions, National Bureau of Standards
  !                 Applied Mathematics Series 55, U. S. Department
  !                 of Commerce, Tenth Printing (1972) or later.
  !               2. D. E. Amos, Computation of Bessel Functions of
  !                 Complex Argument, Report SAND83-0086, Sandia National
  !                 Laboratories, Albuquerque, NM, May 1983.
  !               3. D. E. Amos, Computation of Bessel Functions of
  !                 Complex Argument and Large Order, Report SAND83-0643,
  !                 Sandia National Laboratories, Albuquerque, NM, May
  !                 1983.
  !               4. D. E. Amos, A Subroutine Package for Bessel Functions
  !                 of a Complex Argument and Nonnegative Order, Report
  !                 SAND85-1018, Sandia National Laboratory, Albuquerque,
  !                 NM, May 1985.
  !               5. D. E. Amos, A portable package for Bessel functions
  !                 of a complex argument and nonnegative order, ACM
  !                 Transactions on Mathematical Software, 12 (September
  !                 1986), pp. 265-273.
  !
  !***
  ! **Routines called:**  D1MACH, I1MACH, ZABS, ZACON, ZBKNU, ZBUNK, ZUOIK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890801  REVISION DATE from Version 3.2
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   920128  Category corrected.  (WRB)
  !   920811  Prologue revised.  (DWL)

  !
  !     COMPLEX CY,Z,ZN,ZT,CSGN
  INTEGER i, Ierr, inu, inuh, ir, k, Kode, k1, k2, M, mm, mr, N, nn, nuf, nw, Nz
  REAL(8) :: aa, alim, aln, arg, az, Cyi(N), Cyr(N), dig, elim, fmm, fn, Fnu, &
    fnul, rhpi, rl, r1m5, sgn, str, tol, ufl, Zi, zni, znr, Zr, zti, bb, ascle, &
    rtol, atol, sti, csgnr, csgni
  !
  REAL(8), PARAMETER :: hpi = 1.57079632679489662D0
  !
  !* FIRST EXECUTABLE STATEMENT  ZBESH
  Ierr = 0
  Nz = 0
  IF ( Zr==0.0D0.AND.Zi==0.0D0 ) Ierr = 1
  IF ( Fnu<0.0D0 ) Ierr = 1
  IF ( M<1.OR.M>2 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( N<1 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  nn = N
  !-----------------------------------------------------------------------
  !     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
  !     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
  !     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
  !     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
  !     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
  !     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
  !     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
  !     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
  !     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
  !-----------------------------------------------------------------------
  tol = MAX(D1MACH(4),1.0D-18)
  k1 = I1MACH(15)
  k2 = I1MACH(16)
  r1m5 = D1MACH(5)
  k = MIN(ABS(k1),ABS(k2))
  elim = 2.303D0*(k*r1m5-3.0D0)
  k1 = I1MACH(14) - 1
  aa = r1m5*k1
  dig = MIN(aa,18.0D0)
  aa = aa*2.303D0
  alim = elim + MAX(-aa,-41.45D0)
  fnul = 10.0D0 + 6.0D0*(dig-3.0D0)
  rl = 1.2D0*dig + 3.0D0
  fn = Fnu + (nn-1)
  mm = 3 - M - M
  fmm = mm
  znr = fmm*Zi
  zni = -fmm*Zr
  !-----------------------------------------------------------------------
  !     TEST FOR PROPER RANGE
  !-----------------------------------------------------------------------
  az = ZABS(Zr,Zi)
  aa = 0.5D0/tol
  bb = I1MACH(9)*0.5D0
  aa = MIN(aa,bb)
  IF ( az>aa ) GOTO 300
  IF ( fn>aa ) GOTO 300
  aa = SQRT(aa)
  IF ( az>aa ) Ierr = 3
  IF ( fn>aa ) Ierr = 3
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
  !-----------------------------------------------------------------------
  ufl = D1MACH(1)*1.0D+3
  IF ( az>=ufl ) THEN
    IF ( Fnu>fnul ) THEN
      !-----------------------------------------------------------------------
      !     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
      !-----------------------------------------------------------------------
      mr = 0
      IF ( .NOT.((znr>=0.0D0).AND.(znr/=0.0D0.OR.zni>=0.0D0.OR.M/=2)) ) THEN
        mr = -mm
        IF ( znr==0.0D0.AND.zni<0.0D0 ) THEN
          znr = -znr
          zni = -zni
        END IF
      END IF
      CALL ZBUNK(znr,zni,Fnu,Kode,mr,nn,Cyr,Cyi,nw,tol,elim,alim)
      IF ( nw<0 ) GOTO 200
      Nz = Nz + nw
    ELSE
      IF ( fn>1.0D0 ) THEN
        IF ( fn>2.0D0 ) THEN
          CALL ZUOIK(znr,zni,Fnu,Kode,2,nn,Cyr,Cyi,nuf,tol,elim,alim)
          IF ( nuf<0 ) GOTO 100
          Nz = Nz + nuf
          nn = nn - nuf
          !-----------------------------------------------------------------------
          !     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
          !     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
          !-----------------------------------------------------------------------
          IF ( nn==0 ) THEN
            IF ( znr<0.0D0 ) GOTO 100
            RETURN
          END IF
        ELSEIF ( az<=tol ) THEN
          arg = 0.5D0*az
          aln = -fn*LOG(arg)
          IF ( aln>elim ) GOTO 100
        END IF
      END IF
      IF ( (znr<0.0D0).OR.(znr==0.0D0.AND.zni<0.0D0.AND.M==2) ) THEN
        !-----------------------------------------------------------------------
        !     LEFT HALF PLANE COMPUTATION
        !-----------------------------------------------------------------------
        mr = -mm
        CALL ZACON(znr,zni,Fnu,Kode,mr,nn,Cyr,Cyi,nw,rl,fnul,tol,elim,alim)
        IF ( nw<0 ) GOTO 200
        Nz = nw
      ELSE
        !-----------------------------------------------------------------------
        !     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
        !     YN.GE.0. .OR. M=1)
        !-----------------------------------------------------------------------
        CALL ZBKNU(znr,zni,Fnu,Kode,nn,Cyr,Cyi,Nz,tol,elim,alim)
      END IF
    END IF
    !-----------------------------------------------------------------------
    !     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
    !
    !     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
    !-----------------------------------------------------------------------
    sgn = DSIGN(hpi,-fmm)
    !-----------------------------------------------------------------------
    !     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    !     WHEN FNU IS LARGE
    !-----------------------------------------------------------------------
    inu = INT( Fnu )
    inuh = inu/2
    ir = inu - 2*inuh
    arg = (Fnu-(inu-ir))*sgn
    rhpi = 1.0D0/sgn
    !     ZNI = RHPI*COS(ARG)
    !     ZNR = -RHPI*SIN(ARG)
    csgni = rhpi*COS(arg)
    csgnr = -rhpi*SIN(arg)
    IF ( MOD(inuh,2)/=0 ) THEN
      !     ZNR = -ZNR
      !     ZNI = -ZNI
      csgnr = -csgnr
      csgni = -csgni
    END IF
    zti = -fmm
    rtol = 1.0D0/tol
    ascle = ufl*rtol
    DO i = 1, nn
      !       STR = CYR(I)*ZNR - CYI(I)*ZNI
      !       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
      !       CYR(I) = STR
      !       STR = -ZNI*ZTI
      !       ZNI = ZNR*ZTI
      !       ZNR = STR
      aa = Cyr(i)
      bb = Cyi(i)
      atol = 1.0D0
      IF ( MAX(ABS(aa),ABS(bb))<=ascle ) THEN
        aa = aa*rtol
        bb = bb*rtol
        atol = tol
      END IF
      str = aa*csgnr - bb*csgni
      sti = aa*csgni + bb*csgnr
      Cyr(i) = str*atol
      Cyi(i) = sti*atol
      str = -csgni*zti
      csgni = csgnr*zti
      csgnr = str
    END DO
    RETURN
  END IF
  100  Nz = 0
  Ierr = 2
  RETURN
  200 CONTINUE
  IF ( nw==(-1) ) GOTO 100
  Nz = 0
  Ierr = 5
  RETURN
  300  Nz = 0
  Ierr = 4
END SUBROUTINE ZBESH
