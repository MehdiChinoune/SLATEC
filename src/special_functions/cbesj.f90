!** CBESJ
PURE SUBROUTINE CBESJ(Z,Fnu,Kode,N,Cy,Nz,Ierr)
  !> Compute a sequence of the Bessel functions J(a,z) for complex argument z
  !  and real nonnegative orders a=b,b+1, b+2,... where b>0.
  !  A scaling option is available to help avoid overflow.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10A4
  !***
  ! **Type:**      COMPLEX (CBESJ-C, ZBESJ-C)
  !***
  ! **Keywords:**  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
  !             BESSEL FUNCTIONS OF THE FIRST KIND, J BESSEL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !         On KODE=1, CBESJ computes an N member sequence of complex
  !         Bessel functions CY(L)=J(FNU+L-1,Z) for real nonnegative
  !         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
  !         -pi<arg(Z)<=pi.  On KODE=2, CBESJ returns the scaled functions
  !
  !            CY(L) = exp(-abs(Y))*J(FNU+L-1,Z),  L=1,...,N and Y=Im(Z)
  !
  !         which remove the exponential growth in both the upper and
  !         lower half planes as Z goes to infinity.  Definitions and
  !         notation are found in the NBS Handbook of Mathematical
  !         Functions (Ref. 1).
  !
  !         Input
  !           Z      - Argument of type COMPLEX
  !           FNU    - Initial order of type REAL(SP), FNU>=0
  !           KODE   - A parameter to indicate the scaling option
  !                    KODE=1  returns
  !                            CY(L)=J(FNU+L-1,Z), L=1,...,N
  !                        =2  returns
  !                            CY(L)=J(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
  !                            where Y=Im(Z)
  !           N      - Number of terms in the sequence, N>=1
  !
  !         Output
  !           CY     - Result vector of type COMPLEX
  !           NZ     - Number of underflows set to zero
  !                    NZ=0    Normal return
  !                    NZ>0    CY(L)=0, L=N-NZ+1,...,N
  !           IERR   - Error flag
  !                    IERR=0  Normal return     - COMPUTATION COMPLETED
  !                    IERR=1  Input error       - NO COMPUTATION
  !                    IERR=2  Overflow          - NO COMPUTATION
  !                            (Im(Z) too large on KODE=1)
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
  !         The computation is carried out by the formulae
  !
  !            J(a,z) = exp( a*pi*i/2)*I(a,-i*z),  Im(z)>=0
  !
  !            J(a,z) = exp(-a*pi*i/2)*I(a, i*z),  Im(z)<0
  !
  !         where the I Bessel function is computed as described in the
  !         prologue to CBESI.
  !
  !         For negative orders, the formula
  !
  !            J(-a,z) = J(a,z)*cos(a*pi) - Y(a,z)*sin(a*pi)
  !
  !         can be used.  However, for large orders close to integers, the
  !         the function changes radically.  When a is a large positive
  !         integer, the magnitude of J(-a,z)=J(a,z)*cos(a*pi) is a
  !         large negative power of ten.  But when a is not an integer,
  !         Y(a,z) dominates in magnitude with a large positive power of
  !         ten and the most that the second term can be reduced is by
  !         unit roundoff from the coefficient.  Thus, wide changes can
  !         occur within unit roundoff of a large integer for a.  Here,
  !         large means a>abs(z).
  !
  !         In most complex variable computation, one must evaluate ele-
  !         mentary functions.  When the magnitude of Z or FNU+N-1 is
  !         large, losses of significance by argument reduction occur.
  !         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
  !         losses exceeding half precision are likely and an error flag
  !         IERR=3 is triggered where UR=eps_sp=UNIT ROUNDOFF.  Also,
  !         if either is larger than U2=0.5/UR, then all significance is
  !         lost and IERR=4.  In order to use the INT function, arguments
  !         must be further restricted not to exceed the largest machine
  !         integer, U3=huge_int.  Thus, the magnitude of Z and FNU+N-1
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
  ! **Routines called:**  CBINU, I1MACH, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890801  REVISION DATE from Version 3.2
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   920128  Category corrected.  (WRB)
  !   920811  Prologue revised.  (DWL)
  USE service, ONLY : digits_sp, huge_int, max_exp_sp, min_exp_sp, eps_sp, &
    log10_radix_sp, tiny_sp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Ierr, Nz
  REAL(SP), INTENT(IN) :: Fnu
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Cy(N)
  !
  INTEGER :: i, inu, inuh, ir, k1, k2, nl, k
  COMPLEX(SP) :: ci, csgn, zn
  REAL(SP) :: aa, alim, arg, dig, elim, fnul, rl, r1, r1m5, &
    r2, tol, yy, az, fn, bb, ascle, rtol, atol
  REAL(SP), PARAMETER :: hpi = 1.57079632679489662_SP
  !
  !* FIRST EXECUTABLE STATEMENT  CBESJ
  Ierr = 0
  Nz = 0
  IF( Fnu<0._SP .OR. Kode<1 .OR. Kode>2 .OR. N<1 ) THEN
    Ierr = 1
    RETURN
  END IF
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
  tol = MAX(eps_sp,1.0E-18_SP)
  k1 = min_exp_sp
  k2 = max_exp_sp
  r1m5 = log10_radix_sp
  k = MIN(ABS(k1),ABS(k2))
  elim = 2.303_SP*(k*r1m5-3._SP)
  k1 = digits_sp - 1
  aa = r1m5*k1
  dig = MIN(aa,18._SP)
  aa = aa*2.303_SP
  alim = elim + MAX(-aa,-41.45E0_SP)
  rl = 1.2_SP*dig + 3._SP
  fnul = 10._SP + 6._SP*(dig-3._SP)
  ci = CMPLX(0._SP,1._SP,SP)
  yy = AIMAG(Z)
  az = ABS(Z)
  !-----------------------------------------------------------------------
  !     TEST FOR RANGE
  !-----------------------------------------------------------------------
  aa = 0.5_SP/tol
  bb = huge_int*0.5_SP
  aa = MIN(aa,bb)
  fn = Fnu + (N-1)
  IF( az<=aa ) THEN
    IF( fn<=aa ) THEN
      aa = SQRT(aa)
      IF( az>aa ) Ierr = 3
      IF( fn>aa ) Ierr = 3
      !-----------------------------------------------------------------------
      !     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
      !     WHEN FNU IS LARGE
      !-----------------------------------------------------------------------
      inu = INT( Fnu )
      inuh = inu/2
      ir = inu - 2*inuh
      arg = (Fnu-(inu-ir))*hpi
      r1 = COS(arg)
      r2 = SIN(arg)
      csgn = CMPLX(r1,r2,SP)
      IF( MOD(inuh,2)==1 ) csgn = -csgn
      !-----------------------------------------------------------------------
      !     ZN IS IN THE RIGHT HALF PLANE
      !-----------------------------------------------------------------------
      zn = -Z*ci
      IF( yy<0._SP ) THEN
        zn = -zn
        csgn = CONJG(csgn)
        ci = CONJG(ci)
      END IF
      CALL CBINU(zn,Fnu,Kode,N,Cy,Nz,rl,fnul,tol,elim,alim)
      IF( Nz>=0 ) THEN
        nl = N - Nz
        IF( nl==0 ) RETURN
        rtol = 1._SP/tol
        ascle = tiny_sp*rtol*1.E+3_SP
        DO i = 1, nl
          !       CY(I)=CY(I)*CSGN
          zn = Cy(i)
          aa = REAL(zn,SP)
          bb = AIMAG(zn)
          atol = 1._SP
          IF( MAX(ABS(aa),ABS(bb))<=ascle ) THEN
            zn = zn*CMPLX(rtol,0._SP,SP)
            atol = tol
          END IF
          zn = zn*csgn
          Cy(i) = zn*CMPLX(atol,0._SP,SP)
          csgn = csgn*ci
        END DO
        RETURN
      ELSEIF( Nz==(-2) ) THEN
        Nz = 0
        Ierr = 5
        RETURN
      ELSE
        Nz = 0
        Ierr = 2
        RETURN
      END IF
    END IF
  END IF
  Nz = 0
  Ierr = 4
  !
END SUBROUTINE CBESJ