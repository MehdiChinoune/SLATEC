!*==CBESY.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CBESY
SUBROUTINE CBESY(Z,Fnu,Kode,N,Cy,Nz,Cwrk,Ierr)
  IMPLICIT NONE
  !*--CBESY5
  !***BEGIN PROLOGUE  CBESY
  !***PURPOSE  Compute a sequence of the Bessel functions Y(a,z) for
  !            complex argument z and real nonnegative orders a=b,b+1,
  !            b+2,... where b>0.  A scaling option is available to
  !            help avoid overflow.
  !***LIBRARY   SLATEC
  !***CATEGORY  C10A4
  !***TYPE      COMPLEX (CBESY-C, ZBESY-C)
  !***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
  !             BESSEL FUNCTIONS OF SECOND KIND, WEBER'S FUNCTION,
  !             Y BESSEL FUNCTIONS
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !         On KODE=1, CBESY computes an N member sequence of complex
  !         Bessel functions CY(L)=Y(FNU+L-1,Z) for real nonnegative
  !         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
  !         -pi<arg(Z)<=pi.  On KODE=2, CBESY returns the scaled
  !         functions
  !
  !            CY(L) = exp(-abs(Y))*Y(FNU+L-1,Z),  L=1,...,N, Y=Im(Z)
  !
  !         which remove the exponential growth in both the upper and
  !         lower half planes as Z goes to infinity.  Definitions and
  !         notation are found in the NBS Handbook of Mathematical
  !         Functions (Ref. 1).
  !
  !         Input
  !           Z      - Nonzero argument of type COMPLEX
  !           FNU    - Initial order of type REAL, FNU>=0
  !           KODE   - A parameter to indicate the scaling option
  !                    KODE=1  returns
  !                            CY(L)=Y(FNU+L-1,Z), L=1,...,N
  !                        =2  returns
  !                            CY(L)=Y(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
  !                            where Y=Im(Z)
  !           N      - Number of terms in the sequence, N>=1
  !           CWRK   - A work vector of type COMPLEX and dimension N
  !
  !         Output
  !           CY     - Result vector of type COMPLEX
  !           NZ     - Number of underflows set to zero
  !                    NZ=0    Normal return
  !                    NZ>0    CY(L)=0 for NZ values of L, usually on
  !                            KODE=2 (the underflows may not be in an
  !                            uninterrupted sequence)
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
  ! *Long Description:
  !
  !         The computation is carried out by the formula
  !
  !            Y(a,z) = (H(1,a,z) - H(2,a,z))/(2*i)
  !
  !         where the Hankel functions are computed as described in CBESH.
  !
  !         For negative orders, the formula
  !
  !            Y(-a,z) = Y(a,z)*cos(a*pi) + J(a,z)*sin(a*pi)
  !
  !         can be used.  However, for large orders close to half odd
  !         integers the function changes radically.  When a is a large
  !         positive half odd integer, the magnitude of Y(-a,z)=J(a,z)*
  !         sin(a*pi) is a large negative power of ten.  But when a is
  !         not a half odd integer, Y(a,z) dominates in magnitude with a
  !         large positive power of ten and the most that the second term
  !         can be reduced is by unit roundoff from the coefficient.
  !         Thus,  wide changes can occur within unit roundoff of a large
  !         half odd integer.  Here, large means a>abs(z).
  !
  !         In most complex variable computation, one must evaluate ele-
  !         mentary functions.  When the magnitude of Z or FNU+N-1 is
  !         large, losses of significance by argument reduction occur.
  !         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
  !         losses exceeding half precision are likely and an error flag
  !         IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.  Also,
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
  !***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
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
  !***ROUTINES CALLED  CBESH, I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890801  REVISION DATE from Version 3.2
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   920128  Category corrected.  (WRB)
  !   920811  Prologue revised.  (DWL)
  !***END PROLOGUE  CBESY
  !
  COMPLEX Cwrk, Cy, c1, c2, ex, hci, Z, zu, zv
  REAL elim, ey, Fnu, r1, r2, tay, xx, yy, R1MACH, r1m5, ascle, &
    rtol, atol, tol, aa, bb
  INTEGER i, Ierr, k, Kode, k1, k2, N, Nz, nz1, nz2, I1MACH
  DIMENSION Cy(N), Cwrk(N)
  !***FIRST EXECUTABLE STATEMENT  CBESY
  xx = REAL(Z)
  yy = AIMAG(Z)
  Ierr = 0
  Nz = 0
  IF ( xx==0.0E0.AND.yy==0.0E0 ) Ierr = 1
  IF ( Fnu<0.0E0 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( N<1 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  hci = CMPLX(0.0E0,0.5E0)
  CALL CBESH(Z,Fnu,Kode,1,N,Cy,nz1,Ierr)
  IF ( Ierr/=0.AND.Ierr/=3 ) THEN
    Nz = 0
  ELSE
    CALL CBESH(Z,Fnu,Kode,2,N,Cwrk,nz2,Ierr)
    IF ( Ierr/=0.AND.Ierr/=3 ) THEN
      Nz = 0
    ELSE
      Nz = MIN(nz1,nz2)
      IF ( Kode==2 ) THEN
        tol = MAX(R1MACH(4),1.0E-18)
        k1 = I1MACH(12)
        k2 = I1MACH(13)
        k = MIN(ABS(k1),ABS(k2))
        r1m5 = R1MACH(5)
        !-----------------------------------------------------------------------
        !     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
        !-----------------------------------------------------------------------
        elim = 2.303E0*(k*r1m5-3.0E0)
        r1 = COS(xx)
        r2 = SIN(xx)
        ex = CMPLX(r1,r2)
        ey = 0.0E0
        tay = ABS(yy+yy)
        IF ( tay<elim ) ey = EXP(-tay)
        IF ( yy<0.0E0 ) THEN
          c1 = ex
          c2 = CONJG(ex)*CMPLX(ey,0.0E0)
        ELSE
          c1 = ex*CMPLX(ey,0.0E0)
          c2 = CONJG(ex)
        ENDIF
        Nz = 0
        rtol = 1.0E0/tol
        ascle = R1MACH(1)*rtol*1.0E+3
        DO i = 1, N
          !       CY(I) = HCI*(C2*CWRK(I)-C1*CY(I))
          zv = Cwrk(i)
          aa = REAL(zv)
          bb = AIMAG(zv)
          atol = 1.0E0
          IF ( MAX(ABS(aa),ABS(bb))<=ascle ) THEN
            zv = zv*CMPLX(rtol,0.0E0)
            atol = tol
          ENDIF
          zv = zv*c2*hci
          zv = zv*CMPLX(atol,0.0E0)
          zu = Cy(i)
          aa = REAL(zu)
          bb = AIMAG(zu)
          atol = 1.0E0
          IF ( MAX(ABS(aa),ABS(bb))<=ascle ) THEN
            zu = zu*CMPLX(rtol,0.0E0)
            atol = tol
          ENDIF
          zu = zu*c1*hci
          zu = zu*CMPLX(atol,0.0E0)
          Cy(i) = zv - zu
          IF ( Cy(i)==CMPLX(0.0E0,0.0E0).AND.ey==0.0E0 ) Nz = Nz + 1
        ENDDO
        RETURN
      ELSE
        DO i = 1, N
          Cy(i) = hci*(Cwrk(i)-Cy(i))
        ENDDO
        RETURN
      ENDIF
    ENDIF
  ENDIF
END SUBROUTINE CBESY
