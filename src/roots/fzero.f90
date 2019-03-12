!DECK FZERO
SUBROUTINE FZERO(F,B,C,R,Re,Ae,Iflag)
  IMPLICIT NONE
  REAL F, R1MACH
  !***BEGIN PROLOGUE  FZERO
  !***PURPOSE  Search for a zero of a function F(X) in a given interval
  !            (B,C).  It is designed primarily for problems where F(B)
  !            and F(C) have opposite signs.
  !***LIBRARY   SLATEC
  !***CATEGORY  F1B
  !***TYPE      SINGLE PRECISION (FZERO-S, DFZERO-D)
  !***KEYWORDS  BISECTION, NONLINEAR EQUATIONS, ROOTS, ZEROS
  !***AUTHOR  Shampine, L. F., (SNLA)
  !           Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     FZERO searches for a zero of a REAL function F(X) between the
  !     given REAL values B and C until the width of the interval (B,C)
  !     has collapsed to within a tolerance specified by the stopping
  !     criterion,
  !        ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
  !     The method used is an efficient combination of bisection and the
  !     secant rule and is due to T. J. Dekker.
  !
  !     Description Of Arguments
  !
  !   F     :EXT   - Name of the REAL external function.  This name must
  !                  be in an EXTERNAL statement in the calling program.
  !                  F must be a function of one REAL argument.
  !
  !   B     :INOUT - One end of the REAL interval (B,C).  The value
  !                  returned for B usually is the better approximation
  !                  to a zero of F.
  !
  !   C     :INOUT - The other end of the REAL interval (B,C)
  !
  !   R     :OUT   - A (better) REAL guess of a zero of F which could help
  !                  in speeding up convergence.  If F(B) and F(R) have
  !                  opposite signs, a root will be found in the interval
  !                  (B,R); if not, but F(R) and F(C) have opposite signs,
  !                  a root will be found in the interval (R,C);
  !                  otherwise, the interval (B,C) will be searched for a
  !                  possible root.  When no better guess is known, it is
  !                  recommended that r be set to B or C, since if R is
  !                  not interior to the interval (B,C), it will be
  !                  ignored.
  !
  !   RE    :IN    - Relative error used for RW in the stopping criterion.
  !                  If the requested RE is less than machine precision,
  !                  then RW is set to approximately machine precision.
  !
  !   AE    :IN    - Absolute error used in the stopping criterion.  If
  !                  the given interval (B,C) contains the origin, then a
  !                  nonzero value should be chosen for AE.
  !
  !   IFLAG :OUT   - A status code.  User must check IFLAG after each
  !                  call.  Control returns to the user from FZERO in all
  !                  cases.
  !
  !                1  B is within the requested tolerance of a zero.
  !                   The interval (B,C) collapsed to the requested
  !                   tolerance, the function changes sign in (B,C), and
  !                   F(X) decreased in magnitude as (B,C) collapsed.
  !
  !                2  F(B) = 0.  However, the interval (B,C) may not have
  !                   collapsed to the requested tolerance.
  !
  !                3  B may be near a singular point of F(X).
  !                   The interval (B,C) collapsed to the requested tol-
  !                   erance and the function changes sign in (B,C), but
  !                   F(X) increased in magnitude as (B,C) collapsed, i.e.
  !                     ABS(F(B out)) .GT. MAX(ABS(F(B in)),ABS(F(C in)))
  !
  !                4  No change in sign of F(X) was found although the
  !                   interval (B,C) collapsed to the requested tolerance.
  !                   The user must examine this case and decide whether
  !                   B is near a local minimum of F(X), or B is near a
  !                   zero of even multiplicity, or neither of these.
  !
  !                5  Too many (.GT. 500) function evaluations used.
  !
  !***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
  !                 code, Report SC-TM-70-631, Sandia Laboratories,
  !                 September 1970.
  !               T. J. Dekker, Finding a zero by means of successive
  !                 linear interpolation, Constructive Aspects of the
  !                 Fundamental Theorem of Algebra, edited by B. Dejon
  !                 and P. Henrici, Wiley-Interscience, 1969.
  !***ROUTINES CALLED  R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   700901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  FZERO
  REAL a, acbs, acmb, Ae, aw, B, C, cmb, er, fa, fb, fc, fx, &
    fz, p, q, R, Re, rw, t, tol, z
  INTEGER ic, Iflag, kount
  !***FIRST EXECUTABLE STATEMENT  FZERO
  !
  !   ER is two times the computer unit roundoff value which is defined
  !   here by the function R1MACH.
  !
  er = 2.0E0*R1MACH(4)
  !
  !   Initialize.
  !
  z = R
  IF ( R<=MIN(B,C).OR.R>=MAX(B,C) ) z = C
  rw = MAX(Re,er)
  aw = MAX(Ae,0.E0)
  ic = 0
  t = z
  fz = F(t)
  fc = fz
  t = B
  fb = F(t)
  kount = 2
  IF ( SIGN(1.0E0,fz)/=SIGN(1.0E0,fb) ) THEN
    C = z
  ELSEIF ( z/=C ) THEN
    t = C
    fc = F(t)
    kount = 3
    IF ( SIGN(1.0E0,fz)/=SIGN(1.0E0,fc) ) THEN
      B = z
      fb = fz
    ENDIF
  ENDIF
  a = C
  fa = fc
  acbs = ABS(B-C)
  fx = MAX(ABS(fb),ABS(fc))
  DO
    !
    IF ( ABS(fc)<ABS(fb) ) THEN
      !
      !   Perform interchange.
      !
      a = B
      fa = fb
      B = C
      fb = fc
      C = a
      fc = fa
    ENDIF
    !
    cmb = 0.5E0*(C-B)
    acmb = ABS(cmb)
    tol = rw*ABS(B) + aw
    !
    !   Test stopping criterion and function count.
    !
    IF ( acmb<=tol ) THEN
      !
      !   Finished.  Process results for proper setting of IFLAG.
      !
      IF ( SIGN(1.0E0,fb)==SIGN(1.0E0,fc) ) THEN
        Iflag = 4
        RETURN
      ELSEIF ( ABS(fb)>fx ) THEN
        Iflag = 3
        RETURN
      ELSE
        Iflag = 1
        RETURN
      ENDIF
    ELSEIF ( fb==0.E0 ) THEN
      Iflag = 2
      RETURN
    ELSE
      IF ( kount>=500 ) THEN
        Iflag = 5
        EXIT
      ELSE
        !
        !   Calculate new iterate implicitly as B+P/Q, where we arrange
        !   P .GE. 0.  The implicit form is used to prevent overflow.
        !
        p = (B-a)*fb
        q = fa - fb
        IF ( p<0.E0 ) THEN
          p = -p
          q = -q
        ENDIF
        !
        !   Update A and check for satisfactory reduction in the size of the
        !   bracketing interval.  If not, perform bisection.
        !
        a = B
        fa = fb
        ic = ic + 1
        IF ( ic>=4 ) THEN
          IF ( 8.0E0*acmb>=acbs ) THEN
            !
            !   Use bisection (C+B)/2.
            !
            B = B + cmb
            GOTO 20
          ELSE
            ic = 0
            acbs = acmb
          ENDIF
        ENDIF
        !
        !   Test for too small a change.
        !
        IF ( p<=ABS(q)*tol ) THEN
          !
          !   Increment by TOLerance.
          !
          B = B + SIGN(tol,cmb)
          !
          !   Root ought to be between B and (C+B)/2.
          !
        ELSEIF ( p>=cmb*q ) THEN
          B = B + cmb
        ELSE
          !
          !   Use secant rule.
          !
          B = B + p/q
        ENDIF
      ENDIF
      !
      !   Have completed computation for new iterate B.
      !
      20       t = B
      fb = F(t)
      kount = kount + 1
      !
      !   Decide whether next step is interpolation or extrapolation.
      !
      IF ( SIGN(1.0E0,fb)==SIGN(1.0E0,fc) ) THEN
        C = a
        fc = fa
      ENDIF
    ENDIF
  ENDDO
END SUBROUTINE FZERO
