!** DFZERO
PURE SUBROUTINE DFZERO(F,B,C,R,Re,Ae,Iflag)
  !> Search for a zero of a function F(X) in a given interval (B,C).
  !  It is designed primarily for problems where F(B) and F(C) have opposite signs.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F1B
  !***
  ! **Type:**      DOUBLE PRECISION (FZERO-S, DFZERO-D)
  !***
  ! **Keywords:**  BISECTION, NONLINEAR, ROOTS, ZEROS
  !***
  ! **Author:**  Shampine, L. F., (SNLA)
  !           Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !     DFZERO searches for a zero of a DOUBLE PRECISION function F(X)
  !     between the given DOUBLE PRECISION values B and C until the width
  !     of the interval (B,C) has collapsed to within a tolerance
  !     specified by the stopping criterion,
  !        ABS(B-C) <= 2.*(RW*ABS(B)+AE).
  !     The method used is an efficient combination of bisection and the
  !     secant rule and is due to T. J. Dekker.
  !
  !     Description Of Arguments
  !
  !   F     :EXT   - Name of the DOUBLE PRECISION external function.  This
  !                  name must be in an EXTERNAL statement in the calling
  !                  program.  F must be a function of one DOUBLE
  !                  PRECISION argument.
  !
  !   B     :INOUT - One end of the DOUBLE PRECISION interval (B,C).  The
  !                  value returned for B usually is the better
  !                  approximation to a zero of F.
  !
  !   C     :INOUT - The other end of the DOUBLE PRECISION interval (B,C)
  !
  !   R     :IN    - A (better) DOUBLE PRECISION guess of a zero of F
  !                  which could help in speeding up convergence.  If F(B)
  !                  and F(R) have opposite signs, a root will be found in
  !                  the interval (B,R);  if not, but F(R) and F(C) have
  !                  opposite signs, a root will be found in the interval
  !                  (R,C);  otherwise, the interval (B,C) will be
  !                  searched for a possible root.  When no better guess
  !                  is known, it is recommended that R be set to B or C,
  !                  since if R is not interior to the interval (B,C), it
  !                  will be ignored.
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
  !                  call.  Control returns to the user from DFZERO in all
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
  !                     ABS(F(B out)) > MAX(ABS(F(B in)),ABS(F(C in)))
  !
  !                4  No change in sign of F(X) was found although the
  !                   interval (B,C) collapsed to the requested tolerance.
  !                   The user must examine this case and decide whether
  !                   B is near a local minimum of F(X), or B is near a
  !                   zero of even multiplicity, or neither of these.
  !
  !                5  Too many (> 500) function evaluations used.
  !
  !***
  ! **References:**  L. F. Shampine and H. A. Watts, FZERO, a root-solving
  !                 code, Report SC-TM-70-631, Sandia Laboratories,
  !                 September 1970.
  !               T. J. Dekker, Finding a zero by means of successive
  !                 linear interpolation, Constructive Aspects of the
  !                 Fundamental Theorem of Algebra, edited by B. Dejon
  !                 and P. Henrici, Wiley-Interscience, 1969.
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   700901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : D1MACH

  INTERFACE
    REAL(DP) PURE FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(OUT) :: Iflag
  REAL(DP), INTENT(IN) :: Ae, Re
  REAL(DP), INTENT(INOUT) :: B, C, R
  INTEGER :: ic, kount
  REAL(DP) :: a, acbs, acmb, aw, cmb, er, fa, fb, fc, fx, fz, p, q, rw, t, tol, z
  !
  !* FIRST EXECUTABLE STATEMENT  DFZERO
  !
  !   ER is two times the computer unit roundoff value which is defined
  !   here by the function D1MACH.
  !
  er = 2._DP*D1MACH(4)
  !
  !   Initialize.
  !
  z = R
  IF( R<=MIN(B,C) .OR. R>=MAX(B,C) ) z = C
  rw = MAX(Re,er)
  aw = MAX(Ae,0._DP)
  ic = 0
  t = z
  fz = F(t)
  fc = fz
  t = B
  fb = F(t)
  kount = 2
  IF( SIGN(1._DP,fz)/=SIGN(1._DP,fb) ) THEN
    C = z
  ELSEIF( z/=C ) THEN
    t = C
    fc = F(t)
    kount = 3
    IF( SIGN(1._DP,fz)/=SIGN(1._DP,fc) ) THEN
      B = z
      fb = fz
    END IF
  END IF
  a = C
  fa = fc
  acbs = ABS(B-C)
  fx = MAX(ABS(fb),ABS(fc))
  DO
    !
    IF( ABS(fc)<ABS(fb) ) THEN
      !
      !   Perform interchange.
      !
      a = B
      fa = fb
      B = C
      fb = fc
      C = a
      fc = fa
    END IF
    !
    cmb = 0.5_DP*(C-B)
    acmb = ABS(cmb)
    tol = rw*ABS(B) + aw
    !
    !   Test stopping criterion and function count.
    !
    IF( acmb<=tol ) THEN
      !
      !   Finished.  Process results for proper setting of IFLAG.
      !
      IF( SIGN(1._DP,fb)==SIGN(1._DP,fc) ) THEN
        Iflag = 4
        RETURN
      ELSEIF( ABS(fb)>fx ) THEN
        Iflag = 3
        RETURN
      ELSE
        Iflag = 1
        RETURN
      END IF
    ELSEIF( fb==0._DP ) THEN
      Iflag = 2
      RETURN
    ELSE
      IF( kount>=500 ) THEN
        Iflag = 5
        EXIT
      ELSE
        !
        !   Calculate new iterate implicitly as B+P/Q, where we arrange
        !   P >= 0.  The implicit form is used to prevent overflow.
        !
        p = (B-a)*fb
        q = fa - fb
        IF( p<0._DP ) THEN
          p = -p
          q = -q
        END IF
        !
        !   Update A and check for satisfactory reduction in the size of the
        !   bracketing interval.  If not, perform bisection.
        !
        a = B
        fa = fb
        ic = ic + 1
        IF( ic>=4 ) THEN
          IF( 8._DP*acmb>=acbs ) THEN
            !
            !   Use bisection (C+B)/2.
            !
            B = B + cmb
            GOTO 20
          ELSE
            ic = 0
            acbs = acmb
          END IF
        END IF
        !
        !   Test for too small a change.
        !
        IF( p<=ABS(q)*tol ) THEN
          !
          !   Increment by TOLerance.
          !
          B = B + SIGN(tol,cmb)
          !
          !   Root ought to be between B and (C+B)/2.
          !
        ELSEIF( p>=cmb*q ) THEN
          B = B + cmb
        ELSE
          !
          !   Use secant rule.
          !
          B = B + p/q
        END IF
      END IF
      !
      !   Have completed computation for new iterate B.
      !
      20  t = B
      fb = F(t)
      kount = kount + 1
      !
      !   Decide whether next step is interpolation or extrapolation.
      !
      IF( SIGN(1._DP,fb)==SIGN(1._DP,fc) ) THEN
        C = a
        fc = fa
      END IF
    END IF
  END DO

END SUBROUTINE DFZERO