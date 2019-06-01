!** SDZRO
SUBROUTINE SDZRO(Ae,F,H,N,Nq,Iroot,Re,T,Yh,Uround,B,C,Fb,Fc,Y)
  !>
  !  SDZRO searches for a zero of a function F(N, T, Y, IROOT)
  !            between the given values B and C until the width of the
  !            interval (B, C) has collapsed to within a tolerance
  !            specified by the stopping criterion,
  !              ABS(B - C) .LE. 2.*(RW*ABS(B) + AE).
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      SINGLE PRECISION (SDZRO-S, DDZRO-D, CDZRO-C)
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Description:**
  !
  !     This is a special purpose version of ZEROIN, modified for use with
  !     the SDRIV package.
  !
  !     Sandia Mathematical Program Library
  !     Mathematical Computing Services Division 5422
  !     Sandia Laboratories
  !     P. O. Box 5800
  !     Albuquerque, New Mexico  87115
  !     Control Data 6600 Version 4.5, 1 November 1971
  !
  !     PARAMETERS
  !        F     - Name of the external function, which returns a
  !                real result.  This name must be in an
  !                EXTERNAL statement in the calling program.
  !        B     - One end of the interval (B, C).  The value returned for
  !                B usually is the better approximation to a zero of F.
  !        C     - The other end of the interval (B, C).
  !        RE    - Relative error used for RW in the stopping criterion.
  !                If the requested RE is less than machine precision,
  !                then RW is set to approximately machine precision.
  !        AE    - Absolute error used in the stopping criterion.  If the
  !                given interval (B, C) contains the origin, then a
  !                nonzero value should be chosen for AE.
  !
  !***
  ! **References:**  L. F. Shampine and H. A. Watts, ZEROIN, a root-solving
  !                 routine, SC-TM-70-631, Sept 1970.
  !               T. J. Dekker, Finding a zero by means of successive
  !                 linear interpolation, Constructive Aspects of the
  !                 Fundamental Theorem of Algebra, edited by B. Dejon
  !                 and P. Henrici, 1969.
  !***
  ! **Routines called:**  SDNTP

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.

  INTERFACE
    REAL FUNCTION F(N,T,Y,Iroot)
      INTEGER :: N, Iroot
      REAL :: T, Y(N)
    END FUNCTION F
  END INTERFACE
  INTEGER :: Iroot, N, Nq
  REAL :: Ae, B, C, Fb, Fc, H, Re, T, Uround
  REAL :: Y(N), Yh(N,Nq+1)
  INTEGER :: ic, kount
  REAL :: a, acbs, acmb, cmb, er, fa, p, q, rw, tol
  !* FIRST EXECUTABLE STATEMENT  SDZRO
  er = 4.E0*Uround
  rw = MAX(Re,er)
  ic = 0
  acbs = ABS(B-C)
  a = C
  fa = Fc
  kount = 0
  !                                                    Perform interchange
  100 CONTINUE
  IF ( ABS(Fc)<ABS(Fb) ) THEN
    a = B
    fa = Fb
    B = C
    Fb = Fc
    C = a
    Fc = fa
  END IF
  cmb = 0.5E0*(C-B)
  acmb = ABS(cmb)
  tol = rw*ABS(B) + Ae
  !                                                Test stopping criterion
  IF ( acmb<=tol ) RETURN
  IF ( kount>50 ) RETURN
  !                                    Calculate new iterate implicitly as
  !                                    B + P/Q, where we arrange P .GE. 0.
  !                         The implicit form is used to prevent overflow.
  p = (B-a)*Fb
  q = fa - Fb
  IF ( p<0.E0 ) THEN
    p = -p
    q = -q
  END IF
  !                          Update A and check for satisfactory reduction
  !                          in the size of our bounding interval.
  a = B
  fa = Fb
  ic = ic + 1
  IF ( ic>=4 ) THEN
    IF ( 8.E0*acmb>=acbs ) THEN
      !                                                                 Bisect
      B = 0.5E0*(C+B)
      GOTO 200
    END IF
    ic = 0
  END IF
  acbs = acmb
  !                                            Test for too small a change
  IF ( p<=ABS(q)*tol ) THEN
    !                                                 Increment by tolerance
    B = B + SIGN(tol,cmb)
    !                                               Root ought to be between
    !                                               B and (C + B)/2.
  ELSEIF ( p<cmb*q ) THEN
    !                                                            Interpolate
    B = B + p/q
  ELSE
    !                                                                 Bisect
    B = 0.5E0*(C+B)
  END IF
  !                                             Have completed computation
  !                                             for new iterate B.
  200  CALL SDNTP(H,0,N,Nq,T,B,Yh,Y)
  Fb = F(N,B,Y,Iroot)
  IF ( N==0 ) RETURN
  IF ( Fb==0.E0 ) RETURN
  kount = kount + 1
  !
  !             Decide whether next step is interpolation or extrapolation
  !
  IF ( SIGN(1.0E0,Fb)==SIGN(1.0E0,Fc) ) THEN
    C = a
    Fc = fa
  END IF
  GOTO 100
END SUBROUTINE SDZRO
