!** SROTMG
SUBROUTINE SROTMG(Sd1,Sd2,Sx1,Sy1,Sparam)
  IMPLICIT NONE
  !>
  !***
  !  Construct a modified Givens transformation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B10
  !***
  ! **Type:**      SINGLE PRECISION (SROTMG-S, DROTMG-D)
  !***
  ! **Keywords:**  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
  !***
  ! **Author:**  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***
  ! **Description:**
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !      SD1  single precision scalar
  !      SD2  single precision scalar
  !      SX1  single precision scalar
  !      SY2  single precision scalar
  !   SPARAM  S.P. 5-vector. SPARAM(1)=SFLAG defined below.
  !           Locations 2-5 contain the rotation matrix.
  !
  !     --Output--
  !      SD1  changed to represent the effect of the transformation
  !      SD2  changed to represent the effect of the transformation
  !      SX1  changed to represent the effect of the transformation
  !      SY2  unchanged
  !
  !     Construct the modified Givens transformation matrix H which zeros
  !     the second component of the 2-vector  (SQRT(SD1)*SX1,SQRT(SD2)*
  !     SY2)**T.
  !     With SPARAM(1)=SFLAG, H has one of the following forms:
  !
  !     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
  !
  !       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
  !     H=(          )    (          )    (          )    (          )
  !       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
  !
  !     Locations 2-5 of SPARAM contain SH11, SH21, SH12, and SH22,
  !     respectively.  (Values of 1.E0, -1.E0, or 0.E0 implied by the
  !     value of SPARAM(1) are not stored in SPARAM.)
  !
  !***
  ! **References:**  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   780301  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920316  Prologue corrected.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL Sd1, Sd2, sflag, sh11, sh12, sh21, sh22, sp1, sp2, Sparam(5), sq1, sq2, &
    stemp, su, Sx1, Sy1
  INTEGER igo
  REAL, PARAMETER :: zero = 0.0E0,  one = 1.0E0, two = 2.0E0
  REAL, PARAMETER :: gam = 4096.0E0, gamsq = 1.67772E7, rgamsq = 5.96046E-8
  !* FIRST EXECUTABLE STATEMENT  SROTMG
  IF ( .NOT.Sd1<zero ) THEN
    !     CASE-SD1-NONNEGATIVE
    sp2 = Sd2*Sy1
    IF ( .NOT.sp2==zero ) THEN
      !     REGULAR-CASE..
      sp1 = Sd1*Sx1
      sq2 = sp2*Sy1
      sq1 = sp1*Sx1
      !
      IF ( .NOT.(.NOT.ABS(sq1)>ABS(sq2)) ) THEN
        sh21 = -Sy1/Sx1
        sh12 = sp2/sp1
        !
        su = one - sh12*sh21
        !
        IF ( .NOT.su<=zero ) THEN
          sflag = zero
          Sd1 = Sd1/su
          Sd2 = Sd2/su
          Sx1 = Sx1*su
          !         GO SCALE-CHECK..
          GOTO 200
          !         GO ZERO-H-D-AND-SX1..
        ENDIF
      ELSEIF ( .NOT.sq2<zero ) THEN
        sflag = one
        sh11 = sp1/sp2
        sh22 = Sx1/Sy1
        su = one + sh11*sh22
        stemp = Sd2/su
        Sd2 = Sd1/su
        Sd1 = stemp
        Sx1 = Sy1*su
        !         GO SCALE-CHECK
        GOTO 200
        !         GO ZERO-H-D-AND-SX1..
      ENDIF
    ELSE
      sflag = -two
      Sparam(1) = sflag
      RETURN
    ENDIF
    !       GO ZERO-H-D-AND-SX1..
  ENDIF
  !     PROCEDURE..ZERO-H-D-AND-SX1..
  sflag = -one
  sh11 = zero
  sh12 = zero
  sh21 = zero
  sh22 = zero
  !
  Sd1 = zero
  Sd2 = zero
  Sx1 = zero
  !         RETURN..
  GOTO 1000
  !     PROCEDURE..FIX-H..
  100 CONTINUE
  IF ( .NOT.(.NOT.sflag>=zero) ) THEN
    !
    IF ( .NOT.sflag==zero ) THEN
      sh21 = -one
      sh12 = one
      sflag = -one
    ELSE
      sh11 = one
      sh22 = one
      sflag = -one
    ENDIF
  ENDIF
  SELECT CASE(igo)
    CASE(300)
      GOTO 300
    CASE(500)
      GOTO 500
    CASE(700)
      GOTO 700
    CASE(900)
      GOTO 900
  END SELECT
  !     PROCEDURE..SCALE-CHECK
  200 CONTINUE
  IF ( .NOT.Sd1<=rgamsq ) GOTO 400
  IF ( Sd1==zero ) GOTO 600
  igo = 300
  !              FIX-H..
  GOTO 100
  300  Sd1 = Sd1*gam**2
  Sx1 = Sx1/gam
  sh11 = sh11/gam
  sh12 = sh12/gam
  GOTO 200
  400 CONTINUE
  IF ( .NOT.Sd1>=gamsq ) GOTO 600
  igo = 500
  !              FIX-H..
  GOTO 100
  500  Sd1 = Sd1/gam**2
  Sx1 = Sx1*gam
  sh11 = sh11*gam
  sh12 = sh12*gam
  GOTO 400
  600 CONTINUE
  IF ( .NOT.ABS(Sd2)<=rgamsq ) GOTO 800
  IF ( Sd2==zero ) GOTO 1000
  igo = 700
  !              FIX-H..
  GOTO 100
  700  Sd2 = Sd2*gam**2
  sh21 = sh21/gam
  sh22 = sh22/gam
  GOTO 600
  800 CONTINUE
  IF ( .NOT.ABS(Sd2)>=gamsq ) GOTO 1000
  igo = 900
  !              FIX-H..
  GOTO 100
  900  Sd2 = Sd2/gam**2
  sh21 = sh21*gam
  sh22 = sh22*gam
  GOTO 800
  1000 CONTINUE
  IF ( sflag<0 ) THEN
    Sparam(2) = sh11
    Sparam(3) = sh21
    Sparam(4) = sh12
    Sparam(5) = sh22
  ELSEIF ( sflag==0 ) THEN
    Sparam(3) = sh21
    Sparam(4) = sh12
  ELSE
    Sparam(2) = sh11
    Sparam(5) = sh22
  ENDIF
  Sparam(1) = sflag
  RETURN
END SUBROUTINE SROTMG
