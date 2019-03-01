!DECK DROTMG
SUBROUTINE DROTMG(Dd1,Dd2,Dx1,Dy1,Dparam)
  IMPLICIT NONE
  INTEGER igo
  !***BEGIN PROLOGUE  DROTMG
  !***PURPOSE  Construct a modified Givens transformation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B10
  !***TYPE      DOUBLE PRECISION (SROTMG-S, DROTMG-D)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !      DD1  double precision scalar
  !      DD2  double precision scalar
  !      DX1  double precision scalar
  !      DX2  double precision scalar
  !   DPARAM  D.P. 5-vector. DPARAM(1)=DFLAG defined below.
  !           Locations 2-5 contain the rotation matrix.
  !
  !     --Output--
  !      DD1  changed to represent the effect of the transformation
  !      DD2  changed to represent the effect of the transformation
  !      DX1  changed to represent the effect of the transformation
  !      DX2  unchanged
  !
  !     Construct the modified Givens transformation matrix H which zeros
  !     the second component of the 2-vector  (SQRT(DD1)*DX1,SQRT(DD2)*
  !     DY2)**T.
  !     With DPARAM(1)=DFLAG, H has one of the following forms:
  !
  !     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
  !
  !       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
  !     H=(          )    (          )    (          )    (          )
  !       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
  !
  !     Locations 2-5 of DPARAM contain DH11, DH21, DH12, and DH22,
  !     respectively.  (Values of 1.D0, -1.D0, or 0.D0 implied by the
  !     value of DPARAM(1) are not stored in DPARAM.)
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   780301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920316  Prologue corrected.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DROTMG
  REAL(8) :: gam, one, rgamsq, Dd1, Dd2, dh11, dh12, dh21, &
    dh22, Dparam, dp1, dp2, dq1, dq2, du, Dy1, zero, &
    gamsq, dflag, dtemp, Dx1, two
  DIMENSION Dparam(5)
  SAVE zero, one, two, gam, gamsq, rgamsq
  DATA zero, one, two/0.0D0, 1.0D0, 2.0D0/
  DATA gam, gamsq, rgamsq/4096.0D0, 16777216.D0, 5.9604645D-8/
  !***FIRST EXECUTABLE STATEMENT  DROTMG
  IF ( .NOT.Dd1<zero ) THEN
    !     CASE-DD1-NONNEGATIVE
    dp2 = Dd2*Dy1
    IF ( .NOT.dp2==zero ) THEN
      !     REGULAR-CASE..
      dp1 = Dd1*Dx1
      dq2 = dp2*Dy1
      dq1 = dp1*Dx1
      !
      IF ( .NOT.(.NOT.ABS(dq1)>ABS(dq2)) ) THEN
        dh21 = -Dy1/Dx1
        dh12 = dp2/dp1
        !
        du = one - dh12*dh21
        !
        IF ( .NOT.du<=zero ) THEN
          dflag = zero
          Dd1 = Dd1/du
          Dd2 = Dd2/du
          Dx1 = Dx1*du
          !         GO SCALE-CHECK..
          GOTO 200
          !         GO ZERO-H-D-AND-DX1..
        ENDIF
      ELSEIF ( .NOT.dq2<zero ) THEN
        dflag = one
        dh11 = dp1/dp2
        dh22 = Dx1/Dy1
        du = one + dh11*dh22
        dtemp = Dd2/du
        Dd2 = Dd1/du
        Dd1 = dtemp
        Dx1 = Dy1*du
        !         GO SCALE-CHECK
        GOTO 200
        !         GO ZERO-H-D-AND-DX1..
      ENDIF
    ELSE
      dflag = -two
      Dparam(1) = dflag
      RETURN
    ENDIF
    !       GO ZERO-H-D-AND-DX1..
  ENDIF
  !     PROCEDURE..ZERO-H-D-AND-DX1..
  dflag = -one
  dh11 = zero
  dh12 = zero
  dh21 = zero
  dh22 = zero
  !
  Dd1 = zero
  Dd2 = zero
  Dx1 = zero
  !         RETURN..
  GOTO 1000
  !     PROCEDURE..FIX-H..
  100 CONTINUE
  IF ( .NOT.(.NOT.dflag>=zero) ) THEN
    !
    IF ( .NOT.dflag==zero ) THEN
      dh21 = -one
      dh12 = one
      dflag = -one
    ELSE
      dh11 = one
      dh22 = one
      dflag = -one
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
  IF ( .NOT.Dd1<=rgamsq ) GOTO 400
  IF ( Dd1==zero ) GOTO 600
  igo = 300
  !              FIX-H..
  GOTO 100
  300  Dd1 = Dd1*gam**2
  Dx1 = Dx1/gam
  dh11 = dh11/gam
  dh12 = dh12/gam
  GOTO 200
  400 CONTINUE
  IF ( .NOT.Dd1>=gamsq ) GOTO 600
  igo = 500
  !              FIX-H..
  GOTO 100
  500  Dd1 = Dd1/gam**2
  Dx1 = Dx1*gam
  dh11 = dh11*gam
  dh12 = dh12*gam
  GOTO 400
  600 CONTINUE
  IF ( .NOT.ABS(Dd2)<=rgamsq ) GOTO 800
  IF ( Dd2==zero ) GOTO 1000
  igo = 700
  !              FIX-H..
  GOTO 100
  700  Dd2 = Dd2*gam**2
  dh21 = dh21/gam
  dh22 = dh22/gam
  GOTO 600
  800 CONTINUE
  IF ( .NOT.ABS(Dd2)>=gamsq ) GOTO 1000
  igo = 900
  !              FIX-H..
  GOTO 100
  900  Dd2 = Dd2/gam**2
  dh21 = dh21*gam
  dh22 = dh22*gam
  GOTO 800
  1000 CONTINUE
  IF ( dflag<0 ) THEN
    Dparam(2) = dh11
    Dparam(3) = dh21
    Dparam(4) = dh12
    Dparam(5) = dh22
  ELSEIF ( dflag==0 ) THEN
    Dparam(3) = dh21
    Dparam(4) = dh12
  ELSE
    Dparam(2) = dh11
    Dparam(5) = dh22
  ENDIF
  Dparam(1) = dflag
  RETURN
END SUBROUTINE DROTMG
