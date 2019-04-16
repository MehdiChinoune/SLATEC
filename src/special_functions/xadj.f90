!** XADJ
SUBROUTINE XADJ(X,Ix,Ierror)
  USE XBLK ,ONLY: RADixl, RAD2l, L2, KMAx
  !>
  !***
  !  To provide single-precision floating-point arithmetic
  !            with an extended exponent range.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  A3D
  !***
  ! **Type:**      SINGLE PRECISION (XADJ-S, DXADJ-D)
  !***
  ! **Keywords:**  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !     REAL X
  !     INTEGER IX
  !
  !                  TRANSFORMS (X,IX) SO THAT
  !                  RADIX**(-L) .LE. ABS(X) .LT. RADIX**L.
  !                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
  !                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
  !                  THE NUMBER BASE OF SINGLE-PRECISION ARITHMETIC.
  !
  !***
  ! **See also:**  XSET
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG
  !***
  ! COMMON BLOCKS    XBLK2

  !* REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER Ierror, Ix
  REAL X
  !
  !   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! IS
  !     2*L .LE. KMAX
  !
  ! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE XSET.
  !
  !* FIRST EXECUTABLE STATEMENT  XADJ
  Ierror = 0
  IF ( X==0.0 ) THEN
    Ix = 0
    GOTO 200
  ELSEIF ( ABS(X)>=1.0 ) THEN
    IF ( ABS(X)<RADixl ) GOTO 200
    X = X/RAD2l
    IF ( Ix<=0 ) THEN
      Ix = Ix + L2
      RETURN
    ELSEIF ( Ix<=KMAx-L2 ) THEN
      Ix = Ix + L2
      RETURN
    END IF
  ELSE
    IF ( RADixl*ABS(X)>=1.0 ) GOTO 200
    X = X*RAD2l
    IF ( Ix>=0 ) THEN
      Ix = Ix - L2
      RETURN
    ELSEIF ( Ix>=-KMAx+L2 ) THEN
      Ix = Ix - L2
      RETURN
    END IF
  END IF
  100  CALL XERMSG('SLATEC','XADJ','overflow in auxiliary index',107,1)
  Ierror = 107
  RETURN
  200 CONTINUE
  IF ( ABS(Ix)>KMAx ) GOTO 100
  RETURN
END SUBROUTINE XADJ
