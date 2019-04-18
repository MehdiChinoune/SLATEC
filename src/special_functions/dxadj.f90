!** DXADJ
SUBROUTINE DXADJ(X,Ix,Ierror)
  !>
  !  To provide double-precision floating-point arithmetic
  !            with an extended exponent range.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  A3D
  !***
  ! **Type:**      DOUBLE PRECISION (XADJ-S, DXADJ-D)
  !***
  ! **Keywords:**  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !     DOUBLE PRECISION X
  !     INTEGER IX
  !
  !                  TRANSFORMS (X,IX) SO THAT
  !                  RADIX**(-L) .LE. ABS(X) .LT. RADIX**L.
  !                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
  !                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
  !                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
  !
  !***
  ! **See also:**  DXSET
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG
  !***
  ! COMMON BLOCKS    DXBLK2

  !* REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  USE DXBLK ,ONLY: RADixl, RAD2l, L2, KMAx
  USE service, ONLY : XERMSG
  INTEGER Ierror, Ix
  REAL(8) :: X
  !
  !   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! IS
  !     2*L .LE. KMAX
  !
  ! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE DXSET.
  !
  !* FIRST EXECUTABLE STATEMENT  DXADJ
  Ierror = 0
  IF ( X==0.0D0 ) THEN
    Ix = 0
    GOTO 200
  ELSEIF ( ABS(X)>=1.0D0 ) THEN
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
    IF ( RADixl*ABS(X)>=1.0D0 ) GOTO 200
    X = X*RAD2l
    IF ( Ix>=0 ) THEN
      Ix = Ix - L2
      RETURN
    ELSEIF ( Ix>=-KMAx+L2 ) THEN
      Ix = Ix - L2
      RETURN
    END IF
  END IF
  100  CALL XERMSG('SLATEC','DXADJ','overflow in auxiliary index',207,1)
  Ierror = 207
  RETURN
  200 CONTINUE
  IF ( ABS(Ix)>KMAx ) GOTO 100
  RETURN
END SUBROUTINE DXADJ
