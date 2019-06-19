!** XADJ
SUBROUTINE XADJ(X,Ix,Ierror)
  !> To provide single-precision floating-point arithmetic
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
  !                  RADIX**(-L) <= ABS(X) < RADIX**L.
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
  USE XBLK ,ONLY: radixl_com, rad2l_com, l2_com, kmax_com
  USE service, ONLY : XERMSG
  INTEGER :: Ierror, Ix
  REAL(SP) :: X
  !
  !   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! IS
  !     2*L <= KMAX
  !
  ! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE XSET.
  !
  !* FIRST EXECUTABLE STATEMENT  XADJ
  Ierror = 0
  IF( X==0.0 ) THEN
    Ix = 0
    GOTO 200
  ELSEIF( ABS(X)>=1.0 ) THEN
    IF( ABS(X)<radixl_com ) GOTO 200
    X = X/rad2l_com
    IF( Ix<=0 ) THEN
      Ix = Ix + l2_com
      RETURN
    ELSEIF( Ix<=kmax_com-l2_com ) THEN
      Ix = Ix + l2_com
      RETURN
    END IF
  ELSE
    IF( radixl_com*ABS(X)>=1.0 ) GOTO 200
    X = X*rad2l_com
    IF( Ix>=0 ) THEN
      Ix = Ix - l2_com
      RETURN
    ELSEIF( Ix>=-kmax_com+l2_com ) THEN
      Ix = Ix - l2_com
      RETURN
    END IF
  END IF
  100  CALL XERMSG('XADJ','overflow in auxiliary index',107,1)
  Ierror = 107
  RETURN
  200 CONTINUE
  IF( ABS(Ix)>kmax_com ) GOTO 100
  RETURN
END SUBROUTINE XADJ
