!*==XADJ.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK XADJ
      SUBROUTINE XADJ(X,Ix,Ierror)
      IMPLICIT NONE
!*--XADJ5
!*** Start of declarations inserted by SPAG
      INTEGER Ierror
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  XADJ
!***PURPOSE  To provide single-precision floating-point arithmetic
!            with an extended exponent range.
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      SINGLE PRECISION (XADJ-S, DXADJ-D)
!***KEYWORDS  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!     REAL X
!     INTEGER IX
!
!                  TRANSFORMS (X,IX) SO THAT
!                  RADIX**(-L) .LE. ABS(X) .LT. RADIX**L.
!                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
!                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
!                  THE NUMBER BASE OF SINGLE-PRECISION ARITHMETIC.
!
!***SEE ALSO  XSET
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***COMMON BLOCKS    XBLK2
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  XADJ
      REAL X
      INTEGER Ix
      REAL RADix , RADixl , RAD2l , DLG10r
      INTEGER L , L2 , KMAx
      COMMON /XBLK2 / RADix , RADixl , RAD2l , DLG10r , L , L2 , KMAx
      SAVE /XBLK2 / 
!
!   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
! IS
!     2*L .LE. KMAX
!
! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE XSET.
!
!***FIRST EXECUTABLE STATEMENT  XADJ
      Ierror = 0
      IF ( X==0.0 ) THEN
        Ix = 0
        GOTO 200
      ELSEIF ( ABS(X)>=1.0 ) THEN
        IF ( ABS(X)<RADixl ) GOTO 200
        X = X/RAD2l
        IF ( Ix<=0 ) THEN
          Ix = Ix + L2
          GOTO 99999
        ELSEIF ( Ix<=KMAx-L2 ) THEN
          Ix = Ix + L2
          GOTO 99999
        ENDIF
      ELSE
        IF ( RADixl*ABS(X)>=1.0 ) GOTO 200
        X = X*RAD2l
        IF ( Ix>=0 ) THEN
          Ix = Ix - L2
          GOTO 99999
        ELSEIF ( Ix>=-KMAx+L2 ) THEN
          Ix = Ix - L2
          GOTO 99999
        ENDIF
      ENDIF
 100  CALL XERMSG('SLATEC','XADJ','overflow in auxiliary index',107,1)
      Ierror = 107
      RETURN
 200  IF ( ABS(Ix)>KMAx ) GOTO 100
99999 END SUBROUTINE XADJ
