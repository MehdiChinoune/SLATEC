!*==CATAN2.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK CATAN2
      COMPLEX FUNCTION CATAN2(Csn,Ccs)
      IMPLICIT NONE
!*--CATAN25
!*** Start of declarations inserted by SPAG
      REAL pi
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CATAN2
!***PURPOSE  Compute the complex arc tangent in the proper quadrant.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (CATAN2-C)
!***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, POLAR ANGEL,
!             QUADRANT, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CATAN2(CSN,CCS) calculates the complex trigonometric arc
! tangent of the ratio CSN/CCS and returns a result whose real
! part is in the correct quadrant (within a multiple of 2*PI).  The
! result is in units of radians and the real part is between -PI
! and +PI.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CATAN, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  CATAN2
      COMPLEX Csn , Ccs , CATAN
      SAVE pi
      DATA pi/3.14159265358979323846E0/
!***FIRST EXECUTABLE STATEMENT  CATAN2
      IF ( ABS(Ccs)==0. ) THEN
!
        IF ( ABS(Csn)==0. ) CALL XERMSG('SLATEC','CATAN2',
     &                                  'CALLED WITH BOTH ARGUMENTS ZERO',1,2)
!
        CATAN2 = CMPLX(SIGN(0.5*pi,REAL(Csn)),0.0)
        GOTO 99999
      ENDIF
!
      CATAN2 = CATAN(Csn/Ccs)
      IF ( REAL(Ccs)<0. ) CATAN2 = CATAN2 + pi
      IF ( REAL(CATAN2)>pi ) CATAN2 = CATAN2 - 2.0*pi
      RETURN
!
99999 END FUNCTION CATAN2
