!*==BESKS.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK BESKS
SUBROUTINE BESKS(Xnu,X,Nin,Bk)
  IMPLICIT NONE
  !*--BESKS5
  !*** Start of declarations inserted by SPAG
  REAL Bk , expxi , R1MACH , X , xmax , Xnu
  INTEGER i , n , Nin
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BESKS
  !***PURPOSE  Compute a sequence of modified Bessel functions of the
  !            third kind of fractional order.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B3
  !***TYPE      SINGLE PRECISION (BESKS-S, DBESKS-D)
  !***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
  !             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! BESKS computes a sequence of modified Bessel functions of the third
  ! kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1),
  ! and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... ,
  ! NIN + 1, if NIN is negative.  On return, the vector BK(.) Contains
  ! the results at X for order starting at XNU.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  BESKES, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  BESKS
  DIMENSION Bk(*)
  SAVE xmax
  DATA xmax/0.0/
  !***FIRST EXECUTABLE STATEMENT  BESKS
  IF ( xmax==0.0 ) xmax = -LOG(R1MACH(1))
  !
  IF ( X>xmax ) CALL XERMSG('SLATEC','BESKS','X SO BIG BESSEL K UNDERFLOWS',&
    1,2)
  !
  CALL BESKES(Xnu,X,Nin,Bk)
  !
  expxi = EXP(-X)
  n = ABS(Nin)
  DO i = 1 , n
    Bk(i) = expxi*Bk(i)
  ENDDO
  !
END SUBROUTINE BESKS
