!** BESKS
SUBROUTINE BESKS(Xnu,X,Nin,Bk)
  !>
  !***
  !  Compute a sequence of modified Bessel functions of the
  !            third kind of fractional order.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      SINGLE PRECISION (BESKS-S, DBESKS-D)
  !***
  ! **Keywords:**  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
  !             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESKS computes a sequence of modified Bessel functions of the third
  ! kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1),
  ! and I = 0, 1, ..., NIN - 1, if NIN is positive and I = 0, 1, ... ,
  ! NIN + 1, if NIN is negative.  On return, the vector BK(.) Contains
  ! the results at X for order starting at XNU.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESKES, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL Bk(*), expxi, X, Xnu
  INTEGER i, n, Nin
  REAL :: xmax = 0.0
  !* FIRST EXECUTABLE STATEMENT  BESKS
  IF ( xmax==0.0 ) xmax = -LOG(R1MACH(1))
  !
  IF ( X>xmax ) CALL XERMSG('SLATEC','BESKS','X SO BIG BESSEL K UNDERFLOWS',&
    1,2)
  !
  CALL BESKES(Xnu,X,Nin,Bk)
  !
  expxi = EXP(-X)
  n = ABS(Nin)
  DO i = 1, n
    Bk(i) = expxi*Bk(i)
  END DO
  !
END SUBROUTINE BESKS
